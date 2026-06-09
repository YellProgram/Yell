/*
 Parallelism timing tests for Calculator.h FFT accumulation and Jacobian evaluation.

 Tests forward-calc and derivative-calc at MaxProcessors=1 vs MaxProcessors=10 and
 asserts a minimum speedup.  Run separately to avoid slowing the unit-test suite:

   cmake-build-debug/test-yell-parallel
   cmake-build-debug/test-yell-parallel --gtest_filter="ParallelTiming.*"

 The model has 10 correlations in m-3m symmetry → ~130 pairs at 10 unique r_grids
 (so 10 per-r_grid mutexes, near-zero contention with 10 threads).  The output map
 is 160^3 (64 MB) to exceed typical L3 cache and produce real cache-miss pressure on
 add_pair_to_appropriate_place, making the mutex timing visible.

 TODO: also add a wall-time test that drives Model through a full CeresMinimizer
       round using the parallel.txt file from tests/integration/models/ as a
       reference to validate that the real refinement path is fast end-to-end.
*/

#include <gtest/gtest.h>

#include "basic_classes.h"
#include "model.h"
#include "IntensityMap.h"

#include <chrono>
#include <thread>
#include <atomic>
#include <vector>
#include <cstdio>

// Required global — defined in main.cpp for the binary.
OutputHandler report;

using Clock    = std::chrono::steady_clock;
using Seconds  = std::chrono::duration<double>;

// ─────────────────────────────────────────────────────────────────────────────
// Model with ~130 pairs across 10 unique r_grids.
//
// m-3m symmetry expands each (hkl) shell to its full multiplicity:
//   (1,0,0)x6  (1,1,0)x12  (1,1,1)x8  (2,0,0)x6  (2,1,0)x24
//   (2,1,1)x24 (2,2,0)x12  (2,2,1)x24 (2,2,2)x8  (3,0,0)x6  → 130 pairs
// 10 independent refinable parameters → one per shell.
//
// 160^3 output grid (64 MB) exceeds typical L3 cache so
// add_pair_to_appropriate_place incurs real cache-miss latency.
// ─────────────────────────────────────────────────────────────────────────────
static const char* PARALLEL_MODEL = R"(
Cell 4.0 4.0 4.0  90 90 90
LaueSymmetry m-3m
DiffuseScatteringGrid -8 -8 -8 0.1 0.1 0.1 160 160 160

CalculationMethod fft
FFTGridSize 32 32 32
Refine false

RefinableVariables [
  p1=0.05 p2=0.04 p3=0.03 p4=0.02 p5=0.01
  p6=0.01 p7=0.02 p8=0.03 p9=0.04 p10=0.05
]

UnitCell [
  V = Variant [
    (p=0.5) Au 1  0 0 0  0.003
    (p=0.5) Void
  ]
]

Modes []

Correlations [
  [(1,0,0)   Multiplicity 6  SubstitutionalCorrelation(V,V,0.5+p1)]
  [(1,1,0)   Multiplicity 12 SubstitutionalCorrelation(V,V,0.5+p2)]
  [(1,1,1)   Multiplicity 8  SubstitutionalCorrelation(V,V,0.5+p3)]
  [(2,0,0)   Multiplicity 6  SubstitutionalCorrelation(V,V,0.5+p4)]
  [(2,1,0)   Multiplicity 24 SubstitutionalCorrelation(V,V,0.5+p5)]
  [(2,1,1)   Multiplicity 24 SubstitutionalCorrelation(V,V,0.5+p6)]
  [(2,2,0)   Multiplicity 12 SubstitutionalCorrelation(V,V,0.5+p7)]
  [(2,2,1)   Multiplicity 24 SubstitutionalCorrelation(V,V,0.5+p8)]
  [(2,2,2)   Multiplicity 8  SubstitutionalCorrelation(V,V,0.5+p9)]
  [(3,0,0)   Multiplicity 6  SubstitutionalCorrelation(V,V,0.5+p10)]
]
)";

// params[0]=Scale, params[1..10]=p1..p10
static const std::vector<double> PARAMS = {
    1.0,
    0.05, 0.04, 0.03, 0.02, 0.01,
    0.01, 0.02, 0.03, 0.04, 0.05
};

// ─────────────────────────────────────────────────────────────────────────────
// Helper: time N_RUNS calls of m.calculate() and return average seconds.
// ─────────────────────────────────────────────────────────────────────────────
static double time_forward(Model& m, int n_runs = 3)
{
    auto t0 = Clock::now();
    for (int i = 0; i < n_runs; ++i)
        m.calculate(PARAMS);
    return Seconds(Clock::now() - t0).count() / n_runs;
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper: time computing all structural derivatives sequentially (1 thread per
// derivative call).  Mirrors the single-threaded Jacobian path.
// ─────────────────────────────────────────────────────────────────────────────
static double time_derivatives_sequential(Model& m)
{
    m.calculate(PARAMS);
    int n_params = (int)PARAMS.size();  // includes Scale at index 0

    auto t0 = Clock::now();
    for (int j = 0; j < n_params; ++j)
        m.calculate_derivative(PARAMS, j, /*num_threads=*/1);
    return Seconds(Clock::now() - t0).count();
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper: time computing all structural derivatives in parallel using n_threads
// Model clones and work-stealing (mirrors CeresMinimizer Jacobian path).
// ─────────────────────────────────────────────────────────────────────────────
static double time_derivatives_parallel(Model& m, int n_threads)
{
    m.calculate(PARAMS);
    int n_params = (int)PARAMS.size();

    std::vector<Model*> clones(n_threads);
    for (int t = 0; t < n_threads; ++t)
        clones[t] = m.clone();

    std::atomic<int> next_j(0);
    std::vector<std::thread> workers;

    auto t0 = Clock::now();
    for (int t = 0; t < n_threads; ++t) {
        workers.emplace_back([&, t]() {
            Model* c = clones[t];
            while (true) {
                int j = next_j.fetch_add(1);
                if (j >= n_params) break;
                c->calculate_derivative(PARAMS, j, /*num_threads=*/1);
            }
        });
    }
    for (auto& w : workers) w.join();
    double elapsed = Seconds(Clock::now() - t0).count();

    for (auto* c : clones) delete c;
    return elapsed;
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

TEST(ParallelTiming, ForwardCalcSpeedup)
{
    const int N_THREADS = 10;

    // 1-thread run
    Model m1(PARALLEL_MODEL);
    m1.set_max_processors(1);
    m1.calculate(PARAMS);          // warmup
    double t1 = time_forward(m1);

    // N-thread run
    Model mN(PARALLEL_MODEL);
    mN.set_max_processors(N_THREADS);
    mN.calculate(PARAMS);          // warmup
    double tN = time_forward(mN);

    double speedup = t1 / tN;
    std::printf("\n[ForwardCalc]  1-thread: %.3f s  |  %d-thread: %.3f s  |  speedup: %.2fx\n",
                t1, N_THREADS, tN, speedup);

    // Expect at least 1.3x speedup from 10 threads.
    // Note: debug builds show lower speedup (~1.6x) due to unoptimised FFT
    // and cctbx internal overhead; release builds should reach 5x+.
    // If this fails even at 1.3x, threads are truly serialised (mutex or FFT global state).
    EXPECT_GT(speedup, 1.3)
        << "Forward calc speedup " << speedup << "x < 1.3x with " << N_THREADS
        << " threads.  Likely mutex contention in FFT accumulation.";
}

TEST(ParallelTiming, DerivativeCalcSpeedup)
{
    const int N_THREADS = 10;

    // Sequential: all derivatives on one thread
    Model m1(PARALLEL_MODEL);
    m1.set_max_processors(1);
    double t1 = time_derivatives_sequential(m1);

    // Parallel: N threads, each taking one parameter at a time
    Model mN(PARALLEL_MODEL);
    mN.set_max_processors(N_THREADS);
    double tN = time_derivatives_parallel(mN, N_THREADS);

    double speedup = t1 / tN;
    std::printf("\n[Derivatives]  sequential: %.3f s  |  %d-thread: %.3f s  |  speedup: %.2fx\n",
                t1, N_THREADS, tN, speedup);

    // With 10 parameters and 10 threads the ideal speedup is 10x.
    // Accept >= 3x (30% efficiency accounts for clone overhead + OS scheduling).
    EXPECT_GT(speedup, 3.0)
        << "Derivative speedup " << speedup << "x < 3x with " << N_THREADS
        << " threads.  Check that CeresMinimizer uses the flat work-queue.";
}

TEST(ParallelTiming, PrintHardwareConcurrency)
{
    std::cout << "[Hardware] std::thread::hardware_concurrency() = " << std::thread::hardware_concurrency() << std::endl;
}

// ─────────────────────────────────────────────────────────────────────────────
// Gram–Charlier refinement recovery (Phase 6): generate synthetic diffuse data
// from a known 4th-order coefficient, refine from zero via finite differences,
// and check the value (and Scale) come back. Needs the Ceres minimizer.
// ─────────────────────────────────────────────────────────────────────────────

#include "CeresMinimizer.h"

static std::string gc_refine_model(double d_init)
{
    std::ostringstream oss;
    oss << "Cell 5 5 5  90 90 90\n"
        << "DiffuseScatteringGrid -3 -3 -3  1 1 1  6 6 6\n"   // -1-compatible
        << "CalculationMethod direct\n"
        << "LaueSymmetry -1\n"
        << "RefinableVariables [ d = " << d_init << " ]\n"
        << "UnitCell [\n"
        << "  V = Variant [ (p=0.5) Na = Na 1 0 0 0  0.03 0.03 0.03 0 0 0 "
        << "GramCharlier4[ d 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ] (p=0.5) Void ]\n"
        << "]\n"
        << "Correlations [\n"
        << "  [ (1,0,0) SubstitutionalCorrelation(V,V,0.5) ]\n"
        << "  [ (2,0,0) SubstitutionalCorrelation(V,V,0.4) ]\n"
        << "]\n";
    return oss.str();
}

TEST(GramCharlierRefine, RecoversFourthOrderFromSyntheticData)
{
    const double d_true = 0.004;

    Model truth(gc_refine_model(d_true));
    truth.calculate({1.0, d_true});

    // synthetic diffuse data = full − average (deep copy; IntensityMap copy is shallow)
    IntensityMap Ie(truth.intensity_map.size());
    Ie.set_grid(truth.intensity_map.grid);
    for (int i = 0; i < Ie.size_1d(); ++i)
        Ie.at(i) = truth.intensity_map.at(i) - truth.average_intensity_map.at(i);

    Model fit(gc_refine_model(0.0));   // start the coefficient at zero
    fit.init_asu();                    // populate ASU indices (main.cpp does this pre-refine)
    fit.set_derivatives_mode(FINITE_DIFFERENCE);
    OptionalIntensityMap weights;      // default value 1 everywhere

    CeresMinimizer minimizer;
    vector<double> refined = minimizer.minimize({1.0, 0.0}, &Ie, &fit, &weights,
                                                fit.refinement_options);

    ASSERT_EQ(refined.size(), 2u);
    EXPECT_NEAR(refined[0], 1.0,    1e-3) << "Scale";
    EXPECT_NEAR(refined[1], d_true, 1e-4) << "d1111";
}
