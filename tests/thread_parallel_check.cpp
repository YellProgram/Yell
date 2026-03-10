// Standalone threading sanity check.
//
// Spawns N threads each doing independent heavy floating-point work (volatile
// sin/cos accumulation that cannot be optimised away), measures wall time with
// 1 thread and N threads, and exits with a non-zero status if the speedup is
// below a threshold.
//
// Build (no Yell deps needed):
//   g++ -std=c++20 -O2 -pthread tests/thread_parallel_check.cpp -o thread_parallel_check
//   clang++ -std=c++20 -O2 -pthread tests/thread_parallel_check.cpp -o thread_parallel_check
//
// Usage:
//   ./thread_parallel_check [n_threads]   (default: hardware_concurrency)
//
// Exit code 0 = good speedup, 1 = parallelism broken.
//
// TODO: a "parallel.txt" model with ~100 pairs in tests/integration/ would let
// us run the same timing check through the real Yell code path (forward calc +
// derivative calc), validating that the per-r_grid mutex approach actually
// parallelises on the target machine.  Keep this in mind for an integration
// timing test.

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <thread>
#include <vector>
#include <chrono>
#include <atomic>
#include <numeric>

// Heavy work: sum sin(i * 1e-7) for N_ITER iterations.
// volatile prevents the compiler from eliminating the loop.
static double heavy_work(long n_iter)
{
    volatile double acc = 0.0;
    for (long i = 0; i < n_iter; ++i)
        acc += std::sin(i * 1e-7) * std::cos(i * 1.3e-7);
    return acc;
}

using Clock = std::chrono::steady_clock;

// Run `n_jobs` independent work units across `n_threads` threads.
// Returns elapsed wall-clock seconds.
static double run(int n_threads, int n_jobs, long iter_per_job)
{
    std::atomic<int> next(0);
    std::vector<double> results(n_jobs, 0.0);

    auto worker = [&]() {
        while (true) {
            int j = next.fetch_add(1);
            if (j >= n_jobs) break;
            results[j] = heavy_work(iter_per_job);
        }
    };

    auto t0 = Clock::now();
    if (n_threads <= 1) {
        worker();
    } else {
        std::vector<std::thread> workers;
        workers.reserve(n_threads);
        for (int t = 0; t < n_threads; ++t)
            workers.emplace_back(worker);
        for (auto& w : workers) w.join();
    }
    double elapsed = std::chrono::duration<double>(Clock::now() - t0).count();

    // Use result so the compiler can't eliminate the work.
    volatile double sum = std::accumulate(results.begin(), results.end(), 0.0);
    (void)sum;
    return elapsed;
}

int main(int argc, char** argv)
{
    int n_threads = (int)std::thread::hardware_concurrency();
    if (n_threads <= 0) n_threads = 4;
    if (argc >= 2) n_threads = std::atoi(argv[1]);

    // Calibrate: make single-thread run take ~1 second.
    const int    n_jobs        = n_threads * 4;  // enough jobs to keep all threads busy
    const long   iter_per_job  = 20'000'000L;

    std::printf("Threads available : %d\n", (int)std::thread::hardware_concurrency());
    std::printf("Threads tested    : %d\n", n_threads);
    std::printf("Jobs              : %d  x  %ld iterations\n", n_jobs, iter_per_job);

    // Warm up
    run(1, 1, iter_per_job / 10);

    double t1 = run(1,         n_jobs, iter_per_job);
    double tN = run(n_threads, n_jobs, iter_per_job);
    double speedup = t1 / tN;

    std::printf("Single-thread     : %.3f s\n", t1);
    std::printf("%-2d-thread         : %.3f s\n", n_threads, tN);
    std::printf("Speedup           : %.2fx  (ideal %.1fx)\n", speedup, (double)n_threads);

    // Accept anything >= 50% of ideal (leaves room for OS scheduling noise).
    double min_speedup = 0.5 * n_threads;
    if (n_threads == 1) min_speedup = 0.9;

    if (speedup >= min_speedup) {
        std::printf("PASS: speedup %.2fx >= threshold %.2fx\n", speedup, min_speedup);
        return 0;
    } else {
        std::fprintf(stderr,
            "FAIL: speedup %.2fx < threshold %.2fx -- threads are not running in parallel!\n",
            speedup, min_speedup);
        return 1;
    }
}
