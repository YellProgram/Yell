/*
 Copyright Arkadiy Simonov, Thomas Weber, ETH Zurich 2014

 This file is part of Yell.

 Yell is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Yell is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with Yell.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef YELL_CALCULATOR_H
#define YELL_CALCULATOR_H

#include "AtomicPairs.h"
#include "IntensityMap.h"
#include "diffuser_core.h"

#include <vector>
#include <complex>
#include <thread>
#include <mutex>
#include <atomic>

using namespace std;
using namespace scitbx;

// ─────────────────────────────────────────────────────────────────────────────
// MinimizerCalculator — interface for the minimizer
// ─────────────────────────────────────────────────────────────────────────────

class MinimizerCalculator {
public:
    virtual void calculate(vector<double> p) = 0;
    virtual IntensityMap& data() = 0;
    virtual IntensityMap& get_intensity_map() = 0;
    virtual IntensityMap& get_average_intensity_map() = 0;
    virtual int number_of_observations() = 0;
    virtual bool refine_in_asu() = 0;
    virtual vector<int>& asu_indices() = 0;
};

// ─────────────────────────────────────────────────────────────────────────────
// RefinementOptions
// ─────────────────────────────────────────────────────────────────────────────

struct RefinementOptions {
    int    max_number_of_iterations;
    double tau;
    double thresholds[3];
    double difference;

    // Ceres-specific options
    // num_threads: enables Eigen OpenMP parallelism inside Ceres's QR/Cholesky
    // factorisation.  Safe to set = max_processors because we have a single
    // residual block — Ceres cannot launch parallel block evaluations, so the
    // only effect of num_threads is on Eigen's internal linear algebra threads.
    int    num_threads;          // 0 = use max_processors
    double function_tolerance;   // relative decrease in cost (convergence)
    double gradient_tolerance;   // gradient norm threshold (convergence)
    // LinearSolver: "QR" (default, robust against zero-derivative parameters)
    // or "CHOLESKY" (faster but less stable for ill-conditioned J).
    // QR factorises the augmented [J; sqrt(λ)I] directly, avoiding forming J^T J,
    // so it handles dangling (zero-column) parameters gracefully.
    bool   use_dense_qr;         // true → DENSE_QR (default), false → DENSE_NORMAL_CHOLESKY
    bool   scale_before_refine;  // true → analytically optimise Scale before first Ceres step
    // RefineScale: true → Scale is refined by variable projection (recomputed to its
    // analytic least-squares optimum at every residual evaluation).  false → Scale is
    // held fixed at the value given by the `Scale` keyword (or 1.0) throughout
    // refinement, covariance and output.  When false, scale_before_refine is moot.
    bool   refine_scale;
    // RefineBackground: true → fit an isotropic background B(|q|) = Σ b_k T_k(x),
    // Chebyshev polynomials in normalised |q|, jointly with Scale by linear least
    // squares at every evaluation. The b_k are projected parameters (not Ceres
    // parameters), exactly like Scale. background_degree is the Chebyshev degree;
    // the number of coefficients is background_degree + 1.
    bool   refine_background;
    int    background_degree;
    // Sequential block refinement: cycle through each block separately, running Ceres
    // once per block per supercycle.  0 = off (all active blocks refined together, old behaviour).
    int    num_supercycles;
    // Filename for the iteration-by-iteration trajectory JSON.
    // In supercycle mode this is overridden per block call; set explicitly to control it.
    std::string trajectory_filename;

    static RefinementOptions default_refinement_options() {
        RefinementOptions result;
        result.max_number_of_iterations = 1000;
        result.tau        = 1E-03;
        result.thresholds[0] = result.thresholds[1] = result.thresholds[2] = 1E-17;
        result.difference = 1E-06;
        result.num_threads         = 0;      // 0 → inherit max_processors
        result.function_tolerance  = 1E-6;
        result.gradient_tolerance  = 1E-10;
        result.use_dense_qr        = true;
        result.scale_before_refine = true;
        result.refine_scale        = true;
        result.refine_background   = false;
        result.background_degree   = 2;
        result.num_supercycles     = 0;
        result.trajectory_filename = "refinement_trajectory.json";
        return result;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// IntnsityCalculator
// ─────────────────────────────────────────────────────────────────────────────

class IntnsityCalculator {
public:
    /// FFT-path intensity calculation.
    static void calculate_patterson_map_from_pairs_f(
        const vector<PattersonPeak>& full_peaks,
        const vector<PattersonPeak>& avg_peaks,
        ScattererList&               scatterers,
        IntensityMap&                patterson_map,
        bool                         average_flag,
        vec3<int>                    pair_grid_size,
        vector<bool>                 periodic_directions = vector<bool>(3, false),
        vec3<int>                    border_pixels = vec3<int>(0, 0, 0),
        int                          num_threads = 0)
    {
        const int n_peaks = (int)full_peaks.size();
        double scale = patterson_map.size_1d();

        for (int i = 0; i < patterson_map.size_1d(); ++i)
            patterson_map.at_c(i) = 0;

        Grid grid_for_pairs_p(patterson_map.unit_cell(),
                              patterson_map.grid_steps(),
                              patterson_map.grid_steps().each_mul(-pair_grid_size / 2),
                              patterson_map.grid.reciprocal_flag);
        Grid grid_for_pairs_r = grid_for_pairs_p.reciprocal();

        scatterers.compute_form_factors_on_grid(pair_grid_size, grid_for_pairs_r);

        int n_threads = num_threads;
        if (n_threads <= 0) n_threads = std::thread::hardware_concurrency();
        if (n_threads <= 0) n_threads = 1;

        // Pre-allocate one scratch IntensityMap per thread to avoid per-peak
        // malloc contention (each emplace_back constructs in-place via the
        // IntensityMap(vec3<int>) ctor — safe with af::versa ownership).
        // The fill loop overwrites all elements each iteration, so invert()
        // leaving the buffer in real space is harmless before the next peak.
        vector<IntensityMap> scratch;
        scratch.reserve(n_threads);
        for (int t = 0; t < n_threads; ++t)
            scratch.emplace_back(pair_grid_size);

        // Single global mutex for accumulation.
        // add_pair_to_appropriate_place writes a pair_grid_size-wide region
        // centred at r_grid into the output map.  Nearby peaks (|Δr| < pair_grid_size)
        // write to overlapping regions, so per-r_grid mutexes are insufficient
        // and cause heap corruption on dense models.  The fill + FFT work is
        // done lock-free; only the scatter-add step is serialised.
        std::mutex accum_mutex;

        std::atomic<int> next_peak(0);

        auto worker = [&](int t) {
            IntensityMap& ppm = scratch[t];
            while (true) {
                int k = next_peak.fetch_add(1);
                if (k >= n_peaks) break;

                const PattersonPeak& fpk = full_peaks[k];
                const PattersonPeak& apk = avg_peaks[k];
                const PattersonPeak& pk  = average_flag ? apk : fpk;

                ppm.set_grid(grid_for_pairs_r);

                vec3<int>    r_grid;
                vec3<double> r_res;
                grid_and_residual(apk.r, patterson_map.grid, r_grid, r_res);
                if (!average_flag)
                    r_res += fpk.r - apk.r;

                ppm.init_iterator();
                while (ppm.next()) {
                    complex<double> f1 = scatterers.f_gridded(fpk.type1_idx,
                                                              ppm.current_index());
                    complex<double> f2 = scatterers.f_gridded(fpk.type2_idx,
                                                              ppm.current_index());
                    ppm.current_array_value_c() =
                        scale * calculate_scattering_from_a_pair_in_a_point_c(
                            f1, f2,
                            pk.coefficient,
                            1.0,
                            ppm.current_s(),
                            r_res,
                            pk.U);
                }

                ppm.invert();

                {
                    std::lock_guard<std::mutex> lock(accum_mutex);
                    add_pair_to_appropriate_place(ppm, patterson_map,
                                                  r_grid, periodic_directions, border_pixels);
                }
            }
        };

        if (n_threads <= 1) {
            worker(0);
        } else {
            vector<std::thread> workers;
            for (int t = 0; t < n_threads; ++t) workers.emplace_back(worker, t);
            for (auto& w : workers) w.join();
        }
    }

    /// FFT-path analytical derivative map calculation for a single parameter.
    static void calculate_patterson_map_derivative_from_pairs_f(
        const vector<PattersonPeak>& base_full_peaks,
        const vector<PattersonPeak>& base_avg_peaks,
        const vector<PeakSusceptibility>& full_susc,
        const vector<PeakSusceptibility>& avg_susc,
        ScattererList& scatterers,
        IntensityMap&  deriv_patterson_map,
        bool           average_flag,
        vec3<int>      pair_grid_size,
        vector<bool>   periodic_directions = vector<bool>(3, false),
        vec3<int>      border_pixels = vec3<int>(0, 0, 0),
        int            num_threads = 0)
    {
        const int n_peaks = (int)base_full_peaks.size();
        double scale = deriv_patterson_map.size_1d();

        for (int i = 0; i < deriv_patterson_map.size_1d(); ++i)
            deriv_patterson_map.at_c(i) = 0;

        Grid grid_for_pairs_p(deriv_patterson_map.unit_cell(),
                              deriv_patterson_map.grid_steps(),
                              deriv_patterson_map.grid_steps().each_mul(-pair_grid_size / 2),
                              deriv_patterson_map.grid.reciprocal_flag);
        Grid grid_for_pairs_r = grid_for_pairs_p.reciprocal();

        // Form factors are NOT recomputed here. The caller (Model) is responsible for
        // providing a ScattererList with up-to-date gridded_form_factors_ for this grid.
        // For non-parameterized form factors the Model caches them in thread_scatterer_lists_.
        // When form factors become parameterized, update them in the Model before this call.

        vector<int> active_indices;
        for (int i = 0; i < n_peaks; ++i) {
            const auto& sk = average_flag ? avg_susc[i] : full_susc[i];
            bool has_susc = (sk.d_coefficient != 0.0) || (sk.d_r.length_sq() != 0.0);
            if (!has_susc) {
                for (int m = 0; m < 6; ++m) if (sk.d_U[m] != 0.0) { has_susc = true; break; }
            }
            if (has_susc) active_indices.push_back(i);
        }

        const int n_active = (int)active_indices.size();

        int n_threads = num_threads;
        if (n_threads <= 0) n_threads = std::thread::hardware_concurrency();
        if (n_threads <= 0) n_threads = 1;

        vector<IntensityMap> scratch;
        scratch.reserve(n_threads);
        for (int t = 0; t < n_threads; ++t)
            scratch.emplace_back(pair_grid_size);

        std::mutex accum_mutex;

        std::atomic<int> next_active(0);

        auto worker = [&](int t) {
            IntensityMap& ppm = scratch[t];
            while (true) {
                int idx = next_active.fetch_add(1);
                if (idx >= n_active) break;
                int k = active_indices[idx];

                const PattersonPeak& fpk = base_full_peaks[k];
                const PattersonPeak& apk = base_avg_peaks[k];
                const PattersonPeak& pk  = average_flag ? apk : fpk;
                const PeakSusceptibility& sk = average_flag ? avg_susc[k] : full_susc[k];

                ppm.set_grid(grid_for_pairs_r);

                vec3<int>    r_grid;
                vec3<double> r_res;
                grid_and_residual(apk.r, deriv_patterson_map.grid, r_grid, r_res);
                if (!average_flag)
                    r_res += fpk.r - apk.r;

                ppm.init_iterator();
                while (ppm.next()) {
                    complex<double> f1 = scatterers.f_gridded(fpk.type1_idx, ppm.current_index());
                    complex<double> f2 = scatterers.f_gridded(fpk.type2_idx, ppm.current_index());

                    vec3<double> s = ppm.current_s();
                    double phase_val = M_2PI * (s * r_res);
                    double adp_val   = M2PISQ * (s * pk.U * s);
                    complex<double> E = exp(complex<double>(adp_val, phase_val));

                    double d_phase = M_2PI * (s * sk.d_r);
                    double d_adp   = M2PISQ * (s * sk.d_U * s);

                    complex<double> d_term = sk.d_coefficient * E + pk.coefficient * E * complex<double>(d_adp, d_phase);
                    ppm.current_array_value_c() = scale * conj(f1) * f2 * d_term;
                }

                ppm.invert();

                {
                    std::lock_guard<std::mutex> lock(accum_mutex);
                    add_pair_to_appropriate_place(ppm, deriv_patterson_map,
                                                  r_grid, periodic_directions, border_pixels);
                }
            }
        };

        if (n_threads <= 1) {
            worker(0);
        } else {
            vector<std::thread> workers;
            for (int t = 0; t < n_threads; ++t) workers.emplace_back(worker, t);
            for (auto& w : workers) w.join();
        }
    }

    static void calculate_scattering_from_pairs(vector<AtomicPair> pairs, const Eigen::VectorXd& params, IntensityMap& I, bool average_flag)
    {
        double d_star_square;
        complex<double> f1, f2;
        vec3<double> s;
        AtomicPair* pair;

        I.init_iterator();
        while (I.next()) {
            s            = I.current_s();
            d_star_square = I.current_d_star_square();

            AtomicTypeCollection::update_current_form_factors(s, d_star_square);
            int sz = pairs.size();
            double Intensity = 0;
            for (int i = 0; i < sz; i++) {
                pair = &pairs[i];
                f1   = pair->atomic_type1->current_form_factor;
                f2   = pair->atomic_type2->current_form_factor;
                
                double p_val = pair->p(average_flag)->eval(params);
                vec3<double> r_val(pair->r(average_flag).x->eval(params), pair->r(average_flag).y->eval(params), pair->r(average_flag).z->eval(params));
                sym_mat3<double> U_val(pair->U(average_flag).u11->eval(params), pair->U(average_flag).u22->eval(params), pair->U(average_flag).u33->eval(params),
                                       pair->U(average_flag).u12->eval(params), pair->U(average_flag).u13->eval(params), pair->U(average_flag).u23->eval(params));

                Intensity += real(calculate_scattering_from_a_pair_in_a_point_c(
                    f1, f2,
                    p_val,
                    pair->multiplier,
                    s,
                    r_val,
                    U_val));
            }
            I.current_array_value() = Intensity;
        }
    }

    /// Intensity over a list of PattersonPeaks using integer-indexed ScattererList.
    /// Replaces calculate_scattering_from_pairs for the direct method.
    static void calculate_scattering_from_patterson_peaks(
        const vector<PattersonPeak>& peaks,
        ScattererList&               scatterers,
        IntensityMap&                I)
    {
        I.init_iterator();
        while (I.next()) {
            vec3<double> s       = I.current_s();
            double       d_sq    = I.current_d_star_square();
            scatterers.update(s, d_sq);

            double intensity = 0.0;
            for (const PattersonPeak& pk : peaks) {
                complex<double> f1 = scatterers.f(pk.type1_idx);
                complex<double> f2 = scatterers.f(pk.type2_idx);
                intensity += real(conj(f1) * f2 * pk.coefficient *
                    exp(complex<double>(M2PISQ * (s * pk.U * s),
                                        M_2PI  * (s * pk.r))));
            }
            I.current_array_value() = intensity;
        }
    }

    /// Calculate derivative map dI/dp_j for a single parameter j.
    static void calculate_scattering_derivative_from_patterson_peaks(
        const vector<PattersonPeak>& peaks,
        const vector<PeakSusceptibility>& susceptibilities,
        ScattererList&               scatterers,
        IntensityMap&                dI)
    {
        dI.init_iterator();
        while (dI.next()) {
            vec3<double> s       = dI.current_s();
            double       d_sq    = dI.current_d_star_square();
            scatterers.update(s, d_sq);

            double deriv_val = 0.0;
            for (size_t k = 0; k < peaks.size(); ++k) {
                const PattersonPeak& pk = peaks[k];
                const PeakSusceptibility& sk = susceptibilities[k];

                complex<double> f1f2 = conj(scatterers.f(pk.type1_idx)) * scatterers.f(pk.type2_idx);
                
                double phase_val = M_2PI * (s * pk.r);
                double adp_val   = M2PISQ * (s * pk.U * s);
                complex<double> E = exp(complex<double>(adp_val, phase_val));

                double d_phase = M_2PI * (s * sk.d_r);
                double d_adp   = M2PISQ * (s * sk.d_U * s);

                // d(p * E) = dp * E + p * E * (d_adp + i * d_phase)
                complex<double> term1 = sk.d_coefficient * E;
                complex<double> term2 = pk.coefficient * E * complex<double>(d_adp, d_phase);

                deriv_val += real(f1f2 * (term1 + term2));
            }
            dI.current_array_value() = deriv_val;
        }
    }

    inline static complex<double> calculate_scattering_from_a_pair_in_a_point_c(
        complex<double> const& f1, complex<double> const& f2,
        double const& p,
        double const& N,
        vec3<double> const& s,
        vec3<double> const& r,
        sym_mat3<double> const& U)
    {
        return conj(f1) * f2 * p * N * exp(complex<double>(M2PISQ * (s * U * s), M_2PI * (s * r)));
    }

    class UnknownSymmetry {};
};

#endif // YELL_CALCULATOR_H
