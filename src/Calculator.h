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

using namespace std;
using namespace scitbx;

// ─────────────────────────────────────────────────────────────────────────────
// MinimizerCalculator — interface for the minimizer
// ─────────────────────────────────────────────────────────────────────────────

class MinimizerCalculator {
public:
    virtual void calculate(vector<double> p) = 0;
    virtual IntensityMap& data() = 0;
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

    static RefinementOptions default_refinement_options() {
        RefinementOptions result = {
            1000,                        // max_number_of_iterations
            1E-03,                       // tau
            {1E-17, 1E-17, 1E-17},      // thresholds
            1E-06                        // difference
        };
        return result;
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// IntnsityCalculator
// ─────────────────────────────────────────────────────────────────────────────

class IntnsityCalculator {
public:
    static void calculate_patterson_map_from_pairs_f(
        vector<AtomicPair> pairs,
        IntensityMap& patterson_map,
        bool average_flag,
        vec3<int> pair_grid_size,
        vector<bool> periodic_directions = vector<bool>(3, false))
    {
        double d_star_square;
        complex<double> f1, f2;
        vec3<double> s, r_res;
        vec3<int> r_grid;

        double scale = patterson_map.size_1d();

        for (int i = 0; i < patterson_map.size_1d(); ++i)
            patterson_map.at_c(i) = 0;

        AtomicPair* pair;

        Grid grid_for_pairs_p(patterson_map.unit_cell(),
                              patterson_map.grid_steps(),
                              patterson_map.grid_steps().each_mul(-pair_grid_size / 2),
                              patterson_map.grid.reciprocal_flag);

        Grid grid_for_pairs_r = grid_for_pairs_p.reciprocal();

        AtomicTypeCollection::calculate_form_factors_of_all_atoms_on_grid(pair_grid_size, grid_for_pairs_r);

        omp_set_num_threads(8);
#pragma omp parallel for private(s,r_res,r_grid,d_star_square,f1,f2,pair)
        for (int i = 0; i < pairs.size(); ++i) {
            pair = &pairs[i];
            IntensityMap pair_patterson_map(pair_grid_size);
            pair_patterson_map.set_grid(grid_for_pairs_r);

            grid_and_residual(pair->average_r(), patterson_map.grid, r_grid, r_res);

            if (!average_flag)
                r_res += pair->r() - pair->average_r();

            pair_patterson_map.init_iterator();
            while (pair_patterson_map.next()) {
                s            = pair_patterson_map.current_s();
                d_star_square = pair_patterson_map.current_d_star_square();

                f1 = pair->atomic_type1->gridded_form_factors.at_c(pair_patterson_map.current_index());
                f2 = pair->atomic_type2->gridded_form_factors.at_c(pair_patterson_map.current_index());

                pair_patterson_map.current_array_value_c() =
                    scale * calculate_scattering_from_a_pair_in_a_point_c(
                        f1, f2,
                        pair->p(average_flag),
                        pair->multiplier,
                        s,
                        r_res,
                        pair->U(average_flag));
            }

            pair_patterson_map.invert();

#pragma omp critical
            add_pair_to_appropriate_place(pair_patterson_map, patterson_map, r_grid, periodic_directions);
        }
    }

    static void calculate_scattering_from_pairs(vector<AtomicPair> pairs, IntensityMap& I, bool average_flag)
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
                Intensity += real(calculate_scattering_from_a_pair_in_a_point_c(
                    f1, f2,
                    pair->p(average_flag),
                    pair->multiplier,
                    s,
                    pair->r(average_flag),
                    pair->U(average_flag)));
            }
            I.current_array_value() = Intensity;
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
