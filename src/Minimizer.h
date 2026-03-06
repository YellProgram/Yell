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

#ifndef YELL_MINIMIZER_H
#define YELL_MINIMIZER_H

#include "Calculator.h"
#include "levmar.h"

#include <vector>
#include <cstdlib>

using namespace std;

/**
 * Wraps the levmar library for least-square minimisation of
 * min(I_model(params) - I_experimental).
 * TODO: remove once CeresMinimizer is the sole minimizer.
 */
class Minimizer {
public:
    Minimizer() : covar(NULL) {}

    ~Minimizer() { free(covar); }

    vector<double> minimize(const vector<double> initial_params,
                            IntensityMap* _experimental_data,
                            MinimizerCalculator* _calc,
                            OptionalIntensityMap* _weights,
                            RefinementOptions refinement_options = RefinementOptions::default_refinement_options())
    {
        calc              = _calc;
        experimental_data = _experimental_data;
        weights           = _weights;

        double* p = (double*) malloc(sizeof(double) * initial_params.size());
        for (int i = 0; i < initial_params.size(); i++)
            p[i] = initial_params[i];

        double* x = (double*) malloc(sizeof(double) * experimental_data->size_1d());
        for (int i = 0; i < experimental_data->size_1d(); i++)
            x[i] = 0;

        covar = (double*) malloc(sizeof(double) * initial_params.size() * initial_params.size());

        double opts[5];
        opts[0] = refinement_options.tau;
        for (int i = 0; i < 3; ++i)
            opts[i + 1] = refinement_options.thresholds[i];
        opts[4] = refinement_options.difference;

        double info[10];

        // levmar call commented out — replaced by CeresMinimizer
        // int ret = dlevmar_dif(func_for_levmar, p, x,
        //                       initial_params.size(), experimental_data->size_1d(),
        //                       refinement_options.max_number_of_iterations,
        //                       opts, info, NULL, covar, this);

        static vector<double> result(p, p + initial_params.size());

        delete p;
        delete x;
        return result;
    }

    double* covar;

private:
    static void func_for_levmar(double const* p, double* x,
                                 int parameters_number, int datapoints_number,
                                 void* data)
    {
        Minimizer* _this = reinterpret_cast<Minimizer*>(data);
        vector<double> parameters(p, p + parameters_number);
        _this->calc->calculate(parameters);
        for (int i = 0; i < datapoints_number; i++)
            x[i] = (_this->experimental_data->at(i) - _this->calc->data().at(i))
                   * _this->weights->at(i);
    }

    MinimizerCalculator* calc;
    IntensityMap* experimental_data;
    OptionalIntensityMap* weights;
};

#endif // YELL_MINIMIZER_H
