//
// Created by Arkadiy Simonov on 05.08.24.
//

#ifndef YELL_CERESMINIMIZER_H
#define YELL_CERESMINIMIZER_H

#include "basic_classes.h"

class MinimizerCalculator;
class RefinementOptions;

class CeresMinimizer {
public:
    vector<double> minimize(const vector<double> initial_params,
                            IntensityMap * _experimental_data,
                            MinimizerCalculator * _calc,
                            OptionalIntensityMap * _weights,
                            RefinementOptions & refinement_options);

    bool operator()(double const *const *parameters, double *residuals) const;
    vector<double> covar;
private:
    MinimizerCalculator * calc;
    IntensityMap * experimental_data;
    OptionalIntensityMap *  weights;

    int parameters_number;


    //void func_for_levmar(double *p, double *x, int parameters_number, int datapoints_number, void *data)
    //dlevmar_dif(func_for_levmar, p, x, initial_params.size(),experimental_data->size_1d(), refinement_options.max_number_of_iterations, opts, info, NULL, covar, this);
};


#endif //YELL_CERESMINIMIZER_H
