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
    mutable std::vector<double> last_eval_params_;


};


#endif //YELL_CERESMINIMIZER_H
