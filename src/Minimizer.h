//
// Created by Arkadiy Simonov on 05.08.24.
//

#ifndef YELL_MINIMIZER_H
#define YELL_MINIMIZER_H

//

class CeresMinimizer {
public:
    void minimize(void (*target_function)(double const *p, double *x, int m, int n, void* data),
                  double* p,
                  double* x,
                  int parameter_number, //parameter vector dimension
                  int datapoints_number, //measurement vector dimension
                  int itmax, //maximum number of iterations
                  void* data //extra data. This will be a calculator data to be passed to the calculate function
                  );

    //void func_for_levmar(double *p, double *x, int parameters_number, int datapoints_number, void *data)
    //dlevmar_dif(func_for_levmar, p, x, initial_params.size(),experimental_data->size_1d(), refinement_options.max_number_of_iterations, opts, info, NULL, covar, this);
};


#endif //YELL_MINIMIZER_H
