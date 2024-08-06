//
// Created by Arkadiy Simonov on 05.08.24.
//

#include "Minimizer.h"

//TODO: make code so that ceres could be called as a minimizer to swao-preplace levmar
//DONE: figure out how to define ceres with dynamical number of variables: Use DynamicNumericDiffCostFunction
//LATER:
//TODO: figure out stopping criteria, report it reasonably
//TODO: figure out how to get the covariances out


#include "ceres/ceres.h"
#include "glog/logging.h"

//void func_for_levmar(double *p, double *x, int parameters_number, int datapoints_number, void *data)
class CostFunctor {
public:
    CostFunctor(void (*target_function)(double const *, double *, int, int, void *),
                int parameters_number,
                int datapoints_number,
                void * data): target_function(target_function), parameters_number(parameters_number),
                datapoints_number(datapoints_number), data(data)
                {}

    bool operator()(double const *const *parameters, double *residuals) const {
        target_function(parameters[0], residuals, parameters_number, datapoints_number, data);
        return true;
    }

private:
    void (*target_function)(double const *, double *, int, int, void *);
    void * data;
    int parameters_number;
    int datapoints_number;
};


void CeresMinimizer::minimize(void (*target_function)(double const *, double *, int, int, void *), double *p, double *x, int parameter_number,
                              int datapoints_number, int itmax, void *data) {

    ceres::Problem problem;

    auto cost_function =
            new ceres::DynamicNumericDiffCostFunction<CostFunctor, ceres::CENTRAL>(new CostFunctor(target_function, parameter_number, datapoints_number, data));

    cost_function->AddParameterBlock(parameter_number);
    cost_function->SetNumResiduals(datapoints_number);
    problem.AddResidualBlock(cost_function, nullptr, p);

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = itmax;
    ceres::Solver::Summary summary;
    Solve(options, &problem, &summary);
}

/*
 *
 *
}
Since the sizing of the parameters is done at runtime, you must also specify the sizes after creating the dynamic numeric diff cost function. For example:

auto cost_function = std::make_unique<DynamicNumericDiffCostFunction<MyCostFunctor>>();
cost_function->AddParameterBlock(5);
cost_function->AddParameterBlock(10);
cost_function->SetNumResiduals(21);
 *
 */