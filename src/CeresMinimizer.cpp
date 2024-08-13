//
// Created by Arkadiy Simonov on 05.08.24.
//

#include "CeresMinimizer.h"

//TODO: make code so that ceres could be called as a minimizer to swap-replace levmar
//DONE: figure out how to define ceres with dynamical number of variables: Use DynamicNumericDiffCostFunction
//LATER:
//TODO: figure out stopping criteria, report it reasonably
//TODO: figure out how to get the covariances out


#include "ceres/ceres.h"
#include "glog/logging.h"
#include "IntensityMap.h"


/**
   * Solves the problem of finding parameters which minimize I_model(params)-I_experimental in the Least-square sense.
   * \param _calc - a reference to an object that calculates model diffuse scattering (or PDF). The object should implement MinimizerCalculator interface
   */

vector<double> CeresMinimizer::minimize(const vector<double> initial_params,
                                        IntensityMap * _experimental_data,
                                        MinimizerCalculator * _calc,
                                        OptionalIntensityMap * _weights,
                                        RefinementOptions & refinement_options)
{
    calc = _calc;
    experimental_data = _experimental_data;
    weights = _weights;
    parameters_number = initial_params.size();

    double * p = (double*) malloc(sizeof(double)*initial_params.size());
    for(int i=0; i<initial_params.size(); i++)
        p[i]=initial_params[i];

    ceres::Problem problem;

    auto cost_function =
            new ceres::DynamicNumericDiffCostFunction<CeresMinimizer, ceres::CENTRAL>(this, ceres::DO_NOT_TAKE_OWNERSHIP);

    cost_function->AddParameterBlock(parameters_number);
    cost_function->SetNumResiduals(calc->number_of_observations());
    problem.AddResidualBlock(cost_function, nullptr, p);

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR; //probably denseNormalCholesky
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = refinement_options.max_number_of_iterations;
    ceres::Solver::Summary summary;
    Solve(options, &problem, &summary);

    static vector<double> result(p,p+initial_params.size());



    ceres::Covariance::Options opt;
//    opt.algorithm_type = ceres::DENSE_SVD;
    ceres::Covariance covariance(opt);

    std::vector<std::pair<const double*, const double*> > covariance_blocks;
    covariance_blocks.push_back(make_pair(p, p));

    //CHECK
    covariance.Compute(covariance_blocks, &problem);

    covar = vector<double>(parameters_number*parameters_number);

    covariance.GetCovarianceBlock(p, p, covar.data());

    delete p;
    return result;
}

bool CeresMinimizer::operator()(double const *const *params, double *residuals) const {
    auto p = params[0];
    vector<double> yell_parameters(p, p + parameters_number);

    calc->calculate(yell_parameters);

    auto datapoints_number = calc->number_of_observations();

    if (calc->refine_in_asu()) {
        //copy difference to the *x array
        for(int ii=0; ii<datapoints_number; ii++) {
            auto i = calc->asu_indices()[ii];
            residuals[ii] = (experimental_data->at(i) - calc->data().at(i))*weights->at(i);
        }
    }
    else
    {
        //copy difference to the *x array
        for(int i=0; i<datapoints_number; i++)
            residuals[i] = (experimental_data->at(i) - calc->data().at(i))*weights->at(i);
    }

    return true;
}
