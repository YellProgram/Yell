//
// Created by Arkadiy Simonov on 05.08.24.
//

#include "CeresMinimizer.h"
#include "model.h"

//TODO: make code so that ceres could be called as a minimizer to swap-replace levmar
//DONE: figure out how to define ceres with dynamical number of variables: Use DynamicNumericDiffCostFunction
//LATER:
//TODO: figure out stopping criteria, report it reasonably
//TODO: figure out how to get the covariances out

#include "ceres/ceres.h"
#include "glog/logging.h"
#include "IntensityMap.h"
#include <fstream>

// ─────────────────────────────────────────────────────────────────────────────
// Analytical cost function: computes residuals and Jacobian via ExprPtr trees.
// Used when model->derivatives_mode == ANALYTICAL.
// ─────────────────────────────────────────────────────────────────────────────
class AnalyticalYellCostFunction : public ceres::CostFunction {
public:
    AnalyticalYellCostFunction(Model* model, IntensityMap* exp, OptionalIntensityMap* weights)
        : model_(model), exp_(exp), weights_(weights)
    {
        set_num_residuals(model->number_of_observations());
        mutable_parameter_block_sizes()->push_back(
            (int)model->refinement_parameters.size());
    }

    bool Evaluate(double const* const* parameters,
                  double* residuals,
                  double** jacobians) const override
    {
        int n_params = parameter_block_sizes()[0];
        vector<double> p(parameters[0], parameters[0] + n_params);

        model_->calculate(p);

        int n_obs = model_->number_of_observations();
        if (model_->refine_in_asu()) {
            for (int ii = 0; ii < n_obs; ++ii) {
                int i = model_->asu_indices()[ii];
                residuals[ii] = (exp_->at(i) - model_->data().at(i)) * weights_->at(i);
            }
        } else {
            for (int i = 0; i < n_obs; ++i)
                residuals[i] = (exp_->at(i) - model_->data().at(i)) * weights_->at(i);
        }

        if (jacobians && jacobians[0]) {
            Eigen::MatrixXd J = model_->compute_analytical_jacobian_direct(
                p, *exp_, *weights_);
            // Ceres expects row-major: jacobians[0][ii * n_params + j]
            for (int ii = 0; ii < n_obs; ++ii)
                for (int j = 0; j < n_params; ++j)
                    jacobians[0][ii * n_params + j] = J(ii, j);
        }

        return true;
    }

private:
    Model* model_;
    IntensityMap* exp_;
    OptionalIntensityMap* weights_;
};

class JsonIterationLogger : public ceres::IterationCallback {
public:
    JsonIterationLogger(const double* accepted_params,
                        const std::vector<double>* trial_params,
                        int size, std::string filename)
        : accepted_params_(accepted_params),
          trial_params_(trial_params),
          size_(size),
          filename_(std::move(filename)) {}

    ceres::CallbackReturnType operator()(const ceres::IterationSummary& summary) override {
        const double* src = summary.step_is_successful
            ? accepted_params_
            : trial_params_->data();

        IterationData data;
        data.iteration     = summary.iteration;
        data.cost          = summary.cost;
        data.gradient_norm = summary.gradient_norm;
        data.step_norm     = summary.step_norm;
        data.step_accepted = summary.step_is_successful;
        data.parameters.assign(src, src + size_);
        history_.push_back(std::move(data));
        SaveToFile(filename_);
        return ceres::SOLVER_CONTINUE;
    }

    void SaveToFile(const std::string& filename) const {
        std::ofstream file(filename);
        file << "[\n";
        for (size_t i = 0; i < history_.size(); ++i) {
            const auto& d = history_[i];
            file << "  {\"iteration\": "    << d.iteration
                 << ", \"cost\": "           << d.cost
                 << ", \"gradient_norm\": "  << d.gradient_norm
                 << ", \"step_norm\": "      << d.step_norm
                 << ", \"step_accepted\": "  << (d.step_accepted ? "true" : "false")
                 << ", \"parameters\": [";
            for (int j = 0; j < (int)d.parameters.size(); ++j) {
                if (j) file << ", ";
                file << d.parameters[j];
            }
            file << "]}";
            if (i + 1 < history_.size()) file << ",";
            file << "\n";
        }
        file << "]\n";
    }

private:
    const double* accepted_params_;
    const std::vector<double>* trial_params_;
    int size_;
    std::string filename_;

    struct IterationData {
        int iteration;
        double cost;
        double gradient_norm;
        double step_norm;
        bool step_accepted;
        std::vector<double> parameters;
    };
    std::vector<IterationData> history_;
};

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

    Model* model = dynamic_cast<Model*>(_calc);
    if (model && model->derivatives_mode == ANALYTICAL) {
        auto* cost_function = new AnalyticalYellCostFunction(model, _experimental_data, _weights);
        problem.AddResidualBlock(cost_function, nullptr, p);
    } else {
        auto cost_function =
            new ceres::DynamicNumericDiffCostFunction<CeresMinimizer, ceres::CENTRAL>(
                this, ceres::DO_NOT_TAKE_OWNERSHIP);
        cost_function->AddParameterBlock(parameters_number);
        cost_function->SetNumResiduals(calc->number_of_observations());
        problem.AddResidualBlock(cost_function, nullptr, p);
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR; //
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = refinement_options.max_number_of_iterations;
    last_eval_params_.resize(parameters_number);
    options.update_state_every_iteration = true;
    JsonIterationLogger logger(p, &last_eval_params_, parameters_number, "refinement_trajectory.json");
    options.callbacks.push_back(&logger);
    ceres::Solver::Summary summary;
    Solve(options, &problem, &summary);

    static vector<double> result(p,p+initial_params.size());

    ceres::Covariance::Options opt;
    opt.algorithm_type = ceres::DENSE_SVD;
    opt.null_space_rank = -1;
    ceres::Covariance covariance(opt);

    std::vector<std::pair<const double*, const double*> > covariance_blocks;
    covariance_blocks.push_back(make_pair(p, p));

    covariance.Compute(covariance_blocks, &problem);

    covar = vector<double>(parameters_number*parameters_number);

    covariance.GetCovarianceBlock(p, p, covar.data());

    delete p;
    return result;
}

bool CeresMinimizer::operator()(double const *const *params, double *residuals) const {
    auto p = params[0];
    last_eval_params_.assign(p, p + parameters_number);
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
