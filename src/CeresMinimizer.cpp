//
// Created by Arkadiy Simonov on 05.08.24.
//

#include "CeresMinimizer.h"
#include "model.h"
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "IntensityMap.h"
#include <fstream>
#include <thread>
#include <mutex>
#include <atomic>

// ─────────────────────────────────────────────────────────────────────────────
// Analytical cost function: computes residuals and Jacobian via ExprPtr trees.
// Implements Variable Projection for the global Scale parameter.
// ─────────────────────────────────────────────────────────────────────────────
class AnalyticalYellCostFunction : public ceres::CostFunction {
public:
    AnalyticalYellCostFunction(Model* model, IntensityMap* exp, OptionalIntensityMap* weights)
        : model_(model), exp_(exp), weights_(weights)
    {
        set_num_residuals(model->number_of_observations());
        
        // Structural blocks only (Scale is internal)
        for (auto& block : model->parameter_blocks) {
            mutable_parameter_block_sizes()->push_back((int)block.size());
        }
    }

    bool Evaluate(double const* const* parameters,
                  double* residuals,
                  double** jacobians) const override
    {
        // 1. Flatten structural blocks
        vector<double> p;
        p.push_back(1.0); 
        for (size_t b = 0; b < parameter_block_sizes().size(); ++b) {
            int sz = parameter_block_sizes()[b];
            p.insert(p.end(), parameters[b], parameters[b] + sz);
        }

        // 2. Base model calculation (Unscaled)
        model_->calculate(p);

        int n_obs = model_->number_of_observations();
        const bool use_asu = model_->refine_in_asu();
        const vector<int>& asu = model_->asu_indices();

        // 3. Analytical Scale Optimization
        double num = 0.0, den = 0.0;
        for (int ii = 0; ii < n_obs; ++ii) {
            int i = use_asu ? asu[ii] : ii;
            double w = weights_->at(i);
            double Ic = model_->get_intensity_map().at(i) - model_->get_average_intensity_map().at(i);
            double Ie = exp_->at(i);
            num += w * w * Ie * Ic;
            den += w * w * Ic * Ic;
        }
        double S = (den > 1e-15) ? (num / den) : 1.0;
        if (S < 0) S = 0;
        model_->set_scale(S);
        p[0] = S;

        // 4. Residuals
        for (int ii = 0; ii < n_obs; ++ii) {
            int i = use_asu ? asu[ii] : ii;
            residuals[ii] = (exp_->at(i) - S * (model_->get_intensity_map().at(i) - model_->get_average_intensity_map().at(i))) * weights_->at(i);
        }

        // 5. Parallel Jacobian evaluation
        if (jacobians) {
            unsigned int n_threads = model_->max_processors;
            if (n_threads <= 0) n_threads = std::thread::hardware_concurrency();
            if (n_threads <= 0) n_threads = 1;

            vector<Model*> thread_models(n_threads);
            for (size_t t = 0; t < n_threads; ++t) thread_models[t] = model_->clone();

            int global_offset = 1;
            for (size_t b = 0; b < parameter_block_sizes().size(); ++b) {
                int block_sz = parameter_block_sizes()[b];
                if (jacobians[b]) {
                    std::atomic<int> next_j(0);
                    vector<std::thread> workers;
                    for (size_t t = 0; t < n_threads; ++t) {
                        workers.emplace_back([&, t, p]() {
                            Model* m = thread_models[t];
                            while (true) {
                                int j = next_j.fetch_add(1);
                                if (j >= block_sz) break;

                                // Each thread calculates one derivative serially (num_threads=1)
                                IntensityMap dI_map = m->calculate_derivative(p, global_offset + j, 1);
                                for (int ii = 0; ii < n_obs; ++ii) {
                                    int i = use_asu ? asu[ii] : ii;
                                    jacobians[b][ii * block_sz + j] = -dI_map.at(i) * weights_->at(i);
                                }
                            }
                        });
                    }
                    for (auto& w : workers) w.join();
                }
                global_offset += block_sz;
            }

            for (auto* m : thread_models) delete m;
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
    JsonIterationLogger(Model* model,
                        const vector<double*>& p_pointers,
                        const vector<int>& block_sizes,
                        std::string filename)
        : model_(model),
          p_pointers_(p_pointers),
          block_sizes_(block_sizes),
          filename_(std::move(filename)) {}

    ceres::CallbackReturnType operator()(const ceres::IterationSummary& summary) override {
        IterationData data;
        data.iteration     = summary.iteration;
        data.cost          = summary.cost;
        data.gradient_norm = summary.gradient_norm;
        data.step_norm     = summary.step_norm;
        data.step_accepted = summary.step_is_successful;
        
        data.parameters.push_back(model_->refinement_parameters[0]);

        for (size_t b = 0; b < p_pointers_.size(); ++b) {
            for (int i = 0; i < block_sizes_[b]; ++i) {
                data.parameters.push_back(p_pointers_[b][i]);
            }
        }
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
    Model*          model_;
    vector<double*> p_pointers_;
    vector<int>     block_sizes_;
    std::string     filename_;

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

    Model* model = dynamic_cast<Model*>(_calc);
    ceres::Problem problem;
    vector<double*> p_pointers;
    vector<int> block_sizes;

    if (!model || model->parameter_blocks.empty()) {
        if (initial_params.size() > 1) {
            int rem = (int)initial_params.size() - 1;
            double* p_rest = new double[rem];
            for (int i = 0; i < rem; ++i) p_rest[i] = initial_params[i + 1];
            p_pointers.push_back(p_rest);
            block_sizes.push_back(rem);
        }
    } else {
        int offset = 1;
        for (auto& block : model->parameter_blocks) {
            double* pb = new double[block.size()];
            for (size_t i = 0; i < block.size(); ++i) pb[i] = initial_params[offset + i];
            p_pointers.push_back(pb);
            block_sizes.push_back((int)block.size());
            offset += (int)block.size();
        }
    }

    if (model && (model->derivatives_mode == ANALYTICAL || model->derivatives_mode == MIXED)) {
        auto* cost_function = new AnalyticalYellCostFunction(model, _experimental_data, _weights);
        problem.AddResidualBlock(cost_function, nullptr, p_pointers);
    } else {
        auto cost_function =
            new ceres::DynamicNumericDiffCostFunction<CeresMinimizer, ceres::CENTRAL>(
                this, ceres::DO_NOT_TAKE_OWNERSHIP);
        for (int sz : block_sizes) {
            cost_function->AddParameterBlock(sz);
        }
        cost_function->SetNumResiduals(calc->number_of_observations());
        problem.AddResidualBlock(cost_function, nullptr, p_pointers);
    }

    // Freeze inactive blocks so Ceres neither varies them nor includes them in covariance.
    for (size_t b = 0; b < p_pointers.size(); ++b) {
        if (!model->block_is_active((int)b + 1))
            problem.SetParameterBlockConstant(p_pointers[b]);
    }

    ceres::Solver::Options options;
    // QR factorises the augmented [J; sqrt(λ)I] directly — more robust against
    // near-zero Jacobian columns (dangling parameters).
    // DENSE_NORMAL_CHOLESKY forms J^T J explicitly — faster but less stable.
    // Controlled by "LinearSolver QR|CHOLESKY" in model.txt.
    options.linear_solver_type = refinement_options.use_dense_qr
                                 ? ceres::DENSE_QR
                                 : ceres::DENSE_NORMAL_CHOLESKY;
    REPORT(MAIN) << "Linear solver: "
                 << (refinement_options.use_dense_qr ? "DENSE_QR (robust, factorises J directly)"
                                                     : "DENSE_NORMAL_CHOLESKY (fast, factorises J^T J)")
                 << "\n";
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = refinement_options.max_number_of_iterations;
    options.function_tolerance  = refinement_options.function_tolerance;
    options.gradient_tolerance  = refinement_options.gradient_tolerance;
    // num_threads: controls Eigen/BLAS parallelism in the linear solve and, for
    // DynamicNumericDiff, the parallelism of Jacobian-column evaluations across
    // parameter blocks.  Use model's max_processors unless overridden in model.txt
    // via "CeresThreads N".
    {
        int nt = refinement_options.num_threads > 0
                 ? refinement_options.num_threads
                 : model->max_processors;
        if (nt <= 0) nt = std::thread::hardware_concurrency();
        if (nt <= 0) nt = 1;
        options.num_threads = nt;
    }
    last_eval_params_.resize(parameters_number);
    options.update_state_every_iteration = true;
    
    JsonIterationLogger logger(model, p_pointers, block_sizes, "refinement_trajectory.json");
    options.callbacks.push_back(&logger);
    
    ceres::Solver::Summary summary;
    Solve(options, &problem, &summary);

    vector<double> result;
    result.push_back(model->refinement_parameters[0]); 
    for (size_t b = 0; b < p_pointers.size(); ++b) {
        for (int i = 0; i < block_sizes[b]; ++i) {
            result.push_back(p_pointers[b][i]);
        }
    }
    
    if (model && model->print_covariance_matrix) {
        Eigen::MatrixXd cov = model->compute_full_covariance(result, *experimental_data, *weights);
        covar = vector<double>(cov.data(), cov.data() + cov.size());
    } else if (model) {
        // Always compute covariance to avoid segfault in main
        Eigen::MatrixXd cov = model->compute_full_covariance(result, *experimental_data, *weights);
        covar = vector<double>(cov.data(), cov.data() + cov.size());
    }

    for (double* pb : p_pointers) delete[] pb;
    
    return result;
}

bool CeresMinimizer::operator()(double const *const *params, double *residuals) const {
    vector<double> yell_parameters;
    Model* model = dynamic_cast<Model*>(calc);
    
    yell_parameters.push_back(1.0); 

    if (!model || model->parameter_blocks.empty()) {
        if (parameters_number > 1) {
            for (int i = 0; i < parameters_number - 1; ++i)
                yell_parameters.push_back(params[0][i]);
        }
    } else {
        for (size_t b = 0; b < model->parameter_blocks.size(); ++b) {
            int sz = (int)model->parameter_blocks[b].size();
            for (int i = 0; i < sz; ++i)
                yell_parameters.push_back(params[b][i]);
        }
    }

    last_eval_params_ = yell_parameters;
    calc->calculate(yell_parameters);

    int n_obs = calc->number_of_observations();
    const bool use_asu = calc->refine_in_asu();
    const vector<int>& asu = calc->asu_indices();
    double num = 0.0, den = 0.0;
    for (int ii = 0; ii < n_obs; ++ii) {
        int i = use_asu ? asu[ii] : ii;
        double w = weights->at(i);
        double Ic = calc->get_intensity_map().at(i) - calc->get_average_intensity_map().at(i);
        double Ie = experimental_data->at(i);
        num += w * w * Ie * Ic;
        den += w * w * Ic * Ic;
    }
    double S = (den > 1e-15) ? (num / den) : 1.0;
    if (S < 0) S = 0;
    yell_parameters[0] = S;
    if (model) model->set_scale(S);

    if (use_asu) {
        for(int ii=0; ii<n_obs; ii++) {
            auto i = asu[ii];
            residuals[ii] = (experimental_data->at(i) - S * (calc->get_intensity_map().at(i) - calc->get_average_intensity_map().at(i))) * weights->at(i);
        }
    } else {
        for(int i=0; i<n_obs; i++)
            residuals[i] = (experimental_data->at(i) - S * (calc->get_intensity_map().at(i) - calc->get_average_intensity_map().at(i))) * weights->at(i);
    }

    return true;
}
