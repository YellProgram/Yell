//
// Created by Arkadiy Simonov on 05.08.24.
//

#include "CeresMinimizer.h"
#include "model.h"
#include "ceres/ceres.h"
#include "glog/logging.h"
#include "IntensityMap.h"
#include <fstream>

// ─────────────────────────────────────────────────────────────────────────────
// Analytical cost function: computes residuals and Jacobian via ExprPtr trees.
// Block 0 is always Scale (size 1).
// Blocks 1..N are structural parameters from user input.
// ─────────────────────────────────────────────────────────────────────────────
class AnalyticalYellCostFunction : public ceres::CostFunction {
public:
    AnalyticalYellCostFunction(Model* model, IntensityMap* exp, OptionalIntensityMap* weights)
        : model_(model), exp_(exp), weights_(weights)
    {
        set_num_residuals(model->number_of_observations());
        
        // Block 0: Scale (always size 1)
        mutable_parameter_block_sizes()->push_back(1);
        
        // Structural blocks
        for (auto& block : model->parameter_blocks) {
            mutable_parameter_block_sizes()->push_back((int)block.size());
        }
    }

    bool Evaluate(double const* const* parameters,
                  double* residuals,
                  double** jacobians) const override
    {
        // 1. Flatten blocks into a single parameter vector for model evaluation
        vector<double> p;
        for (size_t b = 0; b < parameter_block_sizes().size(); ++b) {
            int sz = parameter_block_sizes()[b];
            p.insert(p.end(), parameters[b], parameters[b] + sz);
        }

        // 2. Base model calculation
        model_->calculate(p);

        int n_obs = model_->number_of_observations();
        const bool use_asu = model_->refine_in_asu();
        const vector<int>& asu = model_->asu_indices();

        if (use_asu) {
            for (int ii = 0; ii < n_obs; ++ii) {
                int i = asu[ii];
                residuals[ii] = (exp_->at(i) - model_->data().at(i)) * weights_->at(i);
            }
        } else {
            for (int i = 0; i < n_obs; ++i)
                residuals[i] = (exp_->at(i) - model_->data().at(i)) * weights_->at(i);
        }

        // 3. Streaming Jacobian calculation
        if (jacobians) {
            double scale = p[0];
            
            // Block 0: Scale Derivative (dI/dScale = I_full - I_avg)
            if (jacobians[0]) {
                for (int ii = 0; ii < n_obs; ++ii) {
                    int i = use_asu ? asu[ii] : ii;
                    double dI = model_->intensity_map.at(i) - model_->average_intensity_map.at(i);
                    jacobians[0][ii] = -dI * weights_->at(i);
                }
            }

            // Blocks 1..N: Structural Derivatives
            // We reuse memory maps to keep peak usage low.
            int global_offset = 1;
            for (size_t b = 1; b < parameter_block_sizes().size(); ++b) {
                int block_sz = parameter_block_sizes()[b];
                if (jacobians[b]) {
                    for (int j = 0; j < block_sz; ++j) {
                        // calculate_derivative(p, global_idx) handles the full logic (Direct/FFT)
                        IntensityMap dI_map = model_->calculate_derivative(p, global_offset + j);
                        
                        for (int ii = 0; ii < n_obs; ++ii) {
                            int i = use_asu ? asu[ii] : ii;
                            jacobians[b][ii * block_sz + j] = -dI_map.at(i) * weights_->at(i);
                        }
                    }
                }
                global_offset += block_sz;
            }
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
    JsonIterationLogger(const vector<double*>& p_pointers,
                        const vector<int>& block_sizes,
                        int total_size, std::string filename)
        : p_pointers_(p_pointers),
          block_sizes_(block_sizes),
          total_size_(total_size),
          filename_(std::move(filename)) {}

    ceres::CallbackReturnType operator()(const ceres::IterationSummary& summary) override {
        IterationData data;
        data.iteration     = summary.iteration;
        data.cost          = summary.cost;
        data.gradient_norm = summary.gradient_norm;
        data.step_norm     = summary.step_norm;
        data.step_accepted = summary.step_is_successful;
        
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
    vector<double*> p_pointers_;
    vector<int>     block_sizes_;
    int total_size_;
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

    // 1. Scale Block (Size 1)
    double* p_scale = new double[1];
    p_scale[0] = initial_params[0];
    p_pointers.push_back(p_scale);
    block_sizes.push_back(1);

    // 2. Structural Blocks
    if (!model || model->parameter_blocks.empty()) {
        // If no blocks defined, treat all parameters > 0 as one block
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

    // 3. Problem Setup
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

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_QR; 
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = refinement_options.max_number_of_iterations;
    last_eval_params_.resize(parameters_number);
    options.update_state_every_iteration = true;
    
    JsonIterationLogger logger(p_pointers, block_sizes, parameters_number, "refinement_trajectory.json");
    options.callbacks.push_back(&logger);
    
    ceres::Solver::Summary summary;
    Solve(options, &problem, &summary);

    vector<double> result;
    for (size_t b = 0; b < p_pointers.size(); ++b) {
        for (int i = 0; i < block_sizes[b]; ++i) {
            result.push_back(p_pointers[b][i]);
        }
    }
    
    for (double* pb : p_pointers) delete[] pb;
    
    return result;
}

bool CeresMinimizer::operator()(double const *const *params, double *residuals) const {
    vector<double> yell_parameters;
    Model* model = dynamic_cast<Model*>(calc);
    
    // Always Block 0 is Scale
    yell_parameters.push_back(params[0][0]);

    if (!model || model->parameter_blocks.empty()) {
        // Fallback or legacy behavior: all structural params in block 1
        if (parameters_number > 1) {
            for (int i = 0; i < parameters_number - 1; ++i)
                yell_parameters.push_back(params[1][i]);
        }
    } else {
        for (size_t b = 0; b < model->parameter_blocks.size(); ++b) {
            int sz = (int)model->parameter_blocks[b].size();
            for (int i = 0; i < sz; ++i)
                yell_parameters.push_back(params[b + 1][i]);
        }
    }

    last_eval_params_ = yell_parameters;
    calc->calculate(yell_parameters);

    auto datapoints_number = calc->number_of_observations();
    if (calc->refine_in_asu()) {
        for(int ii=0; ii<datapoints_number; ii++) {
            auto i = calc->asu_indices()[ii];
            residuals[ii] = (experimental_data->at(i) - calc->data().at(i))*weights->at(i);
        }
    } else {
        for(int i=0; i<datapoints_number; i++)
            residuals[i] = (experimental_data->at(i) - calc->data().at(i))*weights->at(i);
    }

    return true;
}
