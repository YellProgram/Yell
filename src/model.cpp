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

#include "model.h"
#include "InputFileParser.h"
#include "Calculator.h"
#include <Eigen/LU>
#include <sstream>
#include <unordered_map>
#include <complex>
#include "exceptions.h"

extern OutputHandler report;
typedef iterator_ Iterator;

Iterator line_start(Iterator pos, const Iterator& start) {
  if(pos!=start and *pos=='\n')
    --pos;
  while (pos!=start && *pos != '\n')
    --pos;

  return pos;
}
Iterator line_end(Iterator pos, const Iterator& end) {
  if(pos!=end and *pos=='\n')
    ++pos;
  while (pos!=end && *pos != '\n')
    ++pos;

  return pos;
}
Iterator step_from_newline(const Iterator& pos, const Iterator& end) {
  if(pos!=end and *pos=='\n')
    return pos + 1;
  else
    return pos;
}

void complain_about_error(Iterator start, Iterator current, Iterator end)
{
  //skip empty spaces till the line where the error happened
  InputParser a_parser;
  qi::parse(current,end,a_parser.skipper_no_assignement);

  current = step_from_newline(current,end);

  int line_number = 0;
  for(Iterator it=current; it!=start; it--)
    if(*it=='\n')
      line_number++;

  Iterator lstart = line_start(current, start);

  int error_pos_within_line = 0;
  for(Iterator t = step_from_newline(lstart,end); t!=current; ++t)
    ++error_pos_within_line;

  // Print two lines before the error line
  Iterator context_start = step_from_newline(line_start(line_start(lstart,start),start),end);
  
  Iterator lend = line_end(current,end);

  Iterator next_two_lines_start = step_from_newline(lend,end);
  Iterator next_two_lines_end = line_end(line_end(next_two_lines_start,end),end);

  std::string trouble_line_with_context(context_start,lend);
  std::string next_two_lines(next_two_lines_start,next_two_lines_end);

  std::ostringstream complain;
  
  complain << "Parsing failed with exception around line " << line_number+1 << ":\n\n" << trouble_line_with_context <<"\n";
  for(int i=0; i < error_pos_within_line - 1; ++i)
    complain << ' ';
  complain << "^-----roughly here\n";
  complain << next_two_lines << "\n";
  
  REPORT(ERROR) << complain.str();
}

void Model::parse_model_()
{
  InputParser a_parser;
  a_parser.add_model(this);
  a_parser.formula.initialize_refinable_variables(refined_variable_names, refinement_parameters);
  a_parser.formula.array = refinement_parameters;

  Iterator start = model.begin();
  Iterator end = model.end();
  bool r;
  try {
    r = qi::phrase_parse(start, end, a_parser, a_parser.skipper_no_assignement);
  } catch(const qi::expectation_failure<Iterator>& e) {
    complain_about_error(model.begin(), e.first, model.end());
    throw(TerminateProgram());
  }

  if (!r || start != end) {
    string rest(start, end);
    REPORT(ERROR) << "Parsing failed, the rest of the file is:\n" << rest;
    throw "Parsing failed";
  }

  model_parsed_ = true;
  // Snapshot the scatterer registry into a per-model form-factor cache.
  scatterer_list_ = ScattererList();
}

void Model::calculate(vector<double> params, bool average_flag)
{
  Eigen::VectorXd p;
  if (params.size() == refinement_parameters.size()) {
    refinement_parameters = params;
    p = Eigen::VectorXd::Map(params.data(), params.size());
    yell::EvaluationCache cache;
    for (auto& pad : parameterized_atoms_)
      pad.update(p, &cache);
  } else {
    p = Eigen::VectorXd::Map(refinement_parameters.data(), refinement_parameters.size());
  }
  for (auto* pool : pools)
    pool->pairs.clear();

  vector<AtomicPair> pairs;
  for(auto* pool : pools)
  {
    pool->invoke_correlators(p);
    pairs.insert(pairs.end(), pool->pairs.begin(), pool->pairs.end());
  }
  
  REPORT(FIRST_RUN) << "created " << pairs.size() << " pairs\n\n";
  
  if(dump_pairs)
    for(int i=0; i<pairs.size(); i++)
      REPORT(FIRST_RUN) << pairs[i].to_string(p) <<'\n';

  if(report_pairs_outside_pdf_grid)
  {
    Grid pdf_grid = grid.in_pdf_space();
    for(int i=0; i<pairs.size(); i++)
      if(!pairs[i].pair_is_withing(pdf_grid, p))
        REPORT(FIRST_RUN) << "Warning, pair is outside PDF grid " << pairs[i].to_string(p) << '\n';
   }
  
  pairs = cell.laue_symmetry.apply_patterson_symmetry(pairs, p);
  atomic_pairs = pairs;
  
  vector<PattersonPeak> full_peaks, avg_peaks;
  peaks_from_pairs(pairs, p, scatterer_list_, full_peaks, avg_peaks);
  calculate_from_peaks(full_peaks, avg_peaks);
}

void Model::calculate_from_peaks(const vector<PattersonPeak>& full_peaks,
                                 const vector<PattersonPeak>& avg_peaks)
{
  auto run_calc = [&](const vector<PattersonPeak>& peaks, IntensityMap& out, bool avg) {
    if(direct_diffuse_scattering_calculation) {
      vec3<int> sym_boundary;
      for(int i=0; i<3; ++i) sym_boundary[i] = out.size()[i] > 1;
      IntensityMap padded = out.padded(sym_boundary);
      IntnsityCalculator::calculate_scattering_from_patterson_peaks(peaks, scatterer_list_, padded);
      cell.laue_symmetry.apply_patterson_symmetry(padded);
      out.copy_from_padded(sym_boundary, padded);
    } else {
      vec3<int> sym_boundary;
      for(int i=0; i<3; ++i) sym_boundary[i] = out.size()[i] > 1 && !(periodic_boundaries[i]);
      IntensityMap recipr_padded = out.padded(padding);
      recipr_padded.invert_grid();
      IntensityMap padded = recipr_padded.padded(sym_boundary);
      IntnsityCalculator::calculate_patterson_map_from_pairs_f(full_peaks, avg_peaks, scatterer_list_, padded, avg, fft_grid_size, periodic_boundaries);
      cell.laue_symmetry.apply_patterson_symmetry(padded);
      recipr_padded.copy_from_padded(sym_boundary, padded);
      recipr_padded.invert();
      out.copy_from_padded(padding, recipr_padded);
    }
    apply_resolution_function_if_possible(out);
    apply_reciprocal_space_multipliers_if_possible(out);
  };

  run_calc(full_peaks, intensity_map, false);
  run_calc(avg_peaks,  average_intensity_map, true);

  intensity_map.to_reciprocal();
  average_intensity_map.to_reciprocal();
  report.calculation_is_finished();
}

IntensityMap Model::calculate_derivative_from_peaks(
    const vector<PattersonPeak>& full_peaks,
    const vector<PattersonPeak>& avg_peaks,
    const vector<PeakSusceptibility>& full_susc,
    const vector<PeakSusceptibility>& avg_susc,
    double scale)
{
  IntensityMap res(grid);
  auto run_deriv = [&](const vector<PattersonPeak>& peaks, const vector<PeakSusceptibility>& susc, IntensityMap& out, bool avg) {
    if(direct_diffuse_scattering_calculation) {
      vec3<int> sym_boundary;
      for(int i=0; i<3; ++i) sym_boundary[i] = out.size()[i] > 1;
      IntensityMap padded = out.padded(sym_boundary);
      IntnsityCalculator::calculate_scattering_derivative_from_patterson_peaks(peaks, susc, scatterer_list_, padded);
      cell.laue_symmetry.apply_patterson_symmetry(padded);
      out.copy_from_padded(sym_boundary, padded);
    } else {
      vec3<int> sym_boundary;
      for(int i=0; i<3; ++i) sym_boundary[i] = out.size()[i] > 1 && !(periodic_boundaries[i]);
      IntensityMap recipr_padded = out.padded(padding);
      recipr_padded.invert_grid();
      IntensityMap padded = recipr_padded.padded(sym_boundary);
      IntnsityCalculator::calculate_patterson_map_derivative_from_pairs_f(full_peaks, avg_peaks, full_susc, avg_susc, scatterer_list_, padded, avg, fft_grid_size, periodic_boundaries);
      cell.laue_symmetry.apply_patterson_symmetry(padded);
      recipr_padded.copy_from_padded(sym_boundary, padded);
      recipr_padded.invert();
      out.copy_from_padded(padding, recipr_padded);
    }
  };

  IntensityMap dI_full(grid);
  IntensityMap dI_avg(grid);
  run_deriv(full_peaks, full_susc, dI_full, false);
  run_deriv(avg_peaks,  avg_susc,  dI_avg,  true);

  for (int i = 0; i < res.size_1d(); ++i)
    res.at(i) = scale * (dI_full.at(i) - dI_avg.at(i));

  apply_resolution_function_if_possible(res);
  apply_reciprocal_space_multipliers_if_possible(res);
  res.to_reciprocal();
  return res;
}

IntensityMap Model::calculate_derivative(const vector<double>& params, int param_idx)
{
  if (param_idx == 0) {
    calculate(params);
    IntensityMap res(intensity_map);
    for (int i = 0; i < res.size_1d(); ++i)
      res.at(i) = intensity_map.at(i) - average_intensity_map.at(i);
    return res;
  }

  Eigen::VectorXd q = Eigen::VectorXd::Map(params.data(), params.size());
  double scale = params[0];

  yell::EvaluationCache cache;
  for (auto& pad : parameterized_atoms_) pad.update(q, &cache);
  for (auto* pool : pools) pool->pairs.clear();

  vector<AtomicPair> pairs;
  for (auto* pool : pools) {
    pool->invoke_correlators(q);
    pairs.insert(pairs.end(), pool->pairs.begin(), pool->pairs.end());
  }
  pairs = cell.laue_symmetry.apply_patterson_symmetry(pairs, q);

  vector<PattersonPeak> full_peaks, avg_peaks;
  peaks_from_pairs(pairs, q, scatterer_list_, full_peaks, avg_peaks);

  vector<PeakSusceptibility> full_susc, avg_susc;
  susceptibilities_from_pairs(pairs, q, param_idx, full_susc, avg_susc);

  return calculate_derivative_from_peaks(full_peaks, avg_peaks, full_susc, avg_susc, scale);
}

Eigen::MatrixXd Model::compute_analytical_jacobian_direct(
    const vector<double>& params,
    IntensityMap& exp_map,
    OptionalIntensityMap& wts)
{
  const int n_params = (int)params.size();
  const int n_obs    = number_of_observations();
  const double scale = params[0];

  Eigen::VectorXd q = Eigen::VectorXd::Map(params.data(), n_params);

  for (auto* pool : pools) pool->pairs.clear();

  vector<AtomicPair> pairs;
  for(auto* pool : pools) {
    pool->invoke_correlators(q);
    pairs.insert(pairs.end(), pool->pairs.begin(), pool->pairs.end());
  }
  pairs = cell.laue_symmetry.apply_patterson_symmetry(pairs, q);

  vector<PattersonPeak> full_peaks, avg_peaks;
  peaks_from_pairs(pairs, q, scatterer_list_, full_peaks, avg_peaks);

  Eigen::MatrixXd J(n_obs, n_params);
  const bool use_asu = refine_in_asu();
  const vector<int>& asu = asu_indices();

  calculate_from_peaks(full_peaks, avg_peaks);
  IntensityMap cur_I_full = intensity_map;
  IntensityMap cur_I_avg  = average_intensity_map;

  for (int ii = 0; ii < n_obs; ++ii) {
    int i = use_asu ? asu[ii] : ii;
    double w = wts.at(i);
    J(ii, 0) = -(cur_I_full.at(i) - cur_I_avg.at(i)) * w;
  }

  if (n_params > 1) {
    for (int j = 1; j < n_params; ++j) {
      vector<PeakSusceptibility> full_susc, avg_susc;
      susceptibilities_from_pairs(pairs, q, j, full_susc, avg_susc);

      IntensityMap dI = calculate_derivative_from_peaks(full_peaks, avg_peaks, full_susc, avg_susc, scale);

      for (int ii = 0; ii < n_obs; ++ii) {
        int i = use_asu ? asu[ii] : ii;
        double w = wts.at(i);
        J(ii, j) = -dI.at(i) * w;
      }
    }
  }

  return J;
}

Eigen::MatrixXd Model::compute_jacobian_mixed(
    const vector<double>& params,
    IntensityMap& exp_map,
    OptionalIntensityMap& wts)
{
  const int n_params = (int)params.size();
  const int n_obs    = number_of_observations();
  const double scale = params[0];
  const double eps   = 1e-6;

  Eigen::VectorXd q = Eigen::VectorXd::Map(params.data(), n_params);

  for (auto* pool : pools) pool->pairs.clear();

  vector<AtomicPair> pairs;
  for(auto* pool : pools) {
    pool->invoke_correlators(q);
    pairs.insert(pairs.end(), pool->pairs.begin(), pool->pairs.end());
  }
  pairs = cell.laue_symmetry.apply_patterson_symmetry(pairs, q);
  
  vector<PattersonPeak> full_peaks_base, avg_peaks_base;
  peaks_from_pairs(pairs, q, scatterer_list_, full_peaks_base, avg_peaks_base);

  calculate_from_peaks(full_peaks_base, avg_peaks_base);
  IntensityMap base_I_full = intensity_map;
  IntensityMap base_I_avg  = average_intensity_map;

  Eigen::MatrixXd J(n_obs, n_params);
  const bool use_asu = refine_in_asu();
  const vector<int>& asu = asu_indices();

  for (int ii = 0; ii < n_obs; ++ii) {
    int i = use_asu ? asu[ii] : ii;
    double w = wts.at(i);
    J(ii, 0) = -(base_I_full.at(i) - base_I_avg.at(i)) * w;
  }

  for (int j = 1; j < n_params; ++j) {
    vector<PeakSusceptibility> full_susc, avg_susc;
    susceptibilities_from_pairs(pairs, q, j, full_susc, avg_susc);

    vector<PattersonPeak> full_peaks_pert = full_peaks_base;
    vector<PattersonPeak> avg_peaks_pert  = avg_peaks_base;

    for (size_t k = 0; k < full_peaks_pert.size(); ++k) {
      full_peaks_pert[k].coefficient += eps * full_susc[k].d_coefficient;
      full_peaks_pert[k].r           += eps * full_susc[k].d_r;
      full_peaks_pert[k].U           += eps * full_susc[k].d_U;

      avg_peaks_pert[k].coefficient  += eps * avg_susc[k].d_coefficient;
      avg_peaks_pert[k].r            += eps * avg_susc[k].d_r;
      avg_peaks_pert[k].U            += eps * avg_susc[k].d_U;
    }

    calculate_from_peaks(full_peaks_pert, avg_peaks_pert);
    
    for (int ii = 0; ii < n_obs; ++ii) {
      int i = use_asu ? asu[ii] : ii;
      double w = wts.at(i);
      double dI = (intensity_map.at(i) - average_intensity_map.at(i)) - (base_I_full.at(i) - base_I_avg.at(i));
      J(ii, j) = -scale * (dI / eps) * w;
    }
  }

  intensity_map = base_I_full;
  average_intensity_map = base_I_avg;

  return J;
}

Eigen::MatrixXd Model::compute_full_covariance(
    const vector<double>& params,
    IntensityMap& exp_map,
    OptionalIntensityMap& wts)
{
  const int n_params = (int)params.size(); 
  const int n_obs    = number_of_observations();
  const bool use_asu = refine_in_asu();
  const vector<int>& asu = asu_indices();

  // Update model to current optimum
  calculate(params);
  const double S = refinement_parameters[0];

  // H = J^T * W^2 * J
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n_params, n_params);

  // Streaming approach: process the Jacobian one column at a time to minimize memory.
  // Column 0 is the Scale derivative map: J_S = -(I_full - I_avg)
  vector<double> col0(n_obs);
  for (int ii = 0; ii < n_obs; ++ii) {
      int i = use_asu ? asu[ii] : ii;
      col0[ii] = -(intensity_map.at(i) - average_intensity_map.at(i));
  }

  // Helper to get weight at a residual index
  auto get_w = [&](int ii) {
      return wts.at(use_asu ? asu[ii] : ii);
  };

  // 1. Compute H(0,0) and H(0, j)
  for (int ii = 0; ii < n_obs; ++ii) {
      double w = get_w(ii);
      H(0, 0) += w * w * col0[ii] * col0[ii];
  }

  // 2. Compute other columns
  for (int j = 1; j < n_params; ++j) {
      IntensityMap dI_map = calculate_derivative(params, j);
      for (int ii = 0; ii < n_obs; ++ii) {
          int i = use_asu ? asu[ii] : ii;
          double w = get_w(ii);
          double J_j = -dI_map.at(i);
          
          H(0, j) += w * w * col0[ii] * J_j;
          H(j, 0) = H(0, j);
          
          // Diagonal term
          H(j, j) += w * w * J_j * J_j;
      }
      
      // For cross-terms H(j, k) with k < j, we would need to re-read or store maps.
      // For 5000 params, we'll re-calculate to save memory, though it is slow.
      // Optimization: calculate_derivative is fast (FFT path).
      for (int k = 1; k < j; ++k) {
          IntensityMap dI_map_k = calculate_derivative(params, k);
          for (int ii = 0; ii < n_obs; ++ii) {
              int i = use_asu ? asu[ii] : ii;
              double w = get_w(ii);
              H(k, j) += w * w * (-dI_map_k.at(i)) * (-dI_map.at(i));
          }
          H(j, k) = H(k, j);
      }
  }

  return H.inverse();
}
