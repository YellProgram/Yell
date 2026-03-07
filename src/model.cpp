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
  // Must come after parsing (all Scatterer types are registered by then).
  scatterer_list_ = ScattererList();
}

void Model::calculate(vector<double> params, bool average_flag)
{
  // Update atom parameters from ExprPtr trees when sizes match.
  // (Mismatched size means a legacy timing call — skip update.)
  Eigen::VectorXd p;
  if (params.size() == refinement_parameters.size()) {
    refinement_parameters = params;
    p = Eigen::VectorXd::Map(params.data(), params.size());
    for (auto& pad : parameterized_atoms_)
      pad.update(p);
  } else {
    p = Eigen::VectorXd::Map(refinement_parameters.data(), refinement_parameters.size());
  }
  // Clear cached pairs so invoke_correlators rebuilds them.
  for (auto* pool : pools)
    pool->pairs.clear();

  vector<AtomicPair> pairs;
  vector<AtomicPairPool*>::iterator pool;
  for(pool=pools.begin(); pool!=pools.end(); pool++)
  {
    (*pool)->invoke_correlators(p);
    pairs.insert(pairs.end(),(*pool)->pairs.begin(),(*pool)->pairs.end());
  }
  
//    pairs = cell.laue_symmetry.filter_pairs_from_asymmetric_unit(pairs); //We decided to use another way for multiplicity
  
  REPORT(FIRST_RUN) << "created " << pairs.size() << " pairs\n\n";
  
  // output pairs
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
  atomic_pairs = pairs; //save them to check that everything is fine
  
  IntensityMap* calc_intensity_map;
  if(average_flag==AVERAGE)
    calc_intensity_map=&average_intensity_map;
  else
    calc_intensity_map=&intensity_map;

  if(direct_diffuse_scattering_calculation)
  {
    vec3<int> sym_boundary;
    for(int i=0; i<3; ++i)
      sym_boundary[i]=calc_intensity_map->size()[i]>1;

    IntensityMap padded = calc_intensity_map->padded(sym_boundary);

    vector<PattersonPeak> full_peaks, avg_peaks;
    peaks_from_pairs(pairs, p, scatterer_list_, full_peaks, avg_peaks);
    const vector<PattersonPeak>& active_peaks = average_flag ? avg_peaks : full_peaks;
    IntnsityCalculator::calculate_scattering_from_patterson_peaks(active_peaks, scatterer_list_, padded);

    cell.laue_symmetry.apply_patterson_symmetry(padded);
    calc_intensity_map->copy_from_padded(sym_boundary,padded);
  }else
  {
    vec3<int> sym_boundary;
    for(int i=0; i<3; ++i)
      sym_boundary[i]=calc_intensity_map->size()[i]>1 && !(periodic_boundaries[i]);
    
    IntensityMap recipr_padded = calc_intensity_map->padded(padding); //for better fft accuracy
    recipr_padded.invert_grid();
    IntensityMap padded = recipr_padded.padded(sym_boundary); //for symmetry
    vector<PattersonPeak> full_peaks_fft, avg_peaks_fft;
    peaks_from_pairs(pairs, p, scatterer_list_, full_peaks_fft, avg_peaks_fft);
    IntnsityCalculator::calculate_patterson_map_from_pairs_f(full_peaks_fft,avg_peaks_fft,scatterer_list_,padded,average_flag,fft_grid_size,periodic_boundaries);
    cell.laue_symmetry.apply_patterson_symmetry(padded);
    recipr_padded.copy_from_padded(sym_boundary,padded);
    recipr_padded.invert();
    calc_intensity_map->copy_from_padded(padding,recipr_padded);
  }

  apply_resolution_function_if_possible(*calc_intensity_map);
apply_reciprocal_space_multipliers_if_possible(*calc_intensity_map);
}

IntensityMap Model::calculate_derivative(const vector<double>& params, int param_idx)
{
if (param_idx == 0) { // Scale derivative: ∂I/∂Scale = I_full - I_avg
  calculate(params);
  IntensityMap res(intensity_map);
  for (int i = 0; i < res.size_1d(); ++i)
    res.at(i) = intensity_map.at(i) - average_intensity_map.at(i);
  return res;
}

// Parameter derivative: ∂I/∂p_j = Scale * (∂I_full/∂p_j - ∂I_avg/∂p_j)
Eigen::VectorXd q = Eigen::VectorXd::Map(params.data(), params.size());
double scale = params[0];

for (auto& pad : parameterized_atoms_)
  pad.update(q);

for (auto* pool : pools)
  pool->pairs.clear();

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

IntensityMap dI_full(grid);
IntensityMap dI_avg(grid);

IntnsityCalculator::calculate_scattering_derivative_from_patterson_peaks(full_peaks, full_susc, scatterer_list_, dI_full);
IntnsityCalculator::calculate_scattering_derivative_from_patterson_peaks(avg_peaks,  avg_susc,  scatterer_list_, dI_avg);

IntensityMap res(grid);
for (int i = 0; i < res.size_1d(); ++i)
  res.at(i) = scale * (dI_full.at(i) - dI_avg.at(i));

apply_resolution_function_if_possible(res);
apply_reciprocal_space_multipliers_if_possible(res);

return res;
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

  // Snapshot the model parameters: create peak lists once.
  // Then for each pixel, we evaluate the full expression and get derivatives.
  // Clear cached pairs so invoke_correlators rebuilds them if needed.
  for (auto* pool : pools)
    pool->pairs.clear();

  vector<AtomicPair> pairs;
  for(auto* pool : pools)
  {
    pool->invoke_correlators(q);
    pairs.insert(pairs.end(), pool->pairs.begin(), pool->pairs.end());
  }
  pairs = cell.laue_symmetry.apply_patterson_symmetry(pairs, q);
  
  vector<PattersonPeak> full_peaks, avg_peaks;
  peaks_from_pairs(pairs, q, scatterer_list_, full_peaks, avg_peaks);

  // Output Jacobian (n_obs × n_params).
  Eigen::MatrixXd J(n_obs, n_params);
  const bool use_asu = refine_in_asu();
  const vector<int>& asu = asu_indices();

  // Column 0: ∂r_i/∂Scale = -(I_full_i - I_avg_i) * w_i
  IntensityMap cur_I_full = intensity_map;
  IntensityMap cur_I_avg  = average_intensity_map;
  IntnsityCalculator::calculate_scattering_from_patterson_peaks(full_peaks, scatterer_list_, cur_I_full);
  IntnsityCalculator::calculate_scattering_from_patterson_peaks(avg_peaks,  scatterer_list_, cur_I_avg);

  for (int ii = 0; ii < n_obs; ++ii) {
    int i = use_asu ? asu[ii] : ii;
    double w = wts.at(i);
    J(ii, 0) = -(cur_I_full.at(i) - cur_I_avg.at(i)) * w;
  }

  if (n_params > 1) {
    const int n_pixels = intensity_map.size_1d();
    Eigen::MatrixXd dI_full_all(n_pixels, n_params - 1); dI_full_all.setZero();
    Eigen::MatrixXd dI_avg_all (n_pixels, n_params - 1); dI_avg_all.setZero();

    // ∂I(s)/∂p_j = Σ_pairs conj(f1)*f2*N * ∂(exp(...))/∂p_j
    // ∂(p * exp(arg))/∂p_j = dp/dp_j * exp(arg) + p * exp(arg) * darg/dp_j
    //
    // Since we only want derivatives UNTIL THE LIST OF PAIRS here,
    // and the Intensity stage handles maps, I'll calculate the per-parameter
    // intensity derivative by iterating over pairs and evaluating their Duals.

    cur_I_full.init_iterator();
    int pixel_1d = 0;
    while (cur_I_full.next()) {
      vec3<double> s = cur_I_full.current_s();
      double d_star_sq = cur_I_full.current_d_star_square();
      scatterer_list_.update(s, d_star_sq);

      for (int k = 0; k < (int)pairs.size(); ++k) {
        AtomicPair& pair = pairs[k];
        std::complex<double> f1f2 = std::conj(scatterer_list_.f(full_peaks[k].type1_idx)) * scatterer_list_.f(full_peaks[k].type2_idx);

        auto calc_dI_pair = [&](bool avg) -> Eigen::VectorXd {
            yell::Dual p_dual   = pair.p(avg)->eval_d(q);
            yell::Dual rx_dual  = pair.r(avg).x->eval_d(q);
            yell::Dual ry_dual  = pair.r(avg).y->eval_d(q);
            yell::Dual rz_dual  = pair.r(avg).z->eval_d(q);
            yell::Dual u11_dual = pair.U(avg).u11->eval_d(q);
            yell::Dual u22_dual = pair.U(avg).u22->eval_d(q);
            yell::Dual u33_dual = pair.U(avg).u33->eval_d(q);
            yell::Dual u12_dual = pair.U(avg).u12->eval_d(q);
            yell::Dual u13_dual = pair.U(avg).u13->eval_d(q);
            yell::Dual u23_dual = pair.U(avg).u23->eval_d(q);

            yell::Dual phase = rx_dual * (M_2PI * s[0]) + ry_dual * (M_2PI * s[1]) + rz_dual * (M_2PI * s[2]);
            yell::Dual adp   = u11_dual * (M2PISQ * s[0]*s[0]) + u22_dual * (M2PISQ * s[1]*s[1]) + u33_dual * (M2PISQ * s[2]*s[2]) +
                               u12_dual * (2.0 * M2PISQ * s[0]*s[1]) + u13_dual * (2.0 * M2PISQ * s[0]*s[2]) + u23_dual * (2.0 * M2PISQ * s[1]*s[2]);
            
            double p_v = p_dual.value();
            double adp_v = adp.value();
            double phase_v = phase.value();
            std::complex<double> e_val = std::exp(std::complex<double>(adp_v, phase_v));
            
            Eigen::VectorXd res_g = Eigen::VectorXd::Zero(n_params);
            std::complex<double> term1 = f1f2 * e_val * pair.multiplier;
            std::complex<double> term2 = term1 * p_v;

            res_g += term1.real() * p_dual.derivatives();
            res_g += term2.real() * adp.derivatives();
            res_g -= term2.imag() * phase.derivatives();

            return res_g;
        };

        dI_full_all.row(pixel_1d) += calc_dI_pair(false).tail(n_params - 1);
        dI_avg_all.row(pixel_1d)  += calc_dI_pair(true).tail(n_params - 1);
      }
      ++pixel_1d;
    }

    for (int ii = 0; ii < n_obs; ++ii) {
      int i = use_asu ? asu[ii] : ii;
      double w = wts.at(i);
      for (int j = 1; j < n_params; ++j)
        J(ii, j) = -scale * (dI_full_all(i, j-1) - dI_avg_all(i, j-1)) * w;
    }
  }

  return J;
}