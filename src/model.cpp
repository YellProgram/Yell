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
      REPORT(FIRST_RUN) << pairs[i].to_string() <<'\n';

  
  if(report_pairs_outside_pdf_grid)
  {
    Grid pdf_grid = grid.in_pdf_space();
    for(int i=0; i<pairs.size(); i++)
      if(!pairs[i].pair_is_withing(pdf_grid))
        REPORT(FIRST_RUN) << "Warning, pair is outside PDF grid " << pairs[i].to_string() << '\n';
   }
  
  pairs = cell.laue_symmetry.apply_patterson_symmetry(pairs);
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

    ScattererList scatterers;
    vector<PattersonPeak> full_peaks, avg_peaks;
    peaks_from_pairs(pairs, scatterers, full_peaks, avg_peaks);
    const vector<PattersonPeak>& active_peaks = average_flag ? avg_peaks : full_peaks;
    IntnsityCalculator::calculate_scattering_from_patterson_peaks(active_peaks, scatterers, padded);

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
    IntnsityCalculator::calculate_patterson_map_from_pairs_f(pairs,padded,average_flag,fft_grid_size,periodic_boundaries);
    cell.laue_symmetry.apply_patterson_symmetry(padded);
    recipr_padded.copy_from_padded(sym_boundary,padded);
    recipr_padded.invert();
    calc_intensity_map->copy_from_padded(padding,recipr_padded);
  }

  apply_resolution_function_if_possible(*calc_intensity_map);

  apply_reciprocal_space_multipliers_if_possible(*calc_intensity_map);
  calc_intensity_map->to_reciprocal();
  report.calculation_is_finished();
}

Eigen::MatrixXd Model::compute_analytical_jacobian_direct(
    const vector<double>& params,
    IntensityMap& exp_map,
    OptionalIntensityMap& wts)
{
  // Residuals: r_i = (exp_i - data_i) * w_i,  data_i = Scale*(I_full_i - I_avg_i)
  //
  // ∂r_i/∂Scale = -(I_full_i - I_avg_i) * w_i
  // ∂r_i/∂p_j  = -Scale * (∂I_full_i/∂p_j - ∂I_avg_i/∂p_j) * w_i   (j>0)
  //
  // Direct-method intensity formula:
  //   ∂I(s)/∂p_j = Re[ Σ_pairs conj(f1)*f2*N*exp(M2PISQ*s·U·s + i*M_2PI*s·r)
  //                    * (dp_real_j + p_real*(M2PISQ*s·dU_j·s + i*M_2PI*s·dr_j)) ]
  // For isotropic ADP: s·dU_frac·s = dUiso * d_star_sq  (reciprocal metric).

  const int n_params = (int)params.size();
  const int n_obs    = number_of_observations();
  const double scale = params[0];

  Eigen::VectorXd q = Eigen::VectorXd::Map(params.data(), n_params);

  // Precompute ExprPtr gradient vectors for each parameterized atom.
  struct AtomDuals {
    Eigen::VectorXd dx, dy, dz, dUiso;
    bool isotropic;
  };
  std::unordered_map<Atom*, AtomDuals> atom_duals;
  for (auto& pad : parameterized_atoms_) {
    AtomDuals d;
    d.isotropic = pad.isotropic;
    d.dx = pad.param_exprs[1]->eval_d(q).derivatives();
    d.dy = pad.param_exprs[2]->eval_d(q).derivatives();
    d.dz = pad.param_exprs[3]->eval_d(q).derivatives();
    if (pad.isotropic)
      d.dUiso = pad.param_exprs[4]->eval_d(q).derivatives();
    atom_duals[pad.atom_ptr] = std::move(d);
  }

  // Precompute per-pair gradient vectors (n_pairs × n_params).
  const int n_pairs = (int)atomic_pairs.size();
  Eigen::MatrixXd dp_real_mat(n_pairs, n_params); dp_real_mat.setZero();
  vector<Eigen::VectorXd> dr_x(n_pairs, Eigen::VectorXd::Zero(n_params));
  vector<Eigen::VectorXd> dr_y(n_pairs, Eigen::VectorXd::Zero(n_params));
  vector<Eigen::VectorXd> dr_z(n_pairs, Eigen::VectorXd::Zero(n_params));
  vector<Eigen::VectorXd> dUiso_pair(n_pairs, Eigen::VectorXd::Zero(n_params));

  for (int k = 0; k < n_pairs; ++k) {
    AtomicPair& pair = atomic_pairs[k];
    if (pair.p_real_expr)
      dp_real_mat.row(k) = pair.p_real_expr->eval_d(q).derivatives().transpose();
    auto it1 = atom_duals.find(pair.atom1);
    auto it2 = atom_duals.find(pair.atom2);
    if (it1 != atom_duals.end()) {
      dr_x[k] -= it1->second.dx;
      dr_y[k] -= it1->second.dy;
      dr_z[k] -= it1->second.dz;
      if (it1->second.isotropic) dUiso_pair[k] += it1->second.dUiso;
    }
    if (it2 != atom_duals.end()) {
      dr_x[k] += it2->second.dx;
      dr_y[k] += it2->second.dy;
      dr_z[k] += it2->second.dz;
      if (it2->second.isotropic) dUiso_pair[k] += it2->second.dUiso;
    }
  }

  // Output Jacobian (n_obs × n_params).
  Eigen::MatrixXd J(n_obs, n_params);
  const bool use_asu = refine_in_asu();
  const vector<int>& asu = asu_indices();

  // Column 0: ∂r_i/∂Scale = -(I_full_i - I_avg_i) * w_i
  for (int ii = 0; ii < n_obs; ++ii) {
    int i = use_asu ? asu[ii] : ii;
    double w = wts.at(i);
    J(ii, 0) = -(intensity_map.at(i) - average_intensity_map.at(i)) * w;
  }

  if (n_params > 1) {
    // Accumulate ∂I_full and ∂I_avg into full-size arrays (n_pixels × n_params-1),
    // then select observations by ASU at the end.
    const int n_pixels = intensity_map.size_1d();
    Eigen::MatrixXd dI_full_all(n_pixels, n_params - 1); dI_full_all.setZero();
    Eigen::MatrixXd dI_avg_all (n_pixels, n_params - 1); dI_avg_all.setZero();

    IntensityMap iter_map = intensity_map; // copy to drive iteration
    iter_map.init_iterator();
    int pixel_1d = 0;
    while (iter_map.next()) {
      vec3<double> s = iter_map.current_s();
      double d_star_sq = iter_map.current_d_star_square();
      AtomicTypeCollection::update_current_form_factors(s, d_star_sq);

      for (int k = 0; k < n_pairs; ++k) {
        AtomicPair& pair = atomic_pairs[k];
        std::complex<double> f1 = pair.atomic_type1->current_form_factor;
        std::complex<double> f2 = pair.atomic_type2->current_form_factor;
        double N = pair.multiplier;

        double p_real = pair.p(false);
        std::complex<double> base_real = std::conj(f1) * f2 * N *
            std::exp(std::complex<double>(M2PISQ * (s * pair.U(false) * s),
                                          M_2PI  * (s * pair.r(false))));

        double p_avg = pair.p(true);
        std::complex<double> base_avg = std::conj(f1) * f2 * N *
            std::exp(std::complex<double>(M2PISQ * (s * pair.U(true) * s),
                                          M_2PI  * (s * pair.r(true))));

        for (int j = 1; j < n_params; ++j) {
          double dp_r   = dp_real_mat(k, j);
          double s_dr   = s[0]*dr_x[k][j] + s[1]*dr_y[k][j] + s[2]*dr_z[k][j];
          double s_dU_s = dUiso_pair[k][j] * d_star_sq;

          dI_full_all(pixel_1d, j-1) += std::real(
              base_real * std::complex<double>(dp_r + p_real * M2PISQ * s_dU_s,
                                               p_real * M_2PI * s_dr));
          dI_avg_all(pixel_1d, j-1) += std::real(
              base_avg  * std::complex<double>(p_avg * M2PISQ * s_dU_s,
                                               p_avg * M_2PI * s_dr));
        }
      }
      ++pixel_1d;
    }

    // Fill J columns 1..n_params-1, selecting by ASU.
    for (int ii = 0; ii < n_obs; ++ii) {
      int i = use_asu ? asu[ii] : ii;
      double w = wts.at(i);
      for (int j = 1; j < n_params; ++j)
        J(ii, j) = -scale * (dI_full_all(i, j-1) - dI_avg_all(i, j-1)) * w;
    }
  }

  return J;
}