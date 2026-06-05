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
#include <Eigen/SVD>
#include <limits>
#include <cmath>
#include <sstream>
#include <unordered_map>
#include <complex>
#include "exceptions.h"
#include <thread>
#include <mutex>
#include <atomic>

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
    for (auto* a : all_atoms_)
      a->update_caches(p, &cache);
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
      // Main path: use all available hardware threads
      IntnsityCalculator::calculate_patterson_map_from_pairs_f(full_peaks, avg_peaks, scatterer_list_, padded, avg, fft_grid_size, periodic_boundaries, fft_border_pixels, max_processors);
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
    double scale,
    int num_threads,
    ScattererList& sl)
{
  IntensityMap res(grid);
  auto run_deriv = [&](const vector<PattersonPeak>& peaks, const vector<PeakSusceptibility>& susc, IntensityMap& out, bool avg) {
    if(direct_diffuse_scattering_calculation) {
      vec3<int> sym_boundary;
      for(int i=0; i<3; ++i) sym_boundary[i] = out.size()[i] > 1;
      IntensityMap padded = out.padded(sym_boundary);
      IntnsityCalculator::calculate_scattering_derivative_from_patterson_peaks(peaks, susc, sl, padded);
      cell.laue_symmetry.apply_patterson_symmetry(padded);
      out.copy_from_padded(sym_boundary, padded);
    } else {
      vec3<int> sym_boundary;
      for(int i=0; i<3; ++i) sym_boundary[i] = out.size()[i] > 1 && !(periodic_boundaries[i]);
      IntensityMap recipr_padded = out.padded(padding);
      recipr_padded.invert_grid();
      IntensityMap padded = recipr_padded.padded(sym_boundary);
      IntnsityCalculator::calculate_patterson_map_derivative_from_pairs_f(full_peaks, avg_peaks, full_susc, avg_susc, sl, padded, avg, fft_grid_size, periodic_boundaries, fft_border_pixels, num_threads);
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

IntensityMap Model::calculate_derivative_from_susceptibilities(
    const vector<PattersonPeak>& full_peaks,
    const vector<PattersonPeak>& avg_peaks,
    const vector<AtomicPair>& pairs,
    const Eigen::VectorXd& q,
    int param_idx,
    double scale,
    int num_threads,
    ScattererList& sl)
{
    vector<PeakSusceptibility> full_susc, avg_susc;
    susceptibilities_from_pairs(const_cast<vector<AtomicPair>&>(pairs), q, param_idx, full_susc, avg_susc);
    return calculate_derivative_from_peaks(full_peaks, avg_peaks, full_susc, avg_susc, scale, num_threads, sl);
}

IntensityMap Model::calculate_derivative(const vector<double>& params, int param_idx, int num_threads)
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
  for (auto* a : all_atoms_) a->update_caches(q, &cache);
  for (auto* pool : pools) pool->pairs.clear();

  vector<AtomicPair> pairs;
  for (auto* pool : pools) {
    pool->invoke_correlators(q);
    pairs.insert(pairs.end(), pool->pairs.begin(), pool->pairs.end());
  }
  pairs = cell.laue_symmetry.apply_patterson_symmetry(pairs, q);

  vector<PattersonPeak> full_peaks, avg_peaks;
  peaks_from_pairs(pairs, q, scatterer_list_, full_peaks, avg_peaks);

  return calculate_derivative_from_susceptibilities(full_peaks, avg_peaks, pairs, q, param_idx, scale, num_threads, scatterer_list_);
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
      IntensityMap dI = calculate_derivative_from_susceptibilities(full_peaks, avg_peaks, pairs, q, j, scale, 0, scatterer_list_);
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

double Model::compute_optimal_scale(IntensityMap& exp_map, OptionalIntensityMap& wts)
{
    const int n_obs    = number_of_observations();
    const bool use_asu = refine_in_asu();
    const vector<int>& asu = asu_indices();
    double num = 0.0, den = 0.0;
    for (int ii = 0; ii < n_obs; ++ii) {
        int i = use_asu ? asu[ii] : ii;
        double w  = wts.at(i);
        double Ic = intensity_map.at(i) - average_intensity_map.at(i);
        double Ie = exp_map.at(i);
        num += w * w * Ie * Ic;
        den += w * w * Ic * Ic;
    }
    double S = (den > 1e-15) ? (num / den) : 1.0;
    return (S > 0) ? S : 0.0;
}

void Model::init_background_basis(OptionalIntensityMap& wts)
{
    const int m = background_n_terms();
    // Seed coefficients from the supplied values (Background ...); for a refined
    // background these are the starting values (typically zeros).
    background_coeffs_ = background_input_coeffs_;
    if (m == 0) { background_basis_.resize(0, 0); return; }

    // |q_i| = √(d*²_i) over the full grid, in the same flat order as at(i).
    IntensityMap tmp(grid);
    tmp.init_iterator();
    std::vector<double> s;
    s.reserve(tmp.size_1d());
    while (tmp.next())
        s.push_back(std::sqrt(std::max(0.0, tmp.current_d_star_square())));
    const int n = (int)s.size();

    // A point is "measured" where data exists (weight > 0) and, if a reciprocal-space
    // multiplier is loaded, where it is non-zero (the multiplier marks unmeasured
    // regions of reciprocal space with 0).
    const bool have_mult = reciprocal_space_multiplier.is_loaded;
    auto measured = [&](int i) {
        return wts.at(i) > 0 && (!have_mult || reciprocal_space_multiplier.at(i) > 0);
    };

    // Normalisation domain: |q| range over measured points only.
    double smin = std::numeric_limits<double>::max(), smax = -smin;
    for (int i = 0; i < n; ++i)
        if (measured(i)) { smin = std::min(smin, s[i]); smax = std::max(smax, s[i]); }
    if (!(smax > smin)) smax = smin + 1.0; // degenerate guard

    // x_i = a·s_i + b maps [smin, smax] → [-1, 1].
    const double a = 2.0 / (smax - smin);
    const double b = -(smax + smin) / (smax - smin);

    background_basis_.resize(n, m);
    for (int i = 0; i < n; ++i) {
        double x = a * s[i] + b;
        if (x < -1.0) x = -1.0; else if (x > 1.0) x = 1.0;
        // Chebyshev recurrence: T0=1, T1=x, T_k = 2x·T_{k-1} − T_{k-2}.
        background_basis_(i, 0) = 1.0;
        if (m > 1) background_basis_(i, 1) = x;
        double Tkm1 = 1.0, Tk = x;
        for (int k = 2; k < m; ++k) {
            double Tk1 = 2.0 * x * Tk - Tkm1;
            background_basis_(i, k) = Tk1;
            Tkm1 = Tk; Tk = Tk1;
        }
        // Mask by the reciprocal-space multiplier so the background, like the calculated
        // intensity (which already has it applied), is absent/scaled in unmeasured regions.
        if (have_mult) {
            double mult = reciprocal_space_multiplier.at(i);
            for (int k = 0; k < m; ++k) background_basis_(i, k) *= mult;
        }
    }

    REPORT(MAIN) << "Background: " << m << " Chebyshev term(s) in |q| over ["
                 << smin << ", " << smax << "] A^-1.\n";
}

void Model::compute_optimal_linear_params(IntensityMap& exp_map, OptionalIntensityMap& wts)
{
    const bool fit_scale = refinement_options.refine_scale;
    const int  n_bg      = (refinement_options.refine_background) ? (int)background_coeffs_.size() : 0;
    const bool fit_bg    = n_bg > 0;
    const int  m         = (fit_scale ? 1 : 0) + (fit_bg ? n_bg : 0);
    if (m == 0) return; // nothing linear to fit; scale held fixed

    const int n_obs    = number_of_observations();
    const bool use_asu = refine_in_asu();
    const vector<int>& asu = asu_indices();

    // Weighted normal equations: (Σ w² c cᵀ) θ = Σ w² t c, where c is the design row
    // [Ic (if fit_scale), Φ_i0..Φ_i,n_bg-1 (if fit_bg)] and t is the target
    // (Ie, minus the fixed-scale contribution when Scale is not being fit).
    Eigen::MatrixXd N = Eigen::MatrixXd::Zero(m, m);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(m);
    Eigen::VectorXd c(m);

    for (int ii = 0; ii < n_obs; ++ii) {
        int i = use_asu ? asu[ii] : ii;
        double w  = wts.at(i);
        double w2 = w * w;
        double Ic = intensity_map.at(i) - average_intensity_map.at(i);

        int col = 0;
        if (fit_scale) c[col++] = Ic;
        if (fit_bg) for (int k = 0; k < n_bg; ++k) c[col++] = background_basis_(i, k);

        double t = exp_map.at(i);
        if (!fit_scale) t -= scale_ * Ic;                  // Scale fixed
        if (!fit_bg && !background_coeffs_.empty())
            t -= background_at(i);                          // background fixed (calculated, not refined)

        N.noalias() += w2 * (c * c.transpose());
        rhs.noalias() += (w2 * t) * c;
    }

    // Chebyshev keeps N well-conditioned; LDLT is fine, fall back to SVD if singular.
    Eigen::VectorXd theta = N.ldlt().solve(rhs);
    if (!theta.allFinite())
        theta = N.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(rhs);

    int col = 0;
    if (fit_scale) { double S = theta[col++]; set_scale(S > 0 ? S : 0.0); }
    if (fit_bg) for (int k = 0; k < n_bg; ++k) background_coeffs_[k] = theta[col++];
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

  calculate(params);

  // Background coefficients (if refined) extend the parameter space after the
  // structural params: covariance indices [n_params .. n_params+n_bg).
  // n_bg counts enabled background terms so the covariance has a slot for each (a fixed
  // background gets a zero row/col → ESD 0, like a fixed Scale). fit_bg gates whether the
  // block is actually accumulated (only when the background is being refined).
  const int  n_bg    = (int)background_coeffs_.size();
  const bool fit_bg  = refinement_options.refine_background && n_bg > 0;
  const int  n_total = n_params + n_bg;
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n_total, n_total);

  // Scale column (param index 0).  Only accumulated when Scale is refined; when it
  // is held fixed (RefineScale off) row/col 0 stays zero, so the SVD pseudo-inverse
  // returns zero variance for Scale and the structural block becomes the covariance
  // conditional on the fixed Scale.
  const bool scale_refined = refinement_options.refine_scale;
  vector<double> col0(n_obs);
  for (int ii = 0; ii < n_obs; ++ii) {
      int i = use_asu ? asu[ii] : ii;
      col0[ii] = -(intensity_map.at(i) - average_intensity_map.at(i));
  }
  auto get_w = [&](int ii) { return wts.at(use_asu ? asu[ii] : ii); };
  if (scale_refined) {
      for (int ii = 0; ii < n_obs; ++ii) {
          double w = get_w(ii);
          H(0, 0) += w * w * col0[ii] * col0[ii];
      }
  }

  // Background columns: ∂model/∂b_k = Φ_ik, so the (negated) column is -Φ_ik.
  // Fill the scale-bg and bg-bg blocks now (independent of structural params);
  // struct-bg cross terms are added inside accumulate() below.
  vector<vector<double>> bgcol;
  if (fit_bg) {
      bgcol.assign(n_bg, vector<double>(n_obs));
      for (int k = 0; k < n_bg; ++k)
          for (int ii = 0; ii < n_obs; ++ii) {
              int i = use_asu ? asu[ii] : ii;
              bgcol[k][ii] = -background_basis_(i, k);
          }
      for (int k = 0; k < n_bg; ++k) {
          int K = n_params + k;
          if (scale_refined) {
              double h = 0.0;
              for (int ii = 0; ii < n_obs; ++ii) { double w = get_w(ii); h += w * w * col0[ii] * bgcol[k][ii]; }
              H(0, K) = H(K, 0) = h;
          }
          for (int l = 0; l <= k; ++l) {
              double h = 0.0;
              for (int ii = 0; ii < n_obs; ++ii) { double w = get_w(ii); h += w * w * bgcol[l][ii] * bgcol[k][ii]; }
              H(K, n_params + l) = H(n_params + l, K) = h;
          }
      }
  }

  Eigen::VectorXd q = Eigen::VectorXd::Map(params.data(), n_params);
  const double scale = params[0];

  vector<PattersonPeak> fpeaks, apeaks;
  peaks_from_pairs(atomic_pairs, q, scatterer_list_, fpeaks, apeaks);

  // Collect active parameter indices (skip param 0 = Scale, handled above)
  vector<int> active;
  for (int j = 1; j < n_params; ++j)
      if (param_is_active(j)) active.push_back(j);
  const int n_active = (int)active.size();

  if (n_active == 0) {
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);
      double thr = 1e-12 * svd.singularValues()(0);
      Eigen::VectorXd inv_sv = svd.singularValues();
      for (int i = 0; i < inv_sv.size(); ++i)
          inv_sv[i] = (inv_sv[i] > thr) ? 1.0 / inv_sv[i] : 0.0;
      return svd.matrixV() * inv_sv.asDiagonal() * svd.matrixU().transpose();
  }

  // Memory budget: covariance_batch_size maps at a time.
  // If all active params fit, do one pass. Otherwise split into two halves
  // (3 passes total — no quadratic recomputation).
  int budget = (covariance_batch_size > 0) ? covariance_batch_size : n_active;
  int half   = (n_active + 1) / 2;
  int block_size = (n_active <= budget) ? n_active : half;

  int n_threads = max_processors > 0 ? max_processors : (int)std::thread::hardware_concurrency();
  if (n_threads <= 0) n_threads = 1;

  // Ensure the Model owns one ScattererList per thread, each with its own
  // gridded_form_factors_ cache populated from the forward pass above.
  // This avoids any data race and avoids recomputing form factors per-iteration
  // (they only change when form-factor parameters change, which is the Model's
  // responsibility to invalidate via thread_scatterer_lists_.clear()).
  if ((int)thread_scatterer_lists_.size() != n_threads)
      thread_scatterer_lists_.assign(n_threads, scatterer_list_);

  // Compute derivative columns for active[from .. from+count) in parallel.
  // Each column is a compact n_obs-length vector (ASU-indexed when use_asu),
  // with the negation already applied — full IntensityMap is discarded immediately.
  auto compute_block = [&](int from, int count) {
      vector<vector<double>> cols(count, vector<double>(n_obs));
      std::atomic<int> next(0);
      vector<std::thread> workers;
      for (int t = 0; t < n_threads; ++t) {
          workers.emplace_back([&, from, count, t]() {
              ScattererList& my_sl = thread_scatterer_lists_[t];
              while (true) {
                  int idx = next.fetch_add(1);
                  if (idx >= count) break;
                  IntensityMap dI = calculate_derivative_from_susceptibilities(
                      fpeaks, apeaks, atomic_pairs, q, active[from + idx], scale, 1, my_sl);
                  for (int ii = 0; ii < n_obs; ++ii)
                      cols[idx][ii] = -dI.at(use_asu ? asu[ii] : ii);
              }
          });
      }
      for (auto& w : workers) w.join();
      return cols;
  };

  // Accumulate H for two sets of columns (from_j >= from_k convention so k<=j).
  // Also fills H(0,j) from col0. Columns are already ASU-extracted and negated.
  auto accumulate = [&](int from_j, const vector<vector<double>>& mj,
                         int from_k, const vector<vector<double>>& mk) {
      int nj = (int)mj.size(), nk = (int)mk.size();
      for (int jj = 0; jj < nj; ++jj) {
          int j = active[from_j + jj];
          if (scale_refined) {
              double h0j = 0.0;
              for (int ii = 0; ii < n_obs; ++ii) {
                  double w = get_w(ii);
                  h0j += w * w * col0[ii] * mj[jj][ii];
              }
              H(0, j) += h0j;
              H(j, 0)  = H(0, j);
          }

          // Structural-background cross terms. Only on the diagonal passes
          // (from_j == from_k) so each active j contributes exactly once.
          if (fit_bg && from_j == from_k) {
              for (int k = 0; k < n_bg; ++k) {
                  int K = n_params + k;
                  double h = 0.0;
                  for (int ii = 0; ii < n_obs; ++ii) {
                      double w = get_w(ii);
                      h += w * w * bgcol[k][ii] * mj[jj][ii];
                  }
                  H(K, j) += h;
                  H(j, K)  = H(K, j);
              }
          }

          for (int kk = 0; kk < nk; ++kk) {
              int k = active[from_k + kk];
              if (k > j) continue; // lower triangle only
              double h = 0.0;
              for (int ii = 0; ii < n_obs; ++ii) {
                  double w = get_w(ii);
                  h += w * w * mk[kk][ii] * mj[jj][ii];
              }
              H(k, j) += h;
              if (k != j) H(j, k) = h;
          }
      }
  };

  if (n_active <= budget) {
      // All fit in one pass — compute and accumulate everything at once.
      auto maps = compute_block(0, n_active);
      accumulate(0, maps, 0, maps);
  } else {
      // Two halves, three passes. No quadratic recomputation.
      int n1 = half, n2 = n_active - half;

      auto maps1 = compute_block(0, n1);
      accumulate(0, maps1, 0, maps1);         // H[A,A]

      auto maps2 = compute_block(half, n2);
      accumulate(half, maps2, half, maps2);   // H[B,B]

      // Cross term H[B,A]: keep maps2, recompute maps1 (one extra pass).
      auto maps1b = compute_block(0, n1);
      accumulate(half, maps2, 0, maps1b);     // H[B,A]
  }

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);
  double threshold = 1e-12 * svd.singularValues()(0);
  Eigen::VectorXd inv_sv = svd.singularValues();
  for (int i = 0; i < inv_sv.size(); ++i)
      inv_sv[i] = (inv_sv[i] > threshold) ? 1.0 / inv_sv[i] : 0.0;
  return svd.matrixV() * inv_sv.asDiagonal() * svd.matrixU().transpose();
}
