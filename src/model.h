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

#ifndef MODEL_H
#define MODEL_H

#include "basic_classes.h"
#include "precompiled_header.h"
#include "FormulaParser.h"
#include "ExprFormulaParser.h"
#include <set>
#include <sstream>
#include <unordered_map>
#include <boost/fusion/tuple.hpp>


//COPYPASTE from InputFileParser.h
typedef boost::variant<ChemicalUnit*,ChemicalUnitNode*,ADPMode*> StructurePartRef;
//#include "InputFileParser.h"

enum R_FACTORS {R1, R2};
enum WEIGHTED_OPTIONS {WEIGHTED, UNWEIGHTED};
enum DerivativesMode {FINITE_DIFFERENCE, ANALYTICAL, MIXED};


class Model : public MinimizerCalculator {
private:
  static double sq(double x) {
    return x*x;
  }

  // Flatten all Atom* from a UnitCell in a consistent traversal order.
  // Used by clone() to build the original→clone atom remap.
  static std::vector<Atom*> collect_atoms_(UnitCell& cell) {
    std::vector<Atom*> result;
    for (int i = 0; i < cell.chemical_unit_nodes.size(); ++i) {
      ChemicalUnitNode& node = cell.chemical_unit_nodes[i];
      for (int j = 0; j < node.chemical_units.size(); ++j)
        for (Atom* a : node.chemical_units[j].get_atoms())
          result.push_back(a);
    }
    return result;
  }
public:
  static void register_molecular_scatterers(vector<boost::tuple<string,Scatterer*> > scatterers) {
    for(int i=0; i<scatterers.size();++i)
      AtomicTypeCollection::add(boost::get<0>(scatterers[i]),boost::get<1>(scatterers[i]));
  }



  bool unit_cell_is_initialized() {
    if(!cell_is_initialized)
      REPORT(ERROR) << "Unit cell is not defined\n";
    if(!grid_initialized)
      REPORT(ERROR) << "Diffuse scattering grid is not defined\n";
    
    if(cell.laue_symmetry.label=="")
    {
      REPORT(ERROR) << "Point group is not defined\n";
      return false;
    }
    
    bool symmetry_is_ok = cell.laue_symmetry.is_compatible_with_cell(cell);
    if(!symmetry_is_ok)
      REPORT(ERROR) << "Point group " << cell.laue_symmetry.label << " is incomatible with unit cell " << cell.to_string() << "\n";
    
    return cell_is_initialized && grid_initialized && symmetry_is_ok;
  }
  
  void set_calculate_jacobians(bool inp) {
    calculate_jacobians=inp;
  }
  
  void set_print_covariance_matrix(bool inp){
    print_covariance_matrix=inp;
  }
  
  void set_padding(vector<int> inp) {
    for(int i=0; i<3; ++i)
      padding[i]=inp[i];
  }
  
  void set_report_pairs_outside_pdf_grid(bool inp)
  {
    report_pairs_outside_pdf_grid = inp;
  }
  
  void set_max_number_of_iterations(int inp)
  {
    refinement_options.max_number_of_iterations=inp;
  }
  
  void set_tau(double inp)
  {
    refinement_options.tau = inp;
  }
  
  void set_diff(double inp)
  {
    refinement_options.difference = inp;
  }
  void set_thresholds(vector<double> inp)
  {
    assert(inp.size()==3);
    for(int i=0;i<3;++i)
      refinement_options.thresholds[i]=inp[i];
  }

  void set_ceres_num_threads(int inp)     { refinement_options.num_threads = inp; }
  void set_function_tolerance(double inp) { refinement_options.function_tolerance = inp; }
  void set_gradient_tolerance(double inp) { refinement_options.gradient_tolerance = inp; }
  void set_use_dense_qr(bool use_qr) { refinement_options.use_dense_qr = use_qr; }
  void set_scale_before_refine(bool inp) { refinement_options.scale_before_refine = inp; }
  void set_num_supercycles(int inp)      { refinement_options.num_supercycles = inp; }

  double R_factor(IntensityMap& exp, R_FACTORS r,WEIGHTED_OPTIONS weighted)
  {
    double scale = refinement_parameters[0];
    
    double nom=0,denom=0,result;
    
    if(r==R1) 
      if(weighted==UNWEIGHTED)
      {
        for(int i=0; i<intensity_map.size_1d(); ++i)
        {
          nom += abs(exp.at(i)-scale*(intensity_map.at(i)-average_intensity_map.at(i)));
          denom += abs(exp.at(i));
        }
        result = nom/denom;
      }
      else // WEIGHTED not used
      {
        for(int i=0; i<intensity_map.size_1d(); ++i)
        {
          nom += abs( (exp.at(i)-scale*(intensity_map.at(i)-average_intensity_map.at(i)))*weights.at(i) );
          denom += abs(exp.at(i)*weights.at(i));
        }
        result = nom/denom;
      }
    else //R2
      if(weighted==WEIGHTED)
      {
        for(int i=0; i<intensity_map.size_1d(); ++i)
        {
          nom += sq( (exp.at(i)-scale*(intensity_map.at(i)-average_intensity_map.at(i))) )*weights.at(i) ;
          denom += sq(exp.at(i))*weights.at(i);
        }
        result = sqrt(nom/denom);
      }
      else // UNWEIGHTED not used
      {
        for(int i=0; i<intensity_map.size_1d(); ++i)
        {
          nom += sq(exp.at(i)-scale*(intensity_map.at(i)-average_intensity_map.at(i)));
          denom += sq(exp.at(i));
        }
        result = sqrt(nom/denom);
      }
    
    return result;
    
  }
  
  static void throw_error() {
    throw "error";
  }
  
  void set_periodic_boundaries(vector<bool> inp) {
      periodic_boundaries = inp;
  }
  
  /// Creates an atom with isotropic ADP from ExprPtr parameter trees.
  /// Bakes the ADP conversion (Uiso Å² → fractional) into the ExprPtr at parse time.
  /// param_exprs layout: [mult, x, y, z, Uiso_ang]
  Atom* construct_atom_isotropic_adp(string name, vector<yell::ExprPtr> param_exprs) {
    auto p_eig = refinement_parameters_asEig();
    // TODO: assert cell is initialized before construct_atom is called.
    auto rm = cell.cell.reciprocal_metrical_matrix();
    yell::ExprPtr U_exprs[6];
    for (int i = 0; i < 6; ++i)
        U_exprs[i] = param_exprs[4] * rm[i];  // Uiso * rm[i]
    Atom* atom = new Atom(name, scattering_type,
                          param_exprs[0],  // mult_expr
                          param_exprs[1], param_exprs[2], param_exprs[3],
                          U_exprs[0], U_exprs[1], U_exprs[2],
                          U_exprs[3], U_exprs[4], U_exprs[5]);
    atom->update_caches(p_eig);
    all_atoms_.push_back(atom);
    return atom;
  }

  /// Creates an atom with anisotropic ADP from ExprPtr parameter trees.
  /// Bakes the ADP conversion (U_ij Å² → fractional) into the ExprPtr at parse time.
  /// param_exprs layout: [mult, x, y, z, U11, U22, U33, U12, U13, U23] (ADP in Å²)
  Atom* construct_atom(string name, vector<yell::ExprPtr> param_exprs) {
    auto p_eig = refinement_parameters_asEig();
    // Convert U_ij (Å²) to U_stored[ij] = U_ij * a*_i * a*_j using scalar
    // reciprocal lattice lengths a* = sqrt(G*[i,i]).  Matches the SHELX/CIF
    // DWF convention: T = exp(-2pi^2 (U11 h^2 a*^2 + ... + 2U12 hk a*b* + ...))
    auto G = cell.cell.reciprocal_metrical_matrix();
    const double astar = std::sqrt(G[0]);
    const double bstar = std::sqrt(G[1]);
    const double cstar = std::sqrt(G[2]);
    yell::ExprPtr U_exprs[6] = {
        param_exprs[4] * (astar * astar),
        param_exprs[5] * (bstar * bstar),
        param_exprs[6] * (cstar * cstar),
        param_exprs[7] * (astar * bstar),
        param_exprs[8] * (astar * cstar),
        param_exprs[9] * (bstar * cstar)
    };
    Atom* atom = new Atom(name, scattering_type,
                          param_exprs[0],
                          param_exprs[1], param_exprs[2], param_exprs[3],
                          U_exprs[0], U_exprs[1], U_exprs[2],
                          U_exprs[3], U_exprs[4], U_exprs[5]);
    atom->update_caches(p_eig);
    all_atoms_.push_back(atom);
    return atom;
  }

  void set_derivatives_mode(DerivativesMode m) { derivatives_mode = m; }
  void set_jacobian_multiplier(double v)        { jacobian_multiplier = v; }

  void set_active_blocks_str(const vector<string>& tokens) {
    active_blocks.clear();
    for (const string& s : tokens) {
      size_t dash = s.find('-');
      if (dash != string::npos) {
        int i1 = std::stoi(s.substr(0, dash));
        int i2 = std::stoi(s.substr(dash + 1));
        for (int i = i1; i <= i2; ++i) active_blocks.insert(i);
      } else {
        active_blocks.insert(std::stoi(s));
      }
    }
    std::ostringstream oss;
    oss << "Active blocks set to: ";
    for (int b : active_blocks) oss << b << " ";
    REPORT(MAIN) << oss.str() << "\n";
  }

// TODO: fix this. Users care about 1-based indices, so every time we are talking to users we will add 1, internally keep it 0-based

  // True if the 1-based block index is active (empty set = all active).
  bool block_is_active(int idx_1based) const {
    if (active_blocks.empty()) return true;
    return active_blocks.count(idx_1based) > 0;
  }

  // True if the flat parameter index (0=Scale, 1..N=RefinableVariables) is active.
  bool param_is_active(int param_idx) const {
    if (active_blocks.empty() || param_idx == 0) return true;
    int offset = 1;
    for (size_t b = 0; b < parameter_blocks.size(); ++b) {
      int sz = (int)parameter_blocks[b].size();
      if (param_idx >= offset && param_idx < offset + sz)
        return block_is_active((int)b + 1);
      offset += sz;
    }
    return true;
  }

  void set_scale(double inp) {
    refinement_parameters[0]=inp;
  }
  
  void set_refinable_parameter_blocks(FormulaParser& formula, ExprFormulaParser& expr_formula,
                                     vector<vector<boost::fusion::tuple<string,double>>> blocks) {
    // Flatten for internal expression use, keep blocks for Ceres.
    // Preserve the scale that may have been set by the Scale keyword before this block.
    double saved_scale = refinement_parameters.empty() ? 1.0 : refinement_parameters[0];
    refinement_parameters.clear();
    refined_variable_names.clear();
    parameter_blocks.clear();

    // Scale is global parameter 0, but NOT part of structural parameter_blocks.
    refinement_parameters.push_back(saved_scale);
    refined_variable_names.push_back("Scale");

    set<string> seen_names;
    for (auto& block : blocks) {
      vector<double> b_vals;
      for (auto& p : block) {
        string name = boost::fusion::get<0>(p);
        double val  = boost::fusion::get<1>(p);
        if (!seen_names.insert(name).second) {
          REPORT(ERROR) << "Refinable parameter '" << name << "' is defined more than once in RefinableVariables.\n";
          throw(TerminateProgram());
        }
        refined_variable_names.push_back(name);
        refinement_parameters.push_back(val);
        b_vals.push_back(val);
      }
      if (!b_vals.empty()) parameter_blocks.push_back(b_vals);
    }

    formula.initialize_refinable_variables(refined_variable_names, refinement_parameters);
    expr_formula.initialize_refinable_variables(refined_variable_names, refinement_parameters);
  }

  vector<vector<double>> parameter_blocks;
  set<int> active_blocks; // 1-based indices of blocks to refine; empty = all active
  
  void initialize_unit_cell(vector<double> params)  {
    cell_is_initialized = true;
    cell = UnitCell(params[0],params[1],params[2],params[3],params[4],params[5]);
  }
  
  void set_fft_grid_size(vector<int> params)  {
    fft_grid_size=vec3<int>(params[0],params[1],params[2]);
  }

  /// Extra border pixels added to each PDF peak during FFT accumulation
  /// (see add_pair_to_appropriate_place). Input keyword: FFTBorderPixels.
  void set_fft_border_pixels(vector<int> params)  {
    fft_border_pixels=vec3<int>(params[0],params[1],params[2]);
  }
  
  void set_dump_pairs(bool inp)  {
    dump_pairs=inp;
  }
  
  void set_refinement_flag(bool inp)  {
    refinement_flag=inp;
  }
  
  ///\TODO: rename point group and laue_symmetry to PDF_symmetry consistently thwougout the code
  void set_point_group(string point_goup_symbol)  {
    cell.laue_symmetry = LaueSymmetry(point_goup_symbol,intensity_map.grid);
  }
  
  void set_calculation_method(bool method)  {
    direct_diffuse_scattering_calculation=method;
  }

  void set_use_mixed_derivatives(bool inp) {
    use_mixed_derivatives = inp;
  }

  void set_covariance_batch_size(int inp) {
    covariance_batch_size = inp;
  }

  void set_max_processors(int inp) {
    max_processors = inp;
  }

  void set_scattering_type(ScatteringType t) {
      scattering_type = t;
  }

  void set_recalculation_flag(bool flag)    {
    recalculate_average=flag;
  }
  
  /// Initialize a the grid from ... DiffuseScatteringGrid -6 -6 -6 1 1 1 12 12 12 (lower limits, step sizes,number of pixels)
  void initialize_intensity_grid(vector<double> params)  {
    if(!grid_initialized)
    {
      vec3<double> steps(params[3],params[4],params[5]);
      vec3<double> lower_limits(params[0],params[1],params[2]);
      vec3<int> grid_size(params[6],params[7],params[8]);
      
      grid =Grid(cell.cell,steps,lower_limits,grid_size,true);
      
      intensity_map = IntensityMap(grid);   
      average_intensity_map = IntensityMap(grid);
      data_ = IntensityMap(grid);   
      
      grid_initialized = true;
    }
  }
  
  void add_variant(ChemicalUnitNode* var)  {
    cell.add_node(var);
  }
  static vector<SubstitutionalCorrelation*> correlators_from_cuns_(StructurePartRef _var1,StructurePartRef _var2,vector<yell::ExprPtr> params)  {
    ChemicalUnitNode* var1= boost::get<ChemicalUnitNode*>(_var1);
    ChemicalUnitNode* var2= boost::get<ChemicalUnitNode*>(_var2);
    return correlators_from_cuns(var1,var2,params);
  }
  
  void init_flags()  {
    calculate_jacobians=false;
    direct_diffuse_scattering_calculation = true;
    use_mixed_derivatives = false;
    covariance_batch_size = 0;
    max_processors = 1;
    print_covariance_matrix=false;

    cell_is_initialized=false;
    average_is_calculated = false;
    grid_initialized = false;
    refinement_flag = true;
    recalculate_average=true;
    fft_grid_size=vec3<int>(16,16,16);
    fft_border_pixels=vec3<int>(0,0,0);
    dump_pairs=false;
    scattering_type = XRay;
    refinement_parameters = vector<double>(1,1);
    refined_variable_names = vector<string>(1,"Scale");
    periodic_boundaries = vector<bool>(3,true);
    refinement_options = RefinementOptions::default_refinement_options();
    report_pairs_outside_pdf_grid = false;
    padding = vec3<int>(0,0,0);
    refine_in_asu_val = true;
    model_parsed_ = false;
    derivatives_mode = FINITE_DIFFERENCE;
    jacobian_multiplier = 1.0;
  }
  void parse_model_();

  ///\TODO: test the following part of Model
  Model(string _model) : model(_model)
  {
    init_flags();
    parse_model_();
  }
  ///used for tests
  Model()   { 
    init_flags();
  }
  
  void add_correlations(const vector<AtomicPairPool*>& _pools)  {
    pools = _pools;
  }
  
  ADPMode* create_translational_mode(StructurePartRef chem_unit,int direction)  {
    return translational_mode(boost::get<ChemicalUnit*>(chem_unit),direction,cell.cell.metrical_matrix()); 
  }
  ADPMode* create_rotational_mode(StructurePartRef chem_unit,vector<double> params)  {
    return rot_mode(boost::get<ChemicalUnit*>(chem_unit),vec3<double>(params[0],params[1],params[2]),
                    vec3<double>(params[3],params[4],params[5]),cell.cell.metrical_matrix());
  }
  
  void add_modes(vector<ADPMode*> _modes) {
    modes.concat(_modes);
  }
  
  static DoubleADPMode* create_double_adp_mode(StructurePartRef mode1,StructurePartRef mode2,yell::ExprPtr amplitude)
  {
    return new DoubleADPMode(boost::get<ADPMode*>(mode1),boost::get<ADPMode*>(mode2),amplitude);
  }

  static SizeEffect* create_size_effect(StructurePartRef el1,StructurePartRef el2, yell::ExprPtr amplitude)
  {
    if(boost::get<ChemicalUnit*>(&el1)) //check that first element is ChemicalUnit
      return new SizeEffect(boost::get<ChemicalUnit*>(el1),boost::get<ADPMode*>(el2), amplitude);
    else
      return new SizeEffect(boost::get<ADPMode*>(el1),boost::get<ChemicalUnit*>(el2), amplitude);
  }
  
  static const bool AVERAGE = true;
  static const bool FULL = false;
  
  void apply_resolution_function_if_possible(IntensityMap& map)  {
    if(pdf_multiplier.is_loaded)
    {
      map.to_real();
      
      for(int i=0; i<map.size_1d(); ++i)
        map.at(i) *= pdf_multiplier.at(i);
    }
  }
  
  void apply_reciprocal_space_multipliers_if_possible(IntensityMap& map)  {
    if(reciprocal_space_multiplier.is_loaded)
    {
      map.to_reciprocal();
      
      for(int i=0; i<map.size_1d(); ++i)
        map.at(i) *= reciprocal_space_multiplier.at(i);
    }
  }
  
  /// function for minimizer
  void calculate(vector<double> params)   {
    double scale = params[0];
    
    if(!average_is_calculated || recalculate_average)
    {
      calculate(params, AVERAGE);
      average_is_calculated = true;
    }
    
    calculate(params,FULL);
    
    for(int i=0; i<data_.size_1d(); ++i)
      data_.at(i) = scale*(intensity_map.at(i)-average_intensity_map.at(i));
  }
    
  IntensityMap& model_scaled_to_experiment() {
    return data_;
  }

  int number_of_observations() {
      if(refine_in_asu())
          return asu_indices().size();
      else
          return intensity_map.size_1d();
  }
  
  void calculate(vector<double> params,bool);

  /// Helper: calculate full and average intensity maps from pre-baked peaks.
  void calculate_from_peaks(const vector<PattersonPeak>& full_peaks,
                            const vector<PattersonPeak>& avg_peaks);

  /// Helper: calculate derivative maps from pre-baked peaks and susceptibilities.
  /// sl must have up-to-date gridded form factors for the current FFT grid.
  IntensityMap calculate_derivative_from_peaks(
      const vector<PattersonPeak>& full_peaks,
      const vector<PattersonPeak>& avg_peaks,
      const vector<PeakSusceptibility>& full_susc,
      const vector<PeakSusceptibility>& avg_susc,
      double scale,
      int num_threads,
      ScattererList& sl);

  /// Fast path: calculate derivative map when base peaks are already baked.
  /// sl must have up-to-date gridded form factors for the current FFT grid.
  IntensityMap calculate_derivative_from_susceptibilities(
      const vector<PattersonPeak>& full_peaks,
      const vector<PattersonPeak>& avg_peaks,
      const vector<AtomicPair>& pairs,
      const Eigen::VectorXd& q,
      int param_idx,
      double scale,
      int num_threads,
      ScattererList& sl);

  /// Calculate derivative map dI/dp_j for a single parameter index j.
  /// Corresponds to the index in refinement_parameters / yell::ParameterBlock.
  IntensityMap calculate_derivative(const vector<double>& params, int param_idx, int num_threads = 0);

  UnitCell cell;
  IntensityMap intensity_map, average_intensity_map;
  string model;
  vector<AtomicPairPool*> pools;
  vector<double> refinement_parameters;
  vector<string> refined_variable_names;
  p_vector<ADPMode> modes;
  // Non-owning flat list of all atoms (cell atoms + molecular scatterer atoms).
  // Populated by construct_atom* during parsing; update_caches() called on each in calculate().
  vector<Atom*> all_atoms_;
  // Per-clone private atom copies for molecular-scatterer atoms (not in cell).
  // Empty in the original model; populated and owned by each clone.
  vector<Atom*> mol_owned_atoms_;
  // Per-clone private MolecularScatterer copies with remapped constituent_atoms.
  // Empty in the original model; populated and owned by each clone.
  vector<MolecularScatterer*> mol_owned_scatterers_;
  bool model_parsed_;
  ScattererList scatterer_list_;
  // Per-thread ScattererList copies — each owns its own gridded_form_factors_.
  // Populated lazily in compute_full_covariance(); reused across iterations.
  // When form factors become parameterized, invalidate and repopulate here.
  vector<ScattererList> thread_scatterer_lists_;

  // Destructor: release per-clone private atoms and MolecularScatterer copies.
  // These vectors are populated only in clone(); the original model leaves them empty.
  ~Model() {
      for (auto* a  : mol_owned_atoms_)       delete a;
      for (auto* ms : mol_owned_scatterers_)  delete ms;
  }

  Model* clone() const {
      // Snapshot atom pointers from the original BEFORE copy-constructing, so we
      // can build the remap after p_vector has deep-copied the atom tree.
      auto orig_cell_atoms = collect_atoms_(const_cast<UnitCell&>(cell));

      Model* m = new Model(*this);
      // NOTE: Model(*this) already deep-copies intensity_map / average_intensity_map
      // via IntensityMap's explicit copy constructor.  The assignments below use the
      // compiler-generated IntensityMap::operator= which calls af::versa::operator=
      // (shallow/aliased) — they intentionally re-alias the buffers so clones always
      // start from the same base map.  calculate_derivative() never writes to these
      // maps, so the alias is safe.
      m->intensity_map = intensity_map;
      m->average_intensity_map = average_intensity_map;

      // Deep copy pools (modifiers clone themselves; pairs vector is value-copied).
      m->pools.clear();
      for (size_t i = 0; i < pools.size(); ++i)
          m->pools.push_back(new AtomicPairPool(*pools[i]));

      if (!m->all_atoms_.empty()) {
          // Step 1: remap atoms that live in cell.chemical_unit_nodes (deep-copied
          // by Model(*this) copy-constructor via p_vector).
          auto clone_cell_atoms = collect_atoms_(m->cell);
          std::unordered_map<Atom*, Atom*> remap;
          remap.reserve(orig_cell_atoms.size());
          for (size_t i = 0; i < orig_cell_atoms.size(); ++i)
              remap[orig_cell_atoms[i]] = clone_cell_atoms[i];

          // Step 2: atoms that live in MolecularScatterers (NOT in cell) are NOT
          // deep-copied by the copy-constructor — they remain shared with the global
          // AtomicTypeCollection.  Create per-clone private copies for them.
          for (auto& a_ptr : m->all_atoms_) {
              if (remap.find(a_ptr) == remap.end()) {
                  // Atom not found in cell traversal → molecular scatterer atom.
                  // Allocate a private copy so each clone writes to its own Atom.
                  Atom* private_copy = new Atom(*a_ptr);
                  remap[a_ptr] = private_copy;
                  m->mol_owned_atoms_.push_back(private_copy);
              }
          }

          // Apply full remap to all_atoms_.
          for (auto& a_ptr : m->all_atoms_) {
              auto it = remap.find(a_ptr);
              if (it != remap.end()) a_ptr = it->second;
          }

          // Step 3: fix scatterer_list_.  MolecularScatterer objects in the global
          // AtomicTypeCollection store raw Atom* in constituent_atoms; these are the
          // same pointers we just remapped.  For each affected MolecularScatterer,
          // create a private copy with remapped atoms and register it as an override
          // in the clone's ScattererList.  The original pointer stays in scatterers_
          // so index_of() still works for PattersonPeak type-index assignment.
          for (auto* s_ptr : m->scatterer_list_.scatterers_) {
              if (auto* ms = dynamic_cast<MolecularScatterer*>(s_ptr)) {
                  bool needs_remap = false;
                  for (auto* a : ms->constituent_atoms)
                      if (remap.count(a)) { needs_remap = true; break; }
                  if (needs_remap) {
                      auto* ms_clone = new MolecularScatterer(*ms);
                      for (auto& a : ms_clone->constituent_atoms)
                          if (auto it = remap.find(a); it != remap.end())
                              a = it->second;
                      m->scatterer_list_.add_override(ms, ms_clone);
                      m->mol_owned_scatterers_.push_back(ms_clone);
                  }
              }
          }
      }

      return m;
  }
  
  RefinementOptions refinement_options;

  // TODO: check this is good implementation, remove debug stuff
  Eigen::VectorXd refinement_parameters_asEig() const {
    if (refinement_parameters.empty())
      return {};
    else {
      return Eigen::VectorXd::Map(refinement_parameters.data(), refinement_parameters.size());
    }

  }

  bool refine_in_asu() {
      return refine_in_asu_val;
  }
  bool refine_in_asu_val;

  void init_asu() {
      if(refine_in_asu()) {
          asu_indices_val = cell.laue_symmetry.asymmetric_indices(grid);
      }
  }
    vector<int> & asu_indices() {
        return asu_indices_val;
    }
  vector<int> asu_indices_val;
    
  DerivativesMode derivatives_mode;
  double jacobian_multiplier;

  // Analytical Jacobian for the direct calculation method.
  // Returns matrix of shape (n_observations × n_params).
  // Residuals: r_i = (exp_i - data_i) * w_i, data_i = Scale*(Ifull_i - Iavg_i).
  Eigen::MatrixXd compute_analytical_jacobian_direct(
      const vector<double>& params,
      IntensityMap& exp_map,
      OptionalIntensityMap& wts);

  Eigen::MatrixXd compute_jacobian_mixed(
      const vector<double>& params,
      IntensityMap& exp_map,
      OptionalIntensityMap& wts);

  /// Analytically compute the optimal Scale = (Σ w²·Ie·Ic) / (Σ w²·Ic²).
  /// Requires intensity_map and average_intensity_map to be current (call calculate() first).
  /// Clamps to zero if negative; returns 1.0 if denominator is negligible.
  double compute_optimal_scale(IntensityMap& exp_map, OptionalIntensityMap& wts);

  /// Compute the full covariance matrix (N_params x N_params) including the Scale.
  /// Uses a memory-efficient streaming approach for the Hessian accumulation.
  Eigen::MatrixXd compute_full_covariance(
      const vector<double>& params,
      IntensityMap& exp_map,
      OptionalIntensityMap& wts);

  vec3<int> fft_grid_size;
  vec3<int> fft_border_pixels;
  vector<bool> periodic_boundaries;
//  int number_of_parameters;
  bool cell_is_initialized;
  bool grid_initialized;
  bool average_is_calculated;
  bool recalculate_average;
  bool direct_diffuse_scattering_calculation;
  bool use_mixed_derivatives;
  bool refinement_flag;
  bool dump_pairs;
  bool report_pairs_outside_pdf_grid;
  bool print_covariance_matrix;
  bool calculate_jacobians;
  int  covariance_batch_size;
  int  max_processors;
  ScatteringType scattering_type;
  IntensityMap& data() { return data_; }
  IntensityMap& get_intensity_map() { return intensity_map; }
  IntensityMap& get_average_intensity_map() { return average_intensity_map; }

  IntensityMap data_;
  //TODO: For multithreading I will need to split the class "model" into two parts. All the large arrays like weights, multipliers etc will have to go to another object,
  //this one has to be thread-safe, meaning that calculation operation should be possible to do in an instantiated object
  OptionalIntensityMap weights;
  OptionalIntensityMap reciprocal_space_multiplier;
  OptionalIntensityMap pdf_multiplier;
  vector<AtomicPair> atomic_pairs;
  vec3<int> padding;
  Grid grid;
};

#endif