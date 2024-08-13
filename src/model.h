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
#include <boost/fusion/tuple.hpp>


//COPYPASTE from InputFileParser.h
typedef boost::variant<ChemicalUnit*,ChemicalUnitNode*,ADPMode*> StructurePartRef;
//#include "InputFileParser.h"

enum R_FACTORS {R1, R2};
enum WEIGHTED_OPTIONS {WEIGHTED, UNWEIGHTED};


class Model : public MinimizerCalculator {
private:
  static double sq(double x) {
    return x*x;
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
  
  /**
   * Creates an atom with isotropic ADP
   */
  Atom* construct_atom_isotropic_adp(string name,vector<double> params)  {
    return new Atom(name,params[0],
                    1,//we assign probability to 1 because it will be changed by Variant afterwards anyway
                    params[1],params[2],params[3],params[4],cell.cell.reciprocal_metrical_matrix(),scattering_type);
  }
  
  /**
   * This function creates an atom from atom name and a vector of atomic parameters. The atom should be deleted elsewhere (normally its pointer should be given to AtomicAssembly which will take care of it)
   */
  Atom* construct_atom(string name,vector<double> params)  {
    
    double a=cell.cell.parameters()[0];
    double b=cell.cell.parameters()[1];
    double c=cell.cell.parameters()[2];
    
    //Also transforms input adps from angstroems^2 into fractional values
    return new Atom(name,params[0],
                    1, //we assign probability to 1 because it will be changed by Variant afterwards anyway
                    params[1],params[2],params[3],
                    params[4]/a/a,
                    params[5]/b/b,
                    params[6]/c/c,
                    params[7]/a/b,
                    params[8]/a/c,
                    params[9]/b/c,
                    scattering_type);
  }
  
  void set_scale(double inp) {  
    refinement_parameters[0]=inp;
  }
  
  void set_refinable_parameters(FormulaParser& formula,vector<boost::fusion::tuple<string,double> > inp) {  
    refinement_parameters.resize(inp.size()+1);
    refined_variable_names.resize(inp.size()+1);
    
    for(int i=0; i<inp.size(); ++i)
    {
      
      refined_variable_names[i+1]=boost::fusion::get<0>(inp[i]);
      refinement_parameters[i+1]=boost::fusion::get<1>(inp[i]);
    }
    
    // If this is a first parser run
    formula.initialize_refinable_variables(refined_variable_names,refinement_parameters);
  }
  
  void initialize_unit_cell(vector<double> params)  {
    cell_is_initialized = true;
    cell = UnitCell(params[0],params[1],params[2],params[3],params[4],params[5]);
  }
  
  void set_fft_grid_size(vector<int> params)  {
    fft_grid_size=vec3<int>(params[0],params[1],params[2]);
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
  static vector<SubstitutionalCorrelation*> correlators_from_cuns_(StructurePartRef _var1,StructurePartRef _var2,vector<double> params)  {
    ChemicalUnitNode* var1= boost::get<ChemicalUnitNode*>(_var1);
    ChemicalUnitNode* var2= boost::get<ChemicalUnitNode*>(_var2);
    return correlators_from_cuns(var1,var2,params);
  }
  
  void init_flags()  {
    calculate_jacobians=false;
    direct_diffuse_scattering_calculation = true;
    print_covariance_matrix=false;
    cell_is_initialized=false;
    average_is_calculated = false;
    grid_initialized = false;
    refinement_flag = true;
    recalculate_average=true;
    fft_grid_size=vec3<int>(16,16,16);
    dump_pairs=false;
    scattering_type = XRay;
    refinement_parameters = vector<double>(1,1);
    refined_variable_names = vector<string>(1,"Scale");
    periodic_boundaries = vector<bool>(3,true);
    refinement_options = RefinementOptions::default_refinement_options();
    report_pairs_outside_pdf_grid = false;
    padding = vec3<int>(0,0,0);
    refine_in_asu_val = true;
  }
  ///\TODO: test the following part of Model
  Model(string _model) : model(_model)
  {
    init_flags();
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
  
  static DoubleADPMode* create_double_adp_mode(StructurePartRef mode1,StructurePartRef mode2,double amplitude)
  {
    return new DoubleADPMode(boost::get<ADPMode*> (mode1),boost::get<ADPMode*>(mode2),amplitude);
  }
  
  static SizeEffect* create_size_effect(StructurePartRef el1,StructurePartRef el2, double amplitude)
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

  UnitCell cell;
  IntensityMap intensity_map, average_intensity_map;
  string model;
  vector<AtomicPairPool*> pools;
  vector<double> refinement_parameters;
  vector<string> refined_variable_names;
  p_vector<ADPMode> modes;
  
  RefinementOptions refinement_options;

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
    
  vec3<int> fft_grid_size;
  vector<bool> periodic_boundaries;
//  int number_of_parameters;
  bool cell_is_initialized;
  bool grid_initialized;
  bool average_is_calculated;
  bool recalculate_average;
  bool direct_diffuse_scattering_calculation;
  bool refinement_flag;
  bool dump_pairs;
  bool report_pairs_outside_pdf_grid;
  bool print_covariance_matrix;
  bool calculate_jacobians;
    ScatteringType scattering_type;
  IntensityMap& data() { return data_; }
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