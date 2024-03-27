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

#include "InputFileParser.h"

// templated funcitons somehow do not get binded with phoenix
void report_after_refinement_dbl(double inp) {
  REPORT(AFTER_REFINEMENT) << inp;
}
void report_after_refinement(string inp) {
  REPORT(AFTER_REFINEMENT) << inp;
}

InputParser::InputParser() : InputParser::base_type(start)
{
  using namespace qi;
  using phoenix::bind;
  using phoenix::ref;
  
  //check unit cell is initialized
  
  start = 
      program_options
    > unit_cell [_pass = bind(&Model::unit_cell_is_initialized,*ref(model))]
    > -modes
    > correlations                     [bind(&Model::add_correlations,*ref(model),_1)]
    >> maybe_assignments
    >> *(print_command_parser >> maybe_assignments)
    ;
  
  program_options =
    maybe_assignments
    >> *(program_option >> maybe_assignments)
    ;
  
  modes =
    lit("Modes")
    > "[" >> maybe_assignments
    > (*(mode_assignement >> maybe_assignments))[bind(&Model::add_modes,*ref(model),_1)]
    > "]"
    ;
  
  unit_cell =
    lit("UnitCell")
    > "[" >> maybe_assignments
    > *( (variant | variant_assignement) >> maybe_assignments )[bind(&Model::add_variant,*ref(model),_1)]
    > "]"
  ;
  
  correlations %=
    lit("Correlations")
    > "[" >> maybe_assignments
    > *(atomic_pair_pool >> maybe_assignments)
    > "]"
    ;
  
  maybe_assignments = *(formula.assignment >> ';');
  
  program_option =
      (lit("Cell") > repeat(6)[number])                 [bind(&Model::initialize_unit_cell,*ref(model),_1)]
    | (lit("DiffuseScatteringGrid") > repeat(9)[number])[bind(&Model::initialize_intensity_grid,*ref(model),_1)]
    | (lit("MaxNumberOfIterations") > int_)             [bind(&Model::set_max_number_of_iterations,*ref(model),_1)]
    | (lit("MinimizerTau") > double_)                   [bind(&Model::set_tau,*ref(model),_1)]
    | (lit("MinimizerThresholds") > repeat(3)[double_]) [bind(&Model::set_thresholds,*ref(model),_1)]
    | (lit("MinimizerDiff") > double_)                  [bind(&Model::set_diff,*ref(model),_1)]
    | (lit("FFTGridSize") > repeat(3)[int_])            [bind(&Model::set_fft_grid_size,*ref(model),_1)]
    | (lit("FFTGridPadding") > repeat(3)[int_])         [bind(&Model::set_padding,*ref(model),_1)]
    | (lit("DumpPairs") > bool_)                        [bind(&Model::set_dump_pairs,*ref(model),_1)]
    | (lit("CalculateJacobians") > bool_)               [bind(&Model::set_calculate_jacobians,*ref(model),_1)]
    | program_option1
  ;
  
  program_option1 =
      (lit("Refine") > bool_)                                   [bind(&Model::set_refinement_flag,*ref(model),_1)]
    | (lit("ReportPairsOutsideCalculatedPDF") > bool_)          [bind(&Model::set_report_pairs_outside_pdf_grid,*ref(model),_1)]
    | (lit("PeriodicBoundaries") > repeat(3) [bool_])           [bind(&Model::set_periodic_boundaries,*ref(model),_1)]
    | (lit("RecalculateAverage") > bool_)                       [bind(&Model::set_recalculation_flag,*ref(model),_1)]
    | (lit("CalculationMethod") > calculation_methods)          [bind(&Model::set_calculation_method,*ref(model),_1)]
    | (lit("LaueSymmetry")
       > lexeme[point_group_symbol >> !char_("-:/a-zA-Z0-9")])  [bind(&Model::set_point_group,*ref(model),_1)]
    | (lit("Scale") > double_ > -(omit[lit('(')>int_>lit(')')]))[bind(&Model::set_scale,*ref(model),_1)]
    | (lit("PrintCovarianceMatrix")> bool_)                     [bind(&Model::set_print_covariance_matrix,*ref(model),_1)]
    | refinable_parameters                                      [bind(&Model::set_refinable_parameters,*ref(model),ref(formula),_1)]
    | program_option2
  ;

  /* TODO: fixe possible problem coming from the fact that atoms can be defined in molecular scatterers before the
   * scattering type is defined in such a case some of the atoms will be default XRay, while other will be Neutrons
   * */
  program_option2 =
      molecular_scatterers
    | (lit("Scattering") > scattering_type)                     [bind(&Model::set_scattering_type,*ref(model),_1)]
  ;

  scattering_type.add
          ("x-ray", XRay)
          ("neutron", Neutron)
          ("electron", Electrons)
          ;

  refinable_parameters %=
    lit("RefinableVariables")               
    > "["                                  
    >> *(valid_identifier
         > "="
         > double_
         > -(lit('(') > int_ > lit(')'))
         > -lit(';')                         
         )                                  
    > "]"
    ;
  
  string_in_quotes %=
    '"'
    > *(char_-'"')
    > '"'
  ;
  
  print_command_parser = 
  lit("Print")[bind(&report_after_refinement,"Requested output:")]
    >> *(
         (formula.expr >> !char_(';'))[bind(&report_after_refinement_dbl,_1)]
         | string_in_quotes[bind(&report_after_refinement,_1)]
        )
    >> eps[bind(&report_after_refinement,"\n")]
  ;
  
  calculation_methods.add
    ("direct",true)
    ("exact reciprocal",true)
    ("exact",true)
    ("reciprocal",true)
  
    ("approximate",false)
    ("approximate pdf",false)
    ("pdf",false)
    ("fft",false)
  ;

  
  
  point_group_symbol.add
  ("m-3m",  "m-3m")
  ("m-3",   "m-3")
  ("6/mmm", "6/mmm")
  ("6/m",   "6/m")
  ("-3mH",  "-3mH")
  ("-3mR",  "-3mR")
  ("-3H",   "-3H")
  ("-3R",   "-3R")
  ("4/mmm", "4/mmm")
  ("4/m",   "4/m")
  ("mmm",   "mmm")
  ("2/m",   "2/m")
  ("2/mb",  "2/mb")
  ("-1",    "-1")
  // synonims
  ("-3:R","-3R")
  ("-3:H","-3H")
  ("-3m:H","-3mH")
  ("-3m:R","-3mR")
  ("2/m:b","2/mb")
  ("m3m","m-3m")
  ;

  

  atomic_pair_pool = lit("[")[_val = phoenix::new_<AtomicPairPool>()]
  > *(
      cell_shifter                          [bind(&AtomicPairPool::add_modifier,*_val,_1)]
      | multiplicity_correlation            [bind(&AtomicPairPool::add_modifier,*_val,_1)]
      | substitutional_correlation          [bind(&AtomicPairPool::add_modifiers,*_val,_1)]
      | adp_correlation                     [bind(&AtomicPairPool::add_modifier,*_val,_1)]
      | size_effect                         [bind(&AtomicPairPool::add_modifier,*_val,_1)]
      )
  > "]";
   
  
  translational_mode = 
    lit("TranslationalMode")
    > '(' 
    > ( identifier > ',' > basis_vector )[_val = bind(&Model::create_translational_mode,*ref(model),_1,_2)]
    > ')'
  ;
        
  rotational_mode = 
    lit("RotationalMode")
    > '('
    > ( identifier > repeat(6)[',' > number])[_val = bind(&Model::create_rotational_mode,*ref(model),_1,_2)]
    > ')'
  ;
  

  
  mode_assignement = 
    (
     valid_identifier 
     > "=" 
     > (rotational_mode | translational_mode)[_val=_1] 
    )[bind(InputParser::add_reference,phoenix::ref(references),_1,_2)]
  ;
  

  
  InputParserI();
  InputParserII();
  InputParserIII();
}
