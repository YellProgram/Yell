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

// Mirror a formula assignment (name, evaluated double) into ExprFormulaParser
// as a frozen literal.  Refinable variables have already been registered as
// ParamRef leaves via initialize_refinable_variables(); ordinary assignments
// (ll=-12.8, gs=-ll*2/104, etc.) are constants that do not change during
// refinement, so lit(value) is correct.
static void mirror_to_expr(ExprFormulaParser& efp,
                            FormulaParser::NamedValue nv)
{
    efp.add_assignment(nv.first, yell::lit(nv.second));
}

InputParser::InputParser() : InputParser::base_type(start)
{
  using namespace qi;
  using phoenix::ref;
  
  //check unit cell is initialized
  
  start = 
      program_options
    > unit_cell [_pass = phoenix::bind(&Model::unit_cell_is_initialized,*ref(model))]
    > -modes
    > correlations                     [phoenix::bind(&Model::add_correlations,*ref(model),_1)]
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
    > (*(mode_assignement >> maybe_assignments))[phoenix::bind(&Model::add_modes,*ref(model),_1)]
    > "]"
    ;
  
  unit_cell =
    lit("UnitCell")
    > "[" >> maybe_assignments
    > *( (variant | variant_assignement) >> maybe_assignments )[phoenix::bind(&Model::add_variant,*ref(model),_1)]
    > "]"
  ;
  
  correlations %=
    lit("Correlations")
    > "[" >> maybe_assignments
    > *(atomic_pair_pool >> maybe_assignments)
    > "]"
    ;
  
  maybe_assignments = *(formula.named_assignment[
      phoenix::bind(&mirror_to_expr, phoenix::ref(expr_formula), _1)] >> ';');
  
  program_option =
      (lit("Cell") > repeat(6)[number])                 [phoenix::bind(&Model::initialize_unit_cell,*ref(model),_1)]
    | (lit("DiffuseScatteringGrid") > repeat(9)[number])[phoenix::bind(&Model::initialize_intensity_grid,*ref(model),_1)]
    | (lit("MaxNumberOfIterations") > int_)             [phoenix::bind(&Model::set_max_number_of_iterations,*ref(model),_1)]
    | (lit("MinimizerTau") > double_)                   [phoenix::bind(&Model::set_tau,*ref(model),_1)]
    | (lit("MinimizerThresholds") > repeat(3)[double_]) [phoenix::bind(&Model::set_thresholds,*ref(model),_1)]
    | (lit("MinimizerDiff") > double_)                  [phoenix::bind(&Model::set_diff,*ref(model),_1)]
    | (lit("FFTGridSize") > repeat(3)[int_])            [phoenix::bind(&Model::set_fft_grid_size,*ref(model),_1)]
    | (lit("FFTGridPadding") > repeat(3)[int_])         [phoenix::bind(&Model::set_padding,*ref(model),_1)]
    | (lit("DumpPairs") > bool_)                        [phoenix::bind(&Model::set_dump_pairs,*ref(model),_1)]
    | (lit("CalculateJacobians") > bool_)               [phoenix::bind(&Model::set_calculate_jacobians,*ref(model),_1)]
    | program_option1
  ;
  
  program_option1 =
      (lit("Refine") > bool_)                                   [phoenix::bind(&Model::set_refinement_flag,*ref(model),_1)]
    | (lit("ReportPairsOutsideCalculatedPDF") > bool_)          [phoenix::bind(&Model::set_report_pairs_outside_pdf_grid,*ref(model),_1)]
    | (lit("PeriodicBoundaries") > repeat(3) [bool_])           [phoenix::bind(&Model::set_periodic_boundaries,*ref(model),_1)]
    | (lit("RecalculateAverage") > bool_)                       [phoenix::bind(&Model::set_recalculation_flag,*ref(model),_1)]
    | (lit("CalculationMethod") > calculation_methods)          [phoenix::bind(&Model::set_calculation_method,*ref(model),_1)]
    | (lit("LaueSymmetry")
       > lexeme[point_group_symbol >> !char_("-:/a-zA-Z0-9")])  [phoenix::bind(&Model::set_point_group,*ref(model),_1)]
    | (lit("Scale") > double_ > -(omit[lit('(')>int_>lit(')')]))[phoenix::bind(&Model::set_scale,*ref(model),_1)]
    | (lit("PrintCovarianceMatrix")> bool_)                     [phoenix::bind(&Model::set_print_covariance_matrix,*ref(model),_1)]
    | (lit("CovarianceBatchSize") > int_)                       [phoenix::bind(&Model::set_covariance_batch_size,*ref(model),_1)]
    | refinable_parameters                                      [phoenix::bind(&Model::set_refinable_parameter_blocks,*ref(model),ref(formula),ref(expr_formula),_1)]
    | program_option2
  ;

  /* TODO: fixe possible problem coming from the fact that atoms can be defined in molecular scatterers before the
   * scattering type is defined in such a case some of the atoms will be default XRay, while other will be Neutrons
   * */
  program_option2 =
      molecular_scatterers
    | (lit("Scattering") > scattering_type)                     [phoenix::bind(&Model::set_scattering_type,*ref(model),_1)]
    | (lit("Derivatives") > derivatives_mode_sym)               [phoenix::bind(&Model::set_derivatives_mode,*ref(model),_1)]
    | (lit("JacobianMultiplier") > double_)                     [phoenix::bind(&Model::set_jacobian_multiplier,*ref(model),_1)]
  ;

  derivatives_mode_sym.add
    ("finite_difference", FINITE_DIFFERENCE)
    ("analytical",        ANALYTICAL)
    ("mixed",             MIXED)
    ;

  scattering_type.add
          ("x-ray", XRay)
          ("neutron", Neutron)
          ("electron", Electrons)
          ;

  single_parameter %= 
    valid_identifier
    > "="
    > double_
    > -(lit('(') > int_ > lit(')'))
    > -lit(';')
    ;

  parameter_block %=
    "[" >> *single_parameter >> "]"
    | eps[_val = phoenix::construct<vector<boost::fusion::tuple<std::string,double>>>()] >> +single_parameter
    ;

  refinable_parameters %=
    lit("RefinableVariables")               
    > "["                                  
    >> *parameter_block
    > "]"
    ;
  
  string_in_quotes %=
    '"'
    > *(char_-'"')
    > '"'
  ;
  
  print_command_parser = 
  lit("Print")[phoenix::bind(&report_after_refinement,"Requested output:")]
    >> *(
         (formula.expr >> !char_(';'))[phoenix::bind(&report_after_refinement_dbl,_1)]
         | string_in_quotes[phoenix::bind(&report_after_refinement,_1)]
        )
    >> eps[phoenix::bind(&report_after_refinement,"\n")]
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
      cell_shifter                          [phoenix::bind(&AtomicPairPool::add_modifier,*_val,_1)]
      | multiplicity_correlation            [phoenix::bind(&AtomicPairPool::add_modifier,*_val,_1)]
      | substitutional_correlation          [phoenix::bind(&AtomicPairPool::add_modifiers,*_val,_1)]
      | adp_correlation                     [phoenix::bind(&AtomicPairPool::add_modifier,*_val,_1)]
      | size_effect                         [phoenix::bind(&AtomicPairPool::add_modifier,*_val,_1)]
      )
  > "]";
   
  
  translational_mode = 
    lit("TranslationalMode")
    > '(' 
    > ( identifier > ',' > basis_vector )[_val = phoenix::bind(&Model::create_translational_mode,*ref(model),_1,_2)]
    > ')'
  ;
        
  rotational_mode = 
    lit("RotationalMode")
    > '('
    > ( identifier > repeat(6)[',' > number])[_val = phoenix::bind(&Model::create_rotational_mode,*ref(model),_1,_2)]
    > ')'
  ;
  

  
  mode_assignement = 
    (
     valid_identifier 
     > "=" 
     > (rotational_mode | translational_mode)[_val=_1] 
    )[phoenix::bind(&InputParser::add_reference,phoenix::ref(references),_1,_2)]
  ;
  

  
  InputParserI();
  InputParserII();
  InputParserIII();
}

// ─── InputParserI ─────────────────────────────────────────────────────────────

// Wrapper helpers replacing phoenix::try_/catch_all (removed in Boost 1.88).
// Spirit Qi exposes sequence elements via _1/_2/_3 inside semantic actions;
// _val and _pass are the rule's synthesized attribute and pass flag.
static void do_correlators(
    vector<SubstitutionalCorrelation*>& out,
    StructurePartRef a, StructurePartRef b,
    vector<yell::ExprPtr> c, bool& pass)
{
    try { out = Model::correlators_from_cuns_(a, b, c); }
    catch (...) { pass = false; }
}

static void do_double_adp_mode(
    DoubleADPMode*& out,
    StructurePartRef a, StructurePartRef b,
    double c, bool& pass)
{
    try { out = Model::create_double_adp_mode(a, b, c); }
    catch (...) { pass = false; }
}

static void do_size_effect(
    SizeEffect*& out,
    StructurePartRef a, StructurePartRef b,
    double c, bool& pass)
{
    try { out = Model::create_size_effect(a, b, c); }
    catch (...) { pass = false; }
}

void InputParser::InputParserI()
{

  using namespace qi;
  using phoenix::ref;

  cell_shifter = '(' >
    ( number > ',' > number > ',' > number)[_val =  phoenix::new_<CellShifter>(_1,_2,_3)]
    > ')'
    ;
  multiplicity_correlation = lit("Multiplicity") > number[_val = phoenix::new_<MultiplicityCorrelation>(_1)];

  expr_number = lexeme[expr_formula];

  substitutional_correlation =
    lit("SubstitutionalCorrelation")
    > "("
    > (identifier > "," > identifier > "," > expr_number % ',')[
        phoenix::bind(&do_correlators, _val, _1, _2, _3, _pass)]
    > ")"
    ;

    adp_correlation =
    lit("ADPCorrelation")
    > '('
    > (identifier > ',' > identifier > ',' > number)[
        phoenix::bind(&do_double_adp_mode, _val, _1, _2, _3, _pass)]
    > ')'
    ;

    size_effect =
    lit("SizeEffect")
    > '('
    > (identifier > ',' > identifier > ',' > number)[
        phoenix::bind(&do_size_effect, _val, _1, _2, _3, _pass)]
    > ')'
    ;

  number %= formula | double_;

  skipper_no_assignement = boost::spirit::ascii::space | comment;
  skipper = omit[(formula.named_assignment[
      phoenix::bind(&mirror_to_expr, phoenix::ref(expr_formula), _1)] >> ';')] | skipper_no_assignement;

  comment = lit("#") >> *(char_ - eol) >> eol;

  basis_vector.add("x",0)("y",1)("z",2);

  permutation_component =
  eps[_val = phoenix::construct<vector<double> >(4,0)]
  >> -basis_vector[_val[_1]=1]
  >> *(
         ('+' >> basis_vector)[_val[_1]=1]
       | ('-' >> basis_vector)[_val[_1]=-1]
       | (-lit('+') >> double_ >> lit("/") >> int_)[_val[3]=_1/_2]
       | ('-' >> double_ >> lit("/") >> int_)[_val[3]=-_1/_2]
       | (-lit('+') >> double_)[_val[3]=_1]
       | ('-' >> double_)[_val[3]=-_1]
       )
  ;
}

// ─── InputParserII ────────────────────────────────────────────────────────────

void InputParser::InputParserII()
{

  using namespace qi;
  using phoenix::ref;

  symmetry_element =
  lit("Symmetry")
  > '('
  > permutation_component[_val = phoenix::bind(&update_symmetry_component,_val,_1,0)]
  > ',' > permutation_component[_val = phoenix::bind(&update_symmetry_component,_val,_1,1)]
  > ',' > permutation_component[_val = phoenix::bind(&update_symmetry_component,_val,_1,2)]
  > ')'
  ;

  variant = (
              lit("Variant") [_val = phoenix::new_<ChemicalUnitNode>()]
             > '['
             > *(
                  '(' > omit[ char_("pP") ] > '=' >  number  > ')' >
                  chemical_unit[phoenix::bind(&ChemicalUnitNode::add_chemical_unit,*_val,_1)]
                  ) [phoenix::bind(&ChemicalUnit::set_occupancy,*_2,_1)]
             > ']'
             > eps [_pass = phoenix::bind(&ChemicalUnitNode::complain_if_sum_of_occupancies_is_not_one,_val)]
            );

  variant_assignement = (valid_identifier >> '=' >> variant[_val=_1])[phoenix::bind(&InputParser::add_reference,phoenix::ref(references),_1,phoenix::construct<StructurePartRef>(_val))];

  valid_identifier %= alpha >> *char_("a-zA-Z_0-9") >> !char_("a-zA-Z_0-9");

}

// ─── InputParserIII ───────────────────────────────────────────────────────────

void add_symbol( qi::symbols<char,std::string>& molecular_type, vector<boost::tuple<string,Scatterer*> > inp) {
  for(int i=0; i<inp.size(); ++i) {
    std::string label = boost::get<0>(inp[i]);
    molecular_type.add(label.c_str(),label);
  }
}

void InputParser::InputParserIII()
{
  using namespace qi;
  using phoenix::ref;

  chemical_unit %=
   lit("Void")[_val = phoenix::new_<AtomicAssembly>()]
   | chemical_unit_assignement
   | atomic_assembly
   | symmetric_chemical_unit
   | atom
  ;
  symmetric_chemical_unit = (identifier >> '*' >> symmetry_element)[_val = phoenix::bind(&create_symmetric_chemical_unit,_1,_2)];

  chemical_unit_assignement =
    valid_identifier[_a=_1]
    >> "="
    >> chemical_unit[_val = _1]
    >> eps[phoenix::bind(&InputParser::add_reference,phoenix::ref(references),_a,phoenix::construct<StructurePartRef>(_val))]
    ;

  atomic_assembly = ('[' >> *chemical_unit >> ']')[_val = phoenix::new_<AtomicAssembly>(_1)] ;

  atom =
    (atom_name >> repeat(10)[lexeme[expr_formula]]) [_val = phoenix::bind(&Model::construct_atom,*ref(model),_1,_2)] //Uaniso and p
  | (atom_name > repeat(5)[lexeme[expr_formula]])   [_val = phoenix::bind(&Model::construct_atom_isotropic_adp,*ref(model),_1,_2)]
  ;

  molecular_scatterer =
    chemical_unit[_val = phoenix::new_<MolecularScatterer>(_1)]
    ;

  molecular_scatterers %=
    lit("MolecularScatterers")
    > "["
    >> *(valid_identifier
         >> "="
         >> molecular_scatterer)
    > "]"
    > eps[phoenix::bind(&Model::register_molecular_scatterers,_val)]
    > eps[phoenix::bind(&add_symbol,phoenix::ref(molecular_type),_val)]
    ;

  atom_type.add("H","H")("He","He")("Li","Li")("Be","Be")("B","B")("C","C")("N","N")("O","O")("F","F")("Ne","Ne")("Na","Na")("Mg","Mg")("Al","Al")("Si","Si")("P","P")("S","S")("Cl","Cl")("Ar","Ar")("K","K")("Ca","Ca")("Sc","Sc")("Ti","Ti")("V","V")("Cr","Cr")("Mn","Mn")("Fe","Fe")("Co","Co")("Ni","Ni")("Cu","Cu")("Zn","Zn")("Ga","Ga")("Ge","Ge")("As","As")("Se","Se")("Br","Br")("Kr","Kr")("Rb","Rb")("Sr","Sr")("Y","Y")("Zr","Zr")("Nb","Nb")("Mo","Mo")("Tc","Tc")("Ru","Ru")("Rh","Rh")("Pd","Pd")("Ag","Ag")("Cd","Cd")("In","In")("Sn","Sn")("Sb","Sb")("Te","Te")("I","I")("Xe","Xe")("Cs","Cs")("Ba","Ba")("La","La")("Ce","Ce")("Pr","Pr")("Nd","Nd")("Pm","Pm")("Sm","Sm")("Eu","Eu")("Gd","Gd")("Tb","Tb")("Dy","Dy")("Ho","Ho")("Er","Er")("Tm","Tm")("Yb","Yb")("Lu","Lu")("Hf","Hf")("Ta","Ta")("W","W")("Re","Re")("Os","Os")("Ir","Ir")("Pt","Pt")("Au","Au")("Hg","Hg")("Tl","Tl")("Pb","Pb")("Bi","Bi")("Po","Po")("At","At")("Rn","Rn")("Fr","Fr")("Ra","Ra")("Ac","Ac")("Th","Th")("Pa","Pa")("U","U")("Np","Np")("Pu","Pu")("Am","Am")("Cm","Cm")("Bk","Bk")("Cf","Cf")("Es","Es")("Fm","Fm")("Md","Md")("No","No")("Lr","Lr")("Rf","Rf")("Db","Db")("Sg","Sg")("Bh","Bh")("Hs","Hs")("Mt","Mt")("Ds","Ds")("Rg","Rg")("Cn","Cn")("Uut","Uut")("Fl","Fl")("Uup","Uup")("Lv","Lv")("Uus","Uus")("Uuo","Uuo");

  atom_name %=
    molecular_type
    | atom_type >> *alnum
    ;

  identifier %= references >> !alnum;
  rest_of_the_line = *(qi::char_ - qi::eol) >> qi::eol;
}
