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

#ifndef InputFileParser_H
#define InputFileParser_H

#include "basic_classes.h"
#include "precompiled_header.h"
#include "FormulaParser.h"
#include "model.h"
#include <boost/fusion/tuple.hpp>
#include <boost/fusion/include/boost_tuple.hpp>

typedef boost::variant<ChemicalUnit*,ChemicalUnitNode*,ADPMode*> StructurePartRef;

typedef iterator_ Iterator;
struct InputParser : qi::grammar<Iterator, void(), qi::rule<Iterator,void()> >{
 //public:
  InputParser();
  void InputParserI();
  void InputParserII();
  void InputParserIII();
  //
  //typedef ChemicalUnit* StructurePartRef;
  typedef qi::rule<Iterator,void()> skipper_type;
    
  qi::rule<Iterator, void(), skipper_type> start;
  qi::rule<Iterator, Atom*(), skipper_type> atom;
  qi::rule<Iterator, Atom*(), skipper_type> isotropic_atom;
  qi::rule<Iterator,string()> atom_name;
  qi::rule<Iterator,ChemicalUnit*(), skipper_type > chemical_unit;
  qi::rule<Iterator,qi::locals<string>,ChemicalUnit*(), skipper_type > chemical_unit_assignement;
  qi::rule<Iterator,AtomicAssembly*(), skipper_type> atomic_assembly;
  qi::rule<Iterator,string()> valid_identifier;
  typedef qi::symbols<char,StructurePartRef> ReferenceTable;                                                                                         
  ReferenceTable references;                                       
  qi::rule<Iterator,StructurePartRef()> identifier;
  qi::rule<Iterator,ChemicalUnitNode*(), skipper_type > variant;
  qi::rule<Iterator,ChemicalUnitNode*(), skipper_type > variant_assignement;
  qi::rule<Iterator,SymmetryElement(),skipper_type> symmetry_element; 
  
  qi::symbols<char,int> basis_vector;

  qi::rule<Iterator,vector<double>(),skipper_type> permutation_component;
  qi::rule<Iterator,ChemicalUnit*(), skipper_type> symmetric_chemical_unit;
  qi::rule<Iterator,double()> number;
  FormulaParser formula;
  qi::rule<Iterator,void()> rest_of_the_line;
  qi::rule<Iterator,void(),skipper_type> unit_cell; //this parser does not return anything because it adds the variants directly to model->cell
  qi::rule<Iterator,vector<SubstitutionalCorrelation*>(), skipper_type> substitutional_correlation;
  qi::rule<Iterator,MultiplicityCorrelation*(),skipper_type> multiplicity_correlation;
  qi::rule<Iterator,CellShifter*(),skipper_type> cell_shifter;
  qi::rule<Iterator,AtomicPairPool*(),skipper_type> atomic_pair_pool;
  qi::rule<Iterator,vector<AtomicPairPool*>(),skipper_type> correlations;
  qi::rule<Iterator,vector<boost::fusion::tuple<string,double> >(),skipper_type> refinable_parameters;
  skipper_type skipper;
  skipper_type skipper_no_assignement; //< this is a skipper which does not parse arithmetic assignements. Used in 
  qi::rule<Iterator,void()> comment;
  
  qi::rule<Iterator,ADPMode*(),skipper_type> translational_mode;
  qi::rule<Iterator,ADPMode*(),skipper_type> rotational_mode;
  qi::rule<Iterator,void(),skipper_type> modes;
  qi::rule<Iterator,ADPMode*(),skipper_type> mode_assignement;
  qi::rule<Iterator,DoubleADPMode*(),skipper_type> adp_correlation;
  qi::rule<Iterator,SizeEffect*(),skipper_type> size_effect;
  qi::symbols<char,string> point_group_symbol;
  qi::symbols<char,string> atom_type;
  qi::symbols<char,string> molecular_type;

  
  qi::symbols<char,bool> calculation_methods;
  
  qi::rule<Iterator,void(),skipper_type> print_command_parser;
  qi::rule<Iterator,string()> string_in_quotes;

  
  qi::rule<Iterator,void(),skipper_type> program_option;
  qi::rule<Iterator,void(),skipper_type> program_option1;
  qi::rule<Iterator,void(),skipper_type> program_options;
  qi::rule<Iterator,void(),skipper_type> maybe_assignments;
  qi::rule<Iterator,Scatterer*(),skipper_type> molecular_scatterer;
  qi::rule<Iterator,vector<boost::tuple<string,Scatterer*> >(),skipper_type> molecular_scatterers;
  
  static ChemicalUnit* create_symmetric_chemical_unit(StructurePartRef _unit,SymmetryElement symm)
  {
    ChemicalUnit* res = boost::get<ChemicalUnit*>(_unit)->create_symmetric(symm.permutation_matrix,symm.displacement);
    
    return res;
  }
  
  static void add_reference(ReferenceTable& reference_table, string name,StructurePartRef references)
  {
    reference_table.add(name,references);
  }
  static SymmetryElement update_symmetry_component(SymmetryElement symm,vector<double> update, int component)
  {
    for(int i=0; i<3; i++)
      symm.permutation_matrix[i+component*3] = update[i];
    symm.displacement[component] = update[3];
    
    return symm;
  }
  
  void add_model(Model* _model)
  {
    model = _model;
  }
  Model* model;
};
/* TODO
 Make script for model bind functions and rewrite all functions in model so that they are not static. maybe destroy all inline functions
 May be it should be done in following way: all parsers return the objects. and objects are translated to their references only when they are registered
 */

#endif