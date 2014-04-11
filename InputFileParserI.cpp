/*
 *  InputFileParserI.cpp
 *  diffuser_y
 *
 *  Created by Arkadiy Simonov on 3/23/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include "InputFileParser.h"
#include <boost/spirit/home/phoenix/statement/try_catch.hpp>

void InputParser::InputParserI()
{
  
  using namespace qi;
  using phoenix::bind;
  using phoenix::ref;


  
  cell_shifter = '(' >
    ( number > ',' > number > ',' > number)[_val =  phoenix::new_<CellShifter>(_1,_2,_3)]
    > ')'
    ;
  multiplicity_correlation = lit("Multiplicity") > number[_val = phoenix::new_<MultiplicityCorrelation>(_1)];
  
  substitutional_correlation = 
    lit("SubstitutionalCorrelation")
    > "(" 
  > (identifier > "," > identifier > "," > number % ',')[phoenix::try_
                                                         [
                                                          _val = bind(Model::correlators_from_cuns_,_1,_2,_3)
                                                         ].catch_all
                                                         [
                                                          _pass = false
                                                         ]]
    > ")"
    ;
  adp_correlation =
    lit("ADPCorrelation")
    > '('
  > (identifier > ',' > identifier > ',' > number)[phoenix::try_
                                                   [
                                                    _val = bind(&Model::create_double_adp_mode,_1,_2,_3)
                                                    ].catch_all
                                                   [
                                                    _pass = false
                                                    ]]
    > ')'
    ;
  
  size_effect =
    lit("SizeEffect")
    > '('
  > (identifier > ',' > identifier > ',' > number)[phoenix::try_
                                                   [
                                                    _val = bind(&Model::create_size_effect,_1,_2,_3)
                                                    ].catch_all
                                                   [
                                                    _pass = false
                                                    ]]
    > ')'
    ;
  
  number %= formula | double_;
  
  skipper_no_assignement = boost::spirit::ascii::space | comment;
  skipper = omit[(formula.assignment >> ';')] | skipper_no_assignement;
  
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
};