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
#include <boost/phoenix/statement/try_catch.hpp>
#include <boost/phoenix/core/nothing.hpp>

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
    > "("  //due to some bug, it is necessary to turn try_... into sequence by adding noop in the beginninghttp://stackoverflow.com/questions/32245224/boostphoenix-try-catch-all-construct-fails-to-compile
    > (identifier > "," > identifier > "," > number % ',')[phoenix::nothing,
                                                             phoenix::try_
                                                              [
                                                               _val = bind(Model::correlators_from_cuns_,_1,_2,_3)
                                                               ].catch_all
                                                              [
                                                               _pass = false
                                                               ] ]
    > ")"
    ;
    
    adp_correlation =
    lit("ADPCorrelation")
    > '('
    > (identifier > ',' > identifier > ',' > number)[phoenix::nothing,
                                                     phoenix::try_
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
    > (identifier > ',' > identifier > ',' > number)[phoenix::nothing,
                                                     phoenix::try_
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