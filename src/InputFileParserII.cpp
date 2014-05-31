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

void InputParser::InputParserII()
{
  
  using namespace qi;
  using phoenix::bind;
  using phoenix::ref;
  
  symmetry_element = 
  lit("Symmetry") 
  > '('
  > permutation_component[_val = bind(&update_symmetry_component,_val,_1,0)]
  > ',' > permutation_component[_val = bind(&update_symmetry_component,_val,_1,1)]
  > ',' > permutation_component[_val = bind(&update_symmetry_component,_val,_1,2)]
  > ')'
  ;
  
  variant = (
              lit("Variant") [_val = phoenix::new_<ChemicalUnitNode>()]
             > '['
             > *(
                  '(' > omit[ char_("pP") ] > '=' >  number  > ')' >
                  chemical_unit[bind(&ChemicalUnitNode::add_chemical_unit,*_val,_1)] 
                  ) [bind(&ChemicalUnit::set_occupancy,*_2,_1)]
             > ']'
             > eps [_pass = bind(&ChemicalUnitNode::complain_if_sum_of_occupancies_is_not_one,_val)]
            );
  
  variant_assignement = (valid_identifier >> '=' >> variant[_val=_1])[bind(InputParser::add_reference,phoenix::ref(references),_1,phoenix::construct<StructurePartRef>(_val))];
   
  valid_identifier %= alpha >> *char_("a-zA-Z_0-9") >> !char_("a-zA-Z_0-9");
  
};