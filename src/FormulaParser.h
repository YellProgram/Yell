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
#ifndef H_FORMULA_PARSER
#define H_FORMULA_PARSER


#include <vector>
#include <string>
#include "precompiled_header.h"
#include <math.h>

using namespace std;
namespace phoenix = boost::phoenix;
namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

typedef iterator_ Iterator;
struct FormulaParser : qi::grammar<Iterator, double()>{
  FormulaParser();
  
  typedef qi::rule<Iterator, double()> formula_parser_rule;
  formula_parser_rule start;
  formula_parser_rule term,fact,expr,assignment,identifier,special_function;
  qi::rule<Iterator,string()> valid_identifier;
  typedef qi::symbols<char,double> ReferenceTable;
  ReferenceTable references;
  vector<double> array;
  int current_array_value;
  formula_parser_rule array_element;
  
  void initialize_refinable_variables(vector<string> names,vector<double> values) {
    int number_of_variables=names.size();
    for(int i=1; i<number_of_variables; ++i) // minus scale
    {
      references.add(names.at(i),values.at(i));
    }
  }
  
  static void add_key(ReferenceTable& table,string key,double val)
  {
    if(table.find(key)!=NULL)
      table.remove(key);
    
    table.add(key,val);
  }
  static double extract_array_value(vector<double> arr, int index)
  {
    return arr.at(index);
  }
};



 #endif
