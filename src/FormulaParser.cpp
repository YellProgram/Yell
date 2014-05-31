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

#include "FormulaParser.h"
#include "OutputHandler.h"
extern OutputHandler report;

// for the reasons explained here http://www.boost.org/doc/libs/1_53_0/libs/bind/bind.html it is impossible to bind templated funcitons to spirit semantic actions. Thus, we will have to make adapter to all the arithmetic functions.
double adopted_exp(double x)  { return exp(x);  }
double adopted_log(double x)  { return log(x);  }
double adopted_sin(double x)  { return sin(x);  }
double adopted_cos(double x)  { return cos(x);  }
double adopted_sqrt(double x) { return sqrt(x); }
double adopted_abs(double x)  { return abs(x);  }
double adopted_mod(double x,double y)  { return fmod(x,y);  }
double adopted_pow(double x,double y)  { return pow(x,y);  }

bool log_is_ok(double x) {
  if(x>0)
    return true;
  else
  {
    REPORT(ERROR) << "Calculation of logarithm failed, current argument is " << x << "\n";
    return false;
  }
}
bool sqrt_is_ok(double x) {
  if(x>=0)
    return true;
  else
  {
    REPORT(ERROR) << "Calculation of square root failed, current argument is " << x << "\n";
    return false;
  }
}

FormulaParser::FormulaParser() : FormulaParser::base_type(start),current_array_value(0)
{
  using namespace qi;
  using phoenix::bind;
  using phoenix::ref;
  
  fact %= double_ | '(' >> expr >> ')' | identifier | assignment | array_element | special_function; //
  
  special_function %= 
    ("exp(" >> expr >> ')') [_val = bind(&adopted_exp,_1)]
  | ("log(" >> expr >> ')') [_pass = bind(&log_is_ok,_1)][_val = bind(&adopted_log,_1)]
  | ("sin(" >> expr >> ')') [_val = bind(&adopted_sin,_1)]
  | ("cos(" >> expr >> ')') [_val = bind(&adopted_cos,_1)]
  | ("sqrt(">> expr >> ')') [_pass = bind(&sqrt_is_ok,_1)][_val = bind(&adopted_sqrt,_1)]
  | ("abs(" >> expr >> ')') [_val = bind(&adopted_abs,_1)]
  | ("mod(" >> expr >> ',' >> expr >> ')') [_val = bind(&adopted_mod,_1,_2)]
  | ("pow(" >> expr >> ',' >> expr >> ')') [_val = bind(&adopted_pow,_1,_2)]
  ;
  
  term %= fact[_val = _1]
  >> *( ('*' >> fact)[_val *= _1]
       | ('/' >> fact)[_val/=_1] )
  ;
  expr = ( (-lit('+') >> term)[_val=_1] 
          | ('-' >> term)[_val=- _1] ) 
  >> *( 
       ('+' >> term)[_val+=_1]
       | ('-' >> term)[_val-=_1]
       )
  ;
  // Deprecated
  //array_element =
  //lit("~i")[_val = bind(&extract_array_value,ref(array),ref(current_array_value)++)]
  //;
  
  // copypaste from InputParser better make a separate grammar for this
  // added minus to the charecters
  valid_identifier %= alpha >> *char_("a-zA-Z0-9_");
  identifier %= references >> !(alnum | '=');
  assignment = (valid_identifier >> '=' >> expr[_val = _1])[bind(add_key,ref(references),_1,_2)];//    
  start %= expr;
  

}