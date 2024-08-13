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

#include "diffuser_core.h"
#include <H5Cpp.h> //for reading HDF5
using namespace H5;

//I figured out that i was using complex real() as accessor which it is not. Here is a hack to solve this under windows. ICC and GCC compilers work nicely without this hack
double& real(complex<double>& inp)
{
	return reinterpret_cast<double&> (inp);
}

double round(double d)
{
  return floor(d + 0.5);
}



bool almost_equal(double a, double b)
{
  return abs(a-b)<EPSILON;
}

bool almost_equal(vec3<double> a,vec3<double> b)
{
  return almost_equal<3,vec3<double> >(a,b);
}
bool almost_equal(mat3<double> a,mat3<double> b)
{
  return almost_equal<9,mat3<double> >(a,b);
}
bool almost_equal(scitbx::af::tiny<double,6> a,scitbx::af::tiny<double,6> b)
{
  return almost_equal<6,scitbx::af::tiny<double,6> >(a,b);
}

string format_esd(double val,double esd)
{
  std::ostringstream res;
  
  if(esd==0) //special case, if standard deviation is not set, report the result as-is
  {
    res << val << '(' << esd << ')';
    return res.str();
  }
    
  int leading_power = floor(log(esd)/log((double)10)); //the position of the leading digit. positive are to the right of the zero, negative - to the left
  int leading_digit = floor(esd * pow((double)10,-leading_power)); //the digits of the esd which should be reported
  if(leading_digit==1) {
    leading_power -= 1;
    leading_digit = floor(esd * pow((double)10,-leading_power)); //add one more digit
  }//TODO: else leading_digit = round(esd*pow....
  
  int precision = max(-leading_power,0);

  res.setf(std::ios::fixed);
  res.precision(precision);
  res << val << '(' << leading_digit << ')';
  
  return res.str();
}

