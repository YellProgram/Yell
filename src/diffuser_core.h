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

#ifndef DIFFUSER_CORE
#define DIFFUSER_CORE


#include <complex>
#include <assert.h>
#include "basic_io.h"
#include "math.h"


#include <cctbx/uctbx.h> 

#include "OutputHandler.h"
#include "exceptions.h"

#include <H5Cpp.h>
using namespace H5;

extern OutputHandler report;

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795029L
#endif
const double M_2PI=2*M_PI;
const double M2PISQ=-2*M_PI*M_PI;

using namespace scitbx;
using namespace std;



//I figured out that i was using complex real() as accessor which it is not. Here is a hack to solve this under windows. ICC and GCC compilers work nicely without this hack
double& real(complex<double>& inp);

const double EPSILON = 0.0000000001; //i picked it at random
bool almost_equal(double a, double b);
template<int matrix_size,class FixedSizeMatrix>
bool almost_equal(FixedSizeMatrix a,FixedSizeMatrix b)
{
    bool res=true;
    for(int i=0; i<matrix_size; ++i)
        res = res && almost_equal(a[i],b[i]);
    return res;
}
bool almost_equal(vec3<double> a,vec3<double> b);
bool almost_equal(mat3<double> a,mat3<double> b);
bool almost_equal(scitbx::af::tiny<double,6>,scitbx::af::tiny<double,6>);

string format_esd(double val,double esd);



//element-wise multiplication of vec3-type vectors
template<class vec3type, class another_array>
vec3type ew_mult(vec3type v1,another_array v2)
{
    return vec3type(v1[0]*v2[0],v1[1]*v2[1],v1[2]*v2[2]);
}

#endif