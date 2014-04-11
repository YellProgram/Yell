/*
 *  diffuser_core.cpp
 *  diffuser_y
 *
 *  Created by Arkadiy Simonov on 2/2/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
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

void grid_and_residual(vec3<double> r,Grid grid,vec3<int> & r_grid,vec3<double> & r_res)
{
	r=r-grid.lower_limits;
	
	vec3<double> r_in_steps=r/grid.grid_steps;

	for(int i=0; i<3; i++)
	{
		r_grid[i]=int(  round(r_in_steps[i]) );
		r_res[i]=r_in_steps[i]-r_grid[i];
	}
	
	r_res=r_res.each_mul(grid.grid_steps);
}

IntensityMap ReadHDF5(string filename)
{
  H5File file( filename, H5F_ACC_RDONLY );
  DataSet dataset = file.openDataSet( "data" );
  
  DataSpace dataspace = dataset.getSpace();
  
  hsize_t dims_out[3];
  int no_dimensions = dataspace.getSimpleExtentDims( dims_out, NULL);
  
  //\todo: assert that there is not too much dimensions in input vector
  
  for(int i=no_dimensions; i<3; i++)
    dims_out[i]=1; //fill additional dimensions if input matrix has less then 3 dimensions
    
  
  double * temp_map_holder = (double*) malloc(dims_out[0]*dims_out[1]*dims_out[2] * sizeof(double));
  
  dataset.read( temp_map_holder, PredType::NATIVE_DOUBLE);
  
  IntensityMap result(dims_out[0],dims_out[1],dims_out[2]); ///\todo: copy data directly into versa, not like this, with additional array
  for(int i=0; i<dims_out[0]*dims_out[1]*dims_out[2]; i++)
    result.at(i) = temp_map_holder[i];
  
  free(temp_map_holder);

  return result;
}

void WriteHDF5(string filename,IntensityMap& input)
{
  H5File file( filename, H5F_ACC_TRUNC );
  
  hsize_t dimsf[3]; 
  // fill it from input
  vec3<int> _dims = input.size();
  for(int i=0; i<3; i++)
    dimsf[i] = _dims[i];
  
  DataSpace dataspace( 3, dimsf );
  
  DataSet dataset = file.createDataSet( "data", PredType::NATIVE_DOUBLE, dataspace );
  
  double * temp_map_holder = (double*) malloc(dimsf[0]*dimsf[1]*dimsf[2] * sizeof(double));
  
  for(int i=0; i<dimsf[0]*dimsf[1]*dimsf[2]; i++)
    temp_map_holder[i] = input.at(i);
  
  dataset.write( temp_map_holder, PredType::NATIVE_DOUBLE );
  
  free(temp_map_holder);
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

