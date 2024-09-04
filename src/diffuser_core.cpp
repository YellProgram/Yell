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
  
  //\TODO: assert that there is not too much dimensions in input vector
  
  for(int i=no_dimensions; i<3; i++)
    dims_out[i]=1; //fill additional dimensions if input matrix has less then 3 dimensions

  double * temp_map_holder = (double*) malloc(dims_out[0]*dims_out[1]*dims_out[2] * sizeof(double));
  
  dataset.read(temp_map_holder, PredType::NATIVE_DOUBLE);
  
  IntensityMap result(dims_out[0], dims_out[1], dims_out[2]); ///\TODO: copy data directly into versa, not like this, with additional array
  for(int i=0; i<dims_out[0]*dims_out[1]*dims_out[2]; i++)
    result.at(i) = temp_map_holder[i];
  
  free(temp_map_holder);

  return result;
}


template<typename T>
DataType getH5Type() {}

template<>
DataType getH5Type<int> () {
    return PredType::NATIVE_INT;
}
template<>
DataType getH5Type<bool> () {
    return PredType::NATIVE_HBOOL;
}
template<>
DataType getH5Type<double> () {
    return PredType::NATIVE_DOUBLE;
}

template <typename T>
void creadeAndWriteDataset(H5File& file, string datasetName, T* data, hsize_t n, hsize_t* dims) {
    DataSpace dataspace( n, dims );
    DataSet dataset = file.createDataSet( datasetName, getH5Type<T>(), dataspace );
    dataset.write( data, getH5Type<T>() );
}

template <typename T>
void writeConstant(H5File& file, string datasetName, T data) {
    H5::DataSet ds = file.createDataSet(datasetName, getH5Type<T>(), H5::DataSpace(H5S_SCALAR));
    ds.write(&data, getH5Type<T>());
}

template <typename T, std::size_t N>
void writeVector(H5File& file, string datasetName, af::tiny_plain<T,N> v) {
    T dataAsPlaneArray[N];
    for(int i=0; i<N; ++i)
        dataAsPlaneArray[i]=v[i];

    hsize_t sz[N]={N};
    creadeAndWriteDataset<T>(file, datasetName, dataAsPlaneArray, 1, sz);
}

void writeFormatString(H5File& file) {
    string format = "Yell 1.0";

    H5::StrType h5stringType(H5::PredType::C_S1, H5T_VARIABLE); // + 1 for trailing zero
    H5::DataSet ds = file.createDataSet("format", h5stringType, H5::DataSpace(H5S_SCALAR));
    ds.write(format, h5stringType);
}

void WriteHDF5(string filename, IntensityMap& input)
{
    H5File file( filename, H5F_ACC_TRUNC );

    hsize_t dimsf[3];
    //Fill it from input
    vec3<int> dims = input.size();
    for(int i=0; i<3; i++)
        dimsf[i] = dims[i];

    double * temp_map_holder = (double*) malloc(dimsf[0]*dimsf[1]*dimsf[2] * sizeof(double));

    for(int i=0; i<dimsf[0]*dimsf[1]*dimsf[2]; i++)
        temp_map_holder[i] = input.at(i);

    creadeAndWriteDataset<double>(file, "data", temp_map_holder, 3, dimsf);

    free(temp_map_holder);

    Grid& g=input.grid;
    // Write out the grid settings
    //creadeAndWriteDataset<double>(file, "lower_limits", {g.lower_limits[0], g.lower_limits[1], g.lower_limits[2]}, 1, {3});
    writeVector(file, "lower_limits", g.lower_limits);
    writeConstant(file, "is_direct", !g.reciprocal_flag);
    writeVector(file, "step_sizes", g.grid_steps);
    cctbx::uctbx::unit_cell cell = g.cell;
    if(!g.reciprocal_flag)
        cell=cell.reciprocal();

    writeVector(file, "unit_cell", cell.parameters());
    writeFormatString(file);

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
  
  if(esd==0 || isnan(esd)) //special case, if standard deviation is not set, report the result as-is
  {
    res << val << '(' << esd << ')';
    return res.str();
  }

  int leading_power_parameter=1000;
  if (abs(val)>=1e-32) {
      leading_power_parameter = floor(log(abs(val))/log((double)10));
  }

    //the position of the leading digit. positive are to the right of the zero, negative - to the left
    //we take abs of the esd because in ill-defined cases they can be calculated as negative
  int leading_power_esd = floor(log(abs(esd))/log((double)10));
  int leading_power = min(leading_power_esd, leading_power_parameter); //This way we will make sure that if a parameter is very very small, but not insignificant, it would still be reported
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

