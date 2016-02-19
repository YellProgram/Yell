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

#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/fftpack/complex_to_complex_3d.h>
#include <complex>
#include <assert.h>
#include "basic_io.h"
#include "math.h"


#include <cctbx/uctbx.h> 

#include "OutputHandler.h"
#include "exceptions.h"

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


class Grid {
public:
    Grid()
    {}

    Grid(cctbx::uctbx::unit_cell cell,
         vec3<double> grid_steps,
         vec3<double> lower_limits,
         vec3<int> grid_size = vec3<int>(1, 1, 1),
         bool reciprocal_flag = true) : // by default the grid is defined in reciprocal space. Usually we start working with diffuse scattering
            cell(cell), grid_steps(grid_steps), lower_limits(lower_limits), grid_size(grid_size),
            reciprocal_flag(reciprocal_flag)
    {}

    Grid in_pdf_space()
    {
        if(!reciprocal_flag)
            return *this;
        else
            return reciprocal();
    }

    vec3<double> s_at(af::tiny<int,3> const & index)	{
        return lower_limits+ew_mult(grid_steps,index);
    }
    double d_star_square_at(af::tiny<int,3> const & index)	{
        return cell.d_star_sq(cctbx::miller::index<double>(s_at(index)));
    }


    /// Returns corresponding grid in reciprocal space.
    Grid reciprocal()	{
        vec3<double> res_grid_steps(1, 1, 1);
        vec3<double> res_lower_limits(0, 0, 0);
        for(int i=0; i<3; i++)
        {
            //cout << "lower limit " << i << " " << lower_limits[i];
            if(abs(lower_limits[i])!=0)
            {
                res_grid_steps[i]= -0.5 / lower_limits[i];
                res_lower_limits[i]= -0.5 / grid_steps[i];
            }
        }

        return Grid(cell.reciprocal(), res_grid_steps, res_lower_limits, grid_size, !reciprocal_flag);
    }

    /// Function for unit tests.
    bool operator==(const Grid& inp)  {
        return almost_equal(grid_steps,inp.grid_steps)\
                && almost_equal(lower_limits, inp.lower_limits) \
                && almost_equal(cell.parameters(),inp.cell.parameters()) \
                && grid_size==inp.grid_size;
    }

    vec3<double> upper_limits() {
        return lower_limits+grid_steps.each_mul(grid_size-1);
    }

    bool grid_is_compatible_with_fft() {
        for(int i=0; i<3; ++i)
        {
            if(grid_size[i]==1)
            {
                if(almost_equal(lower_limits[i],0))
                    continue;
                else
                    return false;
            }
            else
            {
                if(grid_size[i]%2==1)
                    return false;

                if(!almost_equal(lower_limits[i]+grid_size[i]/2*grid_steps[i], 0))
                    return false;
            }
        }

        return true;
    }

    /// Returns the volume which the grid covers.
    double volume() {
        return cell.volume()*grid_steps.each_mul(grid_size).product();
    }

    /// Returns the grid with padding by n pixels
    Grid pad(vec3<int> padding) {
        Grid res = *this;
        res.lower_limits-=grid_steps.each_mul(padding);
        res.grid_size+=2*padding;
        return res;
    }

    cctbx::uctbx::unit_cell cell;
    vec3<double> grid_steps; ///< array of steps of the grid in angstroms or reciprical angstroms
    vec3<double> lower_limits; ///< array of lowest bounds of grid
    vec3<int> grid_size;
    bool reciprocal_flag; ///< true if grid is in reciprocal coordinates (\AA^-1) and false if it is in PDF units.
};


/** 
 * Distributes r between vector which is compatible with grid and residual.
 * \param r input vector.
 * \param r_grid component of r which is compatible with grid. Equals number of steps from lowerest bounds to grid point which is closest to r.
 * \param r_res residual component of r.
 */
void grid_and_residual(vec3<double> r,Grid grid,vec3<int> & r_grid,vec3<double> & r_res);

/**
 * Holds 3D array of diffuse scattering intensity. 
 */
class IntensityMap {
public:
    IntensityMap() { }

    IntensityMap(const IntensityMap& inp) {
        Init(inp.data_accessor[0],inp.data_accessor[1],inp.data_accessor[2]);
        if(inp.grid_is_initialized)
            set_grid(inp.grid);

        for(int i=0; i<data_accessor.size_1d(); ++i)
            data[i]=inp.data[i];
    }

    IntensityMap(Grid grid) : grid(grid) {

        Init(grid.grid_size[0],grid.grid_size[1],grid.grid_size[2]);
        grid_is_initialized=true;
    }

    IntensityMap(vec3<int> map_size)
    {
        Init(map_size[0],map_size[1],map_size[2]);
    }
    IntensityMap(int a,int b, int c)
    {
        Init(a,b,c);
    }
    void set_grid(cctbx::uctbx::unit_cell cell,vec3<double> grid_steps, vec3<double> lower_limits)
    {
        set_grid( Grid(cell, grid_steps, lower_limits));
    }
    void set_grid(Grid inp_grid)
    {
        grid = inp_grid;
        grid_is_initialized=true;
    }

    void erase_data()
    {
        for(int i=0; i<data.size(); ++i)
            data[i]=0;
    }

    double& at(int a,int b,int c)
    {
        return real(data[data_accessor(a,b,c)]);
    }
    complex<double>& at_c(af::c_grid<3,int>::index_type _index)
    {
        return data[data_accessor(_index)];
    }
    double& at(af::c_grid<3,int>::index_type _index)
    {
        return real(data[data_accessor(_index)]);
    }
    double& at(int _index)
    {
        return real(data[_index]);
    }
    complex<double>& at_c(int _index)
    {
        return data[_index];
    }

    void subtract(IntensityMap& inp)
    {
        assert(data.size() == inp.data.size());
        for(int i=0; i<data.size(); ++i)
            data[i]-=inp.data[i];
    }

    IntensityMap padded (vec3<int> padding)
    {
        IntensityMap res(size()+2*padding);
        res.set_grid(grid.pad(padding));
        return res;
    }

    void copy_from_padded(vec3<int> padding,IntensityMap& pad)
    {
        for(int i=0; i<size()[0]; ++i)
            for(int j=0; j<size()[1]; ++j)
                for(int k=0; k<size()[2]; ++k)
                    at(i,j,k)=pad.at(i+padding[0],j+padding[1],k+padding[2]);
    }

    void scale_and_fft(double scale)
    {
        assert(grid_is_initialized);

        double pdf_space_volume = grid.volume();

        for(int i=0; i<data.size(); ++i)
            data[i] *= scale/pdf_space_volume/4; //dunno where this 1/4 comes from

        invert();
    }

    /// Performs fourier transform of intensity map and inverts grid
    void invert()
    {
        perform_shifted_fft();
        invert_grid();
    }

    /// Inverts grid, does nothing with intensity map
    void invert_grid()
    {
        if (grid_is_initialized)
            grid=grid.reciprocal();
    }

    /** Applies -1 to every element of array if sum Ni/2 is odd
     *
     * The accurate formula for shifted FFT looks the following
     * X_n-m = FT[ x_n-m * exp( n*2*pi*i*m/N ) ]*exp( (n-m)*2*pi*i*m/N )
     * given that m = N/2 it transforms into
     * X_n-m = FT[ x_n-m * -1^( mod(n,2) ) ]*-1^( mod(n-m,2) )
     * -1^( mod(n,2) ) - is a chessboard and was applied in the following functions
     * -1^( mod(n-m,2) ) = -1^( mod(n,2) ) when N/2 is even. In such cases this funciton does nothing.
     * -1^( mod(n-m,2) ) = -(-1^( mod(n,2) )) when N/2 is odd. In such cases this function applies the minus to all the elements in the array
     */
    void flip_signs_if_nessesary()
    {
        int disp=0;
        for(int i=0; i<3; ++i)
            if(data_accessor[i]>1)
                disp+=data_accessor[i];

        if((disp/2)%2==1)
            for(int i=0; i<data_accessor[0]; i++)
                for(int j=0; j<data_accessor[1]; j++)
                    for(int k=0; k<data_accessor[2]; k++)
                        data[data_accessor(i,j,k)]*=-1;

    }

    /// Multiplies the array to -1 in a chessboard-way. The center N/2,M/2,K/2 keeps positive sign
    void flip_signs_in_chessboard_way()
    {
        for(int i=0; i<data_accessor[0]; i++)
            for(int j=0; j<data_accessor[1]; j++)
                for(int k=0; k<data_accessor[2]; k++)
                    if((i+j+k)%2==1)
                        data[data_accessor(i,j,k)]*=-1;
    }

    /** performs fourier transform, assuming that origin of coordinates is in point grid_size/2.
     */
    void perform_shifted_fft()
    {
        flip_signs_in_chessboard_way();
        perform_fft();
        flip_signs_in_chessboard_way();
        flip_signs_if_nessesary();
    }

    void perform_fft()
    {
        scitbx::fftpack::complex_to_complex_3d<double> fft(data_accessor);
        scitbx::af::ref<std::complex<double>, scitbx::af::c_grid<3,int> >    cmap(data.begin(), data_accessor);

        if (is_in_reciprocal_space())
            fft.backward(cmap);
        else
        {
            fft.forward(cmap);

            // for some reason this forward and backward transforms do not divide by number of pixels in the map. This is why we apply it by hand
            double scale = 1.0/data_accessor.size_1d();
            init_iterator();
            while(next())
                current_array_value_c()*=scale;
        }

    }
    void init_iterator()
    {
        iterator=-1;//because next() will be used before first actions
    }

    bool next()
    {
        return (++iterator<data_accessor.size_1d());
    }

    double& current_array_value()
    {
        return real(data[data_accessor(current_index())]);
    }

    complex<double>& current_array_value_c()
    {
        return data[data_accessor(current_index())];
    }

    vec3<double> current_s()
    {
        return grid.s_at(current_index());
    }

    double current_d_star_square()
    {
        return grid.d_star_square_at(current_index());
    }

    cctbx::uctbx::unit_cell & unit_cell()
    {
        return grid.cell;
    }

    vec3<double> & grid_steps()
    {
        return grid.grid_steps;
    }

    vec3<int> size()
    {
        vec3<int> res;

        for(int i=0; i<3; i++)
            res[i]=data_accessor[i];

        return res;
    }

    int size_1d()
    {
        return data_accessor.size_1d();
    }


    af::c_grid<3,int>::index_type current_index()
    {
        return data_accessor.index_nd(iterator);
    }

    bool is_in_reciprocal_space()
    {
        return grid.reciprocal_flag;
    }

    void to_reciprocal()
    {
        if(!is_in_reciprocal_space())
            invert();
    }
    void to_real()
    {
        if(is_in_reciprocal_space())
            invert();
    }

    bool can_be_fourier_transformed() {
        return grid.grid_is_compatible_with_fft();
    }

    Grid grid;
private:
    void Init(int a,int b, int c)
    {
        data_accessor = af::c_grid<3,int>(a,b,c);
        data = af::versa<complex<double>, af::c_grid<3,int> > (data_accessor);
        grid_is_initialized=false;
    }

    af::versa<complex<double>, af::c_grid<3,int> > data;
    af::c_grid<3,int> data_accessor;

    int iterator;

    bool grid_is_initialized;

};

/**
 * Function for reading in data from hdf5 format. So far just reads in the "/data" dataset in a file.
 * \todo develop the format
 * \todo error handling
 */
IntensityMap ReadHDF5(string filename);

/**
 * function writes intensity map in a hdf5 file. 
 * \todo develop the format
 * \todo error handling
 */
void WriteHDF5(string filename,IntensityMap& input);



/**
 * Class for wrapping all the datasets which are optionaly provided by user, such as weights, measured pixels and resolution function correction.
 */
class OptionalIntensityMap {
public:
    OptionalIntensityMap() : is_loaded(false), default_value(1) //we mainly use these arrays for multipliers, and 1 is neutral there
    { }
    OptionalIntensityMap(double def_val) : is_loaded(false), default_value(def_val)
    { }

    bool is_loaded;
    double default_value;

    void load_data(string filename)
    {
        if(file_exists(filename))
        {
            intensity_map = ReadHDF5(filename);
            is_loaded = true;
        }
    }
    void load_data_and_report(string filename, string dataset_name, vec3<int> expected_size)
    {
        load_data(filename);
        if(is_loaded)
        {
            if(intensity_map.size()!=expected_size)
            {
                REPORT(ERROR) << "Dataset " << filename << " has incorrect size (" << intensity_map.size()[0] << " " << intensity_map.size()[1] << " " << intensity_map.size()[2]
                              << ")  expected (" << expected_size[0] << " " << expected_size[1] << " " << expected_size[2] << ")\n\n";
                cout.flush();
                throw(TerminateProgram());
            }

            for(int i=0; i<intensity_map.size_1d(); ++i)
                if(intensity_map.at(i)!=intensity_map.at(i)) //isnan
                {
                    REPORT(ERROR) << "Dataset " << filename << " contains NaNs\n\n";
                    cout.flush();
                    throw(TerminateProgram());
                }

            REPORT(MAIN) << "Loaded " << dataset_name << "\n";
        }
    }
    double at(int index)
    {
        if(is_loaded)
            return intensity_map.at(index);
        else
            return default_value;
    }
    IntensityMap* get_intensity_map() {
        assert(is_loaded);
        return &intensity_map;
    }
private:
    IntensityMap intensity_map;

};

#endif