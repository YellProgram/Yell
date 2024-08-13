//
// Created by Arkadiy Simonov on 08.08.24.
//

#ifndef YELL_INTENSITYMAP_H
#define YELL_INTENSITYMAP_H

#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/fftpack/complex_to_complex_3d.h>
#include <complex>
#include <assert.h>
#include "basic_io.h"
#include "math.h"



#include "Grid.h"

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

    void copy_from_padded(vec3<int> padding, IntensityMap& pad)
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
        grid_is_initialized = false;
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



#endif //YELL_INTENSITYMAP_H
