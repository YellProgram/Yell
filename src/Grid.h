//
// Created by Arkadiy Simonov on 08.08.24.
//

#ifndef YELL_GRID_H
#define YELL_GRID_H


#include "cctbx/uctbx.h"
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/fftpack/complex_to_complex_3d.h>

#include "diffuser_core.h"

using namespace scitbx;

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

    vec3<double> s_at(scitbx::af::tiny<int,3> const & index)	{
        return lower_limits+ew_mult(grid_steps, index);
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
void grid_and_residual(vec3<double> r, Grid grid, vec3<int> & r_grid,vec3<double> & r_res);



#endif //YELL_GRID_H
