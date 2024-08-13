//
// Created by Arkadiy Simonov on 08.08.24.
//

#include "Grid.h"

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