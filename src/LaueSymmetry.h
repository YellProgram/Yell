//
// Created by Arkadiy Simonov on 09.08.24.
//

#ifndef YELL_LAUESYMMETRY_H
#define YELL_LAUESYMMETRY_H

#include <string>
#include <vector>
#include <cctbx/sgtbx/rt_mx.h>
#include "Grid.h"
#include "IntensityMap.h"


class AtomicPair;
class UnitCell;

using std::string;
using std::vector;
using cctbx::sgtbx::rt_mx;

class LaueSymmetry
{
public:
    string label;

    /** Sometimes grid for diffuse scattering or PDF is incompatible with crystal Laue symmetry. In such cases generators that are incompatible will be applied to AtomicPairs
     */
    vector<string> generators_on_map;

    /// Generators of Laue group that could not be applied to the grid directly
    vector<string> generators_on_vectors;

    LaueSymmetry()
    {}

    /** This function separates the symmetry generators in those which should be applied on the map and those which should be appliet to AtomicPairs.

        Note that if the grid is not provided, Laue class will try to apply all generators on the level of diffuse scattering map and user is responsible to control correctess of the grid steps, and centering.
     */
    LaueSymmetry(string const & sym,
                 Grid const & grid = Grid(cctbx::uctbx::unit_cell(),
                                          vec3<double>(0,0,0),
                                          vec3<double>(0,0,0),
                                          vec3<int>(1,1,1)) ) :
            label(sym)   {
        vector<string> generators = laue_class_generators();

        for (int i=0; i<generators.size(); ++i)
        {
            if (generator_can_be_applied_on_the_map(generators[i], grid))
                generators_on_map.push_back(generators[i]);
            else
            {
                REPORT(FIRST_RUN) << "Diffuse scattering grid is incompatible with generator " << generators[i] << ". Generator will be applied directly to atomic pairs (slow)\n";
                generators_on_vectors.push_back(generators[i]);
            }
        }
    }

    /// used in generator_can_be_applied_on_the_map. Checks if grids in two directions are compatible ie have the same step, nuber of pixels and lower limits
    bool have_same_step_and_ll(Grid const & grid,string const & axes)    {
        int ax1=0,ax2;
        for (int i=1; i<axes.size(); i+=2) {
            //assume we have only xy and xz symbols there
            if(axes[i]=='y')
                ax2=1;
            else
                ax2=2;

            if( !(grid.grid_size[ax1]==grid.grid_size[ax2] &&
                  almost_equal(grid.grid_steps[ax1],grid.grid_steps[ax2]) &&
                  almost_equal(grid.lower_limits[ax1],grid.lower_limits[ax2])) )
                return false;

        }
        return true;
    }

    /// used in generator_can_be_applied_on_the_map. checks if the center of grid in some direction is at pixel number N/2+1 (counting from 1) for N even or in case N=1 that lower limit is equal to 0;
    bool centered_correctly(Grid const & grid,string const & axes)    {
        int ax;
        for(int i=0;i<axes.size();++i)
        {
            // transform symbols xyz to 012
            if (axes[i]=='x')
                ax=0;
            else if(axes[i]=='y')
                ax=1;
            else if(axes[i]=='z')
                ax=2;

            // actual checks
            if (grid.grid_size[ax]==1)
            {
                if(!(almost_equal(grid.lower_limits[ax],0)))
                    return false;
            }
            else if(grid.grid_size[ax]%2==1)
                return false;
            else if(!almost_equal(-grid.lower_limits[ax],0.5*grid.grid_size[ax]*grid.grid_steps[ax]))
                return false;
        }
        return true;
    }

    /// Used to separate generators in two groups. One group is applied to the map (fast) the other is applied to AtomicPairs (slow)
    bool generator_can_be_applied_on_the_map(string generator,Grid const & grid)    {
        if(generator=="-x,-y,-z")
            return centered_correctly(grid,"xyz");
        else if(generator=="-x,y,z")
            return centered_correctly(grid,"x");
        else if(generator=="x,-y,z")
            return centered_correctly(grid,"y");
        else if(generator=="x,y,-z")
            return centered_correctly(grid,"z");
        else if(generator=="z,x,y")
            return have_same_step_and_ll(grid,"xyxz");
        else if(generator=="y,x,z")
            return have_same_step_and_ll(grid,"xy");
        else if(generator=="-y,x-y,z")
            return have_same_step_and_ll(grid,"xy") && centered_correctly(grid,"xy");
        else if(generator=="y,x,-z")
            return have_same_step_and_ll(grid,"xy") && centered_correctly(grid,"z");
        else if(generator=="y,-x,z")
            return have_same_step_and_ll(grid,"xy") && centered_correctly(grid,"xy");
        else
            return true; // or rather throw error
    }

    vector<string> laue_class_generators()  {
        vector<string> generators;
        if(label=="m-3m")
        {
            generators.push_back("z,x,y");
            generators.push_back("-x,-y,-z");
            generators.push_back("-x,y,z");
            generators.push_back("x,-y,z");
            generators.push_back("y,x,z");
        }else if(label=="m-3")
        {
            generators.push_back("z,x,y");
            generators.push_back("-x,-y,-z");
            generators.push_back("-x,y,z");
            generators.push_back("x,-y,z");
        }else if(label=="6/mmm")
        {
            generators.push_back("-x,-y,-z");
            generators.push_back("y,x,z");
            generators.push_back("x,y,-z");
            generators.push_back("-y,x-y,z");
        }
        else if(label=="6/m")
        {
            generators.push_back("-x,-y,-z");
            generators.push_back("x,y,-z");
            generators.push_back("-y,x-y,z");
        }
        else if(label=="-3mH")
        {
            generators.push_back("-x,-y,-z");
            generators.push_back("y,x,-z");
            generators.push_back("-y,x-y,z");
        }
        else if(label=="-3mR")
        {
            generators.push_back("-x,-y,-z");
            generators.push_back("y,x,z");
            generators.push_back("z,x,y");
        }
        else if(label=="-3H")
        {
            generators.push_back("-x,-y,-z");
            generators.push_back("-y,x-y,z");
        }
        else if(label=="-3R")
        {
            generators.push_back("-x,-y,-z");
            generators.push_back("z,x,y");
        }
        else if(label=="4/mmm")
        {
            generators.push_back("-x,-y,-z");
            generators.push_back("-x,y,z");
            generators.push_back("x,-y,z");
            generators.push_back("y,x,z");
        }
        else if(label=="4/m")
        {
            generators.push_back("y,-x,z"); // four fold rotation
            generators.push_back("-x,-y,-z");

        }
        else if(label=="mmm")
        {
            generators.push_back("-x,-y,-z");
            generators.push_back("-x,y,z");
            generators.push_back("x,-y,z");
        }
        else if(label=="2/m")
        {
            generators.push_back("-x,-y,-z");
            generators.push_back("x,y,-z");
        }
        else if(label=="2/mb")
        {
            generators.push_back("-x,-y,-z");
            generators.push_back("x,-y,z");
        }else if(label=="-1")
        {
            generators.push_back("-x,-y,-z");
        }

        return generators;
    }

    vector<AtomicPair> filter_pairs_from_asymmetric_unit(vector<AtomicPair>& pairs);

    void apply_patterson_symmetry(IntensityMap& map)
    {
        for (int i=0; i<generators_on_map.size(); ++i)
            apply_generator(map,generators_on_map[i]);
    }

    vector<AtomicPair> apply_patterson_symmetry(vector<AtomicPair> pairs);

    /// Applies a generator to IntensityMap
    void apply_generator(IntensityMap& map, string generator_symbol);

    bool is_compatible_with_cell(const UnitCell& cell);

    /// Applies transformation matrix to atomic pair
    vector<AtomicPair> multiply_pairs_by_matrix(vector<AtomicPair> pairs,mat3<double> transformation_matrix);

    /// Applies
    vector<AtomicPair> filter_pairs_from_asymmetric_unit_and_scale(vector<AtomicPair>& pairs);

    /// Applies generator which is expressed in matrix form to pairs.
    /**
        \param transformation_matrix in the trasformation matrix to be applied
        \param times the class of generator
     */
    vector<AtomicPair> apply_matrix_generator(vector<AtomicPair> pairs,mat3<double> transformation_matrix,int times);


    //`m-3m` `m-3` `6/mmm` `6/m` `4/mmm` `4/m` `-3:R` `-3:H` `-3m:H` `-3m:R` `mmm` `2/m` `2/m:b` `-1`

    //TODO: split the in asu functions for each particular symmetry, so they are fast. Choose one and use it consistently later.
    /// Checks if a particular h, k, l index is in asu
    template<typename T>
    bool hkl_in_asu(T h, T k, T l){
        if(label=="m-3m") {
            // ASU from here: https://www.globalsino.com/EM/page4547.html
            //0 ≤ x ≤ 1/2 and 0 ≤ y ≤ 1/2 and 0 ≤ z ≤ 1/2 and y ≤ x and z ≤ y
            return h<=k && k<=l && l<=0;
        } else if(label == "m-3") {
            assert(false);
            // 	0≤ x ≤ 2/3 and 0 ≤ y ≤ 1/3 and 0 ≤ z ≤ 1 and x ≤ (1+y)/2 and y ≤ x/2
            return h<=0 && k<=0;
        } else if(label == "6/mmm") {
            //0 ≤ x ≤ 2/3 and 0 ≤ y ≤ 1/3 and 0 ≤ z ≤ 1/2 and x ≤ (1+y)/2 and y ≤ x/2
            assert(false);
            return h<=0 && k <= 0 && l <= 0 && true;
        }


        return true; //in case of doubt, the full space is in ASU. Later die here.
    }

    bool asu_point_multiplicity(double h, double k, double l) {
        return 1;
    }
};



#endif //YELL_LAUESYMMETRY_H
