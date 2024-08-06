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

#ifndef basic_classes_H
#define basic_classes_H

///////#define SCITBX_FFTPACK_COMPLEX_TO_COMPLEX_3D_NO_PRAGMA_OMP
#include <iostream>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/fftpack/complex_to_complex_3d.h>
#include <cctbx/sgtbx/rt_mx.h>
#include <sstream>

#include "diffuser_core.h"

#include <complex>
#include <map>
#include <math.h>

#include "levmar.h" //for minimizer

#include <scitbx/vec3.h>
#include <scitbx/sym_mat3.h>
#include <vector>
#include <string>

#include <cctbx/eltbx/neutron.h> //these three are for AtomicType
#include <cctbx/eltbx/electron_scattering.h>
#include <cctbx/eltbx/xray_scattering.h>
#include <cctbx/eltbx/xray_scattering/gaussian.h>

#include <cctbx/uctbx.h> //These guys are for cell and 
#include <scitbx/array_family/tiny.h>
#include <assert.h>

#include "OutputHandler.h"
using namespace std;
using namespace scitbx;
using namespace cctbx::sgtbx;

extern OutputHandler report;

const int VectorStart=0;
const int VectorEnd=1;

/// Skips check for matrix being symmetric. Just takes upper-triangular part.
sym_mat3<double> trusted_mat_to_sym_mat(mat3<double> inp);

class AtomicPairPool;

class PairModifier {
 public:
  virtual ~PairModifier() { };
  virtual bool generates_pairs()=0;
  virtual void modify_pairs(AtomicPairPool* const)=0;
};




/**	\brief Class of pointer vectors.
 
 Unsafe class which holds a vector of pointers to some objects,
 gives access to them with operator[] and destroys objects that
 it is points to on destruction.
 
 Probably, it should be replaced with boost ptr_vector
 */
template<class Obj>
class p_vector{
 public:
	
  p_vector(p_vector const &in) {
    assert(false); //copy constructur is not implemented
  }
  
  p_vector() : pointer_vector()
  {}
	
  ~p_vector()
  {
    for(int i=0; i<pointer_vector.size(); i++)
	    delete pointer_vector[i];
  }
	
	void push_back(Obj* a)
	{
		pointer_vector.push_back(a);
	}
	
	void concat(vector<Obj*> a)
	{ 
		pointer_vector.insert( pointer_vector.end(), a.begin(), a.end() );
	}
	
	int size()
	{
		return pointer_vector.size();
	}
	
	Obj& operator[](int i) 
	{
		return (*pointer_vector[i]);
	}
	
private:
	vector<Obj*> pointer_vector;
};

class Atom;
class AtomicPair;
class UnitCell;
/// Holds the laue symmetry of crystal and has function to apply the symmetry. Also keeps track of which elements should be applied on the level of 
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
      // Я должен признать поражение
      // в этих лабиринтах я
      // странная манера искать правила и им следовать под наркотой
      // Думать об абстрактных концепциях гораздо круче под этой наркотой. Алгебра, все это решает
      // я заблудился, стоит мне испепелять этот баг, или есть тчо-то более полезное что я все-таки могу сделать?
      for (int i=0; i<generators.size(); ++i)
      {
        if (generator_can_be_applied_on_the_map(generators[i],grid))
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
    // радужная пена - это подарок нам из ЛСД
    // все чувства задействованы
    void apply_patterson_symmetry(IntensityMap& map)
    {
      for (int i=0; i<generators_on_map.size(); ++i)
        apply_generator(map,generators_on_map[i]);
    }
    
    vector<AtomicPair> apply_patterson_symmetry(vector<AtomicPair> pairs)
    {
      
      for(int i=0; i<generators_on_vectors.size(); ++i)
      {
        rt_mx generator(generators_on_vectors[i]);
        pairs = apply_matrix_generator(pairs,generator.r().as_double(),generator.r().order());
      }
      
      return pairs;
    }
    
    /// Applies generator to IntensityMap
    /** 
        Explicitely applies all generators. Hopefully this way is faster.
     */
    void apply_generator(IntensityMap& map,string generator_symbol)
    {
      double t;
      vec3<int> size = map.size();
      if(generator_symbol=="-x,-y,-z") // -1
        for(int i=size[0]/2; i<size[0]; ++i)
          for(int j=1; j<size[1]; ++j)
            for(int k=1; k<size[2]; ++k)
            {
              t=(map.at(i,j,k) + map.at(size[0]-i,size[1]-j,size[2]-k))/2;
              map.at(i,j,k) = t;
              map.at(size[0]-i,size[1]-j,size[2]-k) = t;
            }
      else if(generator_symbol=="-x,y,z") // mx
        for(int i=size[0]/2; i<size[0]; ++i)
          for(int j=0; j<size[1]; ++j)
            for(int k=0; k<size[2]; ++k)
            {
              t=(map.at(i,j,k) + map.at(size[0]-i,j,k))/2;
              map.at(i,j,k) = t;
              map.at(size[0]-i,j,k) = t;
            }
      else if(generator_symbol=="x,-y,z") //my
        for(int i=0; i<size[0]; ++i)
          for(int j=size[1]/2; j<size[1]; ++j)
            for(int k=0; k<size[2]; ++k)
            {
              t=(map.at(i,j,k) + map.at(i,size[1]-j,k))/2;
              map.at(i,j,k) = t;
              map.at(i,size[1]-j,k) = t;
            }
      else if(generator_symbol=="x,y,-z") //mz
        for(int i=0; i<size[0]; ++i)
          for(int j=0; j<size[1]; ++j)
            for(int k=size[2]/2; k<size[2]; ++k)
            {
              t=(map.at(i,j,k) + map.at(i,j,size[2]-k))/2;
              map.at(i,j,k) = t;
              map.at(i,j,size[2]-k) = t;
            }
      else if(generator_symbol=="z,x,y")      //3 along 111
        for(int i=0; i<size[0]; ++i)
          for(int j=0; j<size[1]; ++j)
            for(int k=0; k<size[2]; ++k)
            {
              t=(map.at(i,j,k) + map.at(k,i,j) + map.at(j,k,i))/3;
              map.at(i,j,k) = t;
              map.at(k,i,j) = t;
              map.at(j,k,i) = t;
            }
      else if(generator_symbol=="y,x,z") //mx-y
        for(int i=0; i<size[0]; ++i)
          for(int j=0; j<size[1]; ++j)
            for(int k=0; k<size[2]; ++k)
            {
              t=(map.at(i,j,k) + map.at(j,i,k))/2;
              map.at(i,j,k) = t;
              map.at(j,i,k) = t;
            }
      else if(generator_symbol=="-y,x-y,z") // three-fold rotation axis in hexagonal settings
        if(map.grid.reciprocal_flag) 
        {
          //we need to apply symmetry element reciprocal to -y,x-y,z namely y,-x-y,z and -x-y,x,z
          for(int i=1; i<size[0]; ++i)
            for(int j=max(size[0]/2-i+1,1); j<min(size[0]*3/2-i,size[1]); ++j) //note that here are strange borders
              for(int k=0; k<size[2]; ++k)
              {
                t=(map.at(i,j,k) + map.at(j,size[0]*3/2-i-j,k) + map.at(size[0]*3/2-i-j,i,k))/3;
                map.at(i,j,k) = t;
                map.at(size[0]*3/2-i-j,i,k) = t;
                map.at(j,size[0]*3/2-i-j,k) = t;
              }
        }
        else {
          for(int i=1; i<size[0]; ++i)
            for(int j=max(i-size[0]/2+1,1); j<min(i+size[0]/2,size[1]); ++j) //note that here are strange borders
              for(int k=0; k<size[2]; ++k)
              {
                t=(map.at(i,j,k) + map.at(size[0]-j,size[0]/2+i-j,k) + map.at(size[0]/2+j-i,size[0]-i,k))/3;
                map.at(i,j,k) = t;
                map.at(size[0]-j,size[0]/2+i-j,k) = t;
                map.at(size[0]/2+j-i,size[0]-i,k) = t;
              }
        }
      else if(generator_symbol=="y,x,-z") //2-fold xy direction
        for(int i=0; i<size[0]; ++i)
          for(int j=0; j<size[1]; ++j)
            for(int k=size[2]/2; k<size[2]; ++k)
            {
              t=(map.at(i,j,k) + map.at(j,i,size[2]-k))/2;
              map.at(i,j,k) = t;
              map.at(j,i,size[2]-k) = t;
            }
      else if(generator_symbol=="y,-x,z") // Four-fold rotation applied four times.
        for(int j=size[0]/2; j<size[1]; ++j)
          for(int i=size[1]/2; i<size[0]; ++i)
             for(int k=0; k<size[2]; ++k)
            {
              t=(map.at(i,j,k) + map.at(j,size[1]-i,k) + map.at(size[0]-j,i,k) + map.at(size[0]-i,size[1]-j,k))/4;
              map.at(i,j,k) = t;
              map.at(j,size[1]-i,k) = t;
              map.at(size[0]-j,i,k) = t;
              map.at(size[0]-i,size[1]-j,k) = t;
            }

        // And then a special case for j=0 since there we need to make sure that -y coordinate is again at 0 pixel, or is it????
        int j=0;
        for(int i=1; i<size[0]; ++i)
            for(int k=0; k<size[2]; ++k)
            {
                t=(map.at(i,0,k) + map.at(0,size[1]-i,k) + map.at(0,i,k) + map.at(size[0]-i,0,k))/4;
                map.at(i,0,k) = t;
                map.at(0,size[1]-i,k) = t;
                map.at(0,i,k) = t;
                map.at(size[0]-i,0,k) = t;
            }
    }
    
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
};
class ChemicalUnit
	{
	public:
    virtual ~ChemicalUnit() { };
		virtual double get_occupancy()=0;
		virtual void set_occupancy(double)=0;
		virtual vector<Atom*> get_atoms()  =0;
		virtual ChemicalUnit* create_symmetric(mat3<double>,vec3<double>)  =0;
		virtual ChemicalUnit* operator[](int) =0;
	};

class AtomicAssembly : public ChemicalUnit
	{
	public:
		AtomicAssembly()
		{}
		
    AtomicAssembly(vector<ChemicalUnit*> units)
    {
      for(int i=0; i<units.size(); i++)
        chemical_units.push_back(units[i]);
    }
    
		void add_chemical_unit(ChemicalUnit * unit)
		{
			chemical_units.push_back(unit);
		}
		
		double get_occupancy()
		{	return occupancy;}
		void set_occupancy(double _occ)
		{	
      occupancy=_occ;
      for(int i=0; i<chemical_units.size(); ++i)
        chemical_units[i].set_occupancy(_occ);
    }
		
		vector<Atom*> get_atoms() 
		{
			vector<Atom*> atoms;
			vector<Atom*> t;
			
			for(int i=0; i<chemical_units.size(); i++)
			{
				t=chemical_units[i].get_atoms();
				atoms.insert(atoms.end(),t.begin(),t.end());
			}
			
			return atoms;
		}
		
		AtomicAssembly* create_symmetric(mat3<double> sym_matrix,vec3<double> translation)
		{
			AtomicAssembly* result = new AtomicAssembly();
			
			for(int i=0; i<chemical_units.size(); i++)
				result->add_chemical_unit(chemical_units[i].create_symmetric(sym_matrix,translation));
			
			return result;
		}
		
		ChemicalUnit* operator[](int n)
		{
			return &chemical_units[n];
		}
		
		double occupancy;
		p_vector<ChemicalUnit> chemical_units;
	};

class ChemicalUnitNode{
public:
  ChemicalUnitNode()
  {}
  
  void add_chemical_unit(ChemicalUnit* unit)
  {
    chemical_units.push_back(unit);
  }
  
  p_vector<ChemicalUnit> chemical_units;
  
  bool complain_if_sum_of_occupancies_is_not_one()
  {
    double sum=0;
    for(int i=0; i<chemical_units.size(); ++i)
      sum+=chemical_units[i].get_occupancy();
    
    if(almost_equal(1,sum))
      return true;
    else
      REPORT(ERROR) << "The sum of occupancies should be 1. Here it is " << sum << "\n";
    
    return false;
  }
  
};

class UnitCell
{
public:
  cctbx::uctbx::unit_cell cell;
  LaueSymmetry laue_symmetry;
  
  UnitCell() { }
  
  UnitCell(double a,double b, double c, double alpha, double beta, double gamma) : 
  cell(scitbx::af::tiny<double,6>(a,b,c,alpha,beta,gamma))
  { }
  
  void set_laue_symmetry(string sym)
  {
    laue_symmetry=LaueSymmetry(sym);
  }
  
  void add_node(ChemicalUnitNode* node)
  {
    chemical_unit_nodes.push_back(node);
  }
  
  /**
   * Function for test purposes. For now checks that unit cell parameters are almost equal.
   */
  bool operator==(const UnitCell& inp)
  {
    return almost_equal(cell.parameters(),inp.cell.parameters());
  }
  
  string to_string()
  {
    std::ostringstream res;
    for(int i=0; i<6; ++i)
      res << cell.parameters()[i] << " ";
    
    return res.str();
  }
  
//private:
  p_vector<ChemicalUnitNode> chemical_unit_nodes;
};

class ADPMode;
sym_mat3<double> outer_product(vec3<double> v1, vec3<double> v2);
class SubstitutionalCorrelation;
ADPMode* translational_mode(ChemicalUnit* unit,int direction,sym_mat3<double>);
ADPMode z_rot_mode(ChemicalUnit* unit);
ADPMode* rot_mode(ChemicalUnit* unit,vec3<double> axis , vec3<double> point_on_axis ,sym_mat3<double> metrical_matrix);
ADPMode combine_modes(vector<ADPMode> modes);
vector<SubstitutionalCorrelation*> correlators_from_cuns(ChemicalUnitNode* node1,ChemicalUnitNode* node2,vector<double> corr);

enum ScatteringType {XRay, Neutron, Electrons}; //TODO: after changing to C++11 change to enum class

class Scatterer;

class AtomicTypeCollection
	{
	public:
		AtomicTypeCollection() {	}
		~AtomicTypeCollection();
		static AtomicTypeCollection& get();
		map<string,Scatterer*> types;
    
    static void calculate_form_factors_of_all_atoms_on_grid(vec3<int>map_size, Grid grid);
    
    static void update_current_form_factors(vec3<double> s, double d_star_square);
    
    /** Creates AtomicType and registers it in ~AtomicTypeCollection. if it is already registered returns existing AtType
     */
    static Scatterer* get(string const & label, ScatteringType t);
    static string strip_label(string const & label);
    static void add(string const & label, Scatterer* s);
	};

class Scatterer {
  
public:
  virtual complex<double> form_factor_at_c(vec3<double>s, double d_star_sq)=0;
  virtual double form_factor_at(double d_star_sq)=0;
  virtual ~Scatterer() {};
  void update_current_form_factor(vec3<double>s,double d_star_sq) {
    current_form_factor=form_factor_at_c(s,d_star_sq);
  }
  
  void calculate_form_factors_on_grid(vec3<int>map_size,Grid grid)
	{
		gridded_form_factors = IntensityMap(map_size);
		gridded_form_factors.set_grid(grid);
		
		gridded_form_factors.init_iterator();
		while(gridded_form_factors.next())
		{
      AtomicTypeCollection::update_current_form_factors(gridded_form_factors.current_s(), gridded_form_factors.current_d_star_square());
			gridded_form_factors.current_array_value_c()=current_form_factor;
		}
	}
  
  IntensityMap gridded_form_factors;
	complex<double> current_form_factor; ///< here form factor for direct calculation will be cached
};

//gives access to computation of X-ray atomic form factor
class AtomicType : public Scatterer {
public:
	AtomicType()
	{}
  
	AtomicType(string const & _label) {
		cctbx::eltbx::xray_scattering::wk1995 wk(_label);
		gauss=wk.fetch();
	}
  
	double form_factor_at(double d_star_sq) {
		return gauss.at_d_star_sq(d_star_sq); //TODO: check funny sequence here. The entry with _c seems to be real which is ok for normal scattering though
	}

  complex<double> form_factor_at_c(vec3<double> s, double d_star_sq) {
    return form_factor_at(d_star_sq);
  }
  
private:
	cctbx::eltbx::xray_scattering::gaussian gauss;
};

// TODO: add tests
// TODO: remove copy-pastes in this code
//gives access to computation of Neutron atomic form factor
class NeutronScattererAtom : public Scatterer {
public:
    NeutronScattererAtom()
    {}

    NeutronScattererAtom(string const & _label) :
            scat_table_entry(_label)
    {}

    double form_factor_at(double) {
        return 0; //TODO: check why is this one needed at all
    }

    complex<double> form_factor_at_c(vec3<double>s,double d_star_sq) {
        return scat_table_entry.bound_coh_scatt_length();
    }

private:
    cctbx::eltbx::neutron::neutron_news_1992_table scat_table_entry;
};

class ElectronScattererAtom : public Scatterer {
public:
    ElectronScattererAtom()
    {}

    ElectronScattererAtom(string const & _label) {
        cctbx::eltbx::electron_scattering::peng1996 wk(_label);
        gauss=wk.fetch();
    }

    double form_factor_at(double d_star_sq) {
        return gauss.at_d_star_sq(d_star_sq); //TODO: check funny sequence here. The entry with _c seems to be real which is ok for normal scattering though
    }

    complex<double> form_factor_at_c(vec3<double> s, double d_star_sq) {
        return form_factor_at(d_star_sq);
    }

private:
    cctbx::eltbx::xray_scattering::gaussian gauss;
};


class MolecularScatterer : public Scatterer {
public:
  MolecularScatterer(ChemicalUnit* molecule) {
    constituent_atoms = molecule->get_atoms();
  }
  
  double form_factor_at(double t) {
    return 0;
  }
  
  complex<double> form_factor_at_c(vec3<double>s,double d_star_sq);
  
private:
  vector<Atom*> constituent_atoms;
  
  inline static complex<double> form_factor_in_a_point(complex<double> f,
                                                       double p,
                                                       double N,
                                                       vec3<double> r,
                                                       sym_mat3<double> U,
                                                       vec3<double> s);
};

//This is just a structure which holds occupancy, position and ADP-tensor
//used in atom and AtomicPair
class AtomicParams { 
public:
  double occupancy;
	vec3<double> r;
	sym_mat3<double> U;
	
	//on windows cygwin gets _U defined as numeric constant from somewhere
	AtomicParams(double _occupancy, vec3<double> _r, sym_mat3<double> _Uinp) : occupancy(_occupancy), r(_r), U(_Uinp)
	{}
	AtomicParams()
	{}
  
  /// Operator for test functions. returns true if all variables of both input classes are almost equal
  bool operator== (const AtomicParams& inp) const
  {
    return almost_equal(occupancy,inp.occupancy) && almost_equal(r,inp.r) && almost_equal(U,inp.U);// TODO:check ADPs
  }
};

///TODO: make sure that precessor constructors are not accessible
///TODO: change the name U to beta
class Atom : public ChemicalUnit
{
public:
  string label;
	Scatterer* atomic_type;
  double occupancy;
  double multiplier; //< Variable which holds say a multiplier which comes from symmetry. Like occupancy, but is not affected by the SubstitutionalCorrelation logic
	vec3<double> r;
	sym_mat3<double> U;
  
	Atom()
	{}
  
  /// Atom initialization with Uiso and metric tensor.
  Atom(string const & _label,double _multiplier, double _occupancy, 
       double r1, double r2, double r3,
       double Uiso,
       sym_mat3<double> reciprocal_metric_tensor,
       ScatteringType scattering_type = XRay) :
  multiplier(_multiplier), occupancy(_occupancy),label(_label),r(r1,r2,r3), U(Uiso*reciprocal_metric_tensor)
  {
    atomic_type=AtomicTypeCollection::get(_label, scattering_type);
  }

  /// Function only used for tests.
	Atom(string const & _label, double _occupancy, 
		 double r1, double r2, double r3, //position
		 double U11, double U22, double U33, double U12, double U13, double U23) : //ADP
	occupancy(_occupancy), r(r1,r2,r3), U(U11,U22,U33,U12,U13,U23),label(_label), multiplier(1)
	{  
		atomic_type=AtomicTypeCollection::get(_label, XRay);
	}
  
  /// Atom initialization. Be careful, even though they are called Uij, the input parameters are supposed to be Uij/aiaj.
  Atom(string const & _label, double _multiplier, double _occupancy, 
       double r1, double r2, double r3, //position
       double U11, double U22, double U33, double U12, double U13, double U23,
       ScatteringType scattering_type = XRay) : //ADP
	multiplier(_multiplier), occupancy(_occupancy), r(r1,r2,r3), U(U11,U22,U33,U12,U13,U23),label(_label)
  {
		atomic_type=AtomicTypeCollection::get(_label, scattering_type);
	}
	
	double get_occupancy(){
    return occupancy;
  }
	void set_occupancy(double _occupancy){
    occupancy = _occupancy;
  }
	
	vector<Atom*> get_atoms() {
		vector<Atom*> v;
		v.push_back(this);
		return v;
	}
	
	Atom* create_symmetric(mat3<double> transformation_matrix,vec3<double> translation)
	{
		Atom* result = new Atom( *this );
		
		result->r=transformation_matrix*r+translation;
		result->U=trusted_mat_to_sym_mat(transformation_matrix*U*transformation_matrix.transpose());

		return result;
	}
	
	
	ChemicalUnit* operator[](int n)
	{
		throw "atom is a leaf node";
	}
  
  bool operator==(const Atom & inp)
  {
    return almost_equal(occupancy,inp.occupancy) && almost_equal(r,inp.r) && almost_equal(U,inp.U) && atomic_type==inp.atomic_type;
  }
};

class AtomicPair {
public:
	AtomicPair()
	{}
	
	//Most useful constructor which makes atomic pair form 2 atoms
	AtomicPair(Atom& _atom1, Atom& _atom2) :
	real(_atom1.occupancy * _atom2.occupancy,_atom2.r-_atom1.r, _atom1.U+_atom2.U),
    average(_atom1.occupancy * _atom2.occupancy,_atom2.r-_atom1.r, _atom1.U+_atom2.U),
	atomic_type1(_atom1.atomic_type),
	atomic_type2(_atom2.atomic_type),
	multiplier(_atom1.multiplier * _atom2.multiplier),
	atom1(&_atom1),
	atom2(&_atom2)
	{}
	
	//next three accessors return parameters of either average or real atomic pair
	//depending on the value of average_flag
	double& p(bool average_flag=false)
	{
		return params(average_flag).occupancy;
	}		
	vec3<double>& r(bool average_flag=false)
	{
		return params(average_flag).r;
	}
	sym_mat3<double>& U(bool average_flag=false)
	{
		return params(average_flag).U;
	}		
	
	double& average_p()
	{
		return p(true);
	}	
	double& real_p()
	{
		return p(false);
	}			
	
	vec3<double>& average_r()
	{ return r(true); }
	
	vec3<double>& real_r()
	{ return r(false); }
	
	sym_mat3<double>& average_U()
	{ return U(true); }
	
	sym_mat3<double>& real_U()
	{ return U(false); }
	
	
	Scatterer* atomic_type1;
	Scatterer* atomic_type2;
	double multiplier; //< Takes care of pair multiplicity due to symmetry. Supposed to be equal to 1/multiplicity.
	
	Atom *atom1, *atom2;
  
  bool operator==(const AtomicPair& inp)
  {
    return average==inp.average & real==inp.real & atom1==inp.atom1 & atom2==inp.atom2 & almost_equal(multiplier, inp.multiplier);
  }
  
  bool pair_is_withing(Grid grid)
  {
    for(int i=0;i<3; ++i)
      if(abs(average_r()[i])>abs(grid.lower_limits[i]))
        return false;
    return true;
    
  }
  
  string to_string()
  {
    std::ostringstream oss;
    oss << atom1->label << ' ' << atom2->label << ' ' << multiplier << ' ' << p() << ' ' << r()[0] << ' ' << r()[1] << ' ' << r()[2];
    
    for(int j=0; j<6; ++j)
      oss << ' ' << U()[j];
    
    oss << ' ' << average_p() << ' ' << average_r()[0] << ' ' << average_r()[1] << ' ' << average_r()[2];
    
    for(int j=0; j<6; ++j)
      oss << ' ' << average_U()[j];
    
    return oss.str();
  }
private:
	AtomicParams& params(bool average_flag)
	{
		if(average_flag)
			return average;
		else
			return real;
	}
	AtomicParams real; //parameters of a pair in disordered and
	AtomicParams average; //in average structure
  
  
};

/**
 * Adds values of small_piece array to corresponding vectors of accumulator.
 * \param r vector from lower left corner of accumulator to center of small_piece matrix, center assumed to be in point size/2 + 1
 */
void add_pair_to_appropriate_place(IntensityMap & small_piece,IntensityMap & accumulator,vec3<int> r,vector<bool> periodic);

/**
 * Prints vector in standart output.
 * For test reasons since elements of vec3 are not visible in debugger
 */
template<class NumType>
void print_vector(vec3<NumType> v)
{
	cout << v[0] << " " << v[1] << " " << v[2] << endl;
}

class IntnsityCalculator {
public:
	static void calculate_patterson_map_from_pairs_f(vector<AtomicPair>pairs,IntensityMap& patterson_map,bool average_flag,vec3<int> pair_grid_size,vector<bool> periodic_directions = vector<bool>(3,false))
	{
		double d_star_square;
    complex<double> f1,f2;
		vec3<double> s,r_res;
		vec3<int> r_grid;

    double scale = patterson_map.size_1d(); // we need to multiply the signal by this, to preserve scaling
    
    for(int i=0; i<patterson_map.size_1d(); ++i)
      patterson_map.at_c(i)=0;
    
		//std::vector<AtomicPair>::iterator pair;
		AtomicPair* pair;
		
		//cout << "pmap "<< patterson_map.grid_steps()[0] << patterson_map.grid_steps()[1] << " " << patterson_map.grid_steps()[2] << endl;
		Grid grid_for_pairs_p(patterson_map.unit_cell(),
                          patterson_map.grid_steps(),
                          patterson_map.grid_steps().each_mul(-pair_grid_size/2),//Center is in position pair_grid_size/2+1
                          patterson_map.grid.reciprocal_flag); // going to be false for direct space, i believe
    
		Grid grid_for_pairs_r=grid_for_pairs_p.reciprocal(); // grid in reciprocal (diffse scattering) space
		
		AtomicTypeCollection::calculate_form_factors_of_all_atoms_on_grid(pair_grid_size,grid_for_pairs_r);
omp_set_num_threads(8);
#pragma omp parallel for private(s,r_res,r_grid,d_star_square,f1,f2,pair)
    //shared(patterson_map)
		for (int i=0; i< pairs.size(); ++i)
		{
			pair=&pairs[i];
			IntensityMap pair_patterson_map(pair_grid_size);
			pair_patterson_map.set_grid(grid_for_pairs_r);
			
			grid_and_residual(pair->average_r(),patterson_map.grid,r_grid,r_res);
			
			if (!average_flag)
				r_res+=pair->r()-pair->average_r();
			
			pair_patterson_map.init_iterator();
			while(pair_patterson_map.next())
			{
				s=pair_patterson_map.current_s();
				d_star_square=pair_patterson_map.current_d_star_square();
				
				f1=pair->atomic_type1->gridded_form_factors.at_c(pair_patterson_map.current_index());
				f2=pair->atomic_type2->gridded_form_factors.at_c(pair_patterson_map.current_index());
				
				pair_patterson_map.current_array_value_c()=scale*calculate_scattering_from_a_pair_in_a_point_c(f1,f2,
																										 pair->p(average_flag),
																										 pair->multiplier,
																										 s,
																										 r_res,
																										 pair->U(average_flag));
			}
			
			pair_patterson_map.invert();
			
#pragma omp critical
			add_pair_to_appropriate_place(pair_patterson_map,patterson_map,r_grid,periodic_directions);
		}

	}
	
	static void calculate_scattering_from_pairs(vector<AtomicPair>pairs,IntensityMap& I,bool average_flag)
	{
		double d_star_square;
    complex<double> f1,f2;
		vec3<double> s;
		AtomicPair* pair;
		
		I.init_iterator();
		while(I.next())
		{
			s=I.current_s();
			d_star_square=I.current_d_star_square();
			
			AtomicTypeCollection::update_current_form_factors(s,d_star_square);
			int sz=pairs.size();
			double Intensity_in_current_point = 0;
//#pragma omp parallel for default(shared) private(f1,f2,pair) shared(d_star_square,s,average_flag) reduction(+: Intensity_in_current_point)
			for ( int i=0; i<sz; i++ )
			{
				pair = &pairs[i];
				f1=pair->atomic_type1->current_form_factor;
				f2=pair->atomic_type2->current_form_factor;
				Intensity_in_current_point+=real(calculate_scattering_from_a_pair_in_a_point_c(
																					 f1,
																					 f2,
																					 pair->p(average_flag),
																					 pair->multiplier,
																					 s,
																					 pair->r(average_flag),
																					 pair->U(average_flag)));
			}
			I.current_array_value() = Intensity_in_current_point;
		}
	}
	
	inline static complex<double> calculate_scattering_from_a_pair_in_a_point_c(complex<double> const & f1, complex<double> const & f2, //form factors at the point
																				double const & p, //joint probability
																				double const & N,  //multiplier
																				vec3<double> const & s, //vector of reciprocal coordinates of the point
																				vec3<double> const & r, //vector of distance between atoms in a pair
																				sym_mat3<double> const & U) //tensor of ADP of prevoius vector
	{
		return conj(f1)*f2*p*N*exp( complex<double>(M2PISQ*(s*U*s),M_2PI*(s*r)) );
	}
	
	

	class UnknownSymmetry {};
};



class AtomicPairPool {
public:
	
	AtomicPairPool()
	{}
	
	AtomicPair& get_pair(Atom * atom1,Atom * atom2)
	{
		for(vector<AtomicPair>::iterator pair = pairs.begin(); pair != pairs.end(); pair++)
		{
			if(pair->atom1 == atom1 && pair->atom2 == atom2)
				return *pair;
		}
		
		pairs.push_back(AtomicPair(*atom1,*atom2));
		return pairs.back();
	}
  void add_modifiers(vector<SubstitutionalCorrelation*> _modifiers)
  {
    vector<SubstitutionalCorrelation*>::iterator it;
    for(it=_modifiers.begin(); it!=_modifiers.end(); it++)
      modifiers.push_back((PairModifier*)(*it));
  }
	void add_modifier(PairModifier* modifier)
	{
		modifiers.push_back(modifier);
	}
	
/*	void add_modifiers(vector<PairModifier*> _modifiers)
	{
		modifiers.concat(_modifiers);
	}*/
	
	void invoke_correlators()
	{
		vector<PairModifier*> non_generating_modifiers;
		
		//int s = modifiers.size();
		
		//Invoke all generating modifiers
		for(int i=0; i<modifiers.size(); i++)
		{
			if(modifiers[i].generates_pairs())
			{
				modifiers[i].modify_pairs(this);
			} else {
				non_generating_modifiers.push_back(&modifiers[i]);
			}
		}
		
		vector<PairModifier*>::iterator it;
		for(it=non_generating_modifiers.begin(); it!=non_generating_modifiers.end(); it++)
		{
			(*it)->modify_pairs(this);
		}
		
	}
	

  
	
	vector<AtomicPair> pairs;
	
//private:
	p_vector<PairModifier> modifiers;
};

class CellShifter : public PairModifier
	{
	public:
		CellShifter(double r1, double r2, double r3) : 
		shift(r1,r2,r3)
		{}
		
		bool generates_pairs() 
		{	return false; }
		
		void modify_pairs(AtomicPairPool* const pool)
		{
			for(vector<AtomicPair>::iterator pair = pool->pairs.begin(); pair != pool->pairs.end(); pair++)
			{
				pair->r()+=shift;
				pair->average_r()+=shift;
			}
		}

		bool operator==(const CellShifter& inp)
    {
      return almost_equal(shift,inp.shift);
    }
		vec3<double> shift;
	};

class SubstitutionalCorrelation : public PairModifier{
public:
  SubstitutionalCorrelation(ChemicalUnit* unit1, ChemicalUnit* unit2, double _joint_probability) : 
  joint_probability(_joint_probability)
  {
    chemical_units[0] = unit1;
    chemical_units[1] = unit2;
  }
  
  bool generates_pairs()
  {	return true;	}
  
  void modify_pairs(AtomicPairPool* const pool)
  {
    vector<Atom*> atoms1=chemical_units[0]->get_atoms();
    vector<Atom*> atoms2=chemical_units[1]->get_atoms();
    
    vector<Atom*>::iterator atom1, atom2;
    for(atom1=atoms1.begin(); atom1!=atoms1.end(); atom1++)
    {
      for(atom2=atoms2.begin(); atom2!=atoms2.end(); atom2++)
      {
        pool->get_pair(*atom1,*atom2).p()=joint_probability;
      }
    }
  }
  
  /**
   * Comparison for test purposes. Now only works if the object references to the same chemical units.
   * does not work when the order of chemical units is different, which is correct behavior in current program implementation (because symmetry is not yet attached to chemical untis)
   * \TODO implement comparison when chemical units are not identical but equivalent once the comparison of chemical units is implemented
   */
  bool operator==(const SubstitutionalCorrelation& inp)
  {
    return almost_equal(joint_probability,inp.joint_probability) && chemical_units[0]==inp.chemical_units[0] && chemical_units[1]==inp.chemical_units[1];
  }
  
  
  ChemicalUnit* chemical_units[2];
  double joint_probability;
};

class AtomicDisplacement{
 public:
	Atom* atom;
	vec3<double> displacement_vector;
  
  bool operator==(const AtomicDisplacement& inp)
  {
    return atom==inp.atom && almost_equal(displacement_vector, inp.displacement_vector);
  }
};

class ADPMode
	{
	public:
		
		ADPMode() 
		{}
		
		void add_atom(Atom* _atom, vec3<double> r)
		{
			AtomicDisplacement t;
			t.atom = _atom;
			t.displacement_vector = r;
			atomic_displacements.push_back(t);
		}
		
		vector<AtomicDisplacement> atomic_displacements;
    
    /**
     * Operator for test functions.
     */
    bool operator==(const ADPMode& inp)
    {
      
      if(inp.atomic_displacements.size()!=atomic_displacements.size())
        return false;

      for(int i=0; i<atomic_displacements.size(); i++)
        if(!(atomic_displacements[i]==inp.atomic_displacements[i]))
          return false;
      
      return true;
    }
    
    /**
     * Operator for test functions
     */
    bool operator!=(const ADPMode& inp)
    {
      return !operator==(inp);
    }
	};

class DoubleADPMode : public PairModifier{
 public:
  DoubleADPMode(ADPMode* mode1,ADPMode* mode2,double _amplitude)
  :amplitude(_amplitude)
  {
    modes[0]=mode1;
    modes[1]=mode2;
  }
  
  bool generates_pairs()
  {	return true;	}
  
  void modify_pairs(AtomicPairPool* const pool)
  {
    vector<AtomicDisplacement>::iterator disp1,disp2;
    
    for(disp1=modes[0]->atomic_displacements.begin(); disp1!=modes[0]->atomic_displacements.end(); disp1++)
    {
      for(disp2=modes[1]->atomic_displacements.begin(); disp2!=modes[1]->atomic_displacements.end(); disp2++)
      {
        sym_mat3<double> U=pool->get_pair(disp1->atom,disp2->atom).U();
        pool->get_pair(disp1->atom,disp2->atom).U()+= -amplitude*(outer_product(disp1->displacement_vector,disp2->displacement_vector)+outer_product(disp2->displacement_vector,disp1->displacement_vector));
        // TODO: CHECK WHY THERE IS factor 2 here. I think it is WRONG!. Check some special case with perfect correlation between x and y and see if matrices will become ill-defined
        if(false)
        {
          sym_mat3<double> U1=pool->get_pair(disp1->atom,disp2->atom).U();
          
          sym_mat3<double> O = outer_product(disp1->displacement_vector,disp2->displacement_vector);
          vec3<double> r = pool->get_pair(disp1->atom,disp2->atom).r();
          for(int i=0; i<6; i++)
            cout << " " << U[i];
          cout << endl;
          
          for(int i=0; i<6; i++)
            cout << " " << O[i];
          cout << endl;
          
          for(int i=0; i<3; i++)
            cout << " " << r[i];
          cout << endl;
          
          for(int i=0; i<6; i++)
            cout << " " << U1[i];
          cout << endl;
          cout << "go\n";
        }
      }
    }
  }

  bool operator==(const DoubleADPMode& inp)
  {
    return inp.modes[0]==modes[0] && inp.modes[1]==modes[1] && almost_equal(amplitude, inp.amplitude);
  }
private:
  ADPMode* modes[2];
  double amplitude;
};

template<class T>
vector<T> unique_elements(vector<T> inp)
{
	vector<T> res;
	sort (inp.begin(),inp.end());
	unique_copy( inp.begin(), inp.end(), back_inserter(res) );
	return res;
}

class StaticShift : public PairModifier
	{
	public:
		StaticShift()
		{
		}
		
		bool generates_pairs()
		{	return true;	}
		
		struct ModeAndAmplitude
		{
			ADPMode* mode;
			double amplitude;
		};
		
		void add_displacement(const int& StartOrEnd, ADPMode* mode,double amplitude)
		{
			if (StartOrEnd!=0 && StartOrEnd!=1)
				throw string("StartOrEnd flag should be 0 for start or 1 for end");
			
			ModeAndAmplitude ma={mode,amplitude};
			modes_and_amplitudes[StartOrEnd].push_back(ma);
		}
		
		void modify_pairs(AtomicPairPool* const pool)
		{
			vector<Atom*> atoms[2];
			vector<vec3<double> > displacements[2];
			
			vector<ModeAndAmplitude>::iterator mode_and_amplitude;
			vector<AtomicDisplacement>::iterator at_disp;
			
			for(int st_end=0; st_end<2; st_end++)
				for(mode_and_amplitude=modes_and_amplitudes[st_end].begin(); mode_and_amplitude!=modes_and_amplitudes[st_end].end(); mode_and_amplitude++)
					for(at_disp = mode_and_amplitude->mode->atomic_displacements.begin(); at_disp!=mode_and_amplitude->mode->atomic_displacements.end(); at_disp++)
					{
						atoms[st_end].push_back(at_disp->atom);
						displacements[st_end].push_back(at_disp->displacement_vector*mode_and_amplitude->amplitude);				
					}
			
			vector<Atom*>::iterator atom1, atom2;
			vector<vec3<double> >::iterator displacement1,displacement2;
			vector<Atom*> unique_atoms=unique_elements(atoms[0]);
			
			for(atom1=unique_atoms.begin(); atom1!=unique_atoms.end(); atom1++)
			{
				for(atom2=atoms[1].begin(), displacement2=displacements[1].begin(); atom2!=atoms[1].end(); atom2++,displacement2++)
				{
					pool->get_pair(*atom1,*atom2).r()+=*displacement2;
				}
			}
			
			unique_atoms=unique_elements(atoms[1]);
			
			for(atom1=atoms[0].begin(), displacement1=displacements[0].begin(); atom1!=atoms[0].end(); atom1++,displacement1++)
			{
				for(atom2=unique_atoms.begin(); atom2!=unique_atoms.end(); atom2++)
				{
					pool->get_pair(*atom1,*atom2).r()-=*displacement1;
				}
			}
		}
		
		
		vector<ModeAndAmplitude> modes_and_amplitudes[2]; //< shift modes of start [0] and end[1] PDF vectors
	};

class SizeEffect : public PairModifier {
public:
  bool generates_pairs()
  {	return true;	}
  
  SizeEffect(ADPMode* _mode, ChemicalUnit* _cu, double _amplitude) : 
  cu(_cu),mode(_mode),amplitude(_amplitude),cu_to_adp_mode(false)
  { };
  SizeEffect(ChemicalUnit* _cu,ADPMode* _mode, double _amplitude) : 
  cu(_cu),mode(_mode),amplitude(_amplitude),cu_to_adp_mode(true)
  { };

  void modify_pairs(AtomicPairPool* const pool)
  {
    //Implement
    vector<Atom*> atoms = cu->get_atoms();
    vector<Atom*>::iterator atom;
    vector<AtomicDisplacement>::iterator atomic_displacement;
    
    for(atom = atoms.begin(); atom!=atoms.end(); atom++)
      for(atomic_displacement = mode->atomic_displacements.begin(); atomic_displacement != mode->atomic_displacements.end(); ++atomic_displacement)
        if(cu_to_adp_mode)
          pool->get_pair(*atom,atomic_displacement->atom).r()+=atomic_displacement->displacement_vector*amplitude;
        else
          pool->get_pair(atomic_displacement->atom,*atom).r()-=atomic_displacement->displacement_vector*amplitude;
  }
  
  bool operator==(const SizeEffect& inp)
  {
    return cu==inp.cu && mode==inp.mode && almost_equal(amplitude,inp.amplitude) && cu_to_adp_mode==inp.cu_to_adp_mode;
  }
  
  bool cu_to_adp_mode; ///<Boolean which is true if the modifier is atteched to chemical unit at beginnig of the vector and to adp mode at the end
  ChemicalUnit* cu;
  ADPMode* mode;
  double amplitude;
};

class ZeroVectorCorrelation : public PairModifier
	{
	public:
		ZeroVectorCorrelation() {}
		
		bool generates_pairs()
		{	return false;	}
		
		void modify_pairs(AtomicPairPool* const pool)
		{
			for(vector<AtomicPair>::iterator pair = pool->pairs.begin(); pair != pool->pairs.end(); pair++)
			{
				
				if((pair->r().length()<0.0001)&&(pair->atom1 == pair->atom2))
				{
					pair->U()=sym_mat3<double>(0,0,0,0,0,0);
					pair->p()=pair->atom1->occupancy;
				}
				if((pair->r().length()<0.0001)&&(pair->atom1 != pair->atom2))
				{
					//cout << "missed p = 0" << endl;
					pair->p()=0;
				}
				//cout << "probabilities are " << pair->p() << " " << pair->average_p() << endl;
				//cout << "adps are " << pair->U()[0] << " " << pair->average_U()[0] << endl;
			}
		}
	};

class MultiplicityCorrelation : public PairModifier{
public:
  MultiplicityCorrelation(double _multiplier) :
  multiplier(_multiplier) {}
  
  bool generates_pairs()
  {	return false;	}
  
  void modify_pairs(AtomicPairPool* const pool)
  {
    for(vector<AtomicPair>::iterator pair = pool->pairs.begin(); pair != pool->pairs.end(); pair++)
    {
      pair->multiplier*=multiplier;
    }
  }
  
  
  bool operator==(const MultiplicityCorrelation& inp)
  {
    return almost_equal(multiplier, inp.multiplier);
  }
//private:
  double multiplier;
};

/** 
 * An interface to minima finding routine.
 */
class MinimizerCalculator{
public:
  /**
   * This function shoul fill in the calculated diffuse scattering or PDF map.
   * \param p - is a raw vector of model parameters
   */
  virtual void calculate(vector<double> p) = 0;
  
  /// Accessor to calculated map
  virtual IntensityMap& data() = 0;
};

/**
 * Structure holds parameters to run minimization. For info see levmar documentation
 */
struct RefinementOptions {
  int max_number_of_iterations;
  double tau;
  double thresholds [3];
  double difference;
  
  static RefinementOptions default_refinement_options() {
	  RefinementOptions result = { 1000, //max_number_of_iterations
      1E-03, //tau
      {1E-17,1E-17,1E-17}, //thresholds
      1E-06}; //difference

    return result;

  }
};

/**
 * This class wraps levmar library for least-square solution of a problem min(I_model(params)-I_experimental)
 */
class Minimizer {
 public:
  /**
   * This function makes a setup for minimizer.
   * \todo add initialization of minima finding parameters, like number of iterations and so on
   * \todo make a separate file for the minimizer, or move it to core file
   */
  Minimizer() { 
    covar = NULL; 
  } 
  
  /**
   * Solves the problem of finding parameters which minimize I_model(params)-I_experimental in a least-square sence.
   * \param _calc - a reference to an object that calculates model diffuse scattering (or PDF). The object should implement MinimizerCalculator interface
   */
  vector<double> minimize(const vector<double> initial_params,IntensityMap * _experimental_data, MinimizerCalculator * _calc, OptionalIntensityMap * _weights, RefinementOptions refinement_options = RefinementOptions::default_refinement_options())
  {
    calc = _calc;
    experimental_data = _experimental_data;
    weights = _weights;
    
    double * p = (double*) malloc(sizeof(double)*initial_params.size());
    for(int i=0; i<initial_params.size(); i++)
      p[i]=initial_params[i];
    
    double * x = (double*) malloc(sizeof(double)*experimental_data->size_1d());
    for(int i=0; i<experimental_data->size_1d(); i++)
      x[i]=0;
    
    covar = (double*) malloc(sizeof(double)*initial_params.size()*initial_params.size());

    double opts[5];
    opts[0] = refinement_options.tau;
    for(int i=0; i<3; ++i)
      opts[i+1]=refinement_options.thresholds[i];
    opts[4] = refinement_options.difference;
    
    double info[10];
    
    int ret=dlevmar_dif(func_for_levmar, p, x, initial_params.size(),experimental_data->size_1d(), refinement_options.max_number_of_iterations, opts, info, NULL, covar, this);


   
    if (ret<0)
      REPORT(ERROR) << "Error, least square solution failed\n";

    switch((int)info[6]) {
      case 1: REPORT(MAIN) << "Least squares solution ended normally because of small gradient\n"; break;
      case 2: REPORT(MAIN) << "Least squares solution ended normally because of convergence\n"; break;
      case 3: REPORT(MAIN) << "Least squares solution stopped because maximum iteraction count is exceeded\n"; break;
      case 4: REPORT(ERROR) << "Warning, least square solution stopped because of singular matrix\n"; break;
      case 5: REPORT(MAIN) << "Least squares solution ended normally because no further error reduction is possible.\n"; break;
      case 6: REPORT(MAIN) << "Least squares solution ended normally because experiment and model are almost the same\n"; break;
      case 7: REPORT(ERROR) << "Error, least square solution aborted. Found invalid (i.e. NaN or Inf) values in diffuse scattering\n"; break;
    }
    
    static 
    
    vector<double> result(p,p+initial_params.size());
    
    delete p;
    delete x;
    return result;
  }
  
  ~Minimizer() {
    free(covar);
  }
  
  double * covar;
 private:
  /**
   * the function which is called by levmar. The function has the structure that lev mar supports. the field data is used to store reference to the minimizer object, 
   * so even though this function is technically static, it works like a normal member function
   */
  static void func_for_levmar(double *p, double *x, int parameters_number, int datapoints_number, void *data)
  {
    Minimizer * _this = reinterpret_cast<Minimizer*> (data);
    
    vector<double> parameters(p,p+parameters_number);
    
    _this->calc->calculate(parameters);

    //TODO: save space. Here, use in_asu

    //copy difference to the *x array
    for(int i=0; i<datapoints_number; i++)
      x[i] = (_this->experimental_data->at(i) - _this->calc->data().at(i))*_this->weights->at(i);
  }
  
  MinimizerCalculator * calc; //< model diffuse scattering of PDF calculator.
  IntensityMap * experimental_data;
  OptionalIntensityMap * weights;
  
};

class SymmetryElement{
 public:
  SymmetryElement(mat3<double> _permutation_matrix,vec3<double> _displacement) : 
    permutation_matrix(_permutation_matrix),displacement(_displacement) { }
  
  SymmetryElement() : 
    permutation_matrix(mat3<double>(1,0,0,0,1,0,0,0,1)),displacement(vec3<double>(0,0,0)) { }
    
  bool operator==(SymmetryElement inp)
  {
    return almost_equal(displacement,inp.displacement) && almost_equal(permutation_matrix,inp.permutation_matrix);
  }
  
  vec3<double> displacement;
  mat3<double> permutation_matrix;
    
};

class Error{
public:
  Error(string inp) : message(inp) { }
  string message;
};


#endif