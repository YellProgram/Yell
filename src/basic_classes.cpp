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


#include "basic_classes.h"
#include "diffuser_core.h"
#include <math.h>

double vector_norm(vec3<double> v, sym_mat3<double> M)
{
  return sqrt(v*M*v);
}

sym_mat3<double> outer_product(vec3<double> v1, vec3<double> v2) {
	return sym_mat3<double>(v1[0]*v2[0],v1[1]*v2[1],v1[2]*v2[2],v1[0]*v2[1],v1[0]*v2[2],v1[1]*v2[2]);
}

ADPMode* translational_mode(ChemicalUnit* unit,int direction,sym_mat3<double> metrical_matrix) {
	vector<Atom*> atoms = unit->get_atoms();
	vector<Atom*>::iterator atom;
	
	vec3<double> r(0,0,0);
	r[direction]=1;
  
  r=r/vector_norm(r,metrical_matrix);
	
	ADPMode* mode = new ADPMode();
	
	for(atom=atoms.begin(); atom!=atoms.end(); atom++)
		mode->add_atom(*atom,r);
	
	return mode;
}

ADPMode z_rot_mode(ChemicalUnit* unit) //only for hex system
{
	vector<Atom*> atoms = unit->get_atoms();
	vector<Atom*>::iterator atom;

	ADPMode mode;
	
	mat3<double> rotation(1/sqrt(3.0),-2/sqrt(3.),0,2/sqrt(3.),-1/sqrt(3.),0,0,0,0);
	
	for(atom=atoms.begin(); atom!=atoms.end(); atom++)
		mode.add_atom(*atom,rotation*((*atom)->r));
	
	return mode;
}

ADPMode* rot_mode(ChemicalUnit* unit,vec3<double> axis , vec3<double> point_on_axis ,sym_mat3<double> metrical_matrix) {
	vector<Atom*> atoms = unit->get_atoms();
	vector<Atom*>::iterator atom;
  
  axis = axis/vector_norm(axis, metrical_matrix.inverse());
  
	ADPMode* mode = new ADPMode();
	
	/*
   mat3<double> rotation(1/sqrt(3),-2/sqrt(3),0,2/sqrt(3),-1/sqrt(3),0,0,0,0);
	*/
	
	for(atom=atoms.begin(); atom!=atoms.end(); atom++)    
		mode->add_atom(*atom,metrical_matrix*(axis.cross((*atom)->r - point_on_axis))/sqrt(metrical_matrix.determinant()));
	
	return mode;
}
ADPMode combine_modes(vector<ADPMode> modes) {
	ADPMode result;
	vector<ADPMode>::iterator mode;
	for(mode=modes.begin(); mode!=modes.end(); mode++)
		result.atomic_displacements.insert(result.atomic_displacements.end(),mode->atomic_displacements.begin(),mode->atomic_displacements.end());
	return result;
}


vector<SubstitutionalCorrelation*> correlators_from_cuns(ChemicalUnitNode* node1,ChemicalUnitNode* node2,vector<double> corr) {
	int s1 = node1->chemical_units.size();
	int s2 = node2->chemical_units.size();
	
	vector<double> correlations(s1*s2);
	
  //check there are enough numbers in the matrix
  //fill the correlations
  //if full version, check the last column-row
  
  bool full_version;
  if(corr.size()==(s1-1)*(s2-1)) {
    full_version = false;
    
    //fill upper-right corner of correlation matrix
    for(int i=0; i<s1-1; i++)
      for(int j=0; j<s2-1; j++)
        correlations[i+s1*j] = corr[i+j*(s1-1)];
    
  } else if(corr.size()==s1*s2) {
    full_version = true;
    correlations = corr;
  } else {
    REPORT(ERROR) << "Wrong number of joint probabilities, expected a " << (s1-1)<< "x" << (s2-1) << " or a " << s1 << "x" << s2 << " matrix. Found " << corr.size() << " values.\n";
    throw "fail";
  }
  
	//fill last row
	for(int j=0; j<s2-1; j++){
    correlations[s1*j+s1-1]=node2->chemical_units[j].get_occupancy();
		for(int i=0; i<s1-1; i++)
			correlations[s1*j+s1-1]-=correlations[i+s1*j];

	}
	
	//fill last coloumn
	for(int i=0; i<s1; i++) {
		correlations[(s2-1)*s1+i]=node1->chemical_units[i].get_occupancy();
    
		for(int j=0; j<s2-1; j++)
			correlations[(s2-1)*s1+i] -= correlations[i+s1*j];
	}
	
  if(full_version)
    for(int i=0; i<correlations.size(); ++i)
      if(!(almost_equal(corr[i],correlations[i]))){
        REPORT(ERROR) << "The joint probabilities are incompatible with the average structure.\n";
        throw "fail";
      }
  
	//create correlations
	vector<SubstitutionalCorrelation*> res;
	
	for(int j=0; j<s2; j++)
	{
		for(int i=0; i<s1; i++)
		{
			res.push_back( new SubstitutionalCorrelation(& node1->chemical_units[i],& node2->chemical_units[j], correlations[i+s1*j]));
		}
	}
	
	return res;
}



/// Skips check for matrix being symmetric. Just takes upper-triangular part.
sym_mat3<double> trusted_mat_to_sym_mat(mat3<double> inp) {
return sym_mat3<double>(inp[0],inp[4],inp[8],inp[1],inp[2],inp[5]);
}


// I had to move it here because it depends on AtomicPair class which is defined after LaueSymmetry
vector<AtomicPair> LaueSymmetry::multiply_pairs_by_matrix(vector<AtomicPair> pairs,mat3<double> transformation_matrix) {
  vector<AtomicPair>::iterator pair;
  for(pair=pairs.begin(); pair!=pairs.end(); pair++)
  {
    for(int average=0; average<2; average++)
    {
      pair->r(average)=transformation_matrix*pair->r(average);
      pair->U(average)=trusted_mat_to_sym_mat(transformation_matrix*pair->U(average)*transformation_matrix.transpose());
    }
  }
  
  return pairs;
}


void add_pair_to_appropriate_place(IntensityMap &  small_piece,IntensityMap & accumulator,vec3<int> r,vector<bool> periodic) {
	vec3<int> piece_size=small_piece.size();
	vec3<int> acc_size=accumulator.size();
	
	r=r-(piece_size/2);
	
	vec3<int> borders;
	for(int i=0; i<3; i++)
	{
		borders[i]=0;//(acc_size[i]>1);
	}

	vec3<int> llimits, ulimits;
  for(int i=0; i<3; ++i)
  {
    llimits[i]=r[i];
    ulimits[i]=r[i]+piece_size[i]+borders[i];

    
    if(!periodic.at(i))
    {
      llimits[i]=max(0,llimits[i]);
      ulimits[i]=min(acc_size[i],ulimits[i]);
    }
  }
	
	for(int i=llimits[0];i<ulimits[0]; i++)
	{
		int is=(i-r[0])%piece_size[0];
		for(int j=llimits[1];j<ulimits[1]; j++)
		{
			int js=(j-r[1])%piece_size[1];
			for(int k=llimits[2];k<ulimits[2]; k++)
			{
				int ks=(k-r[2])%piece_size[2];
				accumulator.at((i+acc_size[0])%acc_size[0],
                       (j+acc_size[1])%acc_size[1],
                       (k+acc_size[2])%acc_size[2])+=small_piece.at(is,js,ks);
			}
		}
	}
}

AtomicTypeCollection::~AtomicTypeCollection() {
	map<string,Scatterer*>::iterator atomic_type;
	for(atomic_type=types.begin(); atomic_type!=types.end(); atomic_type++)
		delete atomic_type->second;
}

AtomicTypeCollection& AtomicTypeCollection::get() {
	static AtomicTypeCollection collection;
	return collection;
}

void AtomicTypeCollection::calculate_form_factors_of_all_atoms_on_grid(vec3<int>map_size,Grid grid)
{
  map<string,Scatterer*>::iterator atomic_type;
  for(atomic_type=AtomicTypeCollection::get().types.begin(); atomic_type!=AtomicTypeCollection::get().types.end(); atomic_type++)
  {
    atomic_type->second->calculate_form_factors_on_grid(map_size,grid);
  }
}

void AtomicTypeCollection::update_current_form_factors(vec3<double>s, double d_star_square)
{
  map<string,Scatterer*>::iterator atomic_type;
  for(atomic_type=AtomicTypeCollection::get().types.begin(); atomic_type!=AtomicTypeCollection::get().types.end(); atomic_type++)
  {
    atomic_type->second->update_current_form_factor(s, d_star_square);
  }
}

string AtomicTypeCollection::strip_label(string const & label)
{
  string result;
  result.append(label,0,1);
  if (label[1]>=('A'))
    result.append(label,1,1);
  
  return result;
}

Scatterer* AtomicTypeCollection::get(string const & label, ScatteringType scattering_type)
{
  if(AtomicTypeCollection::get().types.find(label) != AtomicTypeCollection::get().types.end())
     return AtomicTypeCollection::get().types[label];
   else {
     string stripped_label = strip_label(label);
     if(AtomicTypeCollection::get().types.find(stripped_label) != AtomicTypeCollection::get().types.end())
       return AtomicTypeCollection::get().types[stripped_label];
     else
     {
       Scatterer *t;
         switch(scattering_type) {
             case XRay:
                 t = new AtomicType(stripped_label);
                 break;
             case Neutron:
                 t = new NeutronScattererAtom(stripped_label);
                 break;
             case Electrons:
                 t = new ElectronScattererAtom(stripped_label);
                 break;
         }
       AtomicTypeCollection::add(stripped_label, t);
       return t;
     }
   }
}

void AtomicTypeCollection::add(string const & label, Scatterer* s) {
  AtomicTypeCollection::get().types.insert(pair<string,Scatterer*>(label, s));
}

complex<double> MolecularScatterer::form_factor_at_c(vec3<double>s,double d_star_sq) {
  //This method can only be used by AtomicTypeCollection::update_current_form_factor function since it relies on the fact that
  //the form factor of all its underlying atoms is calculated prior to calling this function
  
  complex<double> result;
  vector<Atom*>::iterator it;
  Atom* at;
  for(it=constituent_atoms.begin(); it!=constituent_atoms.end(); ++it) {
    at=*it;
    result += form_factor_in_a_point(at->atomic_type->current_form_factor,
                                     at->occupancy,
                                     at->multiplier,
                                     at->r,
                                     at->U,
                                     s);
  }
  
  return result;
}

complex<double> MolecularScatterer::form_factor_in_a_point(complex<double> f,
                                                     double p,
                                                     double N,
                                                     vec3<double> r,
                                                     sym_mat3<double> U,
                                                     vec3<double> s)
{
  return f*p*N*exp( complex<double>(M2PISQ*(s*U*s),M_2PI*(s*r)) );
}

