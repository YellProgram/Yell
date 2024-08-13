//
// Created by Arkadiy Simonov on 09.08.24.
//

#include "LaueSymmetry.h"
#include "basic_classes.h"

#define Z it->average_r()[2]
#define Y it->average_r()[1]
#define X it->average_r()[0]
#define decrease_multiplicity(m) {result.back().p()/=m; result.back().average_p()/=m;}

///\TODO: move it to appropriate place. Probably to files LaueSymmetry.h and LaueSymmetry.cpp
///\TODO: check that it is in use, maybe delete?
vector<AtomicPair> LaueSymmetry::filter_pairs_from_asymmetric_unit_and_scale(vector<AtomicPair>& pairs)
{
    vector<AtomicPair> result;
    vector<AtomicPair>::iterator it;

    if(label=="6/m" || label=="4/m" || label=="2/m")
    {
        for(it=pairs.begin(); it!=pairs.end(); it++)
        {
            if( Z<0 )
                continue;

            result.push_back(*it);

            if(Z==0)
            decrease_multiplicity(2)
        }

        return result;
    }else if(label=="mmm" || label=="m-3")
    {
        for(it=pairs.begin(); it!=pairs.end(); it++)
        {
            if( Z<0 || X<0 || Y<0 )
                continue;

            result.push_back(*it);

            if(Z==0)
            decrease_multiplicity(2);
            if(X==0)
            decrease_multiplicity(2);
            if(Y==0)
            decrease_multiplicity(2);
        }
        return result;
    }
    else if(label=="4/mmm")
    {
        for(it=pairs.begin(); it!=pairs.end(); it++)
        {
            if(!( Z>=0 && X>=Y && Y>=0 ))
                continue;

            result.push_back(*it);

            if(Z==0)
            decrease_multiplicity(2);

            if(Y==0)
            decrease_multiplicity(2);

            if(X==Y)
            decrease_multiplicity(2);

            if(X==0)
            decrease_multiplicity(2);
        }
        return result;
    }
    else if(label=="m-3m")
    {

        for(it=pairs.begin(); it!=pairs.end(); it++)
        {
            if(!( Z>=-0 && Y >=Z && X>=Y ))
                continue;

            result.push_back(*it);

            if(X==0)// 0 0 0
            decrease_multiplicity(48)
            else if(Y==0 && Z==0) // x 0 0
            decrease_multiplicity(8)
            else if(X==Y && X==Z) // x x x
            decrease_multiplicity(6)

            else if(X==Y && Z==0) // x x 0
            decrease_multiplicity(4)
            else if(Z==0 ||Y==Z || X==Y) // x y 0, x y y, x x z
            decrease_multiplicity(2)
        }
        return result;
    }else if(label=="6/mmm")
    {
        for(it=pairs.begin(); it!=pairs.end(); it++)
        {
            if(!( Z>=0 && Y>=0 && X>=2*Y )) //if y<0 or z<0 or x<2y
                continue;

            result.push_back(*it);

            //DO separately xy and z parts
            if(X==0 && Y==0)// 0 0
            decrease_multiplicity(12)
            else if(Y==0 || X==2*Y) //x 0 and 2y y
            decrease_multiplicity(2);

            if(Z==0)// x y 0
            decrease_multiplicity(2);

        }
        return result;
    }else if(label=="-3mH")
    {
        for(it=pairs.begin(); it!=pairs.end(); it++)
        {
            if(!( X>=-Y && X>=2*Y ))
                continue;

            result.push_back(*it);


            if(X==0 && Y==0)// 0 0
            decrease_multiplicity(6)
            else if(X==-Y || X==2*Y)
            decrease_multiplicity(2);


        }
        return result;
    }else if(label=="-3mR")
    {
        for(it=pairs.begin(); it!=pairs.end(); it++)
        {
            if(!( X>=Y && Y>=Z ))
                continue;

            result.push_back(*it);


            if(X==Y && Y==Z)
            decrease_multiplicity(6)
            else if(X==Y || Y==Z)
            decrease_multiplicity(2);

        }
        return result;
    }
    else //-3R and -3H
        return pairs;

}
#undef Z
#undef X
#undef Y

vector<AtomicPair> LaueSymmetry::filter_pairs_from_asymmetric_unit(vector<AtomicPair>& pairs)
{
    vector<AtomicPair> result = filter_pairs_from_asymmetric_unit_and_scale(pairs);
    vector<AtomicPair>::iterator it;

    double group_multiplicity=0;

    if      (label=="m-3m")     group_multiplicity=48;
    else  if(label=="m-3")    	group_multiplicity=24;
    else  if(label=="6/mmm")    group_multiplicity=24;
    else  if(label=="6/m")    	group_multiplicity=12;
    else  if(label=="-3mH")    	group_multiplicity=12;
    else  if(label=="-3mR")    	group_multiplicity=12;
    else  if(label=="-3H")    	group_multiplicity=6;
    else  if(label=="-3R")    	group_multiplicity=6;
    else  if(label=="4/mmm")    group_multiplicity=16;
    else  if(label=="4/m")    	group_multiplicity=8;
    else  if(label=="mmm")    	group_multiplicity=8;
    else  if(label=="2/m")    	group_multiplicity=4;
    else  if(label=="2/mb")    	group_multiplicity=4/2; //mysterious 1/2 constant
    else  if(label=="-1")    	  group_multiplicity=2;

    assert(!(group_multiplicity==0));

    for(it=result.begin(); it!=result.end(); it++)
        it->multiplier*=group_multiplicity;

    return result;
}

bool LaueSymmetry::is_compatible_with_cell(const UnitCell& cell)
{
    double a=cell.cell.parameters()[0];
    double b=cell.cell.parameters()[1];
    double c=cell.cell.parameters()[2];
    double alpha=cell.cell.parameters()[3];
    double beta=cell.cell.parameters()[4];
    double gamma=cell.cell.parameters()[5];

    //took this from International Tables Vol A p 15
    if (label=="-1")   return true;
    if (label=="2/m")  return alpha==90 && beta==90;
    if (label=="2/mb") return alpha==90 && gamma==90;
    if (label=="mmm")  return alpha==90 && beta==90 && gamma==90;
    if (label=="4/m" || label=="4/mmm")return almost_equal(a, b) && alpha==90 && beta==90 && gamma==90;
    if (label=="-3H" || label=="-3mH") return almost_equal(a, b) && alpha==90 && beta==90 && gamma==120;
    if (label=="-3R" || label=="-3mR") return almost_equal(a, b) && almost_equal(a, c) && almost_equal(alpha, beta) && almost_equal(alpha,gamma);
    if (label=="6/m" || label=="6/mmm")return almost_equal(a, b) && alpha==90 && beta==90 && gamma==120;
    if (label=="m-3" || label=="m-3m") return almost_equal(a, b) && almost_equal(a, c) && alpha==90 && beta==90 && gamma==90;

    return false;
}

vector<AtomicPair> LaueSymmetry::apply_matrix_generator(vector<AtomicPair> pairs,mat3<double> transformation_matrix,int times) {
    vector<AtomicPair> result=pairs;

    for(int repeat=0; repeat<times-1; repeat++)
    {
        pairs = multiply_pairs_by_matrix(pairs,transformation_matrix);
        result.insert(result.end(),pairs.begin(),pairs.end());
    }

    for(vector<AtomicPair>::iterator it=result.begin(); it!=result.end(); ++it)
        it->multiplier/=times;

    return result;
}

vector<AtomicPair> LaueSymmetry::apply_patterson_symmetry(vector<AtomicPair> pairs)
{

    for(int i=0; i<generators_on_vectors.size(); ++i)
    {
        rt_mx generator(generators_on_vectors[i]);
        pairs = apply_matrix_generator(pairs,generator.r().as_double(),generator.r().order());
    }

    return pairs;
}

void LaueSymmetry::apply_generator(IntensityMap& map, string generator_symbol)
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
