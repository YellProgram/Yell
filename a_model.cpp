/*
 *  a_model.cpp
 *  diffuser_z
 *
 *  Created by Arkadiy Simonov on 8/12/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "a_model.h"
#include "basic_classes.h"

#include <scitbx/vec2.h>
#include <scitbx/sym_mat2.h>
#include <scitbx/mat2.h>

#include <fstream>
#include <iostream>


double two_times(double x)
{
	return x*2;
}

void multDoubleArray(double *x,int size)
{
    /* Multiple each element of the array by 3 */
    int i;
    for (i=0;i<size;i++)
	*x++ *= 3;
}

void calculate_scattering_size_effect_fcc(double *scattering,double* inp_params)
{       //size_effect_fcc
		//--------------------------------------------------------------------------------------------------------------------------
		double a=3.5885;
        UnitCell cell(a,a,a,90,90,90);
		
		bool AVERAGE_FLAG = inp_params[0]==1;
		
        int grid_size=40;
		IntensityMap I(grid_size,grid_size,grid_size);
        double ll=-8;
        double gs=0.4;
		I.set_grid(cell.cell,vec3<double>(gs,gs,gs),vec3<double>(ll,ll,ll));
		
		//-------------------------------------------------------------------------------------------------------------------------
        
        double U_Cu=4.7e-04;
        double U_Au=4.3e-04;
        double P_Cu = 0.7000;
		
                
                
		cell.set_laue_symmetry("*432");
			ChemicalUnitNode* node = new ChemicalUnitNode();
				
                Atom* Cu = new Atom("Cu",  P_Cu,0,0,0,U_Cu,U_Cu,U_Cu,0,0,0);
                Atom* Au = new Atom("Au",1-P_Cu,0,0,0,U_Au,U_Au,U_Au,0,0,0);
              
			node->add_chemical_unit(Cu);
            node->add_chemical_unit(Au);

		cell.add_node(node);
		
        //--------------------------------------------------------------------------------------------------------------------------
        
        const int X=0;
		const int Y=1;
        const int Z=2;
		//Create translational modes
		ADPMode cu_x = translational_mode(Cu,X);
		ADPMode cu_y = translational_mode(Cu,Y);
        ADPMode cu_z = translational_mode(Cu,Z);
        ADPMode au_x = translational_mode(Au,X);
		ADPMode au_y = translational_mode(Au,Y);
        ADPMode au_z = translational_mode(Au,Z);
		//--------------------------------------------------------------------------------------------------------------------------
		vector<AtomicPairPool*> pools;
        vector<AtomicPairPool*>::iterator pool;
        
        vector<double> c;
        c.push_back(P_Cu);

        AtomicPairPool zeroth_neighbor;
        zeroth_neighbor.add_modifiers( correlators_from_cuns(node,node,c) );
        zeroth_neighbor.add_modifier( new ZeroVectorCorrelation );
        zeroth_neighbor.add_modifier( new MultiplicityCorrelation(20.0*20) );  
        //pools.push_back(&zeroth_neighbor);
        
         double xx,xy,yx,yy,dr,d,zz,xz,zx,yz,zy,drx,dry,dx,dy;

        vector<ADPMode*> modes_x;
        modes_x.push_back(&cu_x);
        modes_x.push_back(&au_x);
        
        vector<ADPMode*> modes_y;
        modes_y.push_back(&cu_y);
        modes_y.push_back(&au_y);
        
        vector<ADPMode*> modes_z;
        modes_z.push_back(&cu_z);
        modes_z.push_back(&au_z);
        
        c[0]=0.49+inp_params[1];
        xx = inp_params[2];
        yy = xx;
        zz = inp_params[3];
        xy = inp_params[4];
        yx = xy;
        dr = inp_params[5];
		d = dr*0.6*0.5;

        
		AtomicPairPool first_neighbor;
        first_neighbor.add_modifiers( correlators_from_cuns(node,node,c) );
		first_neighbor.add_modifier( new CellShifter(0.5,0.5,0) );
        first_neighbor.add_modifier(new MultiplicityCorrelation(19.0*19));
                 
        StaticShift* shift = new StaticShift();
        shift->add_displacement(0,cu_x,d);
        shift->add_displacement(0,au_x,d-dr);
        shift->add_displacement(1,cu_x,-d);
        shift->add_displacement(1,au_x,dr-d);
        
        shift->add_displacement(0,cu_y,d);
        shift->add_displacement(0,au_y,d-dr);
        shift->add_displacement(1,cu_y,-d);
        shift->add_displacement(1,au_y,dr-d);
        
        first_neighbor.add_modifier(shift);
        
        for(int i=0; i<2; i++)
            for(int j=0; j<2; j++)
            {
                first_neighbor.add_modifier( new DoubleADPMode(modes_x[i],modes_x[j],xx));
                first_neighbor.add_modifier( new DoubleADPMode(modes_y[i],modes_y[j],yy));
                first_neighbor.add_modifier( new DoubleADPMode(modes_x[i],modes_y[j],xy));
                first_neighbor.add_modifier( new DoubleADPMode(modes_y[i],modes_x[j],yx));
                first_neighbor.add_modifier( new DoubleADPMode(modes_z[i],modes_z[j],zz));
            }    
        
		pools.push_back(&first_neighbor);
        
        c[0]=0.49+inp_params[6];
        xx = inp_params[7];
        yy = inp_params[8];
        zz = inp_params[8];
        dr = inp_params[9];
		d = dr*0.6*0.5;
        
        AtomicPairPool second_neighbor;
        second_neighbor.add_modifiers( correlators_from_cuns(node,node,c) );
		second_neighbor.add_modifier( new CellShifter(1,0,0) );
        second_neighbor.add_modifier(new MultiplicityCorrelation(20*18));
        
        shift = new StaticShift();
        shift->add_displacement(0,cu_x,d);
        shift->add_displacement(0,au_x,d-dr);
        shift->add_displacement(1,cu_x,-d);
        shift->add_displacement(1,au_x,dr-d);
        second_neighbor.add_modifier(shift);

        for(int i=0; i<2; i++)
            for(int j=0; j<2; j++)
            {
            second_neighbor.add_modifier( new DoubleADPMode(modes_x[i],modes_x[j],xx));
            second_neighbor.add_modifier( new DoubleADPMode(modes_y[i],modes_y[j],yy));
            second_neighbor.add_modifier( new DoubleADPMode(modes_z[i],modes_z[j],zz));
            }    
        
		pools.push_back(&second_neighbor);
        
        c[0]=0.49+inp_params[10];
        xx = inp_params[11];
        yy = xx;
        zz = inp_params[12];
        xy = inp_params[13];
        yx = xy;
        dr = inp_params[14];

        d = dr*0.6*0.5;
        
		AtomicPairPool third_neighbor;
        third_neighbor.add_modifiers( correlators_from_cuns(node,node,c) );
		third_neighbor.add_modifier( new CellShifter(1,1,0) );
        third_neighbor.add_modifier(new MultiplicityCorrelation(18.0*18));
                 
        shift = new StaticShift();
        shift->add_displacement(0,cu_x,d);
        shift->add_displacement(0,au_x,d-dr);
        shift->add_displacement(1,cu_x,-d);
        shift->add_displacement(1,au_x,dr-d);
        
        shift->add_displacement(0,cu_y,d);
        shift->add_displacement(0,au_y,d-dr);
        shift->add_displacement(1,cu_y,-d);
        shift->add_displacement(1,au_y,dr-d);
        
        third_neighbor.add_modifier(shift);
        
        for(int i=0; i<2; i++)
            for(int j=0; j<2; j++)
            {
                third_neighbor.add_modifier( new DoubleADPMode(modes_x[i],modes_x[j],xx));
                third_neighbor.add_modifier( new DoubleADPMode(modes_y[i],modes_y[j],yy));
                third_neighbor.add_modifier( new DoubleADPMode(modes_x[i],modes_y[j],xy));
                third_neighbor.add_modifier( new DoubleADPMode(modes_y[i],modes_x[j],yx));
                third_neighbor.add_modifier( new DoubleADPMode(modes_z[i],modes_z[j],zz));
            }    
        
		pools.push_back(&third_neighbor);
        
        c[0]=0.49+inp_params[15];
        xx = inp_params[16];
        xy = inp_params[17];
        
        xz = xy;
        yx = xy;
        zx = yx;
        
        yy = inp_params[18];
        zz = yy;
        yz = inp_params[19];
        zy = yz;

        drx = inp_params[20];
        dry = inp_params[21];
        
        dx = drx*0.6*0.5;
        dy = dry*0.6*0.5;
        
		AtomicPairPool fourth_neighbor;
        fourth_neighbor.add_modifiers( correlators_from_cuns(node,node,c) );
		fourth_neighbor.add_modifier( new CellShifter(1,0.5,0.5) );
        fourth_neighbor.add_modifier(new MultiplicityCorrelation(18.0*19*19/20));
                 
        shift = new StaticShift();
        shift->add_displacement(0,cu_x,dx);
        shift->add_displacement(0,au_x,dx-drx);
        shift->add_displacement(1,cu_x,-dx);
        shift->add_displacement(1,au_x,drx-dx);
        
        shift->add_displacement(0,cu_y,dy);
        shift->add_displacement(0,au_y,dy-dry);
        shift->add_displacement(1,cu_y,-dy);
        shift->add_displacement(1,au_y,dry-dy);
        
        shift->add_displacement(0,cu_z,dy);
        shift->add_displacement(0,au_z,dy-dry);
        shift->add_displacement(1,cu_z,-dy);
        shift->add_displacement(1,au_z,dry-dy);
        
        fourth_neighbor.add_modifier(shift);
        
        for(int i=0; i<2; i++)
            for(int j=0; j<2; j++)
            {
                fourth_neighbor.add_modifier( new DoubleADPMode(modes_x[i],modes_x[j],xx));
                fourth_neighbor.add_modifier( new DoubleADPMode(modes_x[i],modes_y[j],xy));
                fourth_neighbor.add_modifier( new DoubleADPMode(modes_x[i],modes_z[j],xz));
                fourth_neighbor.add_modifier( new DoubleADPMode(modes_y[i],modes_x[j],yx));
                fourth_neighbor.add_modifier( new DoubleADPMode(modes_y[i],modes_y[j],yy));
                fourth_neighbor.add_modifier( new DoubleADPMode(modes_y[i],modes_z[j],yz));
                fourth_neighbor.add_modifier( new DoubleADPMode(modes_z[i],modes_x[j],zx));
                fourth_neighbor.add_modifier( new DoubleADPMode(modes_z[i],modes_y[j],zy));
                fourth_neighbor.add_modifier( new DoubleADPMode(modes_z[i],modes_z[j],zz));
            }    
        
		pools.push_back(&fourth_neighbor);
        
        c[0]=0.49+inp_params[22];
        xx = inp_params[23];
        xy = inp_params[24];
        xz = xy;
        yx = xy;
        zx = xy;
        
        yy = xx;
        zz = xx;
        yz = xy;
        zy = xy;

        dr = inp_params[25];

        d = dr*0.6*0.5;
 
		AtomicPairPool fifth_neighbor;
        fifth_neighbor.add_modifiers( correlators_from_cuns(node,node,c) );
		fifth_neighbor.add_modifier( new CellShifter(1,1,1) );
        fifth_neighbor.add_modifier(new MultiplicityCorrelation(18.0*18*18/20));
                 
        shift = new StaticShift();
        shift->add_displacement(0,cu_x,d);
        shift->add_displacement(0,au_x,d-dr);
        shift->add_displacement(1,cu_x,-d);
        shift->add_displacement(1,au_x,dr-d);
        
        shift->add_displacement(0,cu_y,d);
        shift->add_displacement(0,au_y,d-dr);
        shift->add_displacement(1,cu_y,-d);
        shift->add_displacement(1,au_y,dr-d);
        
        shift->add_displacement(0,cu_z,d);
        shift->add_displacement(0,au_z,d-dr);
        shift->add_displacement(1,cu_z,-d);
        shift->add_displacement(1,au_z,dr-d);
        
        fifth_neighbor.add_modifier(shift);
        
        for(int i=0; i<2; i++)
            for(int j=0; j<2; j++)
            {
                fifth_neighbor.add_modifier( new DoubleADPMode(modes_x[i],modes_x[j],xx));
                fifth_neighbor.add_modifier( new DoubleADPMode(modes_x[i],modes_y[j],xy));
                fifth_neighbor.add_modifier( new DoubleADPMode(modes_x[i],modes_z[j],xz));
                fifth_neighbor.add_modifier( new DoubleADPMode(modes_y[i],modes_x[j],yx));
                fifth_neighbor.add_modifier( new DoubleADPMode(modes_y[i],modes_y[j],yy));
                fifth_neighbor.add_modifier( new DoubleADPMode(modes_y[i],modes_z[j],yz));
                fifth_neighbor.add_modifier( new DoubleADPMode(modes_z[i],modes_x[j],zx));
                fifth_neighbor.add_modifier( new DoubleADPMode(modes_z[i],modes_y[j],zy));
                fifth_neighbor.add_modifier( new DoubleADPMode(modes_z[i],modes_z[j],zz));
            }    
        
		pools.push_back(&fifth_neighbor);
    
		//--------------------------------------------------------------------------------------------------------------------------
		vector<AtomicPair> pairs;
		
        for(pool=pools.begin(); pool!=pools.end(); pool++)
        {
            (*pool)->invoke_correlators();
            pairs.insert(pairs.end(),(*pool)->pairs.begin(),(*pool)->pairs.end());
        }
		
		pairs=IntnsityCalculator::filter_pairs_from_asymmetric_unit(pairs,cell.laue_symmetry.label);
		//pairs=IntnsityCalculator::get_half_sphere_of_pairs_from_asymmetric_unit(pairs,cell.laue_symmetry.label);
		
	const bool DIRECT=true;
	if(DIRECT){
		IntnsityCalculator::calculate_scattering_from_pairs(pairs,I,AVERAGE_FLAG);
        I.invert();
	}else{
        I.invert();
		IntnsityCalculator::calculate_patterson_map_from_pairs_f(pairs,I,AVERAGE_FLAG,vec3<int>(16,16,16));
	}


		
		//--------------------------------------------------------------------------------------------------------------------------
		//copy result
		I.init_iterator();
	//int i=0;
		while(I.next())
		{
			*scattering++ = I.current_array_value();
			//cout << "current value is " << I.current_array_value() << endl;
		}
}




void calculate_scattering_tst(double *scattering,double* inp_params)
{      
    //--------------------------------------------------------------------------------------------------------------------------
    UnitCell cell(14.1,14.1,6.93,90,90,90); //Cell for diffuse scattering calculation. It is orthogonal

    bool AVERAGE_FLAG = inp_params[0]==1;

    IntensityMap I(104,104,12);
    double ll=-13;
    double gs=1.0/4;
    I.set_grid(cell.cell,vec3<double>(gs,gs,1),vec3<double>(ll,ll,-6));

    //-------------------------------------------------------------------------------------------------------------------------      
    cell.set_laue_symmetry("*");
        ChemicalUnitNode* node = new ChemicalUnitNode();
            AtomicAssembly* molecule = new AtomicAssembly();
            molecule->set_occupancy(0.5);

                AtomicAssembly* wing = new AtomicAssembly();

                    wing->add_chemical_unit(new Atom("O",  0.5, 0.5,0,0 ,   0.000245,0.000449,0.001837,0.000163,0.000204,0.000449));


               molecule->add_chemical_unit(wing);
               molecule->add_chemical_unit(wing->create_symmetric(mat3<double>(0,-1,0,1,-1,0,0,0,1),vec3<double>(0,0,0))); //3-fold axis once
               molecule->add_chemical_unit(wing->create_symmetric(mat3<double>(-1,1,0,-1,0,0,0,0,1),vec3<double>(0,0,0))); //3-fold axis twice

        node->add_chemical_unit(molecule);

        AtomicAssembly* molecule_refl = molecule->create_symmetric(mat3<double>(1,0,0,0,1,0,0,0,-1),vec3<double>(0,0,0));
        molecule_refl->set_occupancy(0.5);


          node->add_chemical_unit(molecule_refl);


    cell.add_node(node);



    //--------------------------------------------------------------------------------------------------------------------------

    const int X=0;
    const int Y=1;
    const int Z=2;
    //Create translational modes
    //--------------------------------------------------------------------------------------------------------------------------
    vector<AtomicPairPool*> pools;
    vector<AtomicPairPool*>::iterator pool;

    vector<double> c;
    c.push_back(0.5);

    AtomicPairPool zeroth_neighbor;
    zeroth_neighbor.add_modifiers( correlators_from_cuns(node,node,c) );
    //zeroth_neighbor.add_modifier( new ZeroVectorCorrelation );
    pools.push_back(&zeroth_neighbor);
    
		AtomicPairPool first_neighbor;
        first_neighbor.add_modifiers( correlators_from_cuns(node,node,c) );
		first_neighbor.add_modifier( new CellShifter(1,0,0) );

		//pools.push_back(&first_neighbor);

    //--------------------------------------------------------------------------------------------------------------------------
    vector<AtomicPair> pairs;

    for(pool=pools.begin(); pool!=pools.end(); pool++)
    {
        (*pool)->invoke_correlators();
        pairs.insert(pairs.end(),(*pool)->pairs.begin(),(*pool)->pairs.end());
    }

    //pairs=IntnsityCalculator::filter_pairs_from_asymmetric_unit(pairs,cell.laue_symmetry.label);
    //pairs=IntnsityCalculator::get_half_sphere_of_pairs_from_asymmetric_unit(pairs,cell.laue_symmetry.label);
    pairs=multiply_pairs_by_matrix(pairs,mat3<double>(1,-0.5,0,0,sqrt(3)/2,0,0,0,1));

    if(false) {
     for(vector<AtomicPair>::iterator pair=pairs.begin(); pair!=pairs.end(); pair++)
        cout << pair->r()[0] << ' ' << pair->r()[1] << ' ' << pair->r()[2] << endl;
	}
const bool DIRECT=false;
if(DIRECT){
    IntnsityCalculator::calculate_scattering_from_pairs(pairs,I,AVERAGE_FLAG);
    I.invert();
}else{
	const int small_grid_size=64;
    I.invert();
    IntnsityCalculator::calculate_patterson_map_from_pairs_f(pairs,I,AVERAGE_FLAG,vec3<int>(small_grid_size,small_grid_size,small_grid_size));
}



    //--------------------------------------------------------------------------------------------------------------------------
    //copy result
    I.init_iterator();
    while(I.next())
    {
        *scattering++ = I.current_array_value();
    }
}
void add_split_shift(AtomicPairPool& Pool,ADPMode& mode11,ADPMode& mode12,ADPMode& mode21,ADPMode& mode22,double dr, double c)
{
	const int START = 0;
	const int END = 1;
    double dr2=dr*(c+1)/(c-1);
    
    StaticShift* shift = new StaticShift();
    shift->add_displacement(START,mode11,dr2/2);
    shift->add_displacement(END,mode22,-dr2/2);
    Pool.add_modifier(shift);

    shift = new StaticShift();
    shift->add_displacement(START,mode12,dr2/2);
    shift->add_displacement(END,mode21,-dr2/2);
    Pool.add_modifier(shift);

    shift = new StaticShift();
    shift->add_displacement(START,mode11,dr/2);
    shift->add_displacement(END,mode21,-dr/2);
    Pool.add_modifier(shift);

    shift = new StaticShift();
    shift->add_displacement(START,mode12,dr/2);
    shift->add_displacement(END,mode22,-dr/2);
    Pool.add_modifier(shift);
}
#define ADD_NEIGHBOR(corr,x_shift,y_shift,name,name1,name2)\
	c[0] = corr;\
	AtomicPairPool name;\
	name.add_modifiers( correlators_from_cuns(node,node,c) );\
	name.add_modifier( new CellShifter(x_shift,y_shift,0) );\
	pools.push_back(&name);\
	AtomicPairPool name1;\
	name1.add_modifiers( correlators_from_cuns(node,node2,c) );\
	name1.add_modifier( new CellShifter(x_shift,y_shift,0) );\
	name1.add_modifier( new MultiplicityCorrelation(0.5) );\
	pools.push_back(&name1);\
	AtomicPairPool name2;\
	name2.add_modifiers( correlators_from_cuns(node2,node,c) );\
	name2.add_modifier( new CellShifter(x_shift,y_shift,1) );\
	name2.add_modifier( new MultiplicityCorrelation(0.5) );\
	pools.push_back(&name2);


void calculate_scattering(double *scattering,double* inp_params)
{       //Thomas-Philipp model
		//--------------------------------------------------------------------------------------------------------------------------
    //cout << "started calc";
	
double mlt=1;    
UnitCell cell(14.1*mlt,14.1*mlt,6.93*mlt,90,90,90); //Cell for diffuse scattering calculation. It is orthogonal
	
	bool AVERAGE_FLAG = inp_params[0]==1;
	
	const int x_half_size = 607;//607;
	const int x_size = x_half_size*2+1*0;
	
	IntensityMap I(x_size,x_size,12*1); //nado chtobi eta fignya sama schitalas'
	double ll=-12.8677;
	double gs=-ll/x_half_size;
	I.set_grid(cell.cell,vec3<double>(gs,gs,1.0/1),vec3<double>(ll,ll,-6));
	
	//-------------------------------------------------------------------------------------------------------------------------      
	cell.set_laue_symmetry("6*");
		ChemicalUnitNode* node = new ChemicalUnitNode();
			AtomicAssembly* molecule = new AtomicAssembly();
			molecule->set_occupancy(0.5);
				
				AtomicAssembly* wing = new AtomicAssembly();
				
	if(inp_params[1])
				wing->add_chemical_unit(new Atom("C1",  0.5,-0.010700,-0.103000,0.131900   ,0.000287,0.000236,0.000770,0.000136,0.000092,0.000031));
	if(inp_params[2])
				wing->add_chemical_unit(new Atom("C2",  0.5,0.092400,-0.010900,0.132800   ,0.000267,0.000277,0.000708,0.000166,0.000010,0.000072));
	if(inp_params[3])
				wing->add_chemical_unit(new Atom("C3",  0.5,0.195000,-0.017700,0.114000   ,0.000216,0.000277,0.000979,0.000101,0.000010,0.000000));
					
					AtomicAssembly* rot_part = new AtomicAssembly();
					
	if(inp_params[4])
					rot_part->add_chemical_unit(new Atom("O",  0.5,0.270900 ,0.048300,0.011100,   0.000231,0.000423,0.001832,0.000161,0.000205,0.000461));
	if(inp_params[5])
					rot_part->add_chemical_unit(new Atom("N",  0.5,0.198900,-0.096100,0.218000  ,                0.000292,0.000267,0.001978,0.000181,0.000368,0.000368));
	if(inp_params[6])
					rot_part->add_chemical_unit(new Atom("C4",  0.5,0.290700,-0.120200,0.209000   ,0.000262,0.000317,0.002082,0.000186,0.000194,0.000215));
	if(inp_params[7])				
					rot_part->add_chemical_unit(new Atom("C5",  0.5,0.393800,-0.022300,0.286000   , 0.000387,0.000493,0.002915,0.000241,-0.000102,0.000246));
	if(inp_params[8])
					rot_part->add_chemical_unit(new Atom("C6",  0.5,0.255900,-0.218000,0.339400   , 0.000543,0.000553,0.003873,0.000448,0.000737,0.000942));
	 if(inp_params[9])
					rot_part->add_chemical_unit(new Atom("C7",  0.5,0.308100,-0.143700,0.004000   ,  0.000533,0.000765,0.002582,0.000458,0.000368,-0.000133));
	if(inp_params[10])
					rot_part->add_chemical_unit(new Atom("H",  0.5,0.1454,0.8651-1,0.2964,0.078/198.81,0.078/198.81,0.078/48.0249,0.078/-198.81,0,0));
	if(inp_params[11])
					rot_part->add_chemical_unit(new Atom("H",  0.5,-0.019,0.8250-1,0.121,0.056/198.81,0.056/198.81,0.056/48.0249,0.056/-198.81,0,0));
	if(inp_params[12])
					rot_part->add_chemical_unit(new Atom("H",  0.5,0.3815,0.9921-1,0.4164,0.153/198.81,0.153/198.81,0.153/48.0249,0.153/-198.81,0,0));
	if(inp_params[13])
					rot_part->add_chemical_unit(new Atom("H",  0.5,-0.019,0.8250-1,0.121,0.056/198.81,0.056/198.81,0.056/48.0249,0.056/-198.81,0,0));
	if(inp_params[14])
					rot_part->add_chemical_unit(new Atom("H",  0.5,0.3815,0.9921-1,0.4164,0.153/198.81,0.153/198.81,0.153/48.0249,0.153/-198.81,0,0));
	if(inp_params[15])
					rot_part->add_chemical_unit(new Atom("H",  0.5,0.4119,1.0409-1,0.2075,0.153/198.81,0.153/198.81,0.153/48.0249,0.153/-198.81,0,0));
	if(inp_params[16])
					rot_part->add_chemical_unit(new Atom("H",  0.5,0.4530,0.9629-1,0.2829,0.153/198.81,0.153/198.81,0.153/48.0249,0.153/-198.81,0,0));
	if(inp_params[17])
					rot_part->add_chemical_unit(new Atom("H",  0.5,0.2595,0.8043-1,0.4714,0.179/198.81,0.179/198.81,0.179/48.0249,0.179/-198.81,0,0));
	if(inp_params[18])
					rot_part->add_chemical_unit(new Atom("H",  0.5,0.3040,0.7529-1,0.3210,0.179/198.81,0.179/198.81,0.179/48.0249,0.179/-198.81,0,0));
	if(inp_params[19])
					rot_part->add_chemical_unit(new Atom("H",  0.5,0.1822,0.7268-1,0.3089,0.179/198.81,0.179/198.81,0.179/48.0249,0.179/-198.81,0,0));
	if(inp_params[20])
					rot_part->add_chemical_unit(new Atom("H",  0.5,0.2425,    0.7932-1,   -0.0426,   0.174/198.81,    0.174/198.81,   0.174/48.0249,    0.174/-198.81,0,0));
					
			   wing->add_chemical_unit(rot_part);
					
		   molecule->add_chemical_unit(wing->create_symmetric(mat3<double>(0,-1,0,1,-1,0,0,0,1),vec3<double>(0,0,0))); //3-fold axis once
		   molecule->add_chemical_unit(wing->create_symmetric(mat3<double>(-1,1,0,-1,0,0,0,0,1),vec3<double>(0,0,0))); //3-fold axis twice
		   molecule->add_chemical_unit(wing);
		  

		node->add_chemical_unit(molecule);
		
		AtomicAssembly* molecule_refl = molecule->create_symmetric(mat3<double>(1,0,0,0,1,0,0,0,-1),vec3<double>(0,0,0.5));

		molecule_refl->set_occupancy(0.5);
			   
				
		  node->add_chemical_unit(molecule_refl);
				
	cell.add_node(node);
		ChemicalUnitNode* node2 = new ChemicalUnitNode();
			AtomicAssembly* molecule2 = molecule->create_symmetric(mat3<double>(-1,0,0, 0,-1,0, 0,0,1),vec3<double>(0,0,0.5));
			molecule2->set_occupancy(0.5);
		node2->add_chemical_unit(molecule2);
			AtomicAssembly* molecule_refl2 = molecule_refl->create_symmetric(mat3<double>(-1,0,0, 0,-1,0, 0,0,1),vec3<double>(0,0,0.5));
			molecule_refl2->set_occupancy(0.5);
		node2->add_chemical_unit(molecule_refl2);
	cell.add_node(node2);

	//--------------------------------------------------------------------------------------------------------------------------
	
	sym_mat3<double> metrical_matrix(198.81,198.8100,48.0249,-198.8100/2,0,0);
	
	const int X=0;
	const int Y=1;
	const int Z=2;
	//Create translational modes
	ADPMode mol_x = translational_mode(molecule,X);
	ADPMode mol_y = translational_mode(molecule,Y);
	ADPMode mol_z = translational_mode(molecule,Z);
	ADPMode mol_rot = rot_mode(molecule,vec3<double>(0,0,1),vec3<double>(0,0,0),metrical_matrix.inverse());
	
	ADPMode mol_refl_x = translational_mode(molecule_refl,X);
	ADPMode mol_refl_y = translational_mode(molecule_refl,Y);
	ADPMode mol_refl_z = translational_mode(molecule_refl,Z);
	ADPMode mol_refl_rot = rot_mode(molecule_refl,vec3<double>(0,0,1),vec3<double>(0,0,0),metrical_matrix.inverse());
	
	ADPMode mol_x2 = translational_mode(molecule2,X);
	ADPMode mol_y2 = translational_mode(molecule2,Y);
	ADPMode mol_z2 = translational_mode(molecule2,Z);
	ADPMode mol_rot2 = z_rot_mode(molecule2);
	
	ADPMode mol_refl_x2 = translational_mode(molecule_refl2,X);
	ADPMode mol_refl_y2 = translational_mode(molecule_refl2,Y);
	ADPMode mol_refl_z2 = translational_mode(molecule_refl2,Z);
	ADPMode mol_refl_rot2 = z_rot_mode(molecule_refl2);
	
	/*
	vector<ADPMode> rm;
	for(int i=0; i<3; i++)
		rm.push_back(rot_mode((*(*molecule)[i])[3],(*(*molecule)[i])[2]->get_atoms()[0]->r-(*(*molecule)[i])[1]->get_atoms()[0]->r,(*(*molecule)[i])[1]->get_atoms()[0]->r,metrical_matrix.inverse()));
	ADPMode mol_wing_rotations = combine_modes(rm);
	vector<ADPMode> rm_refl;
	for(int i=0; i<3; i++)
		rm_refl.push_back(rot_mode((*(*molecule_refl)[i])[3],(*(*molecule_refl)[i])[2]->get_atoms()[0]->r-(*(*molecule_refl)[i])[1]->get_atoms()[0]->r,(*(*molecule_refl)[i])[1]->get_atoms()[0]->r,metrical_matrix.inverse()));
	ADPMode mol_wing_rotations_refl = combine_modes(rm_refl);

	vector<ADPMode> rm2;
	for(int i=0; i<3; i++)
		rm2.push_back(rot_mode((*(*molecule2)[i])[3],(*(*molecule2)[i])[2]->get_atoms()[0]->r-(*(*molecule2)[i])[1]->get_atoms()[0]->r,(*(*molecule2)[i])[1]->get_atoms()[0]->r,metrical_matrix.inverse()));
	ADPMode mol_wing_rotations2 = combine_modes(rm2);
	vector<ADPMode> rm_refl2;
	for(int i=0; i<3; i++)
		rm_refl2.push_back(rot_mode((*(*molecule_refl2)[i])[3],(*(*molecule_refl2)[i])[2]->get_atoms()[0]->r-(*(*molecule_refl2)[i])[1]->get_atoms()[0]->r,(*(*molecule_refl2)[i])[1]->get_atoms()[0]->r,metrical_matrix.inverse()));
	ADPMode mol_wing_rotations_refl2 = combine_modes(rm_refl2);*/
	//--------------------------------------------------------------------------------------------------------------------------
	vector<AtomicPairPool*> pools;
	vector<AtomicPairPool*>::iterator pool;
	
	double w,xx,xy,yx,yy,rr,dr,dr2,ww,wr,zr,rw,zw,wz,rz,d,zz,xz,zx,yz,zy,drx,dry,dx,dy,dy1,dx1;
	double xx2,yy2,zz2,rr2,xy2,yx2;
	
	vector<double> c;
	c.push_back(0.5);
	c[0] = 0.5;
/*
	xx=inp_params[4];
	zz=inp_params[5];
	yy=xx;
	xy=xx/2;
	yx=xx/2;
	rr=inp_params[6];
	ww=inp_params[22];
	w=inp_params[23];
	wr=inp_params[24];
	zr=inp_params[25];
	zw=inp_params[26];*/

	AtomicPairPool zeroth_neighbor;
	zeroth_neighbor.add_modifiers( correlators_from_cuns(node,node,c) );
	zeroth_neighbor.add_modifier( new MultiplicityCorrelation(0.5) ); // for 2fold symmetry which is already there


	//zeroth_neighbor.add_modifier( new ZeroVectorCorrelation );
	/*
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_x,&mol_x,xx));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_x,&mol_y,xy));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_y,&mol_x,yx));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_y,&mol_y,yy));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_z,&mol_z,zz));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_rot,&mol_rot,rr));

	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_x,&mol_refl_x,xx));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_x,&mol_refl_y,xy));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_y,&mol_refl_x,yx));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_y,&mol_refl_y,yy));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_z,&mol_refl_z,zz));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_rot,&mol_refl_rot,rr));
	 
	for(int i=0; i<3; i++)
	{
		zeroth_neighbor.add_modifier( new DoubleADPMode(&rm[i],&rm[i],w));
		zeroth_neighbor.add_modifier( new DoubleADPMode(&rm_refl[i],&rm_refl[i],w));
	}
	
	
	for(int i=0; i<3; i++)
		for(int j=i; j<3; j++)
		{
			zeroth_neighbor.add_modifier( new DoubleADPMode(&rm[i],&rm[j],ww));
			zeroth_neighbor.add_modifier( new DoubleADPMode(&rm[j],&rm[i],ww));
			zeroth_neighbor.add_modifier( new DoubleADPMode(&rm_refl[i],&rm_refl[j],ww));
			zeroth_neighbor.add_modifier( new DoubleADPMode(&rm_refl[j],&rm_refl[i],ww));
		}
	
	
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_wing_rotations,&mol_rot,wr));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_rot,&mol_wing_rotations,wr));
	
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_z,&mol_rot,zr));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_rot,&mol_z,zr));
	
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_z,&mol_wing_rotations,zw));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_wing_rotations,&mol_z,zw));
	
	//zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_wing_rotations_refl,&mol_wing_rotations_refl,ww));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_wing_rotations_refl,&mol_refl_rot,-wr));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_rot,&mol_wing_rotations_refl,-wr));
	
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_z,&mol_refl_rot,-zr));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_rot,&mol_refl_z,-zr));
	
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_z,&mol_wing_rotations_refl,zw));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_wing_rotations_refl,&mol_refl_z,zw));*/
	

	AtomicPairPool zeroth_neighbor1;
	zeroth_neighbor1.add_modifier( new MultiplicityCorrelation(0.5*0.5) ); 
	zeroth_neighbor1.add_modifiers( correlators_from_cuns(node,node2,c) );
	AtomicPairPool zeroth_neighbor2;
	zeroth_neighbor2.add_modifier( new MultiplicityCorrelation(0.5*0.5) ); 
	zeroth_neighbor2.add_modifiers( correlators_from_cuns(node2,node,c) );
	zeroth_neighbor2.add_modifier( new CellShifter(0,0,1) );
	
	
	/*
	xx=inp_params[27];
	zz=inp_params[28];
	yy=xx;
	xy=xx/2;
	yx=xx/2;
	rr=inp_params[29];
	ww=inp_params[30];
	wr=inp_params[31];
	rw=inp_params[32];
	zr=inp_params[33];
	rz=inp_params[34];
	zw=inp_params[35];
	wz=inp_params[36];
	
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_x,&mol_x2,xx));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_x,&mol_y2,xy));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_y,&mol_x2,yx));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_y,&mol_y2,yy));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_z,&mol_z2,zz));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_rot,&mol_rot2,rr));

	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_x,&mol_refl_x2,xx));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_x,&mol_refl_y2,xy));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_y,&mol_refl_x2,yx));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_y,&mol_refl_y2,yy));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_z,&mol_refl_z2,zz));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_rot,&mol_refl_rot2,rr));
	
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_wing_rotations,&mol_rot2,wr));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_rot,&mol_wing_rotations2,rw));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_z,&mol_rot2,zr));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_rot,&mol_z2,rz));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_z,&mol_wing_rotations2,zw));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_wing_rotations,&mol_z2,wz));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_wing_rotations,&mol_wing_rotations2,ww));
	
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_wing_rotations_refl,&mol_wing_rotations_refl2,ww));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_wing_rotations_refl,&mol_refl_rot2,-rw));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_rot,&mol_wing_rotations_refl2,-wr));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_z,&mol_refl_rot2,-rz));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_rot,&mol_refl_z2,-zr));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_z,&mol_wing_rotations_refl2,wz));
	zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_wing_rotations_refl,&mol_refl_z2,zw));*/
	
	//pools.push_back(&zeroth_neighbor);
	
	if(false)
	{

		c[0]=0.25+inp_params[1];
		
		/*xx=inp_params[7];
		zz=inp_params[8];
		yy=xx;
		xy=xx/2;
		yx=xx/2;
		rr=inp_params[9];
		
		xx2=inp_params[38];
		zz2=inp_params[39];
		yy2=xx2;
		xy2=xx2/2;
		yx2=xx2/2;
		rr2=inp_params[40];*/
		
		dx=inp_params[2];
		dy=inp_params[3];

		
		
		AtomicPairPool first_neighbor;
		first_neighbor.add_modifiers( correlators_from_cuns(node,node,c) );
		first_neighbor.add_modifier( new CellShifter(1,0,0) );
		first_neighbor.add_modifier( new MultiplicityCorrelation(1) );
		/*
		first_neighbor.add_modifier( new DoubleADPMode(&mol_x,&mol_x,xx));
		first_neighbor.add_modifier( new DoubleADPMode(&mol_x,&mol_y,xy));
		first_neighbor.add_modifier( new DoubleADPMode(&mol_y,&mol_x,yx));
		first_neighbor.add_modifier( new DoubleADPMode(&mol_y,&mol_y,yy));
		first_neighbor.add_modifier( new DoubleADPMode(&mol_z,&mol_z,zz));
		first_neighbor.add_modifier( new DoubleADPMode(&mol_rot,&mol_rot,rr));

		first_neighbor.add_modifier( new DoubleADPMode(&mol_refl_x,&mol_refl_x,xx2));
		first_neighbor.add_modifier( new DoubleADPMode(&mol_refl_x,&mol_refl_y,xy2));
		first_neighbor.add_modifier( new DoubleADPMode(&mol_refl_y,&mol_refl_x,yx2));
		first_neighbor.add_modifier( new DoubleADPMode(&mol_refl_y,&mol_refl_y,yy2));
		first_neighbor.add_modifier( new DoubleADPMode(&mol_refl_z,&mol_refl_z,zz2));
		first_neighbor.add_modifier( new DoubleADPMode(&mol_refl_rot,&mol_refl_rot,rr2));*/
		
		add_split_shift( first_neighbor,mol_x,mol_refl_x,mol_x,mol_refl_x,dx+dy/2, c[0]*4-1 );
		add_split_shift( first_neighbor,mol_y,mol_refl_y,mol_y,mol_refl_y,dy, c[0]*4-1 );


		pools.push_back(&first_neighbor);
	

		
		
		AtomicPairPool first_neighbor1;
		first_neighbor1.add_modifiers( correlators_from_cuns(node,node2,c) );
		first_neighbor1.add_modifier( new CellShifter(1,0,0) );
		first_neighbor1.add_modifier( new MultiplicityCorrelation(0.5) );
		
		add_split_shift( first_neighbor1,mol_x,mol_refl_x,mol_x2,mol_refl_x2,dx, c[0]*4-1 );
		add_split_shift( first_neighbor1,mol_y,mol_refl_y,mol_y2,mol_refl_y2,dy, c[0]*4-1 );
		pools.push_back(&first_neighbor1);
		
		
		AtomicPairPool first_neighbor2;
		first_neighbor2.add_modifiers( correlators_from_cuns(node2,node,c) );
		first_neighbor2.add_modifier( new CellShifter(1,0,1) );
		first_neighbor2.add_modifier( new MultiplicityCorrelation(0.5) );
		
		add_split_shift( first_neighbor2,mol_x2,mol_refl_x2,mol_x,mol_refl_x,dx, c[0]*4-1 );
		add_split_shift( first_neighbor2,mol_y2,mol_refl_y2,mol_y,mol_refl_y,dy, c[0]*4-1 );
		pools.push_back(&first_neighbor2);


		
		c[0]=0.25+inp_params[4];
		dx=inp_params[5];
		dy=inp_params[6];
		
		AtomicPairPool second_neighbor;
		second_neighbor.add_modifiers( correlators_from_cuns(node,node,c) );
		second_neighbor.add_modifier( new CellShifter(2,1,0) );
		
		add_split_shift( second_neighbor,mol_x,mol_refl_x,mol_x,mol_refl_x,dx, c[0]*4-1 );
		add_split_shift( second_neighbor,mol_y,mol_refl_y,mol_y,mol_refl_y,dy, c[0]*4-1 );
		pools.push_back(&second_neighbor);
		
		AtomicPairPool second_neighbor1;
		second_neighbor1.add_modifiers( correlators_from_cuns(node,node2,c) );
		second_neighbor1.add_modifier( new CellShifter(2,1,0) );
		second_neighbor1.add_modifier( new MultiplicityCorrelation(0.5) );
		add_split_shift( second_neighbor1,mol_x,mol_refl_x,mol_x2,mol_refl_x2,dx, c[0]*4-1 );
		add_split_shift( second_neighbor1,mol_y,mol_refl_y,mol_y2,mol_refl_y2,dy, c[0]*4-1 );
		pools.push_back(&second_neighbor1);
		
		AtomicPairPool second_neighbor2;
		second_neighbor2.add_modifiers( correlators_from_cuns(node2,node,c) );
		second_neighbor2.add_modifier( new CellShifter(2,1,1) );
		second_neighbor2.add_modifier( new MultiplicityCorrelation(0.5) );
	
		add_split_shift( second_neighbor2,mol_x2,mol_refl_x2,mol_x,mol_refl_x,dx, c[0]*4-1 );
		add_split_shift( second_neighbor2,mol_y2,mol_refl_y2,mol_y,mol_refl_y,dy, c[0]*4-1 );
		pools.push_back(&second_neighbor2);

		c[0]=0.25+inp_params[7];
		
		dx=inp_params[8];
		dy=inp_params[9];
		
		AtomicPairPool third_neighbor;
		third_neighbor.add_modifiers( correlators_from_cuns(node,node,c) );
		third_neighbor.add_modifier( new CellShifter(2,0,0) );
		
		add_split_shift( third_neighbor,mol_x,mol_refl_x,mol_x,mol_refl_x,dx+dy/2, c[0]*4-1 );
		add_split_shift( third_neighbor,mol_y,mol_refl_y,mol_y,mol_refl_y,dx, c[0]*4-1 );
		pools.push_back(&third_neighbor);
			
		AtomicPairPool third_neighbor1;
		third_neighbor1.add_modifiers( correlators_from_cuns(node,node2,c) );
		third_neighbor1.add_modifier( new CellShifter(2,0,0) );
		third_neighbor1.add_modifier( new MultiplicityCorrelation(0.5) );
		add_split_shift( third_neighbor1,mol_x,mol_refl_x,mol_x2,mol_refl_x2,dx+dy/2, c[0]*4-1 );
		add_split_shift( third_neighbor1,mol_y,mol_refl_y,mol_y2,mol_refl_y2,dx, c[0]*4-1 );
		pools.push_back(&third_neighbor1);
		
		AtomicPairPool third_neighbor2;
		third_neighbor2.add_modifiers( correlators_from_cuns(node2,node,c) );
		third_neighbor2.add_modifier( new CellShifter(2,0,1) );
		third_neighbor2.add_modifier( new MultiplicityCorrelation(0.5) );
		add_split_shift( third_neighbor2,mol_x2,mol_refl_x2,mol_x,mol_refl_x,dx+dy/2, c[0]*4-1 );
		add_split_shift( third_neighbor2,mol_y2,mol_refl_y2,mol_y,mol_refl_y,dx, c[0]*4-1 );
		pools.push_back(&third_neighbor2);
		
		ADD_NEIGHBOR(0.25+inp_params[10],3,1,fourth_neighbor,fourth_neighbor1,fourth_neighbor2)
		ADD_NEIGHBOR(0.25+inp_params[11],3-1,-1,fourth_neighbor_prime,fourth_neighbor_prime1,fourth_neighbor_prime2)
	
		ADD_NEIGHBOR(0.25+inp_params[12],3,0,fifth_neighbor,fifth_neighbor1,fifth_neighbor2)
	
		ADD_NEIGHBOR(0.25+inp_params[13],4,2,sixth_neighbor,sixth_neighbor1,sixth_neighbor2)
		ADD_NEIGHBOR(0.25+inp_params[14],4,1,seventh_neighbor,seventh_neighbor1,seventh_neighbor2)
		ADD_NEIGHBOR(0.25+inp_params[15],4-1,-1,seventh_neighbor_prime,seventh_neighbor_prime1,seventh_neighbor_prime2)
		ADD_NEIGHBOR(0.25+inp_params[16],4,0,eighth_neighbor,eighth_neighbor1,eighth_neighbor2)
	/*
		ADD_NEIGHBOR(0.25+inp_params[17],5,2,nineth_neighbor,nineth_neighbor1,nineth_neighbor2)
		ADD_NEIGHBOR(0.25+inp_params[18],5-2,-2,nineth_neighbor_prime,nineth_neighbor_prime1,nineth_neighbor_prime2)
		ADD_NEIGHBOR(0.25+inp_params[19],5,1,tenth_neighbor,tenth_neighbor1,tenth_neighbor2)
		ADD_NEIGHBOR(0.25+inp_params[20],5-1,-1,tenth_neighbor_pime,tenth_neighbor_pime1,tenth_neighbor_pime2)
		ADD_NEIGHBOR(0.25+inp_params[21],5,0,eleventh_neighbor,eleventh_neighbor1,eleventh_neighbor2)
	
		ADD_NEIGHBOR(0.25+inp_params[22],6,3,twelveth_neighbor,twelveth_neighbor1,twelveth_neighbor2)
		ADD_NEIGHBOR(0.25+inp_params[23],6,2,thirteenth_neighbor,thirteenth_neighbor1,thirteenth_neighbor2)
		ADD_NEIGHBOR(0.25+inp_params[24],6-2,-2,thirteenth_neighbor_prime,thirteenth_neighbor_prime1,thirteenth_neighbor_prime2)
		ADD_NEIGHBOR(0.25+inp_params[25],6,1,fourteenth_neighbor,fourteenth_neighbor1,fourteenth_neighbor2)
		ADD_NEIGHBOR(0.25+inp_params[26],6-1,-1,fourteenth_neighbor_prime,fourteenth_neighbor_prime1,fourteenth_neighbor_prime2)
		ADD_NEIGHBOR(0.25+inp_params[27],6,0,fifteenth_neighbor,fifteenth_neighbor1,fifteenth_neighbor2)*/
	}
	//--------------------------------------------------------------------------------------------------------------------------
	vector<AtomicPair> pairs,pairs1;
	
	for(pool=pools.begin(); pool!=pools.end(); pool++)
	{
		(*pool)->invoke_correlators();
		pairs.insert(pairs.end(),(*pool)->pairs.begin(),(*pool)->pairs.end());
	}
	
	//pairs=IntnsityCalculator::filter_pairs_from_asymmetric_unit(pairs,cell.laue_symmetry.label);
	pairs=IntnsityCalculator::get_half_sphere_of_pairs_from_asymmetric_unit(pairs,"3");
	
	//special treatment for zero neighbor
	zeroth_neighbor.invoke_correlators();
	pairs1.insert(pairs.end(),zeroth_neighbor.pairs.begin(),zeroth_neighbor.pairs.end());
	
	zeroth_neighbor1.invoke_correlators();
	pairs.insert(pairs.end(),zeroth_neighbor1.pairs.begin(),zeroth_neighbor1.pairs.end());
	zeroth_neighbor2.invoke_correlators();
	pairs.insert(pairs.end(),zeroth_neighbor2.pairs.begin(),zeroth_neighbor2.pairs.end());


	
	pairs=IntnsityCalculator::filter_pairs_from_asymmetric_unit(pairs,"Pz"); // Cut a part which accounts to mirror plane at z=0.5 
	pairs1=IntnsityCalculator::filter_pairs_from_asymmetric_unit(pairs1,"*"); // Cut a part for mirror at z=0
	pairs.insert(pairs.end(),pairs1.begin(),pairs1.end()); //combine

	//cout << pairs.size();
	/*for(int i=0; i<pairs.size(); i++)
		cout << pairs[i].r()[0] << ' ' << pairs[i].r()[1] << ' ' << pairs[i].r()[2] << ' ' << pairs[i].average_p() << ';' <<endl;*/
	
	pairs=multiply_pairs_by_matrix(pairs,mat3<double>(1,-0.5,0, 0,sqrt(3)/2,0, 0,0,1));
	
	//cout << "number of pairs is " << pairs.size();

	//cout << pairs.size();
	/*for(int i=0; i<pairs.size(); i++)
		cout << pairs[i].r()[0] << ' ' << pairs[i].r()[1] << ' ' << pairs[i].r()[2] << ' ' << pairs[i].average_p() << ';' <<endl;*/
	
	
	const bool DIRECT=true;
	if(DIRECT){
		IntnsityCalculator::calculate_scattering_from_pairs(pairs,I,AVERAGE_FLAG);
		//I.invert();
	}else{
		const int small_grid_size=8;
		I.invert();
		IntnsityCalculator::calculate_patterson_map_from_pairs_f(pairs,I,AVERAGE_FLAG,vec3<int>(small_grid_size,small_grid_size,small_grid_size));
	}

	//-----------------------------------------------------------------------------------------------------------------------------
	//apply symmetry in numerical way
	const int z_half_size = 18;
	
	/*

		for(int i=0; i<x_size; i++)
			for(int j=0; j<x_size; j++)
				for(int k=1; k<z_half_size; k++)
				{
					I.at(i, j, k) += I.at(i, j, z_half_size*2-k); //calc half
					I.at(i, j, z_half_size*2-k) = I.at(i, j, k); //copy
				}
		
		for(int i=0; i<x_size; i++)
			for(int j=0; j<x_size; j++)
			{
				I.at(i, j, 0) *= 2; // double last layer
				I.at(i, j, z_half_size) *= 2; // double middle layer
			}
		
		for(int i=0; i<x_size; i++)
			for(int j=0; j<x_half_size; j++)
				for(int k=0; k<z_half_size*2; k++)
				{
					I.at(i, j, k) += I.at(x_half_size*2-i, x_half_size*2-j, k); //calc half
					I.at(x_half_size*2-i, x_half_size*2-j, k) = I.at(i, j, k);
				}
		
		
		for(int i=0; i<x_half_size; i++)
			for(int j=x_half_size; j<x_half_size+1; j++)
				for(int k=0; k<z_half_size*2; k++)
				{
					I.at(i, j, k) += I.at(x_half_size*2-i, x_half_size*2-j, k); //calc half
					I.at(x_half_size*2-i, x_half_size*2-j, k) = I.at(i, j, k);
				}
	
		
	for(int i=0; i<x_half_size*2; i++)
		for(int j=0; j<x_half_size*2; j++)
			for(int k=12; k<12+12; k++)
			{
				I.at(i,j,k) += I.at(i,j,(k+12)) + I.at(i,j,(k-12)); //wrap periodic boundary in
				*scattering++ = I.at(i,j,k); //and copy to output
			}*/
				
	//--------------------------------------------------------------------------------------------------------------------------
	//copy result
	
	I.init_iterator();
	while(I.next())
	{
		*scattering++ = I.current_array_value();
	}
}

void calculate_scattering_pg(double *scattering,double* inp_params)
{       //Thomas-Philipp model pg
		//--------------------------------------------------------------------------------------------------------------------------
    double mlt=1;    
    UnitCell cell(14.1*mlt,14.1*mlt,6.93*mlt,90,90,90); //Cell for diffuse scattering calculation. It is orthogonal
        
		bool AVERAGE_FLAG = inp_params[0]==1;
		
        
		IntensityMap I(104,104,12); //nado chtobi eta fignya sama schitalas'
        double ll=-12.8749;
        double gs=-ll*2/104;
		I.set_grid(cell.cell,vec3<double>(gs,gs,1),vec3<double>(ll,ll,-6));
		
		//-------------------------------------------------------------------------------------------------------------------------      
		cell.set_laue_symmetry("6*");
			ChemicalUnitNode* node = new ChemicalUnitNode();
                AtomicAssembly* molecule = new AtomicAssembly();
                molecule->set_occupancy(0.5);
                    
                    AtomicAssembly* wing = new AtomicAssembly();
                        
                        wing->add_chemical_unit(new Atom("O",  0.5,0.270900,0.048300,0.011100 ,   0.000245,0.000449,0.001837,0.000163,0.000204,0.000449));
                        wing->add_chemical_unit(new Atom("N",  0.5,0.198900,-0.096100,0.218000  ,   0.000311,0.000286,0.001878,0.000194,0.000316,0.000306));
                        wing->add_chemical_unit(new Atom("C1",  0.5,-0.010700,-0.103000,0.131900   ,0.000311,0.000260,0.000837,0.000153,0.000102,0.000041));
                        wing->add_chemical_unit(new Atom("C2",  0.5,0.092400,-0.010900,0.132800   ,0.000281,0.000291,0.000816,0.000173,0.000000,0.000061));
                        wing->add_chemical_unit(new Atom("C3",  0.5,0.195000,-0.017700,0.114000   ,0.000235,0.000286,0.001102,0.000112,-0.000010,0.000000));
                        wing->add_chemical_unit(new Atom("C4",  0.5,0.290700,-0.120200,0.209000   ,0.000281,0.000347,0.002041,0.000199,0.000184,0.000214));
                        wing->add_chemical_unit(new Atom("C5",  0.5,0.393800,-0.022300,0.286000   , 0.000408,0.000510,0.002857,0.000250,0.000082,-0.000194));
                        wing->add_chemical_unit(new Atom("C6",  0.5,0.256600,-0.217600,0.338800   , 0.000566,0.000566,0.004082,0.000459,-0.000765,-0.000959));
                        wing->add_chemical_unit(new Atom("C7",  0.5,0.308100,-0.143700,0.004000   ,  0.000556,0.000786,0.002510,0.000474,0.000388,-0.000051));

                        
                   molecule->add_chemical_unit(wing);
                   molecule->add_chemical_unit(wing->create_symmetric(mat3<double>(0,-1,0,1,-1,0,0,0,1),vec3<double>(0,0,0))); //3-fold axis once
                   molecule->add_chemical_unit(wing->create_symmetric(mat3<double>(-1,1,0,-1,0,0,0,0,1),vec3<double>(0,0,0))); //3-fold axis twice
                   

			node->add_chemical_unit(molecule);
            
            AtomicAssembly* molecule_refl = molecule->create_symmetric(mat3<double>(1,0,0,0,1,0,0,0,-1),vec3<double>(0,0,0.5));

            molecule_refl->set_occupancy(0.5);
                   
                    
              node->add_chemical_unit(molecule_refl);
                    
            
            
		cell.add_node(node);
            ChemicalUnitNode* node2 = new ChemicalUnitNode();
                AtomicAssembly* molecule2 = molecule->create_symmetric(mat3<double>(-1,0,0, 0,-1,0, 0,0,1),vec3<double>(0,0,-0.5));
                molecule2->set_occupancy(0.5);
            node2->add_chemical_unit(molecule2);
                AtomicAssembly* molecule_refl2 = molecule_refl->create_symmetric(mat3<double>(-1,0,0, 0,-1,0, 0,0,1),vec3<double>(0,0,-0.5));
                molecule_refl2->set_occupancy(0.5);
            node2->add_chemical_unit(molecule_refl2);
        cell.add_node(node2);
        


        //--------------------------------------------------------------------------------------------------------------------------
        
        const int X=0;
		const int Y=1;
        const int Z=2;
		//Create translational modes
        ADPMode mol_x = translational_mode(molecule,X);
        ADPMode mol_y = translational_mode(molecule,Y);
        ADPMode mol_z = translational_mode(molecule,Z);
        ADPMode mol_rot = z_rot_mode(molecule);
        
        ADPMode mol_refl_x = translational_mode(molecule_refl,X);
        ADPMode mol_refl_y = translational_mode(molecule_refl,Y);
        ADPMode mol_refl_z = translational_mode(molecule_refl,Z);
        ADPMode mol_refl_rot = z_rot_mode(molecule_refl);
        
        ADPMode mol_x2 = translational_mode(molecule2,X);
        ADPMode mol_y2 = translational_mode(molecule2,Y);
        ADPMode mol_z2 = translational_mode(molecule2,Z);
        ADPMode mol_rot2 = z_rot_mode(molecule2);
        
        ADPMode mol_refl_x2 = translational_mode(molecule_refl2,X);
        ADPMode mol_refl_y2 = translational_mode(molecule_refl2,Y);
        ADPMode mol_refl_z2 = translational_mode(molecule_refl2,Z);
        ADPMode mol_refl_rot2 = z_rot_mode(molecule_refl2);
		//--------------------------------------------------------------------------------------------------------------------------
		vector<AtomicPairPool*> pools;
        vector<AtomicPairPool*>::iterator pool;
        
        double xx,xy,yx,yy,rr,dr,d,zz,xz,zx,yz,zy,drx,dry,dx,dy;
        
        vector<double> c;
        c.push_back(0.25+inp_params[1]);

        xx=inp_params[2];
        zz=inp_params[3];
        yy=xx;
        xy=xx/2;
        yx=xx/2;
        rr=inp_params[4];

        AtomicPairPool zeroth_neighbor;
        zeroth_neighbor.add_modifiers( correlators_from_cuns(node,node,c) );
        zeroth_neighbor.add_modifiers( correlators_from_cuns(node,node2,c) );
        /*
        zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_x,&mol_x,xx));
        zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_x,&mol_y,xy));
        zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_y,&mol_x,yx));
        zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_y,&mol_y,yy));
        zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_z,&mol_z,zz));
        zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_rot,&mol_rot,rr));

        zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_x,&mol_refl_x,xx));
        zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_x,&mol_refl_y,xy));
        zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_y,&mol_refl_x,yx));
        zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_y,&mol_refl_y,yy));
        zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_z,&mol_refl_z,zz));
        zeroth_neighbor.add_modifier( new DoubleADPMode(&mol_refl_rot,&mol_refl_rot,rr));*/
        
        
        
        //pools.push_back(&zeroth_neighbor);
        

        AtomicPairPool zeroth_neighbor1;
        zeroth_neighbor1.add_modifiers( correlators_from_cuns(node,node,c) );
        zeroth_neighbor1.add_modifiers( correlators_from_cuns(node,node2,c) );
        
        zeroth_neighbor1.add_modifier( new DoubleADPMode(&mol_x,&mol_x,xx));
        zeroth_neighbor1.add_modifier( new DoubleADPMode(&mol_x,&mol_y,xy));
        zeroth_neighbor1.add_modifier( new DoubleADPMode(&mol_y,&mol_x,yx));
        zeroth_neighbor1.add_modifier( new DoubleADPMode(&mol_y,&mol_y,yy));
        zeroth_neighbor1.add_modifier( new DoubleADPMode(&mol_z,&mol_z,zz));
        zeroth_neighbor1.add_modifier( new DoubleADPMode(&mol_rot,&mol_rot,rr));

        zeroth_neighbor1.add_modifier( new DoubleADPMode(&mol_refl_x,&mol_refl_x,xx));
        zeroth_neighbor1.add_modifier( new DoubleADPMode(&mol_refl_x,&mol_refl_y,xy));
        zeroth_neighbor1.add_modifier( new DoubleADPMode(&mol_refl_y,&mol_refl_x,yx));
        zeroth_neighbor1.add_modifier( new DoubleADPMode(&mol_refl_y,&mol_refl_y,yy));
        zeroth_neighbor1.add_modifier( new DoubleADPMode(&mol_refl_z,&mol_refl_z,zz));
        zeroth_neighbor1.add_modifier( new DoubleADPMode(&mol_refl_rot,&mol_refl_rot,rr));
        
        zeroth_neighbor1.add_modifier( new CellShifter(1,0,0) );
        
        //pools.push_back(&zeroth_neighbor1);
        

        
		//--------------------------------------------------------------------------------------------------------------------------
		vector<AtomicPair> pairs;
		
        for(pool=pools.begin(); pool!=pools.end(); pool++)
        {
            (*pool)->invoke_correlators();
            pairs.insert(pairs.end(),(*pool)->pairs.begin(),(*pool)->pairs.end());
        }
		
        
		pairs=IntnsityCalculator::get_half_sphere_of_pairs_from_asymmetric_unit(pairs,cell.laue_symmetry.label);
        
        //special treatment for zero neighbor
        zeroth_neighbor.invoke_correlators();
        pairs.insert(pairs.end(),zeroth_neighbor.pairs.begin(),zeroth_neighbor.pairs.end());
        
        
        pairs=multiply_pairs_by_matrix(pairs,mat3<double>(1,-0.5,0, 0,sqrt(3)/2,0, 0,0,1));
        //pairs=multiply_pairs_by_matrix(pairs,mat3<double>(-1,0,0, 0,1,0, 0,0,1));

	

		const bool DIRECT=false;
		if(DIRECT){
			IntnsityCalculator::calculate_scattering_from_pairs(pairs,I,AVERAGE_FLAG);
			I.invert();
		}else{
			I.invert();
			IntnsityCalculator::calculate_patterson_map_from_pairs_f(pairs,I,AVERAGE_FLAG,vec3<int>(8,8,8));
		}


		
		//--------------------------------------------------------------------------------------------------------------------------
		//copy result
		I.init_iterator();
		while(I.next())
		{
			*scattering++ = I.current_array_value();
		}
}