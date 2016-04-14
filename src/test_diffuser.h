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

#include <cxxtest/TestSuite.h>
#include "basic_classes.h"
#include <boost/fusion/tuple.hpp>

#include "precompiled_header.h"
#include "InputFileParser.h"

#include <cstdlib> 

#define a_skipper a_parser.skipper_no_assignement

#include <iostream>
#include <fstream>

#ifndef PATH_TO_YELL_SRC
#define PATH_TO_YELL_SRC "./src"
#endif

OutputHandler report;

namespace qi = boost::spirit::qi;
typedef std::string::iterator iterator_;

class TestDummy
{
public:
	class Destructed {};
	TestDummy()  {}
	~TestDummy()
	{
		throw Destructed();
	}
};
class TestMinimizerCalculator : public MinimizerCalculator{
 public:
  TestMinimizerCalculator(): data_(10,1,1) { }
  
  void calculate(vector<double> inp_params)
  {
    for(int i=0; i<10; i++)
      data_.at(i)=inp_params[i];
  }
  
  IntensityMap& data() { return data_; }
  
  IntensityMap data_;
};

class MyTestSuite : public CxxTest::TestSuite
{
 public:
	void setUp() {
    unit_matrix=sym_mat3<double>(1,1,1,0,0,0);
		trivial_unit_cell=cctbx::uctbx::unit_cell(scitbx::af::tiny<double,6>(1,1,1,90,90,90));
		atom1 = Atom(string("C"),0.5,0,0,0,1,1,1,0,0,0);
		atom2 = Atom(string("C2"),0.5,.5,.5,.5,1,1,1,0,0,0);
		p_atom1 = new Atom(atom1);
		p_atom2 = new Atom(atom2);
    
	}
	
  InputParser a_parser;
  sym_mat3<double> unit_matrix;
	cctbx::uctbx::unit_cell trivial_unit_cell;
	Atom atom1;
	Atom atom2;
	
	Atom* p_atom1;
	Atom* p_atom2;
  
  void testOutputHandlerSitsShut() {
    REPORT(FIRST_RUN) << "ERROR in testOutputHandlerSitsShut\n";
  }
	
  void testIntensityMapActsLikeArray (void)  {
    IntensityMap I(2,1,1);
    I.at(0,0,0)=1;
    
    TS_ASSERT_EQUALS(I.at(0,0,0),1);
  }
  void testIntensityMap_flip_signs_in_chessboard_way()
  {
    IntensityMap I(2,1,1);
		I.at(0,0,0)=1;
		I.at(1,0,0)=1;
    
    I.flip_signs_in_chessboard_way();
    TS_ASSERT_EQUALS(1,I.at(0,0,0)); // 1 1 -> -1 (1)
    TS_ASSERT_EQUALS(-1,I.at(1,0,0));
    
    IntensityMap I4(4,1,1);
    for(int i=0; i<4; ++i)
      I4.at(i,0,0)=1;
    
    I4.flip_signs_in_chessboard_way();
    TS_ASSERT_EQUALS(1,I4.at(0,0,0)); // 1 1 (1) 1 -> 1 -1 (1) -1
    TS_ASSERT_EQUALS(-1,I4.at(1,0,0)); 
    TS_ASSERT_EQUALS(1,I4.at(2,0,0)); 
    TS_ASSERT_EQUALS(-1,I4.at(3,0,0));     
  }
  void test_grid_gives_grid_in_pdf_space()
  {
    Grid grid(cctbx::uctbx::unit_cell(scitbx::af::tiny<double,6>(1,1,1,90,90,90)), vec3<double>(1,1,1),vec3<double>(-1,0,0));
    TS_ASSERT_EQUALS(false,grid.in_pdf_space().reciprocal_flag);
  }
  void test_grid_checks_pairs_within()
  {
    Grid grid(cctbx::uctbx::unit_cell(scitbx::af::tiny<double,6>(1,1,1,90,90,90)), vec3<double>(1,1,1),vec3<double>(-1,-1,-1),false);
    // we will assume grid is symmetric boundaries
    AtomicPair a_pair(atom1,atom2);
    TS_ASSERT_EQUALS(true, a_pair.pair_is_withing(grid));
    a_pair.average_r()=vec3<double>(2,0,0);
    TS_ASSERT_EQUALS(false, a_pair.pair_is_withing(grid));
  }
  
  
  void testIntensityConditionalInverts()
	{
		IntensityMap I(2,1,1);
    I.set_grid(cctbx::uctbx::unit_cell(scitbx::af::tiny<double,6>(1,1,1,90,90,90)), vec3<double>(1,1,1),vec3<double>(-1,0,0));
    double reciprocal=12;
		I.at(0,0,0)=reciprocal;
		I.at(1,0,0)=17;
    double pdf=(17.0-12);
    
		I.to_reciprocal();
		TS_ASSERT_EQUALS(reciprocal,I.at(0,0,0));
    I.to_real();
		TS_ASSERT_EQUALS(pdf,I.at(0,0,0));
    I.to_real();
		TS_ASSERT_EQUALS(pdf,I.at(0,0,0));
    I.to_reciprocal();
		TS_ASSERT_EQUALS(reciprocal,I.at(0,0,0));
    
	}
	void testIntensityInvert()
	{
		IntensityMap I(2,1,1);
    I.set_grid(cctbx::uctbx::unit_cell(scitbx::af::tiny<double,6>(1,1,1,90,90,90)), vec3<double>(1,1,1),vec3<double>(-1,0,0));
    
    I.at(0,0,0)=1;
		I.at(1,0,0)=1;
		
		I.invert();
		
		TS_ASSERT_EQUALS(I.at(0,0,0),0);
		TS_ASSERT_EQUALS(I.at(1,0,0),2);
		
    
		I.at(0,0,0)=12;
		I.at(1,0,0)=17;
		I.invert();
		
		TS_ASSERT_EQUALS((17.0-12)/2,I.at(0,0,0));
		TS_ASSERT_EQUALS((17.0+12)/2,I.at(1,0,0));

    I.invert();
		
		TS_ASSERT_EQUALS(12,I.at(0,0,0));
		TS_ASSERT_EQUALS(17,I.at(1,0,0));
	}

  
	void test_ew_mult()
	{
		vec3<double> v(1,2,3);

		TS_ASSERT_EQUALS(ew_mult(v,v),vec3<double>(1,4,9));
		TS_ASSERT_EQUALS(ew_mult(v,vec3<double>(2,2,2)),
							vec3<double>(2,4,6));
	}
	
	void testIntensityMapActsLikeIterator()
	{
		IntensityMap I(1,2,3);
		I.init_iterator();
		for(int i=0; i<6; i++)
		{
			TS_ASSERT_EQUALS(I.next(),true);
			I.current_array_value()=i;
		}
		
		TS_ASSERT_EQUALS(I.next(),false);
		TS_ASSERT_EQUALS(I.at(0,1,2),5);
	}

	void testGrid()
	{
		Grid aGrid(cctbx::uctbx::unit_cell(scitbx::af::tiny<double,6>(10,10,10,90,90,90)), //unit cell
					vec3<double>(0.1,0.1,0.1), //grid steps
					vec3<double>(-1,-1,-1)); //lower limits of grid
		TS_ASSERT_EQUALS(aGrid.s_at(af::tiny<int,3>(0,0,0)),vec3<double>(-1,-1,-1));
		TS_ASSERT_EQUALS(aGrid.s_at(af::tiny<int,3>(0,1,2)),vec3<double>(-1,-0.9,-.8));
		
		TS_ASSERT_DELTA(aGrid.d_star_square_at(af::tiny<int,3>(0,0,0)),.03,.0000001);
	}
	
	void testIntensityIsIntegratedWithGrid()
	{
		IntensityMap I(2,2,2);
		I.set_grid(trivial_unit_cell,
					vec3<double>(1,1,1), 
					vec3<double>(-.5,-.5,-.5));
		
		I.set_grid(Grid(trivial_unit_cell,
						vec3<double>(1,1,1), 
						vec3<double>(-.5,-.5,-.5)));
		
		I.init_iterator();
		I.next();
		TS_ASSERT_DELTA(I.current_d_star_square(),0.75,0.00001);
		TS_ASSERT_EQUALS(I.current_s(),vec3<double>(-.5,-.5,-.5));
	}

	void testIntensityCalculatorAtLeastPretendsThatItWorks()
	{
		bool AVERAGE = true;
		vector<AtomicPair> pairs;
		//Atom atom1(string("C"),0.5,0,0,0,1,1,1,0,0,0);
		//Atom atom2(string("C2"),0.5,.5,.5,.5,1,1,1,0,0,0);
		pairs.push_back(AtomicPair(atom1,atom2));
		pairs.push_back(AtomicPair(atom1,atom1));
		pairs.push_back(AtomicPair(atom2,atom2));
		
		IntensityMap I(2,2,2);
		I.set_grid(trivial_unit_cell,vec3<double>(1,1,1),vec3<double>(-.5,-.5,-.5));
		
		TS_ASSERT_THROWS_NOTHING(IntnsityCalculator::calculate_scattering_from_pairs(pairs,I,AVERAGE));
	}
	void testIntensityCalculatorAtom()
	{
		double pi=M_PI;
		TS_ASSERT_DELTA(real(IntnsityCalculator::calculate_scattering_from_a_pair_in_a_point_c(1, 1, 1, 1, vec3<double>(0.5,0.5,0.5), //s
																						vec3<double>(0.25,0.25,0.25), //r
																						sym_mat3<double>(0.01,0.02,0.03,0,0,0))),
						-1/1.41421356237*exp(-2*pi*pi*(0.01+0.02+0.03)*0.5*0.5),0.0001); //1.41421356237 is sqrt(2)
	}
		
	
	void testVectorMultiplication()	{
		vec3<double>a(1,1,1);
		TS_ASSERT_DELTA(a*a,3,0.00001);
	}

	void testIntensityMap_size() 	{
		IntensityMap I(2,12,1);
		TS_ASSERT_EQUALS(12,I.size()[1]);
		TS_ASSERT_EQUALS(vec3<int>(2,12,1),I.size());
	}
  void testIntensityMapCopiesDeeply() {
    IntensityMap I1(1,1,1);
    I1.at(0,0,0)=100;
    IntensityMap I2 = I1;
    I2.at(0,0,0)=200;
    IntensityMap I22 = I2;
    
    TS_ASSERT_EQUALS(100,I1.at(0,0,0));
    TS_ASSERT_EQUALS(200,I2.at(0,0,0));
    TS_ASSERT_EQUALS(200,I22.at(0,0,0));
  }
  
	void testAtomicTypeAndCollection() {
		Scatterer* at=AtomicTypeCollection::get(string("C"));
		Scatterer* at2=AtomicTypeCollection::get(string("C1"));
			
		TS_ASSERT_DELTA(at->form_factor_at(0),5.9972,0.0001);
		TS_ASSERT_EQUALS(at,at2);
		
		TS_ASSERT_THROWS_NOTHING(AtomicTypeCollection::calculate_form_factors_of_all_atoms_on_grid(vec3<int>(2,2,2),Grid(trivial_unit_cell,vec3<double>(1,1,1),vec3<double>(-1,-1,-1))));
	}
  
  void testMolecularType() {
    Atom* an_atom      = new Atom(string("C2"),0.5,.5,.5,.5,0,0,0,0,0,0);
    
    AtomicAssembly a_molecule;
    a_molecule.add_chemical_unit(an_atom);
    
    MolecularScatterer scatterer(&a_molecule);
    
    vec3<double> s(1,0,0);
    
    AtomicTypeCollection::update_current_form_factors(s,0);
  
    //sr = 0.5
    TS_ASSERT_DELTA(-5.9972/2,        scatterer.form_factor_at_c(s,0).real(), 0.0001);
    TS_ASSERT_DELTA(0, scatterer.form_factor_at_c(s,0).imag(), 0.0001);
    // Add interface in the input file
    // check against some old calculated dataset
  }
  
  
	void testStripAtomicLabel() {
		TS_ASSERT_EQUALS("C",AtomicTypeCollection::strip_label("C"));
		TS_ASSERT_EQUALS("Cl",AtomicTypeCollection::strip_label("Cl"));
		TS_ASSERT_EQUALS("Cl",AtomicTypeCollection::strip_label("Cl15690"));
		TS_ASSERT_EQUALS("C",AtomicTypeCollection::strip_label("C25"));
	}
	void testAtomicAssembly()
	{
		AtomicAssembly molecule;
			AtomicAssembly* molecule_fragment = new AtomicAssembly();
			molecule_fragment->add_chemical_unit(p_atom1);
		molecule.add_chemical_unit(molecule_fragment);
		molecule.add_chemical_unit(p_atom2);
		
		vector<Atom*> v;
		v.push_back(p_atom1);
		v.push_back(p_atom2);
		
		TS_ASSERT_EQUALS(v,molecule.get_atoms());
	}
	void testAtom()
	{
		Atom anAtom(string("C12"),1,0,0,0,1,1,1,0,0,0);
		
		TS_ASSERT_EQUALS(anAtom.r,vec3<double>(0,0,0));
		TS_ASSERT_DELTA(anAtom.atomic_type->form_factor_at(0),5.9972,0.0001);
		vector<Atom*> v;
		v.push_back(&anAtom);
		TS_ASSERT_EQUALS(v ,anAtom.get_atoms());
	}
  void testAtomUiso()
	{
    cctbx::uctbx::unit_cell cell(af::double6(1,1,7,90,90,120));
    
		Atom anAtom(string("C12"),1, 1,0,0,0, 0.1,  cell.reciprocal_metrical_matrix());
    
		TS_ASSERT_DELTA(0.4/3,anAtom.U[0],0.00001);
    TS_ASSERT_DELTA(0.4/3,anAtom.U[1],0.00001);
    TS_ASSERT_DELTA(0.1/7/7,anAtom.U[2],0.00001);
    TS_ASSERT_DELTA(0.4/6,anAtom.U[3],0.00001);
    TS_ASSERT_DELTA(0,anAtom.U[4],0.00001);
    TS_ASSERT_DELTA(0,anAtom.U[5],0.00001);

	}
	void testAtomicPair()
	{
		AtomicPair a;
		a.r()=vec3<double>(1,1,1);//real
		a.r(true)=vec3<double>(0,0,0);//average
		
		TS_ASSERT_EQUALS(a.r(false),vec3<double>(1,1,1));
		TS_ASSERT_EQUALS(a.r(true),vec3<double>(0,0,0));
	}
	void testAtomAndPair()
	{
    Atom at1(string("C"),0.12,0.5,0,0,0,1,1,1,0,0,0);
    Atom at2(string("C2"),0.7,0.5,1,1,1,1,1,1,0,0,0);
		AtomicPair pair(at1,at2);
		TS_ASSERT_EQUALS(pair.r(),vec3<double>(1,1,1));
		TS_ASSERT_EQUALS(pair.U(),sym_mat3<double>(2,2,2,0,0,0));
		TS_ASSERT_EQUALS(pair.p(),0.25);
		TS_ASSERT_EQUALS(pair.average_r(),vec3<double>(1,1,1));
		TS_ASSERT_EQUALS(pair.average_U(),sym_mat3<double>(2,2,2,0,0,0));
		TS_ASSERT_EQUALS(pair.average_p(),0.25);
		TS_ASSERT_DELTA(pair.atomic_type1->form_factor_at(0),5.9972,0.0001);
		TS_ASSERT_DELTA(pair.multiplier,0.7*0.12,0.00001);
	}
	
	void helperTestPointerVector()
	{
		p_vector<TestDummy> v;
		v.push_back(new TestDummy());
		TS_ASSERT_THROWS_NOTHING(v[0]);
		TS_ASSERT_EQUALS(v.size(),1);
	}
	void testPointerVector()
	{

		TS_ASSERT_THROWS_ANYTHING(helperTestPointerVector()); //TestDummy will throw exception on destruction
	}
	
	void testItIsPossibleToCreateArray()
	{
		UnitCell cell(1,1,1,90,90,90);
		cell.set_laue_symmetry("4mm");
			ChemicalUnitNode* node = new ChemicalUnitNode();
				AtomicAssembly* cu = new AtomicAssembly();
					Atom* atom = new Atom(string("C"),0.5,0,0,0,1,1,1,0,0,0);
				cu->add_chemical_unit(atom);
					atom = new Atom(string("C2"),0.5,0,0,0,1,1,1,0,0,0);
					
				cu->add_chemical_unit(atom);
			node->add_chemical_unit(cu);
		cell.add_node(node);
		
		
		CreateThomasStructure();
		//maybe check that the construction is valid?
	}
	void CreateThomasStructure()
	{/*
		UnitCell cell(3,3,1,90,90,90);
		cell.set_laue_symmetry("4mm");
			ChemicalUnitNode* node = new ChemicalUnitNode();
				AtomicAssembly* cu1 = new AtomicAssembly();
				cu1->set_occupancy(0.5);
					cu1->add_chemical_unit(new Atom("C",0.5,0.54,0.54,0,0.05,0.05,0,0,0,0));
					cu1->add_chemical_unit(new Atom("C",0.5,-0.54,-0.54,0,0.05,0.05,0,0,0,0));
			node->add_chemical_unit(cu1);
				AtomicAssembly* cu2 = new AtomicAssembly();
				cu2->set_occupancy(0.5);
					cu2->add_chemical_unit(new Atom("C",0.5,-0.54,0.54,0,0.05,0.05,0,0,0,0));
					cu2->add_chemical_unit(new Atom("C",0.5,0.54,-0.54,0,0.05,0.05,0,0,0,0));
			node->add_chemical_unit(cu2);
		cell.add_node(node);
		
		const int X=0;
		const int Y=1;
		//Create translational modes
		ADPMode cu1_x = translational_mode(cu1,X);
		ADPMode cu1_y = translational_mode(cu1,Y);
		
		ADPMode cu2_x = translational_mode(cu2,X);
		ADPMode cu2_y = translational_mode(cu2,Y);
		
		AtomicPairPool zeroth_neighbor;
			//Create manually things that should've been generatrd by node
			
		zeroth_neighbor.add_modifier( new SubstitutionalCorrelation(cu1,cu2,0));
		zeroth_neighbor.add_modifier( new SubstitutionalCorrelation(cu1,cu1,1));
		zeroth_neighbor.add_modifier( new SubstitutionalCorrelation(cu2,cu2,1));
		
			//Now we manually find that amplitude of translational modes of both chemical units equals 0.05. This will be our hard constraint
		zeroth_neighbor.add_modifier( new DoubleADPMode(&cu1_x,&cu1_x,0.05) );
		zeroth_neighbor.add_modifier( new DoubleADPMode(&cu1_y,&cu1_y,0.05) );
		zeroth_neighbor.add_modifier( new DoubleADPMode(&cu2_x,&cu2_x,0.05) );
		zeroth_neighbor.add_modifier( new DoubleADPMode(&cu2_y,&cu2_y,0.05) );
		
		AtomicPairPool first_neighbor;
		
		first_neighbor.add_modifier( new CellShifter(1,0,0) );
		
		vector<double> corr;
		corr.push_back(0.572/2);
		first_neighbor.add_modifiers( correlators_from_cuns(node,node,corr));
		*/
		/*
		идеально было бы конечно иметь нечто вроде
		
		UnitCell 3 3 1 90 90 90
			Symmetry 4mm
			Node node1
				ChemicalUnit cu1 p=0.5
					atom C     0.54  0.54 0   0.05 0.05 0 0 0 0
					atom C    -0.54 -0.54 0   0.05 0.05 0 0 0 0
				ChemicalUnit cu2 p=0.5
					atom C    -0.54  0.54 0   0.05 0.05 0 0 0 0
					atom C     0.54 -0.54 0   0.05 0.05 0 0 0 0
					
		Pool zero
			node1 substitutional correlations
			cu1 x mode
			cu1 y mode
			cu2 x mode
			cu2 y mode
			
		Pool one (1 0 0)
			node1 - node1 substitutional correlation p=0.572
		*/
		
		
		
	}
	void testPoolFindsPairs()
	{
		AtomicPairPool aPool;
		aPool.get_pair(&atom1,&atom1).p()=0;
		aPool.get_pair(&atom1,&atom2);
		
		TS_ASSERT_EQUALS(2,aPool.pairs.size());
		
		TS_ASSERT_EQUALS(0,aPool.get_pair(&atom1,&atom1).p());
		TS_ASSERT_EQUALS(2,aPool.pairs.size());
	}
	void testPoolInvokesCorrelators()
	{
		AtomicPairPool aPool;
		aPool.add_modifier(new CellShifter(1,0,0));
		aPool.add_modifier(new SubstitutionalCorrelation(p_atom1,p_atom1,0));
		
		aPool.invoke_correlators();
		
		TS_ASSERT_EQUALS(0,aPool.get_pair(p_atom1,p_atom1).p());
		TS_ASSERT_EQUALS(vec3<double>(1,0,0),aPool.get_pair(p_atom1,p_atom1).r());
		TS_ASSERT_EQUALS(vec3<double>(1,0,0),aPool.get_pair(p_atom1,p_atom1).average_r());
	}
	
	
	void testCellShifter()
	{
		AtomicPairPool aPool;
		
		aPool.get_pair(&atom1,&atom1);
		
		CellShifter shifter(0,1,2);
		
		shifter.modify_pairs(&aPool);
		
		TS_ASSERT_EQUALS(vec3<double>(0,1,2),aPool.get_pair(&atom1,&atom1).average_r());
		TS_ASSERT_EQUALS(vec3<double>(0,1,2),aPool.get_pair(&atom1,&atom1).r());
	}
	
	void testSubstitutionalCorrelation()
	{
		ChemicalUnitNode node;
		
		Atom* atom3 = new Atom("Si",0.5,2,2,2,0,0,0,0,0,0);
			
		AtomicAssembly unit;
		unit.set_occupancy(0.5);
		unit.add_chemical_unit(p_atom1);
		unit.add_chemical_unit(p_atom2);

		AtomicPairPool aPool;
		
		SubstitutionalCorrelation corr(&unit,atom3,0.5);
		
		corr.modify_pairs(&aPool);
		
		TS_ASSERT_EQUALS(0.5,aPool.get_pair(p_atom1,atom3).p());
		TS_ASSERT_EQUALS(0.25,aPool.get_pair(p_atom1,atom3).average_p());
	}
	
  void testADPMode()
  {
    ADPMode mode1,mode2;
    
    TS_ASSERT_EQUALS(mode1,mode2);
    mode1.add_atom(&atom1, vec3<double>(1,2,3));
    TS_ASSERT(mode1!=mode2);
    mode2.add_atom(&atom1, vec3<double>(1,2,3));
    TS_ASSERT_EQUALS(mode1,mode2);
    mode1.add_atom(&atom2, vec3<double>(1,2,3));
    mode2.add_atom(&atom1, vec3<double>(1,2,3));
    TS_ASSERT(mode1!=mode2);
    mode2.atomic_displacements[1].atom=&atom2;
    TS_ASSERT_EQUALS(mode1,mode2);
    mode2.atomic_displacements[1].displacement_vector = vec3<double>(1,3,2);
    TS_ASSERT(mode1!=mode2);
  }
  
	void testDoubleADPMode()
	{		
		ADPMode mode1,mode2;
		
		mode1.add_atom(p_atom1,vec3<double>(1,0,0));
		
		AtomicPairPool aPool;
		
		DoubleADPMode(&mode1,&mode1,1).modify_pairs(&aPool);
		
		TS_ASSERT_EQUALS(0,aPool.get_pair(p_atom1,p_atom1).U()[0]);
    
    //test equals operator
    DoubleADPMode dmode1(&mode1,&mode1,1);
    DoubleADPMode dmode2(&mode1,&mode2,1);
    DoubleADPMode dmode3(&mode1,&mode1,2);
    TS_ASSERT_EQUALS(dmode1,dmode1);
    TS_ASSERT(!(dmode1==dmode2));
    TS_ASSERT(!(dmode1==dmode3));
	}
	
	void testZeroVectorCorrelation()
	{
		ZeroVectorCorrelation zcor1;
		AtomicPairPool aPool;
		
		aPool.get_pair(p_atom1,p_atom1);
		aPool.get_pair(p_atom1,p_atom2);
		
		zcor1.modify_pairs(&aPool);
		
		sym_mat3<double> zero_ADP(0,0,0,0,0,0);
		sym_mat3<double> nozero_ADP(2,2,2,0,0,0);
		
		TS_ASSERT_EQUALS(zero_ADP,aPool.get_pair(p_atom1,p_atom1).U());
		TS_ASSERT_EQUALS(nozero_ADP,aPool.get_pair(p_atom1,p_atom1).average_U());
		TS_ASSERT_EQUALS(nozero_ADP,aPool.get_pair(p_atom1,p_atom2).U());
		
		aPool = AtomicPairPool();
		
		aPool.get_pair(p_atom1,p_atom1);
		aPool.get_pair(p_atom1,p_atom2);
		
		CellShifter(1,0,0).modify_pairs(&aPool);
		zcor1.modify_pairs(&aPool);
		
		TS_ASSERT_EQUALS(nozero_ADP,aPool.get_pair(p_atom1,p_atom1).U());
	}
	void testMultiplicityCorrelation()
	{
		AtomicPairPool aPool;
		
		double p=aPool.get_pair(p_atom1,p_atom1).p();
		
		double m=15.6;
		
		MultiplicityCorrelation mcor(m);
		
		mcor.modify_pairs(&aPool);
		
		TS_ASSERT_DELTA(p,aPool.get_pair(p_atom1,p_atom1).p(),0.00001);
		TS_ASSERT_DELTA(p,aPool.get_pair(p_atom1,p_atom1).average_p(),0.00001);
    TS_ASSERT_DELTA(m,aPool.get_pair(p_atom1,p_atom1).multiplier,0.00001);

	}
  

	void testTranslationalMode()
	{
		
    UnitCell cell(10,10,10,90,90,90);
    
		AtomicAssembly* cu1 = new AtomicAssembly();
		cu1->add_chemical_unit(p_atom1);
		cu1->add_chemical_unit(p_atom2);
		
		ADPMode* cu1_x = translational_mode(cu1,0,cell.cell.metrical_matrix());
		
		TS_ASSERT(almost_equal(vec3<double>(0.1,0,0),cu1_x->atomic_displacements[0].displacement_vector));
    delete cu1_x;
	}
	void test_rotational_mode()
	{
    UnitCell cell(2,2,3,90,90,90);
    
		AtomicAssembly* cu1 = new AtomicAssembly();
		cu1->add_chemical_unit(p_atom1);//0 0 0
		cu1->add_chemical_unit(p_atom2);//.5 .5 .5
		
		ADPMode* cu1_z = rot_mode(cu1,vec3<double>(0,0,1),vec3<double>(0,0,0),cell.cell.metrical_matrix());//z - rotation
		
		TS_ASSERT_EQUALS(vec3<double>(0,0,0),cu1_z->atomic_displacements[0].displacement_vector);
		TS_ASSERT(almost_equal(vec3<double>(-0.5,0.5,0),cu1_z->atomic_displacements[1].displacement_vector));
		
		
		ADPMode* cu1_hex_z=rot_mode(cu1,vec3<double>(0,0,1),vec3<double>(0,0,0),sym_mat3<double>(1,1,1,0.5,0,0));//z - rotation
		TS_ASSERT_DELTA(-0.2887,cu1_hex_z->atomic_displacements[1].displacement_vector[0],0.0001);
		TS_ASSERT_DELTA(0.2887,cu1_hex_z->atomic_displacements[1].displacement_vector[1],0.0001);
	}
	
	void testCorrelationFromCUNS()
	{
		p_atom2->occupancy=0.3;
		Atom* p_atom3 = new Atom("Si",0.2,2,2,2,0,0,0,0,0,0);
		
		Atom* p_atom21 = new Atom(*p_atom3);
		Atom* p_atom22 = new Atom(*p_atom2);
		Atom* p_atom23 = new Atom(*p_atom1);
		
		ChemicalUnitNode node1;
		
		node1.add_chemical_unit(p_atom1);
		node1.add_chemical_unit(p_atom2);
		node1.add_chemical_unit(p_atom3);
		
		ChemicalUnitNode node2;
		
		node2.add_chemical_unit(p_atom21);
		node2.add_chemical_unit(p_atom22);
		node2.add_chemical_unit(p_atom23);
		
		
		vector<double> corr;
		corr.push_back(0.1);
		corr.push_back(0.1);
		corr.push_back(0.2);
		corr.push_back(0.0);
		
		// Probabilities will be:
		//     0.5 0.3 0.2
		// 0.2 0.1 0.1 0.0
		// 0.3 0.2 0.0 0.1
		// 0.5 0.2 0.2 0.1
		//
		// in an array in row-coloumn order
		
		
		vector<SubstitutionalCorrelation*> correlations = correlators_from_cuns(&node1,&node2,corr);
	
		TS_ASSERT_EQUALS(9,correlations.size());
		SubstitutionalCorrelation corr32(p_atom3, p_atom22, 0.3);
		SubstitutionalCorrelation corr33(p_atom3, p_atom23, 0.0);
		
		TS_ASSERT_EQUALS(0.0,correlations[2]->joint_probability);
		TS_ASSERT_EQUALS(p_atom3,correlations[2]->chemical_units[0]);
		TS_ASSERT_EQUALS(p_atom21,correlations[2]->chemical_units[1]);
		
		TS_ASSERT_DELTA(0.1,correlations[5]->joint_probability,0.000001);
		TS_ASSERT_EQUALS(p_atom3,correlations[5]->chemical_units[0]);
		TS_ASSERT_EQUALS(p_atom22,correlations[5]->chemical_units[1]);		
		
		TS_ASSERT_DELTA(0.1,correlations[8]->joint_probability,0.000001);
		TS_ASSERT_EQUALS(p_atom3,correlations[8]->chemical_units[0]);
		TS_ASSERT_EQUALS(p_atom23,correlations[8]->chemical_units[1]);	
		
		for(int i=0; i<correlations.size(); i++)
			delete correlations[i];
	}
	
	void testOuterProduct()
	{
		TS_ASSERT_EQUALS(sym_mat3<double>(0,1,6,0,0,3), outer_product(vec3<double>(0,1,2),vec3<double>(0,1,3)));
	}	
	void testStaticShift(){

    
    
		StaticShift shift;
		ADPMode* atom1_x=translational_mode(&atom1,0,unit_matrix);
		ADPMode* atom2_x=translational_mode(&atom2,0,unit_matrix);
		ADPMode* atom1_y=translational_mode(&atom1,1,unit_matrix);
		ADPMode* atom2_y=translational_mode(&atom2,1,unit_matrix);

		shift.add_displacement(VectorStart,atom1_x,0.25);
		shift.add_displacement(VectorEnd,atom2_x,-0.25);
		shift.add_displacement(VectorStart,atom1_y,0.25);
		shift.add_displacement(VectorEnd,atom2_y,-0.25);

		AtomicPairPool pool;
		
		shift.modify_pairs(&pool);
		
		TS_ASSERT_DELTA(0,pool.get_pair(&atom1,&atom2).r()[0],0.0001);
		TS_ASSERT_DELTA(0,pool.get_pair(&atom1,&atom2).r()[1],0.0001);
	}
  void testSizeEffect()
  {
    ADPMode* atom2_x=translational_mode(&atom2,0,unit_matrix);
    SizeEffect size_effect_cu_adp(&atom1,atom2_x,0.5);//initially distance is (0.5,0.5,0.5) -> (1,0.5,0.5)
    SizeEffect size_effect_adp_cu(atom2_x,&atom1,0.5);// dist is -(0.5,0.5,0.5) -> -(1,0.5,0.5)
    
    AtomicPairPool pool;
    
    
		size_effect_cu_adp.modify_pairs(&pool);
    size_effect_adp_cu.modify_pairs(&pool);
    
    
    TS_ASSERT(almost_equal( vec3<double>(1,0.5,0.5) , pool.get_pair(&atom1,&atom2).r() ));
    TS_ASSERT(almost_equal( vec3<double>(-1,-0.5,-0.5) , pool.get_pair(&atom2,&atom1).r() ));
  }
    
  

  ///TODO: make this test
	void testGridReciprocal()
	{
		//still too lazy to write this test
		//Grid grid(trivial_unit_cell,vec3<double>(1,1,1),vec3<double>(-.5,-.5,-.5)));
	}
	void test_grid_and_residual()
	{
		vec3<double> r_res;
		vec3<int> r_grid;
		
		
		grid_and_residual(vec3<double>(2.1,3.8,-0.51),
						  Grid(trivial_unit_cell,vec3<double>(1,1,0.5),vec3<double>(0,0,0)),
						  r_grid,r_res);
		
		TS_ASSERT_DELTA(0.1,r_res[0],0.00001);
		TS_ASSERT_DELTA(-0.2,r_res[1],0.00001);
		TS_ASSERT_DELTA(-0.01,r_res[2],0.00001);
		TS_ASSERT_EQUALS(vec3<int>(2,4,-1),r_grid);
		
		grid_and_residual(vec3<double>(2.1,3.8,-0.51),
						  Grid(trivial_unit_cell,vec3<double>(1,1,0.5),vec3<double>(-1,-1,-1)),
						  r_grid,r_res);
		
		TS_ASSERT_EQUALS(vec3<int>(3,5,1),r_grid);
	}
	
	void test_add_pair_to_appropriate_place()
	{
		IntensityMap P(20,20,1);
		P.set_grid(trivial_unit_cell,vec3<double>(1,1,1),vec3<double>(-10,-10,0));
		
		P.at(14,12,0)=0;
		
		IntensityMap part_of_P(4,4,1);
		part_of_P.set_grid(trivial_unit_cell,vec3<double>(1,1,1),vec3<double>(-2,-2,0));
		part_of_P.at(0,0,0)=150;
		
		add_pair_to_appropriate_place(part_of_P,P,vec3<int>(16,14,0),vector<bool> (3,false));
		
		TS_ASSERT_EQUALS(part_of_P.at(0,0,0),P.at(14,12,0));
		
		TS_ASSERT_THROWS_NOTHING(add_pair_to_appropriate_place(part_of_P,P,vec3<int>(20000,20000,0),vector<bool> (3,false)));
	}
  
  void test_add_pair_to_appropriate_place_with_periodic_boundary_conditions()
	{
		IntensityMap P(2,2,1);
		P.set_grid(trivial_unit_cell,vec3<double>(1,1,1),vec3<double>(-1,-1,0));
		
		P.erase_data();
		
		IntensityMap part_of_P(2,2,1);
		part_of_P.set_grid(trivial_unit_cell,vec3<double>(1,1,1),vec3<double>(-1,-1,0));
    // 1 2 
    // 3 4
		part_of_P.at(0,0,0)=1.0;
    part_of_P.at(0,1,0)=2.0;
    part_of_P.at(1,0,0)=3.0;
    part_of_P.at(1,1,0)=4;
		
    // 0 0
    // 0 1
    P.erase_data();
		add_pair_to_appropriate_place(part_of_P,P,vec3<int>(2,2,0),vector<bool>(3,false));
		TS_ASSERT_EQUALS(0,P.at(0,0,0));
		TS_ASSERT_EQUALS(0,P.at(0,1,0));
		TS_ASSERT_EQUALS(0,P.at(1,0,0));
		TS_ASSERT_EQUALS(1.0,P.at(1,1,0));
    P.erase_data();

    
    // 4 3
    // 2 1
		add_pair_to_appropriate_place(part_of_P,P,vec3<int>(2,2,0),vector<bool>(3,true));
		TS_ASSERT_EQUALS(4,P.at(0,0,0));
		TS_ASSERT_EQUALS(3,P.at(0,1,0));
		TS_ASSERT_EQUALS(2,P.at(1,0,0));
		TS_ASSERT_EQUALS(1,P.at(1,1,0));
    
    vector<bool> only_y(3,false);
    only_y.at(1)=true;
    
    P.erase_data();
    
    // 0 0
    // 2 1
		add_pair_to_appropriate_place(part_of_P,P,vec3<int>(2,2,0),only_y);
		TS_ASSERT_EQUALS(0,P.at(0,0,0));
		TS_ASSERT_EQUALS(0,P.at(0,1,0));
		TS_ASSERT_EQUALS(2.0,P.at(1,0,0));
		TS_ASSERT_EQUALS(1.0,P.at(1,1,0));
		}
	
	void test_misc()
	{
		TS_ASSERT_EQUALS(0,1/2);
		TS_ASSERT_EQUALS(vec3<int>(0,1,1),vec3<int>(1,2,3)/2);
		TS_ASSERT_EQUALS(0.5-0.5,0);
	}

	void test_unique()
	{
		const int SIZE = 10;
		int a1[ SIZE ] = { 1, 3, 5, 7, 7, 9, 1, 3, 5, 7 };
		std::vector< int > inp( a1, a1 + SIZE ); // copy of a
		std::vector< int > res( a1, a1 + 4 ); // Unique values of a
		res.push_back(9);
		
		TS_ASSERT_EQUALS(res,unique_elements(inp));
	}
	
	void test_apply_symmetry_element_to_atomic_assembly()
	{
		AtomicAssembly group;
		AtomicAssembly* subgroup =  new AtomicAssembly();
		subgroup->add_chemical_unit(p_atom1);
		subgroup->add_chemical_unit(p_atom2);
		group.add_chemical_unit(subgroup);
		
		ChemicalUnit* symmetry_equivalent = group.create_symmetric(mat3<double>(1,0,0,0,1,0,0,0,1),vec3<double>(0,0,0));
		TS_ASSERT_EQUALS(2,symmetry_equivalent->get_atoms().size());
		delete symmetry_equivalent;
		
    //TODO: make the test
		//TS_FAIL("check occupancy of AtomicAssembly");
		
		Atom atom = Atom(string("C"),0.5,1,2,3,1,2,3,4,5,6);
		Atom* sym_atom = atom.create_symmetric(mat3<double>(-1,0,0,0,1,0,0,0,1), vec3<double>(0,0,100));
		TS_ASSERT_EQUALS(vec3<double>(-1,2,103),sym_atom->r);
		TS_ASSERT_EQUALS(sym_mat3<double>(1,2,3,-4,-5,6),sym_atom->U);
		delete sym_atom;
		
	}
	
	void test_minimizer()
	{
		Minimizer a_minimizer;
		
		vector<double> params(10,0);
		for(int i=0; i<10; i++)
			params[i]=(double)(10-i)/10;

    TestMinimizerCalculator calc;
    
    IntensityMap exp_data(10,1,1);
    
    for(int i=0; i<10; i++)
      exp_data.at(i)=7;
    
    OptionalIntensityMap weights;
    
    report.shut_up();
		params = a_minimizer.minimize(params, &exp_data, &calc, &weights);
    
    TS_ASSERT_DELTA(params[0],7,0.001);
	}
	
  void test_hdf5_reader()
  {
    IntensityMap map_123 = ReadHDF5(PATH_TO_YELL_SRC "/test_files/123.h5");
    
    TS_ASSERT_DELTA(1,map_123.at(0),0.01);
    TS_ASSERT_DELTA(2,map_123.at(1),0.01);
    TS_ASSERT_DELTA(3,map_123.at(2),0.01);
  }
  
  void testOptionalIntensityMap()
  {
    OptionalIntensityMap map(138690);
    
    TS_ASSERT_DELTA(138690, map.at(31), 0.01);
    
    map.load_data(PATH_TO_YELL_SRC "/test_files/123.h5");
    
    TS_ASSERT_DELTA(1,map.at(0),0.01);
    TS_ASSERT_DELTA(2,map.at(1),0.01);
    TS_ASSERT_DELTA(3,map.at(2),0.01);
    
    //test that unloaded functions do not break
    OptionalIntensityMap unloaded_map(0);
    
    TS_ASSERT_DELTA(0, unloaded_map.at(312), 0.01);
    
    unloaded_map.load_data("file_does_not_exist.h5");
    
    TS_ASSERT_DELTA(0,unloaded_map.at(0),0.01);
    TS_ASSERT_DELTA(0,unloaded_map.at(1),0.01);
  }
  
  void test_file_exists()
  {
    TS_ASSERT_EQUALS(true, file_exists(PATH_TO_YELL_SRC "/test_files/123.h5"));
    TS_ASSERT_EQUALS(false, file_exists("does.not.exist"));
  }
  
  void test_R_factor()
  {
    Model a_model;
    IntensityMap I(2,1,1);
    I.at(0)=2;
    I.at(1)=10;
    
    IntensityMap I0(2,1,1);
    I0.at(0)=0;
    I0.at(1)=0;

    IntensityMap I_half(2,1,1);
    I_half.at(0)=1;
    I_half.at(1)=5;

    a_model.intensity_map = I;
    a_model.average_intensity_map = I0;
    
    TS_ASSERT_DELTA(0,a_model.R_factor(I,R1,UNWEIGHTED), 0.001);
    TS_ASSERT_DELTA(1, a_model.R_factor(I_half,R1,UNWEIGHTED),0.001); 

  }
  
  void test_padded_intensity_map()
  {
    IntensityMap I(2,1,1);
    
    IntensityMap I_padded(4,1,1);
    for(int i=0;i<4;++i)
      I_padded.at(i)=i;

    I.copy_from_padded(vec3<int> (1,0,0),I_padded);
    TS_ASSERT_EQUALS(1,I.at(0));
    TS_ASSERT_EQUALS(2,I.at(1)); 
  }
  
  void test_create_padded_map()
  {
    Grid a_grid (trivial_unit_cell,vec3<double>(0.1,0.1,0.1), vec3<double>(-1,-2,-3),vec3<int>(10,10,10));
    Grid padded_grid = a_grid.pad(vec3<int>(3,2,1));
    TS_ASSERT_EQUALS(vec3<double>(-1.3,-2.2,-3.1),padded_grid.lower_limits);
    TS_ASSERT_EQUALS(vec3<int>(16,14,12),padded_grid.grid_size);
  }
  
  // the following functions are used to test qi parsers
  template<class Attr,class parser_type,class skipper_type>
  bool test_parser(Attr expected, string input,parser_type& a_parser,skipper_type skipper)
  {
    iterator_ start = input.begin();
    iterator_ end = input.end();
    Attr result;
    return qi::phrase_parse(start,end,a_parser,skipper,result) && start==end && result == expected ;
  }

  template<class parser_type,class Attr>
  bool test_parser(Attr expected,string input,parser_type& a_parser)
  {
    iterator_ start = input.begin();
    iterator_ end = input.end();
    Attr result;
    return qi::parse(start,end,a_parser,result) && start==end && result == expected ;
  }
  template<class parser_type> 
  bool test_parser(double expected,string input,parser_type& a_parser)
  {
    iterator_ start = input.begin();
    iterator_ end = input.end();
    double result;
    return qi::parse(start,end,a_parser,result) && start==end && almost_equal(result, expected) ;
  }
  
  template<class parser_type,class skipper_type>
  bool test_parser_nores(string input,parser_type& a_parser,skipper_type skipper)
  {
    iterator_ start = input.begin();
    iterator_ end = input.end();

    return qi::phrase_parse(start,end,a_parser,skipper) && start==end;
  }
  template<class parser_type>
  bool test_parser_nores(string input,parser_type& a_parser)
  {
    iterator_ start = input.begin();
    iterator_ end = input.end();
    return qi::parse(start,end,a_parser) && start==end ;
  }
  
  
  template<class parser_type,class skipper_type,class result_class>
  bool run_parser(string input,parser_type a_parser, skipper_type skipper,result_class& result )
  {
    iterator_ start = input.begin();
    iterator_ end = input.end();
    
    bool r = qi::phrase_parse(start,end,a_parser,skipper,result);
    
    TS_ASSERT(r);
    TS_ASSERT(start==end);
    
    return r && (start==end);
  }
  template<class parser_type,class result_class>
  bool run_parser(string input,parser_type a_parser, result_class& result )
  {
    iterator_ start = input.begin();
    iterator_ end = input.end();
    
    bool r = qi::parse(start,end,a_parser,result);
    
    TS_ASSERT(r);
    TS_ASSERT(start==end);
    return r && (start==end);
  }
  
  void test_atom_parser()    {
    Atom * atom;

    a_parser.add_model(new Model());
    test_parser_nores("Cell 1 1 1  90. 90. 90.",a_parser.program_options,a_skipper);
    run_parser("C1  0.3 -0.0107 0.8970 0.1319 0.066 0.066 0.066 -0.033 0 0",a_parser.atom,a_skipper,atom);
    
    TS_ASSERT_DELTA(0.3,atom->multiplier,0.01);
    TS_ASSERT_DELTA(-0.0107,atom->r[0],0.01);
    TS_ASSERT_DELTA(0.066,atom->U[0],0.01);
    delete atom;
    
    //check how Uiso parser works
    run_parser("C1 0.5 0 0 0  0.1",a_parser.atom,a_skipper,atom);
    
    TS_ASSERT_EQUALS(Atom("C1",0.5,1,0,0,0, 0.1,0.1,0.1,0,0,0),*atom);

    delete atom;
  }
  
  void test_molecular_form_factor() {
    Scatterer *scatterer, *not_initialized_scatterer;
    
    a_parser.add_model(new Model());
    test_parser_nores("Cell 1 1 1  90. 90. 90.",a_parser.program_options,a_skipper);
    
    run_parser("[C1 0.5 0 0 0  0.1\n C1 0.5 0 0 0  0.1]",a_parser.molecular_scatterer,a_skipper,scatterer);
    
    TS_ASSERT_DIFFERS(scatterer, not_initialized_scatterer);
    
    vector<boost::tuple<string,Scatterer*> > scatterers;
    run_parser("MolecularScatterers[ C_molecule=[C1 0.5 0 0 0  0.1\n C1 0.5 0 0 0  0.1] Nop=[] ]",a_parser.molecular_scatterers,a_skipper,scatterers);
    
    TS_ASSERT_EQUALS(2,scatterers.size());
    TS_ASSERT_EQUALS(boost::get<0>(scatterers[0]), "C_molecule");
    TS_ASSERT_DIFFERS(boost::get<1>(scatterers[0]), not_initialized_scatterer);
    TS_ASSERT_EQUALS(boost::get<1>(scatterers[0]), AtomicTypeCollection::get("C_molecule")); //it is done by molecular_scatterers do register it
  }
  
  void test_atomic_assembly_parser()  {
    
    AtomicAssembly* assembly;
    run_parser("[ C1  21.0 0 0 0 0 0 0 0 0 0 C2  0 10 0 0 0 0 0 0 0 0 ]",a_parser.atomic_assembly,a_skipper,assembly);
    
    TS_ASSERT_DELTA(21.0,(*assembly)[0]->get_atoms()[0]->multiplier,0.01);
    TS_ASSERT_DELTA(10.0,(*assembly)[1]->get_atoms()[0]->r[0],0.01);
    delete assembly;
  }
  void test_chemical_unit_parser()  {
    
    ChemicalUnit * unit;
    
    run_parser("C1  0.1238 -0.0107 0.8970 0.1319 0.066 0.066 0.066 -0.033 0 0",a_parser.chemical_unit,a_skipper,unit);
    
    TS_ASSERT_DELTA(0.1238,unit->get_atoms()[0]->multiplier,0.0001);//should wrap the atom
    
    delete unit;
    
    run_parser("[ C1  0.3215 -0.0107 0.8970 0.1319 0.066 0.066 0.066 -0.033 0 0 ]",a_parser.chemical_unit,a_skipper,unit);
    
    TS_ASSERT_DELTA(0.3215,unit->get_atoms()[0]->multiplier,0.0001);//should wrap the atom atomic assembly
    
    delete unit;
  }
  
  void test_identifier_assignement()  {
    
    ChemicalUnit* unit;
    StructurePartRef assigned_unit;
    run_parser("an_atom1 =[ C1  0.3215 -0.0107 0.8970 0.1319 0.066 0.066 0.066 -0.033 0 0 ]",a_parser.chemical_unit,a_skipper,unit);
    run_parser("an_atom1",a_parser.identifier,assigned_unit);
    
    TS_ASSERT_EQUALS(StructurePartRef(unit),assigned_unit);
    
    delete unit;
  }
 
  void test_variant()  {
    
    ChemicalUnitNode * a_variant;
    
   run_parser("Variant[ (p=0.9876) an_atom1 =[ C1  0.3215 -0.0107 0.8970 0.1319 0.066 0.066 0.066 -0.033 0 0 ] (p=1-0.9876)Void]",a_parser.variant,a_skipper,a_variant);
    TS_ASSERT_DELTA(0.9876,a_variant->chemical_units[0].get_atoms()[0]->occupancy,0.0001);
      
  }
  
  void test_variant_assignement()  {
    
    ChemicalUnitNode * a_variant;
    StructurePartRef assigned_unit;
    
    run_parser("our_variant = Variant[ (p=0.9876) an_atom1 =[ C1  0.3215 -0.0107 0.8970 0.1319 0.066 0.066 0.066 -0.033 0 0 ] (p=1-0.9876)[]]",a_parser.variant_assignement,a_skipper,a_variant);
    run_parser("our_variant",a_parser.identifier,assigned_unit);
    TS_ASSERT_EQUALS(StructurePartRef(a_variant),assigned_unit);
  }
  
  void test_symmetry_parser()  {
    SymmetryElement unity(mat3<double>(1,0,0,0,1,0,0,0,1),vec3<double>(0,0,0));
    TS_ASSERT_EQUALS(SymmetryElement(),unity); //check that we implemented operator==
    TS_ASSERT_DIFFERS(SymmetryElement(mat3<double>(.9,0,0,0,1,0,0,0,1),vec3<double>(0,0,0)),unity);
     
    TS_ASSERT(test_parser(0,"x",a_parser.basis_vector));
    
    double r[] = {1,-1,1,0.5};
    vector<double> correct;
    correct.insert(correct.end(),r,r+4);
    
    vector<double> res;
    
    run_parser("x-y+z+1/2",a_parser.permutation_component,a_skipper,res);
    
    TS_ASSERT_EQUALS(correct,res);
    // -x; +x; x-y-1/2;x+1/2;

    
    TS_ASSERT(test_parser(unity, "Symmetry(x,y,z)",a_parser.symmetry_element,a_skipper));//unity
    mat3<double> xyz(1,0,0,0,1,0,0,0,1);
    TS_ASSERT(test_parser(SymmetryElement(xyz,vec3<double>(1,0,0)), "Symmetry(x+1,y,z)",a_parser.symmetry_element,a_skipper));
    TS_ASSERT(test_parser(SymmetryElement(xyz,vec3<double>(3,0,0)), "Symmetry(3+x,y,z)",a_parser.symmetry_element,a_skipper));
    TS_ASSERT(test_parser(SymmetryElement(xyz,vec3<double>(3,125.0/6,0)), "Symmetry(3+x,125/6+y,z)",a_parser.symmetry_element,a_skipper));
  }
  
  void test_symmetric_atomic_assembly()  {
    ChemicalUnit * unchanged_atom;
    run_parser("some_atom = C 0.123 1 2 3 11 22 33 12 13 23",a_parser.chemical_unit,a_skipper,unchanged_atom);
    ChemicalUnit * changed_atom;
    run_parser("some_atom*Symmetry(y,x,z)",a_parser.chemical_unit,a_skipper,changed_atom);
    
    // some_atom*Symmetry[y,x,z] == C 0.123 2 1 3 22 11 33 12 23 13
    TS_ASSERT_EQUALS(Atom("C",0.123,1,2,1,3,22,11,33,12,23,13),*(changed_atom->get_atoms()[0]));
    

    TS_ASSERT(test_parser_nores("Void",a_parser.chemical_unit,a_skipper));
              
    
    
    delete unchanged_atom;
    delete changed_atom;
  }
  
  void test_expression_calculator()  {
    FormulaParser formula_parser;

    TS_ASSERT(test_parser(-12.0,"-12",formula_parser));
    TS_ASSERT(test_parser(4.0,"2*2",formula_parser));
    TS_ASSERT(test_parser(4.0,"2+2",formula_parser));
    TS_ASSERT(test_parser(5.0,"1+2*2",formula_parser));
    TS_ASSERT(test_parser(6.0,"(1+2)*2",formula_parser));
    TS_ASSERT(test_parser(-12.8749,"ll=-12.8749",formula_parser));
    TS_ASSERT(test_parser(-12.8749,"ll",formula_parser)); 
    TS_ASSERT(test_parser(12.8749,"-ll",formula_parser)); 
    TS_ASSERT(test_parser(-12.8749*2,"2*ll",formula_parser));
    
    TS_ASSERT(!test_parser_nores("",formula_parser));
    
  }
  void test_expression_calculator_variable_redefinition() {
    FormulaParser formula_parser;
    
    TS_ASSERT(test_parser(2,"a=2",formula_parser));
    TS_ASSERT(test_parser(2,"a",formula_parser));
    TS_ASSERT(test_parser(5.2,"a=5.2",formula_parser));
    TS_ASSERT(!test_parser(2,"a",formula_parser));
    TS_ASSERT(test_parser(5.2,"a",formula_parser));
    
    
  }
  void test_special_funcitons_in_expression_calculator() {
    FormulaParser formula_parser;
    
    //Special funcitons. Tests generated with the following matlab code:
    /*functions = {'log','sin','cos','exp','sqrt','abs'};
     for i=1:length(functions)
     res = eval([functions{i} '(2)']);
     fprintf('TS_ASSERT(test_parser(%.10f,"%s(2)",formula_parser));\n',res,functions{i});
     end*/
    TS_ASSERT(test_parser( 0.6931471806,"log(2)",formula_parser));
    TS_ASSERT(test_parser( 0.9092974268,"sin(2)",formula_parser));
    TS_ASSERT(test_parser(-0.4161468365,"cos(2)",formula_parser));
    TS_ASSERT(test_parser( 7.3890560989,"exp(2)",formula_parser));
    TS_ASSERT(test_parser( 1.4142135624,"sqrt(2)",formula_parser));
    TS_ASSERT(test_parser( 2.0000000000,"abs(2)",formula_parser));
    
    TS_ASSERT(test_parser(1.5,"mod(1.5,2)",formula_parser));
    TS_ASSERT(test_parser(9,"pow(3,2)",formula_parser));
  }
  

  void test_crystal_parameters_parser()  {
    Model a_model;
    a_parser.add_model(&a_model);
    
    /*
     Cell 14.115 14.115 6.93 90 90 90
     ll=-12.8749
     gs=-ll*2/104
     DiffuseScatteringGrid ll ll -6 gs gs 1 12 12 12
     */
    
    TS_ASSERT(test_parser_nores("Cell 14.115 14.115 6.93 90 90 90\n  \nDiffuseScatteringGrid -6 -6 -6 1 1 1 12 12 12",a_parser.program_options,a_skipper)); 
    TS_ASSERT(UnitCell(14.115,14.115,6.93,90,90,90) == a_model.cell);
    TS_ASSERT_EQUALS(Grid(UnitCell(14.115,14.115,6.93,90,90,90).cell,vec3<double>(1,1,1),vec3<double>(-6,-6,-6),vec3<int>(12,12,12)),a_model.intensity_map.grid);
    TS_ASSERT(test_parser_nores("Cell 14.115 14.115 6.93 90 90 90 gs=1; ll=-6; DiffuseScatteringGrid ll ll ll gs gs gs 12 12 12",a_parser.program_options,a_skipper)); //gs=1
    TS_ASSERT(UnitCell(14.115,14.115,6.93,90,90,90) == a_model.cell);
    TS_ASSERT_EQUALS(Grid(UnitCell(14.115,14.115,6.93,90,90,90).cell,vec3<double>(1,1,1),vec3<double>(-6,-6,-6),vec3<int>(12,12,12)),a_model.intensity_map.grid); 
  }
  void test_unit_cell_parser()  {
    Model a_model;
    a_model.cell = UnitCell(1,2,3,90,90,90);
    a_parser.add_model(&a_model);
    
    /*
     UnitCell
     [
     Var = Variant[
     (p=0.5)
     C1  0.5 -0.0107 0.8970 0.1319 0.066 0.066 0.066 -0.033 0 0
     ]
     ]
     */
    
    TS_ASSERT_EQUALS(0,a_model.cell.chemical_unit_nodes.size()); //start with empty model
    TS_ASSERT(test_parser_nores("     UnitCell \
                                [ \
                                Var = Variant[ \
                                (p=0.5) \
                                C1  0.5 -0.0107 0.8970 0.1319 0.066 0.066 0.066 -0.033 0 0 \
                                (p=0.5) \
                                [] \
                                ] \
                                ]",
                                a_parser.unit_cell,a_skipper));
    TS_ASSERT_EQUALS(1,a_model.cell.chemical_unit_nodes.size()); //check that we now have one unit node
    TS_ASSERT_EQUALS(Atom("C1", 0.5,-0.0107,0.8970,0.1319,0.066,0.066/4,0.066/9,-0.033/2,0,0),*a_model.cell.chemical_unit_nodes[0].chemical_units[0].get_atoms()[0]);
   }

  void test_substitutional_correlation_parser()  {
    ChemicalUnitNode* var;
    vector<SubstitutionalCorrelation*> corrs,expected_corrs;
    run_parser("a_variant = Variant[(p=0.5) C1  0.5 -0.0107 0.8970 0.1319 0.066 0.066 0.066 -0.033 0 0 (p=0.5)[] ]",a_parser.variant_assignement,a_skipper,var);
    run_parser("SubstitutionalCorrelation(a_variant,a_variant,0.5)",a_parser.substitutional_correlation,a_skipper,corrs);
    expected_corrs = correlators_from_cuns(var,var,vector<double>(1,0.5));
    for(int i=0; i<expected_corrs.size(); i++)
    {
      TS_ASSERT_EQUALS(*expected_corrs[i],*corrs[i]);
      delete expected_corrs[i];
      delete corrs[i];
    }
    delete var;
  }
  
  void test_multiplicity_correlation_parser()  {
    MultiplicityCorrelation* mlt;
    run_parser("Multiplicity 0.5",a_parser.multiplicity_correlation,a_skipper,mlt);
    TS_ASSERT_EQUALS(MultiplicityCorrelation(0.5),*mlt);
    delete mlt;
  }
  
  void test_cell_shifter_parser()  {
    CellShifter* cell_shifter;
    run_parser("(0,1,2)",a_parser.cell_shifter,a_skipper,cell_shifter);
    TS_ASSERT_EQUALS(CellShifter(0,1,2),*cell_shifter);
    delete cell_shifter;
  }
  
  void test_pool_parser()  {
    AtomicPairPool* pool;
    run_parser("[(0,1,2) Multiplicity 0.5]",a_parser.atomic_pair_pool,a_skipper,pool);
    // todo test that multiplicity correlations also work
    TS_ASSERT_EQUALS(2,pool->modifiers.size());
    delete pool;
  }
  
  void test_correlations_parser()  {
    //todo change behavior of run_parser so that it returns true if parsed and false if not, warnings showing what's wrong
    vector<AtomicPairPool*> pools;
    run_parser("Correlations[ [(1,2,3)] [Multiplicity 0.5 (3,2,1)] ]",a_parser.correlations,a_skipper,pools);
    TS_ASSERT_EQUALS(2,pools.size());
    for(int i=0; i<pools.size(); i++)
      delete pools[i];
  }
  
  void test_input_file_parser_grammar()  {
    Model a_model;
    a_parser.add_model(&a_model);

    const char* crystal_params =     
    "Cell 14.115 14.115 6.93 90 90 90 \n\
    ll=-12.8749; \n\
    gs=-ll*2/104; \n\
    DiffuseScatteringGrid ll ll -6 gs gs 1 104 104 12\n\
    LaueSymmetry mmm\n";
        
    const char* unit_cell =
    "UnitCell \
    [ \
     Var = Variant[ \
                   (p=1) \
                   C1   0.5 -0.010700 -0.103000 0.131900    0.000287 0.000236 0.000770 0.000136 0.000092 0.000031 \
                   ] \
     ]\n";
    
    const char* modes= "Modes [ ] \n";
    
    const char* correlations = 
    "Correlations \
    [ \
    [(0,0,0) \
    Multiplicity 0.5] \
    ]";
    
    string cryst =     crystal_params;
    cryst +=    unit_cell ;
    
    cryst +=    modes;
    cryst +=    correlations;
    

    TS_ASSERT(test_parser_nores(crystal_params,a_parser.program_options,a_skipper));
    TS_ASSERT(test_parser_nores(unit_cell,a_parser.unit_cell,a_skipper));
    TS_ASSERT(test_parser_nores(correlations,a_parser.correlations,a_skipper));
    
    TS_ASSERT(test_parser_nores(cryst,a_parser,a_skipper));
  }
  
  void test_translational_mode_parser()  {
    Model a_model;
    a_parser.add_model(&a_model);
    
    ChemicalUnit* one_atom;
    const char cell [] = "Cell 1 1 1 90 90 90 \n\
                   DiffuseScatteringGrid 0 0 0 0 0 0 0 0 0\n\
                   LaueSymmetry -1";
    
    test_parser_nores(cell,a_parser.program_options,a_skipper);
    run_parser("one_atom = [C1   0.5 -0.010700 -0.103000 0.131900    0.000287 0.000236 0.000770 0.000136 0.000092 0.000031]",a_parser.chemical_unit_assignement,a_skipper,one_atom);
    ADPMode* answer = translational_mode(one_atom, 0,unit_matrix);
    ADPMode* res;
                                                                  
    run_parser("TranslationalMode(one_atom,x)",a_parser.translational_mode,a_skipper,res);
    TS_ASSERT_EQUALS(*answer,*res);
  }
  
  void test_rotational_mode_parser()  {
    Model a_model;
    
    double prm[]={1,1,1,90,90,120};
    vector<double> unit_cell_params(prm,prm+6);
    
    a_model.initialize_unit_cell(unit_cell_params);
    
    a_parser.add_model(&a_model);
    
    ChemicalUnit* one_atom;
    run_parser("atom_rot = [C1   0.5 -0.010700 -0.103000 0.131900    0.000287 0.000236 0.000770 0.000136 0.000092 0.000031]",a_parser.chemical_unit_assignement,a_skipper,one_atom);
    ADPMode* answer = rot_mode(one_atom, vec3<double>(0,0,1), vec3<double>(0,0,0),a_model.cell.cell.metrical_matrix());
    ADPMode* res;
    
    run_parser("RotationalMode(atom_rot, 0,0,1, 0,0,0)",a_parser.rotational_mode,a_skipper,res);
    TS_ASSERT_EQUALS(*answer,*res);
  }
  
  void test_modes_parser()  {
    Model a_model;
    a_model.cell = UnitCell(1,1,1,90,90,90);
    a_parser.add_model(&a_model);
    
    ChemicalUnit* one_atom;
    run_parser("modes_atom = [C1   0.5 -0.010700 -0.103000 0.131900    0.000287 0.000236 0.000770 0.000136 0.000092 0.000031]",a_parser.chemical_unit_assignement,a_skipper,one_atom);
    ADPMode* answer = translational_mode(one_atom, 1,unit_matrix);
    StructurePartRef res;
    
    TS_ASSERT(test_parser_nores("Modes[ mode1 = TranslationalMode(modes_atom,y) ]",a_parser.modes,a_skipper));
    run_parser("mode1",a_parser.identifier,a_skipper,res);
    ADPMode* result=boost::get<ADPMode*>(res);

    TS_ASSERT_EQUALS(*answer,*result);
    TS_ASSERT_EQUALS(1,a_model.modes.size());
    TS_ASSERT_EQUALS(boost::get<ADPMode*>(res),&a_model.modes[0]);
  }
  
  void test_adp_correlation_parser()  {
    Model a_model;
    a_parser.add_model(&a_model);
    
    test_parser_nores("modes_atom2 = [C1   0.5 -0.010700 -0.103000 0.131900    0.000287 0.000236 0.000770 0.000136 0.000092 0.000031]",a_parser.chemical_unit_assignement,a_skipper);
    test_parser_nores("Modes[ mode_x = TranslationalMode(modes_atom,x) mode_y = TranslationalMode(modes_atom,y)]",a_parser.modes,a_skipper);
    
    ADPMode *mol_x = &a_model.modes[0];
    ADPMode *mol_y = &a_model.modes[1];
    DoubleADPMode correct(mol_x,mol_y,0.123);
    DoubleADPMode *res;
    
    TS_ASSERT(run_parser("ADPCorrelation(mode_x,mode_y,0.123)",a_parser.adp_correlation,a_skipper,res));
    TS_ASSERT_EQUALS(*res,correct);
  }

  void test_initialization_of_refinable_variables() {
    Model a_model;
    a_parser.add_model(&a_model);
    
    vector<boost::fusion::tuple<string,double> > res;
    run_parser("RefinableVariables[ a= 1 b=0; c=0 d=0.12(19);]",a_parser.refinable_parameters,a_skipper,res);
    TS_ASSERT_EQUALS(4, res.size());
    TS_ASSERT_EQUALS("a",boost::fusion::get<0>(res[0]));
    TS_ASSERT_EQUALS(1,boost::fusion::get<1>(res[0]));
    TS_ASSERT_EQUALS("b",boost::fusion::get<0>(res[1]));
    TS_ASSERT_EQUALS(0,boost::fusion::get<1>(res[1]));
    TS_ASSERT_EQUALS(0.12,boost::fusion::get<1>(res[3]));
  }
  
  void test_size_effect_corrlation_parser()  {
    Model a_model;
    a_parser.add_model(&a_model);

    ChemicalUnit* se_atom;
    run_parser("se_atom = [C1   0.5 -0.010700 -0.103000 0.131900    0.000287 0.000236 0.000770 0.000136 0.000092 0.000031]",a_parser.chemical_unit_assignement,a_skipper,se_atom);

    ADPMode* mode;
    run_parser("se_atom_x = TranslationalMode(se_atom,x)",a_parser.mode_assignement,a_skipper,mode);

    SizeEffect* res;
    TS_ASSERT(run_parser("SizeEffect(se_atom,se_atom_x,0.123)",a_parser.size_effect,a_skipper,res));
    TS_ASSERT_EQUALS(SizeEffect(se_atom,mode,0.123),*res);
    TS_ASSERT(run_parser("SizeEffect(se_atom_x,se_atom,0.123)",a_parser.size_effect,a_skipper,res)); //other order of operands
    TS_ASSERT_EQUALS(SizeEffect(mode,se_atom,0.123),*res);
  }
  
  void notest_integrated_minimizer()  {
    // create model that should be refined
    // unit cell 1 1 1 90 90 90
    // atoms [0.5 0;-0.5 0;][ 0 0.5;0 -0.5]
    string  model =     "Cell 1 1 1 90 90 90 \n\
    DiffuseScatteringGrid -1 -1 -1   0.1 0.1 0.1   20 20 20 \n\
    InitialParams 1 0.0 \
    UnitCell \
    [ \
    Var = Variant[ \
    (p=0.5) \
    [    C1   0.5    0.5 0 0   0 0 0 0 0 0 \
         C2   0.5   -0.5 0 0   0 0 0 0 0 0 \
    ]\
    (p=0.5)\
    [    C3   0.5  0  0.5 0   0 0 0 0 0 0 \
         C4   0.5  0 -0.5 0   0 0 0 0 0 0 \
    ] \
    ] \
    ] \
    Correlations \
    [ \
    [(0,0,0) \
      Multiplicity 0.5 \
      SubstitutionalCorrelation(Var,Var,0.25+~i)] \
    ]";

    // calculate diffuse scattering from the model with known correlation
    Model a_model(model); 
    
    //hook to raw internal function of diffuse scattering caclulation
    vector<double> params(1,1);
    params.push_back(0.25);
    
    a_model.calculate(params);
    IntensityMap diffuse_map(a_model.data().size());
    for(int i=0; i<diffuse_map.size_1d(); ++i)
      diffuse_map.at(i) = a_model.data().at(i);

    
    // run minimizer and compare results
    Minimizer a_minimizer;
    //vector<double> initial_params(1,0);
    vector<double> refined_params(2,0);
    
    OptionalIntensityMap weights;
    refined_params = a_minimizer.minimize(a_model.refinement_parameters, &diffuse_map, &a_model, &weights);
    
    TS_ASSERT_DELTA(1.0,refined_params[0],0.01);
    TS_ASSERT_DELTA(0.25,refined_params[1],0.01);
  }
  
  void test_laue_symmetry_is_compatible_with_lattice() {

    TS_ASSERT(LaueSymmetry("-1").is_compatible_with_cell(UnitCell(1,2,3,70,80,99)));
    TS_ASSERT_EQUALS(false,LaueSymmetry("2/m").is_compatible_with_cell(UnitCell(1,2,3,70,80,99)));
    TS_ASSERT_EQUALS(false,LaueSymmetry("6/m").is_compatible_with_cell(UnitCell(1,1,1,90,90,90)));
    TS_ASSERT_EQUALS(true,LaueSymmetry("-3mR").is_compatible_with_cell(UnitCell(1,1,1,90,90,90)));
  }
  
  ///TODO: write this test, since I destrroyed the old one
  void test_laue_symmetry_apply_symmetry_to_pairs()
	{
    Grid incompatible_grid(trivial_unit_cell, vec3<double> (1,2,3), vec3<double> (1,2,3),vec3<int> (1,2,3));
    LaueSymmetry laue("4/m",incompatible_grid);
    
		vector<AtomicPair> pairs;
		Atom atom3(string("C2"),1,  0.5,0.25,0,  1,0,0,-1,0,0);
		AtomicPair aPair(atom1,atom3);
		pairs.push_back(aPair);
		
		pairs = laue.apply_patterson_symmetry(pairs);
		
		TS_ASSERT_EQUALS(8,pairs.size());
    
    pairs.clear();
    pairs.push_back(aPair);
    
    incompatible_grid = Grid(trivial_unit_cell, vec3<double> (0.1,0.1,0.1), vec3<double> (-9,-9,0), vec3<int>(180,180,1));
    laue = LaueSymmetry("m-3m",incompatible_grid);
    pairs = laue.apply_patterson_symmetry(pairs);
		
		TS_ASSERT_EQUALS(3,pairs.size());
    
    //with the following symmetry and grid there is one symmetry element which should be applied to pairs
    // namely -y,x-y,z - three fold rotation
    incompatible_grid = Grid(trivial_unit_cell, 
                             vec3<double> (1,0.5,1), //centering is all right
                             vec3<double> (-9,-4.5,0), // but grid steps are different. 
                             vec3<int>(18,18,1));
    laue = LaueSymmetry("-3H",incompatible_grid);
    
    pairs.clear();
    pairs.push_back(aPair); // coordinates are 0.5 0.25 0
    // symmetry equivalent are -0.25,0.25,0   and -0.25,-0.5,0
    
    pairs = laue.apply_patterson_symmetry(pairs);
		
		TS_ASSERT_EQUALS(3,pairs.size());
    TS_ASSERT(almost_equal(vec3<double>(-0.25,0.25,0),pairs[1].r()));
    TS_ASSERT(almost_equal(vec3<double>(-0.25,-0.5,0),pairs[2].r()));
	}
  
  void test_laue_symmetry_apply_matrix_generator()
	{
    LaueSymmetry laue;
    
		mat3<double> T(-1,0,0,0,1,0,0,0,1); //mirror plane
		vector<AtomicPair> pairs;
		pairs.push_back(AtomicPair(atom1,atom2));// 0 0 0->0.5 0.5 0
		pairs[0].U()=sym_mat3<double>(1,2,3,1,0,0);
		
		pairs=laue.apply_matrix_generator(pairs,T,2);
		TS_ASSERT_EQUALS(1,(T*sym_mat3<double>(1,2,3,-1,0,0)*T.transpose())[3]);
		TS_ASSERT_EQUALS(2,pairs.size());
		TS_ASSERT_EQUALS(vec3<double>(-0.5,0.5,0.5),pairs[1].r());
		TS_ASSERT_EQUALS(sym_mat3<double>(1,2,3,-1,0,0),pairs[1].U());
	}
  
  void test_Laue_have_same_step()  {
    LaueSymmetry laue_symmetry;
    Grid grid_xy (trivial_unit_cell,vec3<double>(1,1,0.5),vec3<double>(-10,-10,-10),vec3<int>(20,20,20));
    TS_ASSERT_EQUALS(true, laue_symmetry.have_same_step_and_ll(grid_xy,"xy"));
    TS_ASSERT_EQUALS(false, laue_symmetry.have_same_step_and_ll(grid_xy,"xyxz"));
    
    Grid grid_xy_again (trivial_unit_cell,vec3<double>(1,1,1),vec3<double>(-10,-10,-10),vec3<int>(20,20,30));
    TS_ASSERT_EQUALS(true, laue_symmetry.have_same_step_and_ll(grid_xy_again,"xy"));
    TS_ASSERT_EQUALS(false, laue_symmetry.have_same_step_and_ll(grid_xy_again,"xyxz"));
                                      
    Grid grid_xyz (trivial_unit_cell,vec3<double>(1,1,1),vec3<double>(-10,-10,-10),vec3<int>(20,20,20));
    TS_ASSERT_EQUALS(true, laue_symmetry.have_same_step_and_ll(grid_xyz,"xyxz"));
  }
  void test_Laue_centered_correctly()  {
    LaueSymmetry laue_symmetry;
    Grid grid_xy (trivial_unit_cell,vec3<double>(1,1,0.5),vec3<double>(-10,-10,-10),vec3<int>(20,20,20));
    
    TS_ASSERT_EQUALS(true, laue_symmetry.centered_correctly(grid_xy,"xy"));
    TS_ASSERT_EQUALS(false, laue_symmetry.centered_correctly(grid_xy,"z"));
    TS_ASSERT_EQUALS(false, laue_symmetry.centered_correctly(grid_xy,"xyz"));
    
    // we do not support grids with uneven lengths
    Grid grid_xy_different_step (trivial_unit_cell,vec3<double>(1,0.5,1),vec3<double>(-10,-5,-10),vec3<int>(20,20,21));
    TS_ASSERT_EQUALS(true, laue_symmetry.centered_correctly(grid_xy_different_step,"xy"));
    TS_ASSERT_EQUALS(false, laue_symmetry.centered_correctly(grid_xy_different_step,"z"));
    
    Grid grid_z (trivial_unit_cell,vec3<double>(1,1,1),vec3<double>(-12,-1,-5),vec3<int>(10,10,10));
    TS_ASSERT_EQUALS(true, laue_symmetry.centered_correctly(grid_z,"z"));    
  }
  void test_Laue_generator_initialization()  {
    Grid grid(trivial_unit_cell,vec3<double>(0.1,0.1,0.1),vec3<double>(-9,-9,0),vec3<int>(180,180,1));
    LaueSymmetry laue_symmetry("m-3m",grid);
    TS_ASSERT_EQUALS(4, laue_symmetry.generators_on_map.size())
    
    
  }
  //The following block was generated with cctbx
  void test_Laue_Symmetry_1bar()  {
    IntensityMap int_map(10,10,10);
    
    int_map.at(5+3,5+2,5+1)=2; // the group class is 2

    LaueSymmetry a_symmetry("-1");
    
    a_symmetry.apply_patterson_symmetry(int_map);
    
    TS_ASSERT_EQUALS(int_map.at(5+3,5+2,5+1),1); //xyz
    int_map.at(5+3,5+2,5+1)=0;
    TS_ASSERT_EQUALS(int_map.at(5-3,5-2,5-1),1); //-x,-y,-z
    int_map.at(5-3,5-2,5-1)=0;
    
    for(int k=1;k<10;++k)
      for(int j=1;j<10;++j)
        for(int i=1;i<10;++i)
          TS_ASSERT_EQUALS(int_map.at(i,j,k),0); //some unexpected symmetry
  }
  void test_Laue_Symmetry_m3m()  {
    //hall symbol -P 4 2 3
    

    IntensityMap int_map(10,10,10);
    
    int x=3,y=2,z=1;
    int_map.at(5+x,5+y,5+z)=48; // the group class
    
    LaueSymmetry a_symmetry("m-3m");
    
    a_symmetry.apply_patterson_symmetry(int_map);
    
    
    /* this tests are generated with the following cctbx.python script:
     
     from cctbx import sgtbx
     import re
     sg = sgtbx.space_group("-P 4 2 3")
     for i in range(0,sg.order_z()):
      line = re.sub(',',',5+',sg(i).as_xyz())
      print "TS_ASSERT_EQUALS(int_map.at(5+" + line + "),1); //" + sg(i).as_xyz()
      print "int_map.at(5+" + line + ")=0; "
     */
  
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+z),1); //x,y,z
    int_map.at(5+x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+x,5+z),1); //-y,x,z
    int_map.at(5+-y,5+x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+z),1); //-x,-y,z
    int_map.at(5+-x,5+-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+-x,5+z),1); //y,-x,z
    int_map.at(5+y,5+-x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+-y,5+-z),1); //x,-y,-z
    int_map.at(5+x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+x,5+-z),1); //y,x,-z
    int_map.at(5+y,5+x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+y,5+-z),1); //-x,y,-z
    int_map.at(5+-x,5+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+-x,5+-z),1); //-y,-x,-z
    int_map.at(5+-y,5+-x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+z,5+x,5+y),1); //z,x,y
    int_map.at(5+z,5+x,5+y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+z,5+y),1); //-x,z,y
    int_map.at(5+-x,5+z,5+y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-z,5+-x,5+y),1); //-z,-x,y
    int_map.at(5+-z,5+-x,5+y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+-z,5+y),1); //x,-z,y
    int_map.at(5+x,5+-z,5+y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+z,5+-x,5+-y),1); //z,-x,-y
    int_map.at(5+z,5+-x,5+-y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+z,5+-y),1); //x,z,-y
    int_map.at(5+x,5+z,5+-y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-z,5+x,5+-y),1); //-z,x,-y
    int_map.at(5+-z,5+x,5+-y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-z,5+-y),1); //-x,-z,-y
    int_map.at(5+-x,5+-z,5+-y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+z,5+x),1); //y,z,x
    int_map.at(5+y,5+z,5+x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+-z,5+-x),1); //y,-z,-x
    int_map.at(5+y,5+-z,5+-x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+z,5+y,5+-x),1); //z,y,-x
    int_map.at(5+z,5+y,5+-x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+z,5+-x),1); //-y,z,-x
    int_map.at(5+-y,5+z,5+-x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-z,5+-y,5+-x),1); //-z,-y,-x
    int_map.at(5+-z,5+-y,5+-x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+-z,5+x),1); //-y,-z,x
    int_map.at(5+-y,5+-z,5+x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+z,5+-y,5+x),1); //z,-y,x
    int_map.at(5+z,5+-y,5+x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-z,5+y,5+x),1); //-z,y,x
    int_map.at(5+-z,5+y,5+x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+-z),1); //-x,-y,-z
    int_map.at(5+-x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+-x,5+-z),1); //y,-x,-z
    int_map.at(5+y,5+-x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+-z),1); //x,y,-z
    int_map.at(5+x,5+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+x,5+-z),1); //-y,x,-z
    int_map.at(5+-y,5+x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+y,5+z),1); //-x,y,z
    int_map.at(5+-x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+-x,5+z),1); //-y,-x,z
    int_map.at(5+-y,5+-x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+-y,5+z),1); //x,-y,z
    int_map.at(5+x,5+-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+x,5+z),1); //y,x,z
    int_map.at(5+y,5+x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-z,5+-x,5+-y),1); //-z,-x,-y
    int_map.at(5+-z,5+-x,5+-y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+-z,5+-y),1); //x,-z,-y
    int_map.at(5+x,5+-z,5+-y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+z,5+x,5+-y),1); //z,x,-y
    int_map.at(5+z,5+x,5+-y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+z,5+-y),1); //-x,z,-y
    int_map.at(5+-x,5+z,5+-y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-z,5+x,5+y),1); //-z,x,y
    int_map.at(5+-z,5+x,5+y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-z,5+y),1); //-x,-z,y
    int_map.at(5+-x,5+-z,5+y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+z,5+-x,5+y),1); //z,-x,y
    int_map.at(5+z,5+-x,5+y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+z,5+y),1); //x,z,y
    int_map.at(5+x,5+z,5+y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+-z,5+-x),1); //-y,-z,-x
    int_map.at(5+-y,5+-z,5+-x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+z,5+x),1); //-y,z,x
    int_map.at(5+-y,5+z,5+x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-z,5+-y,5+x),1); //-z,-y,x
    int_map.at(5+-z,5+-y,5+x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+-z,5+x),1); //y,-z,x
    int_map.at(5+y,5+-z,5+x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+z,5+y,5+x),1); //z,y,x
    int_map.at(5+z,5+y,5+x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+z,5+-x),1); //y,z,-x
    int_map.at(5+y,5+z,5+-x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-z,5+y,5+-x),1); //-z,y,-x
    int_map.at(5+-z,5+y,5+-x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+z,5+-y,5+-x),1); //z,-y,-x
    int_map.at(5+z,5+-y,5+-x)=0; 
    
    for(int k=1;k<10;++k)
      for(int j=1;j<10;++j)
        for(int i=1;i<10;++i)
          TS_ASSERT_EQUALS(int_map.at(i,j,k),0); //some unexpected symmetry
                                                                        
  }
  void test_Laue_symmetry_m3bar()  {
    // hall symbol is -P 2 2 3
    
    IntensityMap int_map(10,10,10);
    
    int x=3,y=2,z=1;
    int_map.at(5+x,5+y,5+z)=24; // the group class
    
    LaueSymmetry a_symmetry("m-3");
    
    a_symmetry.apply_patterson_symmetry(int_map);
        
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+z),1); //x,y,z
    int_map.at(5+x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+z),1); //-x,-y,z
    int_map.at(5+-x,5+-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+-y,5+-z),1); //x,-y,-z
    int_map.at(5+x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+y,5+-z),1); //-x,y,-z
    int_map.at(5+-x,5+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+z,5+x,5+y),1); //z,x,y
    int_map.at(5+z,5+x,5+y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-z,5+-x,5+y),1); //-z,-x,y
    int_map.at(5+-z,5+-x,5+y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+z,5+-x,5+-y),1); //z,-x,-y
    int_map.at(5+z,5+-x,5+-y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-z,5+x,5+-y),1); //-z,x,-y
    int_map.at(5+-z,5+x,5+-y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+z,5+x),1); //y,z,x
    int_map.at(5+y,5+z,5+x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+-z,5+-x),1); //y,-z,-x
    int_map.at(5+y,5+-z,5+-x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+z,5+-x),1); //-y,z,-x
    int_map.at(5+-y,5+z,5+-x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+-z,5+x),1); //-y,-z,x
    int_map.at(5+-y,5+-z,5+x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+-z),1); //-x,-y,-z
    int_map.at(5+-x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+-z),1); //x,y,-z
    int_map.at(5+x,5+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+y,5+z),1); //-x,y,z
    int_map.at(5+-x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+-y,5+z),1); //x,-y,z
    int_map.at(5+x,5+-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-z,5+-x,5+-y),1); //-z,-x,-y
    int_map.at(5+-z,5+-x,5+-y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+z,5+x,5+-y),1); //z,x,-y
    int_map.at(5+z,5+x,5+-y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-z,5+x,5+y),1); //-z,x,y
    int_map.at(5+-z,5+x,5+y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+z,5+-x,5+y),1); //z,-x,y
    int_map.at(5+z,5+-x,5+y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+-z,5+-x),1); //-y,-z,-x
    int_map.at(5+-y,5+-z,5+-x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+z,5+x),1); //-y,z,x
    int_map.at(5+-y,5+z,5+x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+-z,5+x),1); //y,-z,x
    int_map.at(5+y,5+-z,5+x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+z,5+-x),1); //y,z,-x
    int_map.at(5+y,5+z,5+-x)=0; 
    
    for(int k=1;k<10;++k)
      for(int j=1;j<10;++j)
        for(int i=1;i<10;++i)
          TS_ASSERT_EQUALS(int_map.at(i,j,k),0); //some unexpected symmetry
    
  }
  void test_Laue_symmetry_6mmm()  {
    // hall symbol is -P 6 2
    
    IntensityMap int_map(10,10,10);
	  int_map.grid.reciprocal_flag=false;
    
    int x=3,y=2,z=1;
    int_map.at(5+x,5+y,5+z)=24; // the group class
    
    LaueSymmetry a_symmetry("6/mmm");

    TS_ASSERT_EQUALS(4,a_symmetry.generators_on_map.size());
    
    a_symmetry.apply_patterson_symmetry(int_map);
    
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+z),1); //x,y,z
    int_map.at(5+x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x-y,5+x,5+z),1); //x-y,x,z
    int_map.at(5+x-y,5+x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+x-y,5+z),1); //-y,x-y,z
    int_map.at(5+-y,5+x-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+z),1); //-x,-y,z
    int_map.at(5+-x,5+-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x+y,5+-x,5+z),1); //-x+y,-x,z
    int_map.at(5+-x+y,5+-x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+-x+y,5+z),1); //y,-x+y,z
    int_map.at(5+y,5+-x+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+-x,5+-z),1); //-y,-x,-z
    int_map.at(5+-y,5+-x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x-y,5+-y,5+-z),1); //x-y,-y,-z
    int_map.at(5+x-y,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+x-y,5+-z),1); //x,x-y,-z
    int_map.at(5+x,5+x-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+x,5+-z),1); //y,x,-z
    int_map.at(5+y,5+x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x+y,5+y,5+-z),1); //-x+y,y,-z
    int_map.at(5+-x+y,5+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-x+y,5+-z),1); //-x,-x+y,-z
    int_map.at(5+-x,5+-x+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+-z),1); //-x,-y,-z
    int_map.at(5+-x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x+y,5+-x,5+-z),1); //-x+y,-x,-z
    int_map.at(5+-x+y,5+-x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+-x+y,5+-z),1); //y,-x+y,-z
    int_map.at(5+y,5+-x+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+-z),1); //x,y,-z  //Это как высоченные стены таинственных письмен
    int_map.at(5+x,5+y,5+-z)=0; //В комнате стало жарко, как будто одного лаптома достаточно чтобы нагреть всю комнату
    TS_ASSERT_EQUALS(int_map.at(5+x-y,5+x,5+-z),1); //x-y,x,-z
    int_map.at(5+x-y,5+x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+x-y,5+-z),1); //-y,x-y,-z
    int_map.at(5+-y,5+x-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+x,5+z),1); //y,x,z
    int_map.at(5+y,5+x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x+y,5+y,5+z),1); //-x+y,y,z
    int_map.at(5+-x+y,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-x+y,5+z),1); //-x,-x+y,z
    int_map.at(5+-x,5+-x+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+-x,5+z),1); //-y,-x,z
    int_map.at(5+-y,5+-x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x-y,5+-y,5+z),1); //x-y,-y,z
    int_map.at(5+x-y,5+-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+x-y,5+z),1); //x,x-y,z
    int_map.at(5+x,5+x-y,5+z)=0; 
    
    for(int k=1;k<10;++k)
      for(int j=1;j<10;++j)
        for(int i=1;i<10;++i)
          TS_ASSERT_EQUALS(int_map.at(i,j,k),0); //some unexpected symmetry

  }
  void test_Laue_symmetry_6m()  {
    // hall symbol is -P 6
    
    IntensityMap int_map(10,10,10);
	  int_map.grid.reciprocal_flag=false;

    for(int k=1;k<10;++k)
      for(int j=1;j<10;++j)
        for(int i=1;i<10;++i)
          TS_ASSERT_EQUALS(int_map.at(i,j,k),0); //some unexpected symmetry
    
    int x=3,y=2,z=1;
    int_map.at(5+x,5+y,5+z)=12; // the group class
    
    LaueSymmetry a_symmetry("6/m");
    
    a_symmetry.apply_patterson_symmetry(int_map);
    
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+z),1); //x,y,z
    int_map.at(5+x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x-y,5+x,5+z),1); //x-y,x,z
    int_map.at(5+x-y,5+x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+x-y,5+z),1); //-y,x-y,z
    int_map.at(5+-y,5+x-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+z),1); //-x,-y,z
    int_map.at(5+-x,5+-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x+y,5+-x,5+z),1); //-x+y,-x,z
    int_map.at(5+-x+y,5+-x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+-x+y,5+z),1); //y,-x+y,z
    int_map.at(5+y,5+-x+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+-z),1); //-x,-y,-z
    int_map.at(5+-x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x+y,5+-x,5+-z),1); //-x+y,-x,-z
    int_map.at(5+-x+y,5+-x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+-x+y,5+-z),1); //y,-x+y,-z
    int_map.at(5+y,5+-x+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+-z),1); //x,y,-z
    int_map.at(5+x,5+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x-y,5+x,5+-z),1); //x-y,x,-z
    int_map.at(5+x-y,5+x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+x-y,5+-z),1); //-y,x-y,-z
    int_map.at(5+-y,5+x-y,5+-z)=0; 
    
    for(int k=1;k<10;++k)
      for(int j=1;j<10;++j)
        for(int i=1;i<10;++i)
          TS_ASSERT_EQUALS(int_map.at(i,j,k),0); //some unexpected symmetry
    
  }
  void test_Laue_symmetry_3barmH()  {
    // hall symbol is -R 3 2"
    
    IntensityMap int_map(10,10,10);
	  int_map.grid.reciprocal_flag=false;
    
    int x=3,y=2,z=1;
    int_map.at(5+x,5+y,5+z)=12; // the group class
    
    LaueSymmetry a_symmetry("-3mH");
    
    a_symmetry.apply_patterson_symmetry(int_map);
    
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+z),1); //x,y,z
    int_map.at(5+x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+x-y,5+z),1); //-y,x-y,z
    int_map.at(5+-y,5+x-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x+y,5+-x,5+z),1); //-x+y,-x,z
    int_map.at(5+-x+y,5+-x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+x,5+-z),1); //y,x,-z
    int_map.at(5+y,5+x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-x+y,5+-z),1); //-x,-x+y,-z
    int_map.at(5+-x,5+-x+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x-y,5+-y,5+-z),1); //x-y,-y,-z
    int_map.at(5+x-y,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+-z),1); //-x,-y,-z
    int_map.at(5+-x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+-x+y,5+-z),1); //y,-x+y,-z
    int_map.at(5+y,5+-x+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x-y,5+x,5+-z),1); //x-y,x,-z
    int_map.at(5+x-y,5+x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+-x,5+z),1); //-y,-x,z
    int_map.at(5+-y,5+-x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+x-y,5+z),1); //x,x-y,z
    int_map.at(5+x,5+x-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x+y,5+y,5+z),1); //-x+y,y,z
    int_map.at(5+-x+y,5+y,5+z)=0; 
    
    for(int k=1;k<10;++k)
      for(int j=1;j<10;++j)
        for(int i=1;i<10;++i)
          TS_ASSERT_EQUALS(int_map.at(i,j,k),0); //some unexpected symmetry
    
  }
  void test_Laue_symmetry_3barmR()  {
    // hall symbol is -P 3* 2
    
    IntensityMap int_map(10,10,10);
    
    int x=3,y=2,z=1;
    int_map.at(5+x,5+y,5+z)=12; // the group class
    
    LaueSymmetry a_symmetry("-3mR");
    
    a_symmetry.apply_patterson_symmetry(int_map);
    
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+z),1); //x,y,z
    int_map.at(5+x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+z,5+x,5+y),1); //z,x,y
    int_map.at(5+z,5+x,5+y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+z,5+x),1); //y,z,x
    int_map.at(5+y,5+z,5+x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+-x,5+-z),1); //-y,-x,-z
    int_map.at(5+-y,5+-x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-z,5+-y,5+-x),1); //-z,-y,-x
    int_map.at(5+-z,5+-y,5+-x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-z,5+-y),1); //-x,-z,-y
    int_map.at(5+-x,5+-z,5+-y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+-z),1); //-x,-y,-z
    int_map.at(5+-x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-z,5+-x,5+-y),1); //-z,-x,-y
    int_map.at(5+-z,5+-x,5+-y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+-z,5+-x),1); //-y,-z,-x
    int_map.at(5+-y,5+-z,5+-x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+x,5+z),1); //y,x,z
    int_map.at(5+y,5+x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+z,5+y,5+x),1); //z,y,x
    int_map.at(5+z,5+y,5+x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+z,5+y),1); //x,z,y
    int_map.at(5+x,5+z,5+y)=0; //Я помню это показалось мне смешной мыслью - если хочешь чтобы все запомнили как тебя зовут, заставь людей вводить твои имя и фамилию чтобы добраться до интернета
    
    for(int k=1;k<10;++k)
      for(int j=1;j<10;++j)
        for(int i=1;i<10;++i) //Блин, эту задачу решить совершенно невозможно
          TS_ASSERT_EQUALS(int_map.at(i,j,k),0); //some unexpected symmetry
    //Так, надо доделать и уметываться отсюда
	  //тут слишком много чего отвлекает
	  //wdtnf цвета конечно крутые
  }
  void test_Laue_symmetry_3barH()  {
    // hall symbol is -P 3
    
    IntensityMap int_map(10,10,10);
	  int_map.grid.reciprocal_flag=false;
    
    int x=3,y=2,z=1;
    int_map.at(5+x,5+y,5+z)=6; // the group class
    
    LaueSymmetry a_symmetry("-3H");
    
    a_symmetry.apply_patterson_symmetry(int_map);
    
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+z),1); //x,y,z
    int_map.at(5+x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+x-y,5+z),1); //-y,x-y,z
    int_map.at(5+-y,5+x-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x+y,5+-x,5+z),1); //-x+y,-x,z
    int_map.at(5+-x+y,5+-x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+-z),1); //-x,-y,-z
    int_map.at(5+-x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+-x+y,5+-z),1); //y,-x+y,-z
    int_map.at(5+y,5+-x+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x-y,5+x,5+-z),1); //x-y,x,-z
    int_map.at(5+x-y,5+x,5+-z)=0; 
    
    for(int k=1;k<10;++k)
      for(int j=1;j<10;++j)
        for(int i=1;i<10;++i)
          TS_ASSERT_EQUALS(int_map.at(i,j,k),0); //some unexpected symmetry
    
  }
  void test_Laue_symmetry_m3barR()  {
    // hall symbol is -P 3*
    
    IntensityMap int_map(10,10,10);
    
    int x=3,y=2,z=1;
    int_map.at(5+x,5+y,5+z)=6; // the group class
    
    LaueSymmetry a_symmetry("-3R");
    
    a_symmetry.apply_patterson_symmetry(int_map);
    
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+z),1); //x,y,z
    int_map.at(5+x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+z,5+x,5+y),1); //z,x,y
    int_map.at(5+z,5+x,5+y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+z,5+x),1); //y,z,x
    int_map.at(5+y,5+z,5+x)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+-z),1); //-x,-y,-z
    int_map.at(5+-x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-z,5+-x,5+-y),1); //-z,-x,-y
    int_map.at(5+-z,5+-x,5+-y)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+-z,5+-x),1); //-y,-z,-x
    int_map.at(5+-y,5+-z,5+-x)=0; 
    
    for(int k=1;k<10;++k)
      for(int j=1;j<10;++j)
        for(int i=1;i<10;++i)
          TS_ASSERT_EQUALS(int_map.at(i,j,k),0); //some unexpected symmetry
    
  }
  void test_Laue_symmetry_4mmm()  {
    // hall symbol is -P 4 2
    
    IntensityMap int_map(10,10,10);
    
    int x=3,y=2,z=1;
    int_map.at(5+x,5+y,5+z)=16; // the group class
    
    LaueSymmetry a_symmetry("4/mmm");
    
    a_symmetry.apply_patterson_symmetry(int_map);
    
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+z),1); //x,y,z
    int_map.at(5+x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+x,5+z),1); //-y,x,z
    int_map.at(5+-y,5+x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+z),1); //-x,-y,z
    int_map.at(5+-x,5+-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+-x,5+z),1); //y,-x,z
    int_map.at(5+y,5+-x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+-y,5+-z),1); //x,-y,-z
    int_map.at(5+x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+x,5+-z),1); //y,x,-z
    int_map.at(5+y,5+x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+y,5+-z),1); //-x,y,-z
    int_map.at(5+-x,5+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+-x,5+-z),1); //-y,-x,-z
    int_map.at(5+-y,5+-x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+-z),1); //-x,-y,-z
    int_map.at(5+-x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+-x,5+-z),1); //y,-x,-z
    int_map.at(5+y,5+-x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+-z),1); //x,y,-z
    int_map.at(5+x,5+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+x,5+-z),1); //-y,x,-z
    int_map.at(5+-y,5+x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+y,5+z),1); //-x,y,z
    int_map.at(5+-x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+-x,5+z),1); //-y,-x,z
    int_map.at(5+-y,5+-x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+-y,5+z),1); //x,-y,z
    int_map.at(5+x,5+-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+x,5+z),1); //y,x,z
    int_map.at(5+y,5+x,5+z)=0; 
    
    
    for(int k=1;k<10;++k)
      for(int j=1;j<10;++j)
        for(int i=1;i<10;++i)
          TS_ASSERT_EQUALS(int_map.at(i,j,k),0); //some unexpected symmetry
    
  }
  void test_Laue_symmetry_4m()  {
    // hall symbol is -P 4
    
    IntensityMap int_map(10,10,10);
    
    int x=3,y=2,z=1;
    int_map.at(5+x,5+y,5+z)=8; // the group class
    
    LaueSymmetry a_symmetry("4/m");
    
    a_symmetry.apply_patterson_symmetry(int_map);
    
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+z),1); //x,y,z
    int_map.at(5+x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+x,5+z),1); //-y,x,z
    int_map.at(5+-y,5+x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+z),1); //-x,-y,z
    int_map.at(5+-x,5+-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+-x,5+z),1); //y,-x,z
    int_map.at(5+y,5+-x,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+-z),1); //-x,-y,-z
    int_map.at(5+-x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+y,5+-x,5+-z),1); //y,-x,-z
    int_map.at(5+y,5+-x,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+-z),1); //x,y,-z
    int_map.at(5+x,5+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-y,5+x,5+-z),1); //-y,x,-z
    int_map.at(5+-y,5+x,5+-z)=0; 
    
    for(int k=1;k<10;++k)
      for(int j=1;j<10;++j)
        for(int i=1;i<10;++i)
          TS_ASSERT_EQUALS(int_map.at(i,j,k),0); //some unexpected symmetry
    
  }
  void test_Laue_symmetry_mmm()  {
    // hall symbol is -P 2 2
    
    IntensityMap int_map(10,10,10);
    
    int x=3,y=2,z=1;
    int_map.at(5+x,5+y,5+z)=8; // the group class
    
    LaueSymmetry a_symmetry("mmm");
    
    a_symmetry.apply_patterson_symmetry(int_map);
    
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+z),1); //x,y,z
    int_map.at(5+x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+z),1); //-x,-y,z
    int_map.at(5+-x,5+-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+-y,5+-z),1); //x,-y,-z
    int_map.at(5+x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+y,5+-z),1); //-x,y,-z
    int_map.at(5+-x,5+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+-z),1); //-x,-y,-z
    int_map.at(5+-x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+-z),1); //x,y,-z
    int_map.at(5+x,5+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+y,5+z),1); //-x,y,z
    int_map.at(5+-x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+-y,5+z),1); //x,-y,z
    int_map.at(5+x,5+-y,5+z)=0; 
    
    for(int k=1;k<10;++k)
      for(int j=1;j<10;++j)
        for(int i=1;i<10;++i)
          TS_ASSERT_EQUALS(int_map.at(i,j,k),0); //some unexpected symmetry
    
  }
  void test_Laue_symmetry_2m()  {
    // hall symbol is -P 2
    
    IntensityMap int_map(10,10,10);
    
    int x=3,y=2,z=1;
    int_map.at(5+x,5+y,5+z)=4; // the group class
    
    LaueSymmetry a_symmetry("2/m");
    
    a_symmetry.apply_patterson_symmetry(int_map);
    
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+z),1); //x,y,z
    int_map.at(5+x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+z),1); //-x,-y,z
    int_map.at(5+-x,5+-y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+-z),1); //-x,-y,-z
    int_map.at(5+-x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+-z),1); //x,y,-z
    int_map.at(5+x,5+y,5+-z)=0; 
    
    for(int k=1;k<10;++k)
      for(int j=1;j<10;++j)
        for(int i=1;i<10;++i)
          TS_ASSERT_EQUALS(int_map.at(i,j,k),0); //some unexpected symmetry
    
  }
  void test_Laue_symmetry_2mb()  {
    // hall symbol is -P 2y
    
    IntensityMap int_map(10,10,10);
    
    int x=3,y=2,z=1;
    int_map.at(5+x,5+y,5+z)=4; // the group class
    
    LaueSymmetry a_symmetry("2/mb");
    
    a_symmetry.apply_patterson_symmetry(int_map);
    
    TS_ASSERT_EQUALS(int_map.at(5+x,5+y,5+z),1); //x,y,z
    int_map.at(5+x,5+y,5+z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+y,5+-z),1); //-x,y,-z
    int_map.at(5+-x,5+y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+-x,5+-y,5+-z),1); //-x,-y,-z
    int_map.at(5+-x,5+-y,5+-z)=0; 
    TS_ASSERT_EQUALS(int_map.at(5+x,5+-y,5+z),1); //x,-y,z
    int_map.at(5+x,5+-y,5+z)=0; 
    
    for(int k=1;k<10;++k)
      for(int j=1;j<10;++j)
        for(int i=1;i<10;++i)
          TS_ASSERT_EQUALS(int_map.at(i,j,k),0); //some unexpected symmetry
    
  }
  
  void test_comment_parser()  {
    TS_ASSERT(test_parser_nores("#Ruby-type comment\n",a_parser.skipper));
  }
  
  void test_filter_asymmetric_unit_m()  {
    // The symmetry m is a dummy symmetry for the z mirror part of 6/m 4/m 2/m 
    vector<AtomicPair> pairs;
    AtomicPair ethalon(atom1,atom2);
    for(int i=0; i<3; ++i)
      pairs.push_back(AtomicPair(atom1,atom2));
    
    pairs[1].average_r()=vec3<double>(10,20,0); //occupancy should be 1/2
    pairs[2].average_r()=-pairs[2].average_r(); //pair should be deleted
    
    LaueSymmetry a_symmetry("6/m");
    
    pairs=a_symmetry.filter_pairs_from_asymmetric_unit(pairs);
    
    TS_ASSERT_DIFFERS(ethalon,pairs[1]);
    TS_ASSERT_EQUALS(ethalon.average_p()/2,pairs[1].average_p()); // x y 0 should decrease occupancy
    TS_ASSERT_EQUALS(2,pairs.size()); 
  }
  void test_filter_asymmetric_unit_mmm()  {
    vector<AtomicPair> pairs;
    for(int i=0; i<12; ++i)
      pairs.push_back(AtomicPair(atom1,atom2));
    
    pairs[0].average_r()=vec3<double>(0,0,0); //occupancy should be 1/8
    pairs[1].average_r()=vec3<double>(10,20,0); //occupancy should be 1/2
    pairs[2].average_r()=vec3<double>(10,0,13); //occupancy should be 1/2
    pairs[3].average_r()=vec3<double>(0,10,20); //occupancy should be 1/2
    
    pairs[4].average_r()=vec3<double>(1,2,3); //This one stays
    AtomicPair ethalon=pairs[4];
    pairs[5].average_r()=vec3<double>(1,2,-3); //pair should be deleted
    pairs[6].average_r()=vec3<double>(1,-2,3); //pair should be deleted
    pairs[7].average_r()=vec3<double>(1,-2,-3); //pair should be deleted
    
    pairs[8].average_r()=vec3<double>(-1,2,3); //pair should be deleted
    pairs[9].average_r()=vec3<double>(-1,2,-3); //pair should be deleted
    pairs[10].average_r()=vec3<double>(-1,-2,3); //pair should be deleted
    pairs[11].average_r()=vec3<double>(-1,-2,-3); //pair should be deleted
    
    LaueSymmetry a_symmetry("mmm");
    
    pairs=a_symmetry.filter_pairs_from_asymmetric_unit(pairs);
    
    TS_ASSERT_EQUALS(ethalon.average_p()/8,pairs[0].average_p()); 
    
    for(int i=1; i<4; ++i)
      TS_ASSERT_EQUALS(ethalon.average_p()/2,pairs[i].average_p()); 
    
    TS_ASSERT_EQUALS(5,pairs.size()); //  all the rest are deleted
  }
  void test_filter_asymmetric_unit_4mmm()  {
    vector<AtomicPair> pairs;
    for(int i=0; i<4; ++i)
      pairs.push_back(AtomicPair(atom1,atom2));
    
    pairs[0].average_r()=vec3<double>(0,0,0); //occupancy should be 1/16
    pairs[1].average_r()=vec3<double>(10,10,12); //occupancy should be 1/2

    
    pairs[2].average_r()=vec3<double>(3,2,1); //This one stays
    AtomicPair ethalon=pairs[2];
    pairs[3].average_r()=vec3<double>(2,3,1); //pair should be deleted

    LaueSymmetry a_symmetry("4/mmm");
    
    pairs=a_symmetry.filter_pairs_from_asymmetric_unit(pairs);
    
    TS_ASSERT_EQUALS(ethalon.average_p()/16,pairs[0].average_p()); 
    
    
    TS_ASSERT_EQUALS(ethalon.average_p()/2,pairs[1].average_p()); 
    
//    TS_ASSERT_EQUALS(ethalon,pairs[2]); // does not change
    TS_ASSERT_EQUALS(3,pairs.size()); 
  }
  void test_filter_asymmetric_unit_m3m()  {
    vector<AtomicPair> pairs;
    for(int i=0; i<5; ++i)
      pairs.push_back(AtomicPair(atom1,atom2));
    
    pairs[0].average_r()=vec3<double>(0,0,0); //occupancy should be 1/48
    pairs[1].average_r()=vec3<double>(10,0,0); //occupancy should be 1/8
    
    pairs[2].average_r()=vec3<double>(3,2,1); //This one does not change
    AtomicPair ethalon=pairs[2];
    pairs[3].average_r()=vec3<double>(2,3,1); //pair should be deleted
    pairs[4].average_r()=vec3<double>(2,1,3); //pair should also be deleted
    
    LaueSymmetry a_symmetry("m-3m");
    
    pairs=a_symmetry.filter_pairs_from_asymmetric_unit(pairs);
    
    TS_ASSERT_EQUALS(ethalon.average_p()/48,pairs[0].average_p()); 
    TS_ASSERT_EQUALS(ethalon.average_p()/8,pairs[1].average_p()); 
    
//    TS_ASSERT_EQUALS(ethalon,pairs[2]); // does not change
    TS_ASSERT_EQUALS(3,pairs.size()); 
  }
  void test_filter_asymmetric_unit_6mmm()  {
    vector<AtomicPair> pairs;
    for(int i=0; i<4; ++i)
      pairs.push_back(AtomicPair(atom1,atom2));
    
    pairs[0].average_r()=vec3<double>(0,0,0); //occupancy should be 1/24
    pairs[1].average_r()=vec3<double>(10,0,0); //occupancy should be 1/4
    
    pairs[2].average_r()=vec3<double>(5,2,1); //This one does not change
    AtomicPair ethalon=pairs[2];
    pairs[3].average_r()=vec3<double>(3,2,1); //pair should be deleted
    
    LaueSymmetry a_symmetry("6/mmm");
    
    pairs=a_symmetry.filter_pairs_from_asymmetric_unit(pairs);
    
    TS_ASSERT_EQUALS(ethalon.average_p()/24,pairs[0].average_p()); 
    TS_ASSERT_EQUALS(ethalon.average_p()/4,pairs[1].average_p()); 
    
//    TS_ASSERT_EQUALS(ethalon,pairs[2]); // does not change
    TS_ASSERT_EQUALS(3,pairs.size()); 
  }
  void test_filter_asymmetric_unit_3barmH()  {
    vector<AtomicPair> pairs;
    for(int i=0; i<7; ++i)
      pairs.push_back(AtomicPair(atom1,atom2));
    
    pairs[0].average_r()=vec3<double>(0,0,0); //occupancy should be 1/6! Even though the group class is 12
    pairs[1].average_r()=vec3<double>(0,0,7); //occupancy should be 1/6
    pairs[2].average_r()=vec3<double>(5,-5,1); //occupancy should be 1/2
    pairs[3].average_r()=vec3<double>(6,3,2); //occupancy should be 1/2
    
    pairs[4].average_r()=vec3<double>(4,0,0); //This one does not change even though this is a two fold rotation axis
    AtomicPair ethalon=pairs[4];
    pairs[5].average_r()=vec3<double>(4,4,0); //pair should be deleted
    pairs[6].average_r()=vec3<double>(0,4,0); //pair should be deleted
    
    LaueSymmetry a_symmetry("-3mH");
    
    pairs=a_symmetry.filter_pairs_from_asymmetric_unit(pairs);
    
    TS_ASSERT_EQUALS(ethalon.average_p()/6,pairs[0].average_p()); 
    TS_ASSERT_EQUALS(ethalon.average_p()/6,pairs[1].average_p()); 
    TS_ASSERT_EQUALS(ethalon.average_p()/2,pairs[2].average_p()); 
    TS_ASSERT_EQUALS(ethalon.average_p()/2,pairs[3].average_p()); 
    
        
//    TS_ASSERT_EQUALS(ethalon,pairs[4]); // does not change
    TS_ASSERT_EQUALS(5,pairs.size()); 
  }
  void test_filter_asymmetric_unit_3barmR()  {
    vector<AtomicPair> pairs;
    for(int i=0; i<7; ++i)
      pairs.push_back(AtomicPair(atom1,atom2));
    
    pairs[0].average_r()=vec3<double>(0,0,0); //occupancy should be 1/6! Even though the group class is 12
    pairs[1].average_r()=vec3<double>(1,1,1); //occupancy should be 1/6
    pairs[2].average_r()=vec3<double>(5,5,1); //occupancy should be 1/2
    pairs[3].average_r()=vec3<double>(6,3,3); //occupancy should be 1/2
    
    pairs[4].average_r()=vec3<double>(4,0,-4); //This one does not change even though this is a two fold rotation axis
    AtomicPair ethalon=pairs[4];
    pairs[5].average_r()=vec3<double>(-4,4,0); //pair should be deleted
    pairs[6].average_r()=vec3<double>(0,-4,4); //pair should be deleted
    
    LaueSymmetry a_symmetry("-3mR");
    
    pairs=a_symmetry.filter_pairs_from_asymmetric_unit(pairs);
    
    TS_ASSERT_EQUALS(ethalon.average_p()/6,pairs[0].average_p()); 
    TS_ASSERT_EQUALS(ethalon.average_p()/6,pairs[1].average_p()); 
    TS_ASSERT_EQUALS(ethalon.average_p()/2,pairs[2].average_p()); 
    TS_ASSERT_EQUALS(ethalon.average_p()/2,pairs[3].average_p()); 
    
//    TS_ASSERT_EQUALS(ethalon,pairs[4]); // does not change
    TS_ASSERT_EQUALS(5,pairs.size()); 
  }
  
  void test_format_esd()
  {
    TS_ASSERT_EQUALS("1(2)", format_esd(1,2));
    TS_ASSERT_EQUALS("0.1(2)", format_esd(0.1,0.2));
    TS_ASSERT_EQUALS("0.1(2)", format_esd(0.133162,0.2123154));
    TS_ASSERT_EQUALS("0.13(12)", format_esd(0.133162,0.12123154));
    TS_ASSERT_EQUALS("0.0003(7)", format_esd(0.00033175,0.000717));
    TS_ASSERT_EQUALS("0.0000003(7)", format_esd(0.00000033175,0.000000717));
    TS_ASSERT_EQUALS("0.0000000003(7)", format_esd(0.00000000033175,0.000000000717));
    TS_ASSERT_EQUALS("0.123456(0)", format_esd(0.123456,0));
  }
  
	void nottestTODO()
	{
		TS_FAIL("Test SetUp contains memory leak");
		TS_FAIL("Rewrite grid with hkl-prefix for all vals which are in HKL - coordinates to separate them from reciprocal coordinates");
		TS_FAIL("Rewrite symmetry check and multiplication");
		TS_FAIL("Rewrite modifiers so they will have priority to make possible to create zero vector modifier");
	}
	
	
};
