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

#include "InputFileParser.h"

void add_symbol( qi::symbols<char,std::string>& molecular_type, vector<boost::tuple<string,Scatterer*> > inp) {
  for(int i=0; i<inp.size(); ++i) {
    std::string label = boost::get<0>(inp[i]);
    molecular_type.add(label.c_str(),label);
  }
}

void InputParser::InputParserIII()
{
  using namespace qi;
  using phoenix::bind;
  using phoenix::ref;
  
  chemical_unit %=
   lit("Void")[_val = phoenix::new_<AtomicAssembly>()]
   | chemical_unit_assignement
   | atomic_assembly
   | symmetric_chemical_unit
   | atom
  ;
  symmetric_chemical_unit = (identifier >> '*' >> symmetry_element)[_val = bind(create_symmetric_chemical_unit,_1,_2)];
  
  chemical_unit_assignement = 
    valid_identifier[_a=_1] 
    >> "="
    >> chemical_unit[_val = _1]
    >> eps[bind(InputParser::add_reference,phoenix::ref(references),_a,phoenix::construct<StructurePartRef>(_val))]
    ;
  
  atomic_assembly = ('[' >> *chemical_unit >> ']')[_val = phoenix::new_<AtomicAssembly>(_1)] ;
  
  atom = 
    (atom_name >> repeat(10)[lexeme[formula]]) [_val = bind(&Model::construct_atom,*ref(model),_1,_2)] //Uaniso and p
  | (atom_name > repeat(5)[lexeme[formula]])   [_val = bind(&Model::construct_atom_isotropic_adp,*ref(model),_1,_2)]
  ;
  
  molecular_scatterer =
    chemical_unit[_val = phoenix::new_<MolecularScatterer>(_1)]
    ;
  
  molecular_scatterers %=
    lit("MolecularScatterers")
    > "["
    >> *(valid_identifier
         >> "="
         >> molecular_scatterer)
    > "]"
    > eps[bind(&Model::register_molecular_scatterers,_val)]
    > eps[bind(&add_symbol,ref(molecular_type),_val)]
    ;
  
  atom_type.add("H","H")("He","He")("Li","Li")("Be","Be")("B","B")("C","C")("N","N")("O","O")("F","F")("Ne","Ne")("Na","Na")("Mg","Mg")("Al","Al")("Si","Si")("P","P")("S","S")("Cl","Cl")("Ar","Ar")("K","K")("Ca","Ca")("Sc","Sc")("Ti","Ti")("V","V")("Cr","Cr")("Mn","Mn")("Fe","Fe")("Co","Co")("Ni","Ni")("Cu","Cu")("Zn","Zn")("Ga","Ga")("Ge","Ge")("As","As")("Se","Se")("Br","Br")("Kr","Kr")("Rb","Rb")("Sr","Sr")("Y","Y")("Zr","Zr")("Nb","Nb")("Mo","Mo")("Tc","Tc")("Ru","Ru")("Rh","Rh")("Pd","Pd")("Ag","Ag")("Cd","Cd")("In","In")("Sn","Sn")("Sb","Sb")("Te","Te")("I","I")("Xe","Xe")("Cs","Cs")("Ba","Ba")("La","La")("Ce","Ce")("Pr","Pr")("Nd","Nd")("Pm","Pm")("Sm","Sm")("Eu","Eu")("Gd","Gd")("Tb","Tb")("Dy","Dy")("Ho","Ho")("Er","Er")("Tm","Tm")("Yb","Yb")("Lu","Lu")("Hf","Hf")("Ta","Ta")("W","W")("Re","Re")("Os","Os")("Ir","Ir")("Pt","Pt")("Au","Au")("Hg","Hg")("Tl","Tl")("Pb","Pb")("Bi","Bi")("Po","Po")("At","At")("Rn","Rn")("Fr","Fr")("Ra","Ra")("Ac","Ac")("Th","Th")("Pa","Pa")("U","U")("Np","Np")("Pu","Pu")("Am","Am")("Cm","Cm")("Bk","Bk")("Cf","Cf")("Es","Es")("Fm","Fm")("Md","Md")("No","No")("Lr","Lr")("Rf","Rf")("Db","Db")("Sg","Sg")("Bh","Bh")("Hs","Hs")("Mt","Mt")("Ds","Ds")("Rg","Rg")("Cn","Cn")("Uut","Uut")("Fl","Fl")("Uup","Uup")("Lv","Lv")("Uus","Uus")("Uuo","Uuo");
  
  atom_name %=
    molecular_type
    | atom_type >> *alnum
    ;
   
  identifier %= references >> !alnum;
  rest_of_the_line = *(qi::char_ - qi::eol) >> qi::eol;
};