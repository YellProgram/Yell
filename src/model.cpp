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

#include "model.h"
#include "InputFileParser.h"
#include <sstream>
#include "exceptions.h"
extern OutputHandler report;
typedef iterator_ Iterator;

Iterator line_start(Iterator pos, const Iterator& start) {
  if(pos!=start and *pos=='\n')
    --pos;
  while (pos!=start && *pos != '\n')
    --pos;

  return pos;
}
Iterator line_end(Iterator pos, const Iterator& end) {
  if(pos!=end and *pos=='\n')
    ++pos;
  while (pos!=end && *pos != '\n')
    ++pos;

  return pos;
}
Iterator step_from_newline(const Iterator& pos, const Iterator& end) {
  if(pos!=end and *pos=='\n')
    return pos + 1;
  else
    return pos;
}

void complain_about_error(Iterator start, Iterator current, Iterator end)
{
  //skip empty spaces till the line where the error happened
  InputParser a_parser;
  qi::parse(current,end,a_parser.skipper_no_assignement);

  current = step_from_newline(current,end);

  int line_number = 0;
  for(Iterator it=current; it!=start; it--)
    if(*it=='\n')
      line_number++;

  Iterator lstart = line_start(current, start);

  int error_pos_within_line = 0;
  for(Iterator t = step_from_newline(lstart,end); t!=current; ++t)
    ++error_pos_within_line;

  // Print two lines before the error line
  Iterator context_start = step_from_newline(line_start(line_start(lstart,start),start),end);
  
  Iterator lend = line_end(current,end);

  Iterator next_two_lines_start = step_from_newline(lend,end);
  Iterator next_two_lines_end = line_end(line_end(next_two_lines_start,end),end);

  std::string trouble_line_with_context(context_start,lend);
  std::string next_two_lines(next_two_lines_start,next_two_lines_end);

  std::ostringstream complain;
  
  complain << "Parsing failed with exception around line " << line_number+1 << ":\n\n" << trouble_line_with_context <<"\n";
  for(int i=0; i < error_pos_within_line - 1; ++i)
    complain << ' ';
  complain << "^-----roughly here\n";
  complain << next_two_lines << "\n";
  
  REPORT(ERROR) << complain.str();
}

void Model::calculate(vector<double> params,bool average_flag)
{
  ///\TODO: since calculate function funs several times, move parser creation to initialization stage
  InputParser a_parser;
  a_parser.add_model(this);
  a_parser.formula.initialize_refinable_variables(refined_variable_names,params);
  a_parser.formula.array = params;

  Iterator start = model.begin();
  Iterator end = model.end();
  bool r;
  try{
    r = qi::phrase_parse(start,end,a_parser,a_parser.skipper_no_assignement);
  }catch(const qi::expectation_failure<Iterator>& e)
  {
    complain_about_error(model.begin(), e.first, model.end());
    throw(TerminateProgram());
  }
  
  if(!r || start!=end)
  {
    string rest(start,end);
    REPORT(ERROR) << "Parsing failed, the rest of the file is:\n" << rest;
    throw "Parsing failed";
  }
  
  vector<AtomicPair> pairs;
  vector<AtomicPairPool*>::iterator pool;
  for(pool=pools.begin(); pool!=pools.end(); pool++)
  {
    (*pool)->invoke_correlators();
    pairs.insert(pairs.end(),(*pool)->pairs.begin(),(*pool)->pairs.end());
  }
  
//    pairs = cell.laue_symmetry.filter_pairs_from_asymmetric_unit(pairs); //We decided to use another way for multiplicity
  
  REPORT(FIRST_RUN) << "created " << pairs.size() << " pairs\n\n";
  
  // output pairs
  if(dump_pairs)
    for(int i=0; i<pairs.size(); i++)
      REPORT(FIRST_RUN) << pairs[i].to_string() <<'\n';

  
  if(report_pairs_outside_pdf_grid)
  {
    Grid pdf_grid = grid.in_pdf_space();
    for(int i=0; i<pairs.size(); i++)
      if(!pairs[i].pair_is_withing(pdf_grid))
        REPORT(FIRST_RUN) << "Warning, pair is outside PDF grid " << pairs[i].to_string() << '\n';
   }
  
  pairs = cell.laue_symmetry.apply_patterson_symmetry(pairs);
  atomic_pairs = pairs; //save them to check that everything is fine
  
  IntensityMap* calc_intensity_map;
  if(average_flag==AVERAGE)
    calc_intensity_map=&average_intensity_map;
  else
    calc_intensity_map=&intensity_map;
  

  
  if(direct_diffuse_scattering_calculation)
  {
    vec3<int> sym_boundary;
    for(int i=0; i<3; ++i)
      sym_boundary[i]=calc_intensity_map->size()[i]>1;
    
    IntensityMap padded = calc_intensity_map->padded(sym_boundary); //this will avoid the problem with symmetry.
    IntnsityCalculator::calculate_scattering_from_pairs(pairs,padded,average_flag);
    cell.laue_symmetry.apply_patterson_symmetry(padded);
    calc_intensity_map->copy_from_padded(sym_boundary,padded);
  }else
  {
    vec3<int> sym_boundary;
    for(int i=0; i<3; ++i)
      sym_boundary[i]=calc_intensity_map->size()[i]>1 && !(periodic_boundaries[i]);
    
    IntensityMap recipr_padded = calc_intensity_map->padded(padding); //for better fft accuracy
    recipr_padded.invert_grid();
    IntensityMap padded = recipr_padded.padded(sym_boundary); //for symmetry
    IntnsityCalculator::calculate_patterson_map_from_pairs_f(pairs,padded,average_flag,fft_grid_size,periodic_boundaries);
    cell.laue_symmetry.apply_patterson_symmetry(padded);
    recipr_padded.copy_from_padded(sym_boundary,padded);
    recipr_padded.invert();
    calc_intensity_map->copy_from_padded(padding,recipr_padded);
  }

  apply_resolution_function_if_possible(*calc_intensity_map);

  apply_reciprocal_space_multipliers_if_possible(*calc_intensity_map);
  calc_intensity_map->to_reciprocal();
  report.calculation_is_finished();
}