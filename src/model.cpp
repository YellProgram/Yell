/*
 *  model.cpp
 *  diffuser_y
 *
 *  Created by Arkadiy Simonov on 3/15/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#include "model.h"
#include "InputFileParser.h"
#include <sstream>
extern OutputHandler report;
typedef iterator_ Iterator;

void complain_about_error(Iterator start, Iterator current, Iterator end)
{
  //skip empty spaces till the line where the error happened
  InputParser a_parser;
  qi::parse(current,end,a_parser.skipper_no_assignement);

  int line_number=0;
  int line_pos=0;
  
  Iterator line_start = current;
  while (line_start!=start && *line_start != '\n')
  {
    line_start--;
    line_pos++;
  }
  
  Iterator line_end = current;
  while (line_end!=end && *line_end != '\n')
    line_end++;
  
  for(Iterator it=current; it!=start; it--)
    if(*it=='\n')
      line_number++;

  std::string trouble_line(line_start,line_end);

  std::ostringstream complain;
  
  complain << "Parsing failed with exception around line " << line_number+1 << ":\n" << trouble_line <<"\n";
  for(int i=0; i<line_pos-1; ++i)
    complain << ' ';
  complain << "^-----roughly here\n\n";
  
  REPORT(ERROR) << complain.str();
  std::cout.flush();
  
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
    terminate();
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