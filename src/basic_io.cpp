//
//  basic_io.cpp
//  diffuser_y
//
//  Created by Arkadiy Simonov on 3/27/13.
//
//

#include "basic_io.h"
#include <string>
#include <fstream>

bool file_exists(string fname)
{
  ifstream ifile(fname.c_str());
  return ifile.good();
}

string read_input_file(string filename)
{
  string buf;
  string line;
  ifstream in(filename.c_str());
  while(getline(in,line))
    buf += line + "\n";
  
  return buf;
}