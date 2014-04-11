//
//  basic_io.h
//  diffuser_y
//
//  Created by Arkadiy Simonov on 3/27/13.
//
//

#ifndef __diffuser_y__basic_io__
#define __diffuser_y__basic_io__

#include <iostream>
#include <string>
using std::string;
using std::ifstream;
/// checks if file file exists
bool file_exists(string filename);
/// reads the whole file in a string
string read_input_file(string filename);


#endif /* defined(__diffuser_y__basic_io__) */
