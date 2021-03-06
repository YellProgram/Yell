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
