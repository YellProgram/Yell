/*
 *  ReferenceTable.h
 *  diffuser_y
 *
 *  Created by Arkadiy Simonov on 3/22/11.
 *  Copyright 2011 ETH Zurich. All rights reserved.
 *
 */

#ifndef REFERENCE_TABLE_H
#define REFERENCE_TABLE_H

#include "InputFileParser.h"

struct ReferenceTable : qi::grammar<Iterator, StructurePartRef()>{
  
  ReferenceTable();
  
  typedef qi::symbols<char,StructurePartRef> _table;                                                                                         
  _table references; 
  
  qi::rule<Iterator,StructurePartRef()> identifier;
  
  static void add_reference(_table& table, string name,StructurePartRef references)
  {
    table.add(name,references);
  }
};

#endif