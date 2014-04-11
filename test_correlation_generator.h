//
//  test_symmetry_generator.h
//  diffuser_y
//
//  Created by Arkadiy Simonov on 5/3/12.
//  Copyright (c) 2012 ETH Zurich. All rights reserved.
//



#ifndef diffuser_y_test_correlation_generator_h
#define diffuser_y_test_correlation_generator_h

#include "CorrelationGenerator.cpp"
#include <cxxtest/TestSuite.h>

class SymmetryGeneratorTestSuite : public CxxTest::TestSuite
{
public:
    void test_generate_pairs()
    {
        //vector<Pair> pairs = generate_pairs(
        //TS_ASSERT_EQUALS(<#x#>, <#y#>)
    }
    
    void test_fail()
    {
        TS_FAIL("test successefully failed!");    
    }

};

#endif
