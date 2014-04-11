/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/ErrorPrinter.h>

int main() {
 return CxxTest::ErrorPrinter().run();
}
#include "test_correlation_generator.h"

static SymmetryGeneratorTestSuite suite_SymmetryGeneratorTestSuite;

static CxxTest::List Tests_SymmetryGeneratorTestSuite = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_SymmetryGeneratorTestSuite( "test_correlation_generator.h", 17, "SymmetryGeneratorTestSuite", suite_SymmetryGeneratorTestSuite, Tests_SymmetryGeneratorTestSuite );

static class TestDescription_SymmetryGeneratorTestSuite_test_generate_pairs : public CxxTest::RealTestDescription {
public:
 TestDescription_SymmetryGeneratorTestSuite_test_generate_pairs() : CxxTest::RealTestDescription( Tests_SymmetryGeneratorTestSuite, suiteDescription_SymmetryGeneratorTestSuite, 20, "test_generate_pairs" ) {}
 void runTest() { suite_SymmetryGeneratorTestSuite.test_generate_pairs(); }
} testDescription_SymmetryGeneratorTestSuite_test_generate_pairs;

static class TestDescription_SymmetryGeneratorTestSuite_test_fail : public CxxTest::RealTestDescription {
public:
 TestDescription_SymmetryGeneratorTestSuite_test_fail() : CxxTest::RealTestDescription( Tests_SymmetryGeneratorTestSuite, suiteDescription_SymmetryGeneratorTestSuite, 27, "test_fail" ) {}
 void runTest() { suite_SymmetryGeneratorTestSuite.test_fail(); }
} testDescription_SymmetryGeneratorTestSuite_test_fail;

#include <cxxtest/Root.cpp>
