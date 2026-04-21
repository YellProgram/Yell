/*
 Tests for expr.hpp (ParameterBlock, ExprPtr, eval, eval_d),
 ExprFormulaParser, ParameterizedAtom, and Model parse-once behaviour.

 Added alongside test_diffuser.h for the parameterized_model branch.
*/

#pragma once
#include <cxxtest/TestSuite.h>

// — expr.hpp ——————————————————————————————————————————————————————————————————
#include "expr.hpp"

// — ExprFormulaParser —————————————————————————————————————————————————————————
#include "ExprFormulaParser.h"

// — integration (Model parse-once) ————————————————————————————————————————————
#include "basic_classes.h"
#include "InputFileParser.h"

#include <cmath>
#include <stdexcept>

namespace qi = boost::spirit::qi;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

static Eigen::VectorXd make_p(std::initializer_list<double> vals)
{
    Eigen::VectorXd p(vals.size());
    int i = 0;
    for (double v : vals) p[i++] = v;
    return p;
}

// Parse a formula string with ExprFormulaParser and return the ExprPtr.
static yell::ExprPtr parse_expr(ExprFormulaParser& efp, const std::string& s)
{
    std::string input = s;
    std::string::iterator first = input.begin(), last = input.end();
    yell::ExprPtr result;
    bool ok = qi::parse(first, last, efp, result);
    if (!ok || first != last)
        throw std::runtime_error("ExprFormulaParser failed on: " + s);
    return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// 1. ParameterBlock
// ─────────────────────────────────────────────────────────────────────────────

class ParameterBlockTests : public CxxTest::TestSuite
{
public:
    void testAddAndEval()
    {
        yell::ParameterBlock b;
        yell::ExprPtr x = b.add("x", 3.0);
        TS_ASSERT_DELTA(3.0, x->eval(b.values()), 1e-12);
    }

    void testMultipleParams()
    {
        yell::ParameterBlock b;
        yell::ExprPtr x = b.add("x", 1.0);
        yell::ExprPtr y = b.add("y", 2.0);
        TS_ASSERT_EQUALS(2, b.size());
        TS_ASSERT_DELTA(1.0, x->eval(b.values()), 1e-12);
        TS_ASSERT_DELTA(2.0, y->eval(b.values()), 1e-12);
    }

    void testDuplicateThrows()
    {
        yell::ParameterBlock b;
        b.add("x", 1.0);
        TS_ASSERT_THROWS(b.add("x", 2.0), const std::invalid_argument&);
    }

    void testUnknownThrows()
    {
        yell::ParameterBlock b;
        TS_ASSERT_THROWS(b["no_such"], const std::out_of_range&);
    }

    void testIndexOf()
    {
        yell::ParameterBlock b;
        b.add("Scale", 1.0);
        b.add("x",     0.5);
        TS_ASSERT_EQUALS(0, b.index_of("Scale"));
        TS_ASSERT_EQUALS(1, b.index_of("x"));
    }

    void testSetValues()
    {
        yell::ParameterBlock b;
        yell::ExprPtr x = b.add("x", 0.0);
        yell::ExprPtr y = b.add("y", 0.0);
        Eigen::VectorXd v(2); v << 7.0, 9.0;
        b.set_values(v);
        TS_ASSERT_DELTA(7.0, x->eval(b.values()), 1e-12);
        TS_ASSERT_DELTA(9.0, y->eval(b.values()), 1e-12);
    }

    void testSetValuesSizeMismatch()
    {
        yell::ParameterBlock b;
        b.add("x", 1.0);
        Eigen::VectorXd v(3); v << 1.0, 2.0, 3.0;
        TS_ASSERT_THROWS(b.set_values(v), const std::invalid_argument&);
    }

    void testSetByName()
    {
        yell::ParameterBlock b;
        yell::ExprPtr x = b.add("x", 0.0);
        b.set("x", 5.5);
        TS_ASSERT_DELTA(5.5, x->eval(b.values()), 1e-12);
    }

    void testNames()
    {
        yell::ParameterBlock b;
        b.add("a", 1.0);
        b.add("b", 2.0);
        TS_ASSERT_EQUALS(2u, b.names().size());
        TS_ASSERT_EQUALS("a", b.names()[0]);
        TS_ASSERT_EQUALS("b", b.names()[1]);
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// 2. ExprPtr — eval (plain double)
// ─────────────────────────────────────────────────────────────────────────────

class ExprEvalTests : public CxxTest::TestSuite
{
public:
    void testLiteral()
    {
        auto p = make_p({});
        TS_ASSERT_DELTA(3.14, yell::lit(3.14)->eval(p), 1e-12);
    }

    void testLiteralZero()
    {
        auto p = make_p({});
        TS_ASSERT_DELTA(0.0, yell::lit(0.0)->eval(p), 1e-12);
    }

    void testParamRef()
    {
        yell::ParameterBlock b;
        yell::ExprPtr x = b.add("x", 0.5);
        TS_ASSERT_DELTA(0.5, x->eval(b.values()), 1e-12);
    }

    void testAdd()
    {
        auto p = make_p({2.0, 3.0});
        auto e = std::make_shared<yell::ParamRef>(0) + std::make_shared<yell::ParamRef>(1);
        TS_ASSERT_DELTA(5.0, e->eval(p), 1e-12);
    }

    void testSub()
    {
        auto p = make_p({5.0, 2.0});
        auto e = std::make_shared<yell::ParamRef>(0) - std::make_shared<yell::ParamRef>(1);
        TS_ASSERT_DELTA(3.0, e->eval(p), 1e-12);
    }

    void testMul()
    {
        auto p = make_p({4.0, 3.0});
        auto e = std::make_shared<yell::ParamRef>(0) * std::make_shared<yell::ParamRef>(1);
        TS_ASSERT_DELTA(12.0, e->eval(p), 1e-12);
    }

    void testDiv()
    {
        auto p = make_p({9.0, 3.0});
        auto e = std::make_shared<yell::ParamRef>(0) / std::make_shared<yell::ParamRef>(1);
        TS_ASSERT_DELTA(3.0, e->eval(p), 1e-12);
    }

    void testUnaryNeg()
    {
        auto p = make_p({4.0});
        auto e = -std::make_shared<yell::ParamRef>(0);
        TS_ASSERT_DELTA(-4.0, e->eval(p), 1e-12);
    }

    void testLiteralMixedArithmetic()
    {
        auto p = make_p({2.0});
        auto x = std::make_shared<yell::ParamRef>(0);
        auto e = x * yell::lit(3.0) + yell::lit(1.0);
        TS_ASSERT_DELTA(7.0, e->eval(p), 1e-12);
    }

    void testSin()
    {
        auto p = make_p({M_PI / 6.0});
        auto e = yell::sin(std::make_shared<yell::ParamRef>(0));
        TS_ASSERT_DELTA(0.5, e->eval(p), 1e-9);
    }

    void testCos()
    {
        auto p = make_p({0.0});
        auto e = yell::cos(std::make_shared<yell::ParamRef>(0));
        TS_ASSERT_DELTA(1.0, e->eval(p), 1e-12);
    }

    void testExp()
    {
        auto p = make_p({1.0});
        auto e = yell::exp(std::make_shared<yell::ParamRef>(0));
        TS_ASSERT_DELTA(std::exp(1.0), e->eval(p), 1e-12);
    }

    void testLog()
    {
        auto p = make_p({std::exp(2.0)});
        auto e = yell::log(std::make_shared<yell::ParamRef>(0));
        TS_ASSERT_DELTA(2.0, e->eval(p), 1e-12);
    }

    void testSqrt()
    {
        auto p = make_p({4.0});
        auto e = yell::sqrt(std::make_shared<yell::ParamRef>(0));
        TS_ASSERT_DELTA(2.0, e->eval(p), 1e-12);
    }

    void testAbsPositive()
    {
        auto p = make_p({3.0});
        auto e = yell::abs(std::make_shared<yell::ParamRef>(0));
        TS_ASSERT_DELTA(3.0, e->eval(p), 1e-12);
    }

    void testAbsNegative()
    {
        auto p = make_p({-3.0});
        auto e = yell::abs(std::make_shared<yell::ParamRef>(0));
        TS_ASSERT_DELTA(3.0, e->eval(p), 1e-12);
    }

    void testPowConstantExponent()
    {
        auto p = make_p({3.0});
        auto e = yell::pow(std::make_shared<yell::ParamRef>(0), 2.0);
        TS_ASSERT_DELTA(9.0, e->eval(p), 1e-12);
    }

    void testPowGeneralExponent()
    {
        auto p = make_p({2.0, 3.0});
        auto e = yell::pow(std::make_shared<yell::ParamRef>(0),
                           std::make_shared<yell::ParamRef>(1));
        TS_ASSERT_DELTA(8.0, e->eval(p), 1e-12);
    }

    void testDoubleOperatorMix()
    {
        // 2.0 * x + 1.0 / x  at x = 2 → 4 + 0.5 = 4.5
        auto p = make_p({2.0});
        auto x = std::make_shared<yell::ParamRef>(0);
        auto e = 2.0 * x + 1.0 / x;
        TS_ASSERT_DELTA(4.5, e->eval(p), 1e-12);
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// 3. ExprPtr — eval_d (gradient)
// ─────────────────────────────────────────────────────────────────────────────

class ExprGradTests : public CxxTest::TestSuite
{
public:
    void testLiteralGradIsZero()
    {
        auto p = make_p({1.0, 2.0});
        auto d = yell::lit(5.0)->eval_d(p);
        TS_ASSERT_DELTA(5.0, d.value(), 1e-12);
        TS_ASSERT_DELTA(0.0, d.derivatives()[0], 1e-12);
        TS_ASSERT_DELTA(0.0, d.derivatives()[1], 1e-12);
    }

    void testParamRefGrad()
    {
        auto p = make_p({3.0, 7.0});
        auto d0 = std::make_shared<yell::ParamRef>(0)->eval_d(p);
        TS_ASSERT_DELTA(3.0, d0.value(), 1e-12);
        TS_ASSERT_DELTA(1.0, d0.derivatives()[0], 1e-12);
        TS_ASSERT_DELTA(0.0, d0.derivatives()[1], 1e-12);

        auto d1 = std::make_shared<yell::ParamRef>(1)->eval_d(p);
        TS_ASSERT_DELTA(7.0, d1.value(), 1e-12);
        TS_ASSERT_DELTA(0.0, d1.derivatives()[0], 1e-12);
        TS_ASSERT_DELTA(1.0, d1.derivatives()[1], 1e-12);
    }

    void testAddGrad()
    {
        // f = x + y  → df/dx=1, df/dy=1
        auto p = make_p({2.0, 3.0});
        auto e = std::make_shared<yell::ParamRef>(0) + std::make_shared<yell::ParamRef>(1);
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(5.0, d.value(), 1e-12);
        TS_ASSERT_DELTA(1.0, d.derivatives()[0], 1e-12);
        TS_ASSERT_DELTA(1.0, d.derivatives()[1], 1e-12);
    }

    void testMulGrad()
    {
        // f = x * y  → df/dx=y, df/dy=x  at (2,3)
        auto p = make_p({2.0, 3.0});
        auto e = std::make_shared<yell::ParamRef>(0) * std::make_shared<yell::ParamRef>(1);
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(6.0, d.value(), 1e-12);
        TS_ASSERT_DELTA(3.0, d.derivatives()[0], 1e-12);
        TS_ASSERT_DELTA(2.0, d.derivatives()[1], 1e-12);
    }

    void testDivGrad()
    {
        // f = x / y  → df/dx=1/y, df/dy=-x/y²  at (4,2)
        auto p = make_p({4.0, 2.0});
        auto e = std::make_shared<yell::ParamRef>(0) / std::make_shared<yell::ParamRef>(1);
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(2.0,  d.value(), 1e-12);
        TS_ASSERT_DELTA(0.5,  d.derivatives()[0], 1e-12);
        TS_ASSERT_DELTA(-1.0, d.derivatives()[1], 1e-12);
    }

    void testNegGrad()
    {
        // f = -x  → df/dx = -1
        auto p = make_p({5.0});
        auto e = -std::make_shared<yell::ParamRef>(0);
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(-5.0, d.value(), 1e-12);
        TS_ASSERT_DELTA(-1.0, d.derivatives()[0], 1e-12);
    }

    void testSinGrad()
    {
        // f = sin(x)  → df/dx = cos(x)  at x = pi/4
        double xv = M_PI / 4.0;
        auto p = make_p({xv});
        auto e = yell::sin(std::make_shared<yell::ParamRef>(0));
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(std::sin(xv), d.value(), 1e-12);
        TS_ASSERT_DELTA(std::cos(xv), d.derivatives()[0], 1e-12);
    }

    void testCosGrad()
    {
        // f = cos(x)  → df/dx = -sin(x)  at x = 1
        double xv = 1.0;
        auto p = make_p({xv});
        auto e = yell::cos(std::make_shared<yell::ParamRef>(0));
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(std::cos(xv),  d.value(), 1e-12);
        TS_ASSERT_DELTA(-std::sin(xv), d.derivatives()[0], 1e-12);
    }

    void testExpGrad()
    {
        // f = exp(x)  → df/dx = exp(x)  at x = 2
        double xv = 2.0;
        auto p = make_p({xv});
        auto e = yell::exp(std::make_shared<yell::ParamRef>(0));
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(std::exp(xv), d.value(), 1e-12);
        TS_ASSERT_DELTA(std::exp(xv), d.derivatives()[0], 1e-12);
    }

    void testSqrtGrad()
    {
        // f = sqrt(x)  → df/dx = 1/(2*sqrt(x))  at x = 4
        double xv = 4.0;
        auto p = make_p({xv});
        auto e = yell::sqrt(std::make_shared<yell::ParamRef>(0));
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(2.0,  d.value(), 1e-12);
        TS_ASSERT_DELTA(0.25, d.derivatives()[0], 1e-12);
    }

    void testSqrtAtZeroGradIsNaN()
    {
        // sqrt'(0) is infinity; Eigen produces NaN
        auto p = make_p({0.0});
        auto e = yell::sqrt(std::make_shared<yell::ParamRef>(0));
        auto d = e->eval_d(p);
        TS_ASSERT(std::isnan(d.derivatives()[0]) || std::isinf(d.derivatives()[0]));
    }

    void testAbsPositiveGrad()
    {
        // abs'(x) = +1 for x > 0
        auto p = make_p({3.0});
        auto e = yell::abs(std::make_shared<yell::ParamRef>(0));
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(1.0, d.derivatives()[0], 1e-12);
    }

    void testAbsNegativeGrad()
    {
        // abs'(x) = -1 for x < 0
        auto p = make_p({-3.0});
        auto e = yell::abs(std::make_shared<yell::ParamRef>(0));
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(-1.0, d.derivatives()[0], 1e-12);
    }

    void testLogGrad()
    {
        // f = log(x)  → df/dx = 1/x  at x = e²
        double xv = std::exp(2.0);
        auto p = make_p({xv});
        auto e = yell::log(std::make_shared<yell::ParamRef>(0));
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(2.0,    d.value(), 1e-12);
        TS_ASSERT_DELTA(1.0/xv, d.derivatives()[0], 1e-12);
    }

    void testPowGrad()
    {
        // f = x^3  → df/dx = 3x²  at x = 2
        auto p = make_p({2.0});
        auto e = yell::pow(std::make_shared<yell::ParamRef>(0), 3.0);
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(8.0,  d.value(), 1e-12);
        TS_ASSERT_DELTA(12.0, d.derivatives()[0], 1e-12);
    }

    void testComplexExprGrad()
    {
        // f = 2*x*x + 3*y - 1  → df/dx = 4x, df/dy = 3  at (2, 5)
        auto p = make_p({2.0, 5.0});
        auto x = std::make_shared<yell::ParamRef>(0);
        auto y = std::make_shared<yell::ParamRef>(1);
        auto e = 2.0 * x * x + 3.0 * y - yell::lit(1.0);
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(2*4 + 3*5 - 1, d.value(), 1e-12);  // 8 + 15 - 1 = 22
        TS_ASSERT_DELTA(8.0, d.derivatives()[0], 1e-12);
        TS_ASSERT_DELTA(3.0, d.derivatives()[1], 1e-12);
    }

    void testSharedRefUsedTwice()
    {
        // f = x * x  → df/dx = 2x  at x = 3
        auto p = make_p({3.0});
        auto x = std::make_shared<yell::ParamRef>(0);
        auto e = x * x;  // same ExprPtr used twice
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(9.0, d.value(), 1e-12);
        TS_ASSERT_DELTA(6.0, d.derivatives()[0], 1e-12);
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// 4. ExprFormulaParser
// ─────────────────────────────────────────────────────────────────────────────

class ExprFormulaParserTests : public CxxTest::TestSuite
{
    ExprFormulaParser efp;
    Eigen::VectorXd   p;

public:
    void setUp()
    {
        // Set up two refinable variables: Scale(idx0)=1, x(idx1)=0.5, y(idx2)=0.3
        std::vector<std::string> names = {"Scale", "x", "y"};
        std::vector<double>      vals  = {1.0, 0.5, 0.3};
        efp.initialize_refinable_variables(names, vals);
        p = Eigen::VectorXd(3); p << 1.0, 0.5, 0.3;
    }

    void testParseConstant()
    {
        auto e = parse_expr(efp, "3.14");
        TS_ASSERT_DELTA(3.14, e->eval(p), 1e-10);
    }

    void testParseZero()
    {
        auto e = parse_expr(efp, "0");
        TS_ASSERT_DELTA(0.0, e->eval(p), 1e-12);
    }

    void testParseNegativeConstant()
    {
        auto e = parse_expr(efp, "-2.5");
        TS_ASSERT_DELTA(-2.5, e->eval(p), 1e-12);
    }

    void testParseVariable()
    {
        auto e = parse_expr(efp, "x");
        TS_ASSERT_DELTA(0.5, e->eval(p), 1e-12);
        // update p and re-eval — should track the new value
        Eigen::VectorXd p2(3); p2 << 1.0, 0.9, 0.3;
        TS_ASSERT_DELTA(0.9, e->eval(p2), 1e-12);
    }

    void testParseAdd()
    {
        auto e = parse_expr(efp, "x+1.0");
        TS_ASSERT_DELTA(1.5, e->eval(p), 1e-12);
    }

    void testParseSub()
    {
        auto e = parse_expr(efp, "x-y");
        TS_ASSERT_DELTA(0.2, e->eval(p), 1e-9);
    }

    void testParseMul()
    {
        auto e = parse_expr(efp, "x*2.0");
        TS_ASSERT_DELTA(1.0, e->eval(p), 1e-12);
    }

    void testParseDiv()
    {
        auto e = parse_expr(efp, "x/y");
        TS_ASSERT_DELTA(0.5/0.3, e->eval(p), 1e-9);
    }

    void testParseParentheses()
    {
        auto e = parse_expr(efp, "(x+y)*2.0");
        TS_ASSERT_DELTA((0.5+0.3)*2.0, e->eval(p), 1e-12);
    }

    void testParseSin()
    {
        auto e = parse_expr(efp, "sin(x)");
        TS_ASSERT_DELTA(std::sin(0.5), e->eval(p), 1e-12);
    }

    void testParseCos()
    {
        auto e = parse_expr(efp, "cos(x)");
        TS_ASSERT_DELTA(std::cos(0.5), e->eval(p), 1e-12);
    }

    void testParseExp()
    {
        auto e = parse_expr(efp, "exp(x)");
        TS_ASSERT_DELTA(std::exp(0.5), e->eval(p), 1e-12);
    }

    void testParseLog()
    {
        auto e = parse_expr(efp, "log(x)");
        TS_ASSERT_DELTA(std::log(0.5), e->eval(p), 1e-12);
    }

    void testParseSqrt()
    {
        auto e = parse_expr(efp, "sqrt(x)");
        TS_ASSERT_DELTA(std::sqrt(0.5), e->eval(p), 1e-12);
    }

    void testParseAbs()
    {
        auto e = parse_expr(efp, "abs(x)");
        TS_ASSERT_DELTA(0.5, e->eval(p), 1e-12);
    }

    void testParseCompound()
    {
        // 2.0 * x * x + y
        auto e = parse_expr(efp, "2.0*x*x+y");
        TS_ASSERT_DELTA(2.0*0.25 + 0.3, e->eval(p), 1e-12);
    }

    // Gradient tests via eval_d

    void testVariableGradient()
    {
        auto e = parse_expr(efp, "x");
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(0.5, d.value(), 1e-12);
        TS_ASSERT_DELTA(0.0, d.derivatives()[0], 1e-12); // Scale
        TS_ASSERT_DELTA(1.0, d.derivatives()[1], 1e-12); // x
        TS_ASSERT_DELTA(0.0, d.derivatives()[2], 1e-12); // y
    }

    void testAddGradient()
    {
        auto e = parse_expr(efp, "x+y");
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(0.8, d.value(), 1e-9);
        TS_ASSERT_DELTA(0.0, d.derivatives()[0], 1e-12);
        TS_ASSERT_DELTA(1.0, d.derivatives()[1], 1e-12);
        TS_ASSERT_DELTA(1.0, d.derivatives()[2], 1e-12);
    }

    void testMulGradient()
    {
        // f = x * y  → df/dx = y, df/dy = x
        auto e = parse_expr(efp, "x*y");
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(0.5*0.3, d.value(), 1e-12);
        TS_ASSERT_DELTA(0.3, d.derivatives()[1], 1e-12); // df/dx
        TS_ASSERT_DELTA(0.5, d.derivatives()[2], 1e-12); // df/dy
    }

    void testUnaryMinusGradient()
    {
        auto e = parse_expr(efp, "-x");
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(-0.5,  d.value(), 1e-12);
        TS_ASSERT_DELTA(-1.0,  d.derivatives()[1], 1e-12);
    }

    void testSinGradient()
    {
        auto e = parse_expr(efp, "sin(x)");
        auto d = e->eval_d(p);
        TS_ASSERT_DELTA(std::sin(0.5),  d.value(), 1e-12);
        TS_ASSERT_DELTA(std::cos(0.5),  d.derivatives()[1], 1e-12);
    }

    // Assignment: defines a derived constant; subsequent refs to it should work.
    void testAssignmentConstant()
    {
        // Mirror an assignment as if the skipper had fired
        efp.add_assignment("a", yell::lit(0.7));
        auto e = parse_expr(efp, "a");
        TS_ASSERT_DELTA(0.7, e->eval(p), 1e-12);
    }

    void testAssignmentDerived()
    {
        // a = x * 2;  then use a in another expr
        yell::ExprPtr x_ref = efp.get_expr("x");
        efp.add_assignment("a", x_ref * yell::lit(2.0));
        auto e = parse_expr(efp, "a");
        TS_ASSERT_DELTA(1.0, e->eval(p), 1e-12);  // x=0.5 → a=1.0

        // check gradient flows through
        Eigen::VectorXd p2(3); p2 << 1.0, 0.8, 0.3;
        TS_ASSERT_DELTA(1.6, e->eval(p2), 1e-12);
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// 5. Atom::update_caches (replaces removed ParameterizedAtomData tests)
// ─────────────────────────────────────────────────────────────────────────────

class ParameterizedAtomTests : public CxxTest::TestSuite
{
    cctbx::uctbx::unit_cell cubic_cell;
    cctbx::uctbx::unit_cell ortho_cell;

public:
    void setUp()
    {
        cubic_cell = cctbx::uctbx::unit_cell(scitbx::af::tiny<double,6>(5,5,5,90,90,90));
        ortho_cell = cctbx::uctbx::unit_cell(scitbx::af::tiny<double,6>(2,3,4,90,90,90));
    }

    // ── isotropic ──────────────────────────────────────────────────────────

    void testIsotropicConstantAtom()
    {
        // Build atom with literal ExprPtrs — same as construct_atom_isotropic_adp.
        auto rm = cubic_cell.reciprocal_metrical_matrix();
        auto uiso = yell::lit(0.02);
        yell::ExprPtr U[6];
        for (int i = 0; i < 6; ++i) U[i] = uiso * rm[i];
        Atom a("C", XRay, yell::lit(1.0),
               yell::lit(0.1), yell::lit(0.2), yell::lit(0.3),
               U[0], U[1], U[2], U[3], U[4], U[5]);

        Eigen::VectorXd p;
        a.update_caches(p);

        TS_ASSERT_DELTA(1.0, a.occ_cache, 1e-12);
        TS_ASSERT_DELTA(0.1, a.r_cache[0], 1e-12);
        TS_ASSERT_DELTA(0.2, a.r_cache[1], 1e-12);
        TS_ASSERT_DELTA(0.3, a.r_cache[2], 1e-12);
    }

    void testIsotropicADPConversionCubic()
    {
        // Cubic a=5Å, Uiso=0.025 Å²: U_frac = Uiso * rm[i]
        // For cubic: rm[0..2] = 1/25 = 0.04, rm[3..5] = 0
        auto rm = cubic_cell.reciprocal_metrical_matrix();
        auto uiso = yell::lit(0.025);
        yell::ExprPtr U[6];
        for (int i = 0; i < 6; ++i) U[i] = uiso * rm[i];
        Atom a("C", XRay, yell::lit(1.0),
               yell::lit(0.0), yell::lit(0.0), yell::lit(0.0),
               U[0], U[1], U[2], U[3], U[4], U[5]);

        Eigen::VectorXd p;
        a.update_caches(p);

        TS_ASSERT_DELTA(0.025 / (5.0*5.0), a.U_cache[0], 1e-10);
        TS_ASSERT_DELTA(0.025 / (5.0*5.0), a.U_cache[1], 1e-10);
        TS_ASSERT_DELTA(0.025 / (5.0*5.0), a.U_cache[2], 1e-10);
        TS_ASSERT_DELTA(0.0, a.U_cache[3], 1e-10);
        TS_ASSERT_DELTA(0.0, a.U_cache[4], 1e-10);
        TS_ASSERT_DELTA(0.0, a.U_cache[5], 1e-10);
    }

    void testIsotropicParameterizedPosition()
    {
        yell::ParameterBlock block;
        block.add("Scale", 1.0);
        yell::ExprPtr x_expr = block.add("x", 0.25);

        auto rm = cubic_cell.reciprocal_metrical_matrix();
        auto uiso = yell::lit(0.01);
        yell::ExprPtr U[6];
        for (int i = 0; i < 6; ++i) U[i] = uiso * rm[i];
        Atom a("C", XRay, yell::lit(1.0),
               x_expr, yell::lit(0.0), yell::lit(0.0),
               U[0], U[1], U[2], U[3], U[4], U[5]);

        a.update_caches(block.values());
        TS_ASSERT_DELTA(0.25, a.r_cache[0], 1e-12);

        block.set("x", 0.75);
        a.update_caches(block.values());
        TS_ASSERT_DELTA(0.75, a.r_cache[0], 1e-12);
    }

    // ── anisotropic ────────────────────────────────────────────────────────

    void testAnisotropicConstantAtom()
    {
        // params in Å²; ADP conversion baked into ExprPtrs as in construct_atom
        auto rm = cubic_cell.reciprocal_metrical_matrix();
        double astar = std::sqrt(rm[0]), bstar = std::sqrt(rm[1]), cstar = std::sqrt(rm[2]);
        yell::ExprPtr U[6] = {
            yell::lit(0.04) * (astar*astar), yell::lit(0.04) * (bstar*bstar),
            yell::lit(0.04) * (cstar*cstar), yell::lit(0.0)  * (astar*bstar),
            yell::lit(0.0)  * (astar*cstar), yell::lit(0.0)  * (bstar*cstar)
        };
        Atom a("C", XRay, yell::lit(1.0),
               yell::lit(0.1), yell::lit(0.2), yell::lit(0.3),
               U[0], U[1], U[2], U[3], U[4], U[5]);

        Eigen::VectorXd p;
        a.update_caches(p);

        TS_ASSERT_DELTA(0.1, a.r_cache[0], 1e-12);
        TS_ASSERT_DELTA(0.04 / 25.0, a.U_cache[0], 1e-12);
        TS_ASSERT_DELTA(0.04 / 25.0, a.U_cache[1], 1e-12);
        TS_ASSERT_DELTA(0.04 / 25.0, a.U_cache[2], 1e-12);
        TS_ASSERT_DELTA(0.0, a.U_cache[3], 1e-12);
    }

    void testAnisotropicADPConversionOrthorhombic()
    {
        auto rm = ortho_cell.reciprocal_metrical_matrix();
        double astar = std::sqrt(rm[0]), bstar = std::sqrt(rm[1]), cstar = std::sqrt(rm[2]);
        yell::ExprPtr U[6] = {
            yell::lit(1.0) * (astar*astar), yell::lit(1.0) * (bstar*bstar),
            yell::lit(1.0) * (cstar*cstar), yell::lit(1.0) * (astar*bstar),
            yell::lit(1.0) * (astar*cstar), yell::lit(1.0) * (bstar*cstar)
        };
        Atom a("C", XRay, yell::lit(1.0),
               yell::lit(0.0), yell::lit(0.0), yell::lit(0.0),
               U[0], U[1], U[2], U[3], U[4], U[5]);

        Eigen::VectorXd p;
        a.update_caches(p);

        TS_ASSERT_DELTA(1.0/(2.0*2.0), a.U_cache[0], 1e-12); // U11/a²
        TS_ASSERT_DELTA(1.0/(3.0*3.0), a.U_cache[1], 1e-12); // U22/b²
        TS_ASSERT_DELTA(1.0/(4.0*4.0), a.U_cache[2], 1e-12); // U33/c²
        TS_ASSERT_DELTA(1.0/(2.0*3.0), a.U_cache[3], 1e-12); // U12/(a*b)
        TS_ASSERT_DELTA(1.0/(2.0*4.0), a.U_cache[4], 1e-12); // U13/(a*c)
        TS_ASSERT_DELTA(1.0/(3.0*4.0), a.U_cache[5], 1e-12); // U23/(b*c)
    }

    void testAnisotropicParameterizedU()
    {
        yell::ParameterBlock block;
        block.add("Scale", 1.0);
        yell::ExprPtr u11_expr = block.add("U11", 0.04);

        auto rm = cubic_cell.reciprocal_metrical_matrix();
        double astar = std::sqrt(rm[0]), bstar = std::sqrt(rm[1]), cstar = std::sqrt(rm[2]);
        yell::ExprPtr U[6] = {
            u11_expr * (astar*astar), yell::lit(0.04) * (bstar*bstar),
            yell::lit(0.04) * (cstar*cstar), yell::lit(0.0) * (astar*bstar),
            yell::lit(0.0) * (astar*cstar),  yell::lit(0.0) * (bstar*cstar)
        };
        Atom a("C", XRay, yell::lit(1.0),
               yell::lit(0.0), yell::lit(0.0), yell::lit(0.0),
               U[0], U[1], U[2], U[3], U[4], U[5]);

        a.update_caches(block.values());
        TS_ASSERT_DELTA(0.04/25.0, a.U_cache[0], 1e-12);

        block.set("U11", 0.09);
        a.update_caches(block.values());
        TS_ASSERT_DELTA(0.09/25.0, a.U_cache[0], 1e-12);
    }

    void testUpdateChangesOccupancy()
    {
        yell::ParameterBlock block;
        block.add("Scale", 1.0);
        yell::ExprPtr mult = block.add("mult", 0.6);

        auto rm = cubic_cell.reciprocal_metrical_matrix();
        auto uiso = yell::lit(0.01);
        yell::ExprPtr U[6];
        for (int i = 0; i < 6; ++i) U[i] = uiso * rm[i];
        Atom a("C", XRay, mult,
               yell::lit(0.0), yell::lit(0.0), yell::lit(0.0),
               U[0], U[1], U[2], U[3], U[4], U[5]);

        a.update_caches(block.values());
        TS_ASSERT_DELTA(0.6, a.occ_cache, 1e-12);

        block.set("mult", 0.4);
        a.update_caches(block.values());
        TS_ASSERT_DELTA(0.4, a.occ_cache, 1e-12);
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// 6. Model: parse-once integration tests
// ─────────────────────────────────────────────────────────────────────────────

class ModelParseOnceTests : public CxxTest::TestSuite
{
    // Minimal valid Yell input for a simple model.
    // Uses the direct method so we don't need FFT grid.
    static std::string simple_model_str(double x_val, double Uiso_val)
    {
        std::ostringstream oss;
        // RefinableVariables: index 0=Scale, 1=x, 2=Uiso
        // Atom (isotropic, 5 params): mult x_coord y z Uiso
        // mult=1 (literal), x_coord=x (refinable), y=0, z=0, Uiso=Uiso (refinable)
        oss << "Cell 4 4 4  90 90 90\n"
            << "DiffuseScatteringGrid -1 -1 -1  1 1 1  3 3 3\n"
            << "CalculationMethod direct\n"
            << "LaueSymmetry -1\n"
            << "RefinableVariables [ x = " << x_val << " Uiso = " << Uiso_val << " ]\n"
            << "UnitCell [\n"
            << "  V = Variant [ (p=1) C 1 x 0 0 Uiso ]\n"
            << "]\n"
            << "Correlations [\n"
            << "  [ (0,0,0) Multiplicity 1 ]\n"
            << "]\n";
        return oss.str();
    }

public:
    void testModelCalculatesWithoutCrash()
    {
        Model m(simple_model_str(0.25, 0.01));
        TS_ASSERT_THROWS_NOTHING(m.calculate({1.0, 0.25, 0.01}));
    }

    void testModelParsesOnlyOnce()
    {
        // Call calculate() twice and verify model_parsed_ flag is set.
        Model m(simple_model_str(0.25, 0.01));
        m.calculate({1.0, 0.25, 0.01});
        TS_ASSERT(m.model_parsed_);
        m.calculate({1.0, 0.30, 0.02});
        // Still true — should not have been reset or re-parsed
        TS_ASSERT(m.model_parsed_);
    }

    void testModelParameterChangeAffectsResult()
    {
        // Two separate models with different params should give different intensities.
        Model m1(simple_model_str(0.25, 0.01));
        Model m2(simple_model_str(0.45, 0.05));

        m1.calculate({1.0, 0.25, 0.01});
        m2.calculate({1.0, 0.45, 0.05});

        // They should differ at some point.
        bool differs = false;
        for (int i = 0; i < m1.intensity_map.size_1d(); i++)
            if (std::abs(m1.intensity_map.at(i) - m2.intensity_map.at(i)) > 1e-10)
                { differs = true; break; }
        TS_ASSERT(differs);
    }

    void testModelSecondCallDifferentFromFirst()
    {
        // Single model called with two different param vectors → different results.
        Model m(simple_model_str(0.25, 0.01));
        m.calculate({1.0, 0.25, 0.01});
        IntensityMap first = m.intensity_map;

        m.calculate({1.0, 0.45, 0.05});
        IntensityMap second = m.intensity_map;

        bool differs = false;
        for (int i = 0; i < first.size_1d(); i++)
            if (std::abs(first.at(i) - second.at(i)) > 1e-10)
                { differs = true; break; }
        TS_ASSERT(differs);
    }
};
