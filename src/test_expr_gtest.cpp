/*
 Tests for expr.hpp (ParameterBlock, ExprPtr, eval, eval_d),
 ExprFormulaParser, ParameterizedAtom, and Model parse-once behaviour.

 Google Test version — replaces the CxxTest test_expr.h.
*/

#include <gtest/gtest.h>

#include "expr.hpp"
#include "ExprFormulaParser.h"
#include "basic_classes.h"
#include "InputFileParser.h"

#include <cmath>
#include <stdexcept>

// Required global — defined in main.cpp for the binary, here for the test binary.
OutputHandler report;

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

TEST(ParameterBlock, AddAndEval)
{
    yell::ParameterBlock b;
    yell::ExprPtr x = b.add("x", 3.0);
    EXPECT_NEAR(3.0, x->eval(b.values()), 1e-12);
}

TEST(ParameterBlock, MultipleParams)
{
    yell::ParameterBlock b;
    yell::ExprPtr x = b.add("x", 1.0);
    yell::ExprPtr y = b.add("y", 2.0);
    EXPECT_EQ(2, b.size());
    EXPECT_NEAR(1.0, x->eval(b.values()), 1e-12);
    EXPECT_NEAR(2.0, y->eval(b.values()), 1e-12);
}

TEST(ParameterBlock, DuplicateThrows)
{
    yell::ParameterBlock b;
    b.add("x", 1.0);
    EXPECT_THROW(b.add("x", 2.0), std::invalid_argument);
}

TEST(ParameterBlock, UnknownThrows)
{
    yell::ParameterBlock b;
    EXPECT_THROW(b["no_such"], std::out_of_range);
}

TEST(ParameterBlock, IndexOf)
{
    yell::ParameterBlock b;
    b.add("Scale", 1.0);
    b.add("x", 0.5);
    EXPECT_EQ(0, b.index_of("Scale"));
    EXPECT_EQ(1, b.index_of("x"));
}

TEST(ParameterBlock, SetValues)
{
    yell::ParameterBlock b;
    yell::ExprPtr x = b.add("x", 0.0);
    yell::ExprPtr y = b.add("y", 0.0);
    Eigen::VectorXd v(2); v << 7.0, 9.0;
    b.set_values(v);
    EXPECT_NEAR(7.0, x->eval(b.values()), 1e-12);
    EXPECT_NEAR(9.0, y->eval(b.values()), 1e-12);
}

TEST(ParameterBlock, SetValuesSizeMismatch)
{
    yell::ParameterBlock b;
    b.add("x", 1.0);
    Eigen::VectorXd v(3); v << 1.0, 2.0, 3.0;
    EXPECT_THROW(b.set_values(v), std::invalid_argument);
}

TEST(ParameterBlock, SetByName)
{
    yell::ParameterBlock b;
    yell::ExprPtr x = b.add("x", 0.0);
    b.set("x", 5.5);
    EXPECT_NEAR(5.5, x->eval(b.values()), 1e-12);
}

TEST(ParameterBlock, Names)
{
    yell::ParameterBlock b;
    b.add("a", 1.0);
    b.add("b", 2.0);
    EXPECT_EQ(2u, b.names().size());
    EXPECT_EQ("a", b.names()[0]);
    EXPECT_EQ("b", b.names()[1]);
}

// ─────────────────────────────────────────────────────────────────────────────
// 2. ExprPtr — eval (plain double)
// ─────────────────────────────────────────────────────────────────────────────

TEST(ExprEval, Literal)
{
    auto p = make_p({});
    EXPECT_NEAR(3.14, yell::lit(3.14)->eval(p), 1e-12);
}

TEST(ExprEval, LiteralZero)
{
    auto p = make_p({});
    EXPECT_NEAR(0.0, yell::lit(0.0)->eval(p), 1e-12);
}

TEST(ExprEval, ParamRef)
{
    yell::ParameterBlock b;
    yell::ExprPtr x = b.add("x", 0.5);
    EXPECT_NEAR(0.5, x->eval(b.values()), 1e-12);
}

TEST(ExprEval, Add)
{
    auto p = make_p({2.0, 3.0});
    auto e = std::make_shared<yell::ParamRef>(0) + std::make_shared<yell::ParamRef>(1);
    EXPECT_NEAR(5.0, e->eval(p), 1e-12);
}

TEST(ExprEval, Sub)
{
    auto p = make_p({5.0, 2.0});
    auto e = std::make_shared<yell::ParamRef>(0) - std::make_shared<yell::ParamRef>(1);
    EXPECT_NEAR(3.0, e->eval(p), 1e-12);
}

TEST(ExprEval, Mul)
{
    auto p = make_p({4.0, 3.0});
    auto e = std::make_shared<yell::ParamRef>(0) * std::make_shared<yell::ParamRef>(1);
    EXPECT_NEAR(12.0, e->eval(p), 1e-12);
}

TEST(ExprEval, Div)
{
    auto p = make_p({9.0, 3.0});
    auto e = std::make_shared<yell::ParamRef>(0) / std::make_shared<yell::ParamRef>(1);
    EXPECT_NEAR(3.0, e->eval(p), 1e-12);
}

TEST(ExprEval, UnaryNeg)
{
    auto p = make_p({4.0});
    auto e = -std::make_shared<yell::ParamRef>(0);
    EXPECT_NEAR(-4.0, e->eval(p), 1e-12);
}

TEST(ExprEval, LiteralMixedArithmetic)
{
    auto p = make_p({2.0});
    auto x = std::make_shared<yell::ParamRef>(0);
    auto e = x * yell::lit(3.0) + yell::lit(1.0);
    EXPECT_NEAR(7.0, e->eval(p), 1e-12);
}

TEST(ExprEval, Sin)
{
    auto p = make_p({M_PI / 6.0});
    auto e = yell::sin(std::make_shared<yell::ParamRef>(0));
    EXPECT_NEAR(0.5, e->eval(p), 1e-9);
}

TEST(ExprEval, Cos)
{
    auto p = make_p({0.0});
    auto e = yell::cos(std::make_shared<yell::ParamRef>(0));
    EXPECT_NEAR(1.0, e->eval(p), 1e-12);
}

TEST(ExprEval, Exp)
{
    auto p = make_p({1.0});
    auto e = yell::exp(std::make_shared<yell::ParamRef>(0));
    EXPECT_NEAR(std::exp(1.0), e->eval(p), 1e-12);
}

TEST(ExprEval, Log)
{
    auto p = make_p({std::exp(2.0)});
    auto e = yell::log(std::make_shared<yell::ParamRef>(0));
    EXPECT_NEAR(2.0, e->eval(p), 1e-12);
}

TEST(ExprEval, Sqrt)
{
    auto p = make_p({4.0});
    auto e = yell::sqrt(std::make_shared<yell::ParamRef>(0));
    EXPECT_NEAR(2.0, e->eval(p), 1e-12);
}

TEST(ExprEval, AbsPositive)
{
    auto p = make_p({3.0});
    auto e = yell::abs(std::make_shared<yell::ParamRef>(0));
    EXPECT_NEAR(3.0, e->eval(p), 1e-12);
}

TEST(ExprEval, AbsNegative)
{
    auto p = make_p({-3.0});
    auto e = yell::abs(std::make_shared<yell::ParamRef>(0));
    EXPECT_NEAR(3.0, e->eval(p), 1e-12);
}

TEST(ExprEval, PowConstantExponent)
{
    auto p = make_p({3.0});
    auto e = yell::pow(std::make_shared<yell::ParamRef>(0), 2.0);
    EXPECT_NEAR(9.0, e->eval(p), 1e-12);
}

TEST(ExprEval, PowGeneralExponent)
{
    auto p = make_p({2.0, 3.0});
    auto e = yell::pow(std::make_shared<yell::ParamRef>(0),
                       std::make_shared<yell::ParamRef>(1));
    EXPECT_NEAR(8.0, e->eval(p), 1e-12);
}

TEST(ExprEval, DoubleOperatorMix)
{
    // 2.0 * x + 1.0 / x  at x=2 → 4 + 0.5 = 4.5
    auto p = make_p({2.0});
    auto x = std::make_shared<yell::ParamRef>(0);
    auto e = 2.0 * x + 1.0 / x;
    EXPECT_NEAR(4.5, e->eval(p), 1e-12);
}

// ─────────────────────────────────────────────────────────────────────────────
// 3. ExprPtr — eval_d (gradient)
// ─────────────────────────────────────────────────────────────────────────────

TEST(ExprGrad, LiteralGradIsZero)
{
    auto p = make_p({1.0, 2.0});
    auto d = yell::lit(5.0)->eval_d(p);
    EXPECT_NEAR(5.0, d.value(), 1e-12);
    EXPECT_NEAR(0.0, d.derivatives()[0], 1e-12);
    EXPECT_NEAR(0.0, d.derivatives()[1], 1e-12);
}

TEST(ExprGrad, ParamRefGrad)
{
    auto p = make_p({3.0, 7.0});
    auto d0 = std::make_shared<yell::ParamRef>(0)->eval_d(p);
    EXPECT_NEAR(3.0, d0.value(), 1e-12);
    EXPECT_NEAR(1.0, d0.derivatives()[0], 1e-12);
    EXPECT_NEAR(0.0, d0.derivatives()[1], 1e-12);

    auto d1 = std::make_shared<yell::ParamRef>(1)->eval_d(p);
    EXPECT_NEAR(7.0, d1.value(), 1e-12);
    EXPECT_NEAR(0.0, d1.derivatives()[0], 1e-12);
    EXPECT_NEAR(1.0, d1.derivatives()[1], 1e-12);
}

TEST(ExprGrad, AddGrad)
{
    // f = x + y  → df/dx=1, df/dy=1
    auto p = make_p({2.0, 3.0});
    auto e = std::make_shared<yell::ParamRef>(0) + std::make_shared<yell::ParamRef>(1);
    auto d = e->eval_d(p);
    EXPECT_NEAR(5.0, d.value(), 1e-12);
    EXPECT_NEAR(1.0, d.derivatives()[0], 1e-12);
    EXPECT_NEAR(1.0, d.derivatives()[1], 1e-12);
}

TEST(ExprGrad, MulGrad)
{
    // f = x * y  → df/dx=y, df/dy=x  at (2,3)
    auto p = make_p({2.0, 3.0});
    auto e = std::make_shared<yell::ParamRef>(0) * std::make_shared<yell::ParamRef>(1);
    auto d = e->eval_d(p);
    EXPECT_NEAR(6.0, d.value(), 1e-12);
    EXPECT_NEAR(3.0, d.derivatives()[0], 1e-12);
    EXPECT_NEAR(2.0, d.derivatives()[1], 1e-12);
}

TEST(ExprGrad, DivGrad)
{
    // f = x / y  → df/dx=1/y, df/dy=-x/y²  at (4,2)
    auto p = make_p({4.0, 2.0});
    auto e = std::make_shared<yell::ParamRef>(0) / std::make_shared<yell::ParamRef>(1);
    auto d = e->eval_d(p);
    EXPECT_NEAR(2.0,  d.value(), 1e-12);
    EXPECT_NEAR(0.5,  d.derivatives()[0], 1e-12);
    EXPECT_NEAR(-1.0, d.derivatives()[1], 1e-12);
}

TEST(ExprGrad, NegGrad)
{
    // f = -x  → df/dx = -1
    auto p = make_p({5.0});
    auto e = -std::make_shared<yell::ParamRef>(0);
    auto d = e->eval_d(p);
    EXPECT_NEAR(-5.0, d.value(), 1e-12);
    EXPECT_NEAR(-1.0, d.derivatives()[0], 1e-12);
}

TEST(ExprGrad, SinGrad)
{
    // f = sin(x)  → df/dx = cos(x)  at x = pi/4
    double xv = M_PI / 4.0;
    auto p = make_p({xv});
    auto e = yell::sin(std::make_shared<yell::ParamRef>(0));
    auto d = e->eval_d(p);
    EXPECT_NEAR(std::sin(xv), d.value(), 1e-12);
    EXPECT_NEAR(std::cos(xv), d.derivatives()[0], 1e-12);
}

TEST(ExprGrad, CosGrad)
{
    // f = cos(x)  → df/dx = -sin(x)  at x = 1
    double xv = 1.0;
    auto p = make_p({xv});
    auto e = yell::cos(std::make_shared<yell::ParamRef>(0));
    auto d = e->eval_d(p);
    EXPECT_NEAR(std::cos(xv),  d.value(), 1e-12);
    EXPECT_NEAR(-std::sin(xv), d.derivatives()[0], 1e-12);
}

TEST(ExprGrad, ExpGrad)
{
    // f = exp(x)  → df/dx = exp(x)  at x = 2
    double xv = 2.0;
    auto p = make_p({xv});
    auto e = yell::exp(std::make_shared<yell::ParamRef>(0));
    auto d = e->eval_d(p);
    EXPECT_NEAR(std::exp(xv), d.value(), 1e-12);
    EXPECT_NEAR(std::exp(xv), d.derivatives()[0], 1e-12);
}

TEST(ExprGrad, SqrtGrad)
{
    // f = sqrt(x)  → df/dx = 1/(2*sqrt(x))  at x = 4
    double xv = 4.0;
    auto p = make_p({xv});
    auto e = yell::sqrt(std::make_shared<yell::ParamRef>(0));
    auto d = e->eval_d(p);
    EXPECT_NEAR(2.0,  d.value(), 1e-12);
    EXPECT_NEAR(0.25, d.derivatives()[0], 1e-12);
}

TEST(ExprGrad, SqrtAtZeroGradIsNaNOrInf)
{
    auto p = make_p({0.0});
    auto e = yell::sqrt(std::make_shared<yell::ParamRef>(0));
    auto d = e->eval_d(p);
    EXPECT_TRUE(std::isnan(d.derivatives()[0]) || std::isinf(d.derivatives()[0]));
}

TEST(ExprGrad, AbsPositiveGrad)
{
    auto p = make_p({3.0});
    auto e = yell::abs(std::make_shared<yell::ParamRef>(0));
    auto d = e->eval_d(p);
    EXPECT_NEAR(1.0, d.derivatives()[0], 1e-12);
}

TEST(ExprGrad, AbsNegativeGrad)
{
    auto p = make_p({-3.0});
    auto e = yell::abs(std::make_shared<yell::ParamRef>(0));
    auto d = e->eval_d(p);
    EXPECT_NEAR(-1.0, d.derivatives()[0], 1e-12);
}

TEST(ExprGrad, LogGrad)
{
    // f = log(x)  → df/dx = 1/x  at x = e²
    double xv = std::exp(2.0);
    auto p = make_p({xv});
    auto e = yell::log(std::make_shared<yell::ParamRef>(0));
    auto d = e->eval_d(p);
    EXPECT_NEAR(2.0,    d.value(), 1e-12);
    EXPECT_NEAR(1.0/xv, d.derivatives()[0], 1e-12);
}

TEST(ExprGrad, PowGrad)
{
    // f = x^3  → df/dx = 3x²  at x = 2
    auto p = make_p({2.0});
    auto e = yell::pow(std::make_shared<yell::ParamRef>(0), 3.0);
    auto d = e->eval_d(p);
    EXPECT_NEAR(8.0,  d.value(), 1e-12);
    EXPECT_NEAR(12.0, d.derivatives()[0], 1e-12);
}

TEST(ExprGrad, ComplexExprGrad)
{
    // f = 2*x*x + 3*y - 1  → df/dx = 4x, df/dy = 3  at (2, 5)
    auto p = make_p({2.0, 5.0});
    auto x = std::make_shared<yell::ParamRef>(0);
    auto y = std::make_shared<yell::ParamRef>(1);
    auto e = 2.0 * x * x + 3.0 * y - yell::lit(1.0);
    auto d = e->eval_d(p);
    EXPECT_NEAR(2*4 + 3*5 - 1, d.value(), 1e-12);  // 8 + 15 - 1 = 22
    EXPECT_NEAR(8.0, d.derivatives()[0], 1e-12);
    EXPECT_NEAR(3.0, d.derivatives()[1], 1e-12);
}

TEST(ExprGrad, SharedRefUsedTwice)
{
    // f = x * x  → df/dx = 2x  at x = 3
    auto p = make_p({3.0});
    auto x = std::make_shared<yell::ParamRef>(0);
    auto e = x * x;
    auto d = e->eval_d(p);
    EXPECT_NEAR(9.0, d.value(), 1e-12);
    EXPECT_NEAR(6.0, d.derivatives()[0], 1e-12);
}

// ─────────────────────────────────────────────────────────────────────────────
// 4. ExprFormulaParser
// ─────────────────────────────────────────────────────────────────────────────

class ExprFormulaFixture : public ::testing::Test
{
protected:
    ExprFormulaParser efp;
    Eigen::VectorXd   p;

    void SetUp() override
    {
        std::vector<std::string> names = {"Scale", "x", "y"};
        std::vector<double>      vals  = {1.0, 0.5, 0.3};
        efp.initialize_refinable_variables(names, vals);
        p = Eigen::VectorXd(3); p << 1.0, 0.5, 0.3;
    }
};

TEST_F(ExprFormulaFixture, ParseConstant)
{
    auto e = parse_expr(efp, "3.14");
    EXPECT_NEAR(3.14, e->eval(p), 1e-10);
}

TEST_F(ExprFormulaFixture, ParseZero)
{
    auto e = parse_expr(efp, "0");
    EXPECT_NEAR(0.0, e->eval(p), 1e-12);
}

TEST_F(ExprFormulaFixture, ParseNegativeConstant)
{
    auto e = parse_expr(efp, "-2.5");
    EXPECT_NEAR(-2.5, e->eval(p), 1e-12);
}

TEST_F(ExprFormulaFixture, ParseVariable)
{
    auto e = parse_expr(efp, "x");
    EXPECT_NEAR(0.5, e->eval(p), 1e-12);
    Eigen::VectorXd p2(3); p2 << 1.0, 0.9, 0.3;
    EXPECT_NEAR(0.9, e->eval(p2), 1e-12);
}

TEST_F(ExprFormulaFixture, ParseAdd)
{
    auto e = parse_expr(efp, "x+1.0");
    EXPECT_NEAR(1.5, e->eval(p), 1e-12);
}

TEST_F(ExprFormulaFixture, ParseSub)
{
    auto e = parse_expr(efp, "x-y");
    EXPECT_NEAR(0.2, e->eval(p), 1e-9);
}

TEST_F(ExprFormulaFixture, ParseMul)
{
    auto e = parse_expr(efp, "x*2.0");
    EXPECT_NEAR(1.0, e->eval(p), 1e-12);
}

TEST_F(ExprFormulaFixture, ParseDiv)
{
    auto e = parse_expr(efp, "x/y");
    EXPECT_NEAR(0.5/0.3, e->eval(p), 1e-9);
}

TEST_F(ExprFormulaFixture, ParseParentheses)
{
    auto e = parse_expr(efp, "(x+y)*2.0");
    EXPECT_NEAR((0.5+0.3)*2.0, e->eval(p), 1e-12);
}

TEST_F(ExprFormulaFixture, ParseSin)
{
    auto e = parse_expr(efp, "sin(x)");
    EXPECT_NEAR(std::sin(0.5), e->eval(p), 1e-12);
}

TEST_F(ExprFormulaFixture, ParseCos)
{
    auto e = parse_expr(efp, "cos(x)");
    EXPECT_NEAR(std::cos(0.5), e->eval(p), 1e-12);
}

TEST_F(ExprFormulaFixture, ParseExp)
{
    auto e = parse_expr(efp, "exp(x)");
    EXPECT_NEAR(std::exp(0.5), e->eval(p), 1e-12);
}

TEST_F(ExprFormulaFixture, ParseLog)
{
    auto e = parse_expr(efp, "log(x)");
    EXPECT_NEAR(std::log(0.5), e->eval(p), 1e-12);
}

TEST_F(ExprFormulaFixture, ParseSqrt)
{
    auto e = parse_expr(efp, "sqrt(x)");
    EXPECT_NEAR(std::sqrt(0.5), e->eval(p), 1e-12);
}

TEST_F(ExprFormulaFixture, ParseAbs)
{
    auto e = parse_expr(efp, "abs(x)");
    EXPECT_NEAR(0.5, e->eval(p), 1e-12);
}

TEST_F(ExprFormulaFixture, ParseCompound)
{
    auto e = parse_expr(efp, "2.0*x*x+y");
    EXPECT_NEAR(2.0*0.25 + 0.3, e->eval(p), 1e-12);
}

TEST_F(ExprFormulaFixture, VariableGradient)
{
    auto e = parse_expr(efp, "x");
    auto d = e->eval_d(p);
    EXPECT_NEAR(0.5, d.value(), 1e-12);
    EXPECT_NEAR(0.0, d.derivatives()[0], 1e-12); // Scale
    EXPECT_NEAR(1.0, d.derivatives()[1], 1e-12); // x
    EXPECT_NEAR(0.0, d.derivatives()[2], 1e-12); // y
}

TEST_F(ExprFormulaFixture, AddGradient)
{
    auto e = parse_expr(efp, "x+y");
    auto d = e->eval_d(p);
    EXPECT_NEAR(0.8, d.value(), 1e-9);
    EXPECT_NEAR(0.0, d.derivatives()[0], 1e-12);
    EXPECT_NEAR(1.0, d.derivatives()[1], 1e-12);
    EXPECT_NEAR(1.0, d.derivatives()[2], 1e-12);
}

TEST_F(ExprFormulaFixture, MulGradient)
{
    auto e = parse_expr(efp, "x*y");
    auto d = e->eval_d(p);
    EXPECT_NEAR(0.5*0.3, d.value(), 1e-12);
    EXPECT_NEAR(0.3, d.derivatives()[1], 1e-12); // df/dx
    EXPECT_NEAR(0.5, d.derivatives()[2], 1e-12); // df/dy
}

TEST_F(ExprFormulaFixture, UnaryMinusGradient)
{
    auto e = parse_expr(efp, "-x");
    auto d = e->eval_d(p);
    EXPECT_NEAR(-0.5, d.value(), 1e-12);
    EXPECT_NEAR(-1.0, d.derivatives()[1], 1e-12);
}

TEST_F(ExprFormulaFixture, SinGradient)
{
    auto e = parse_expr(efp, "sin(x)");
    auto d = e->eval_d(p);
    EXPECT_NEAR(std::sin(0.5), d.value(), 1e-12);
    EXPECT_NEAR(std::cos(0.5), d.derivatives()[1], 1e-12);
}

TEST_F(ExprFormulaFixture, AssignmentConstant)
{
    efp.add_assignment("a", yell::lit(0.7));
    auto e = parse_expr(efp, "a");
    EXPECT_NEAR(0.7, e->eval(p), 1e-12);
}

TEST_F(ExprFormulaFixture, AssignmentDerived)
{
    yell::ExprPtr x_ref = efp.get_expr("x");
    efp.add_assignment("a", x_ref * yell::lit(2.0));
    auto e = parse_expr(efp, "a");
    EXPECT_NEAR(1.0, e->eval(p), 1e-12);  // x=0.5 → a=1.0

    Eigen::VectorXd p2(3); p2 << 1.0, 0.8, 0.3;
    EXPECT_NEAR(1.6, e->eval(p2), 1e-12);
}

// ─────────────────────────────────────────────────────────────────────────────
// 5. Atom::update_caches (replaces ParameterizedAtomData)
// ─────────────────────────────────────────────────────────────────────────────

class ParameterizedAtomFixture : public ::testing::Test
{
protected:
    cctbx::uctbx::unit_cell cubic_cell;
    cctbx::uctbx::unit_cell ortho_cell;

    void SetUp() override
    {
        cubic_cell = cctbx::uctbx::unit_cell(scitbx::af::tiny<double,6>(5,5,5,90,90,90));
        ortho_cell = cctbx::uctbx::unit_cell(scitbx::af::tiny<double,6>(2,3,4,90,90,90));
    }
};

TEST_F(ParameterizedAtomFixture, IsotropicConstantAtom)
{
    auto rm = cubic_cell.reciprocal_metrical_matrix();
    auto uiso = yell::lit(0.02);
    yell::ExprPtr U[6];
    for (int i = 0; i < 6; ++i) U[i] = uiso * yell::lit(rm[i]);
    Atom a("C", XRay, yell::lit(1.0),
           yell::lit(0.1), yell::lit(0.2), yell::lit(0.3),
           U[0], U[1], U[2], U[3], U[4], U[5]);

    Eigen::VectorXd p;
    a.update_caches(p);

    EXPECT_NEAR(1.0, a.occ_cache, 1e-12);
    EXPECT_NEAR(0.1, a.r_cache[0], 1e-12);
    EXPECT_NEAR(0.2, a.r_cache[1], 1e-12);
    EXPECT_NEAR(0.3, a.r_cache[2], 1e-12);
}

TEST_F(ParameterizedAtomFixture, IsotropicADPConversionCubic)
{
    // For cubic a=5Å, Uiso=0.025 Å²: U11 = Uiso / a² = 0.025/25 = 0.001
    auto rm = cubic_cell.reciprocal_metrical_matrix();
    auto uiso = yell::lit(0.025);
    yell::ExprPtr U[6];
    for (int i = 0; i < 6; ++i) U[i] = uiso * yell::lit(rm[i]);
    Atom a("C", XRay, yell::lit(1.0),
           yell::lit(0.0), yell::lit(0.0), yell::lit(0.0),
           U[0], U[1], U[2], U[3], U[4], U[5]);

    Eigen::VectorXd p;
    a.update_caches(p);

    EXPECT_NEAR(0.025 / (5.0*5.0), a.U_cache[0], 1e-10);
    EXPECT_NEAR(0.025 / (5.0*5.0), a.U_cache[1], 1e-10);
    EXPECT_NEAR(0.025 / (5.0*5.0), a.U_cache[2], 1e-10);
    EXPECT_NEAR(0.0, a.U_cache[3], 1e-10);
    EXPECT_NEAR(0.0, a.U_cache[4], 1e-10);
    EXPECT_NEAR(0.0, a.U_cache[5], 1e-10);
}

TEST_F(ParameterizedAtomFixture, IsotropicParameterizedPosition)
{
    yell::ParameterBlock block;
    block.add("Scale", 1.0);
    yell::ExprPtr x_expr = block.add("x", 0.25);

    auto rm = cubic_cell.reciprocal_metrical_matrix();
    auto uiso = yell::lit(0.01);
    yell::ExprPtr U[6];
    for (int i = 0; i < 6; ++i) U[i] = uiso * yell::lit(rm[i]);
    Atom a("C", XRay, yell::lit(1.0),
           x_expr, yell::lit(0.0), yell::lit(0.0),
           U[0], U[1], U[2], U[3], U[4], U[5]);

    a.update_caches(block.values());
    EXPECT_NEAR(0.25, a.r_cache[0], 1e-12);

    block.set("x", 0.75);
    a.update_caches(block.values());
    EXPECT_NEAR(0.75, a.r_cache[0], 1e-12);
}

TEST_F(ParameterizedAtomFixture, AnisotropicConstantAtom)
{
    // cubic a=5: astar = 1/5, so U_cache[i] = Uaniso[i] * astar_i * astar_j
    auto rm = cubic_cell.reciprocal_metrical_matrix();
    yell::ExprPtr U[6];
    for (int i = 0; i < 6; ++i) U[i] = yell::lit(0.04) * yell::lit(rm[i]);
    Atom a("C", XRay, yell::lit(1.0),
           yell::lit(0.1), yell::lit(0.2), yell::lit(0.3),
           U[0], U[1], U[2], U[3], U[4], U[5]);

    Eigen::VectorXd p;
    a.update_caches(p);

    EXPECT_NEAR(0.1, a.r_cache[0], 1e-12);
    EXPECT_NEAR(0.04 / 25.0, a.U_cache[0], 1e-12);
    EXPECT_NEAR(0.04 / 25.0, a.U_cache[1], 1e-12);
    EXPECT_NEAR(0.04 / 25.0, a.U_cache[2], 1e-12);
    EXPECT_NEAR(0.0, a.U_cache[3], 1e-12);
}

TEST_F(ParameterizedAtomFixture, AnisotropicADPConversionOrthorhombic)
{
    // ortho a=2, b=3, c=4: astar=1/2, bstar=1/3, cstar=1/4
    // Off-diagonals use astar*bstar etc (not rm[3..5] which are 0 for ortho cells)
    auto rm = ortho_cell.reciprocal_metrical_matrix();
    const double astar = std::sqrt(rm[0]);
    const double bstar = std::sqrt(rm[1]);
    const double cstar = std::sqrt(rm[2]);
    yell::ExprPtr U[6] = {
        yell::lit(1.0) * (astar * astar),
        yell::lit(1.0) * (bstar * bstar),
        yell::lit(1.0) * (cstar * cstar),
        yell::lit(1.0) * (astar * bstar),
        yell::lit(1.0) * (astar * cstar),
        yell::lit(1.0) * (bstar * cstar)
    };
    Atom a("C", XRay, yell::lit(1.0),
           yell::lit(0.0), yell::lit(0.0), yell::lit(0.0),
           U[0], U[1], U[2], U[3], U[4], U[5]);

    Eigen::VectorXd p;
    a.update_caches(p);

    EXPECT_NEAR(1.0/(2.0*2.0), a.U_cache[0], 1e-12); // U11/a²
    EXPECT_NEAR(1.0/(3.0*3.0), a.U_cache[1], 1e-12); // U22/b²
    EXPECT_NEAR(1.0/(4.0*4.0), a.U_cache[2], 1e-12); // U33/c²
    EXPECT_NEAR(1.0/(2.0*3.0), a.U_cache[3], 1e-12); // U12/(a*b)
    EXPECT_NEAR(1.0/(2.0*4.0), a.U_cache[4], 1e-12); // U13/(a*c)
    EXPECT_NEAR(1.0/(3.0*4.0), a.U_cache[5], 1e-12); // U23/(b*c)
}

TEST_F(ParameterizedAtomFixture, AnisotropicParameterizedU)
{
    yell::ParameterBlock block;
    block.add("Scale", 1.0);
    yell::ExprPtr u11_expr = block.add("U11", 0.04);

    auto rm = cubic_cell.reciprocal_metrical_matrix();
    yell::ExprPtr U[6];
    U[0] = u11_expr * yell::lit(rm[0]);
    for (int i = 1; i < 6; ++i) U[i] = yell::lit(0.04) * yell::lit(rm[i]);
    Atom a("C", XRay, yell::lit(1.0),
           yell::lit(0.0), yell::lit(0.0), yell::lit(0.0),
           U[0], U[1], U[2], U[3], U[4], U[5]);

    a.update_caches(block.values());
    EXPECT_NEAR(0.04/25.0, a.U_cache[0], 1e-12);

    block.set("U11", 0.09);
    a.update_caches(block.values());
    EXPECT_NEAR(0.09/25.0, a.U_cache[0], 1e-12);
}

TEST_F(ParameterizedAtomFixture, UpdateChangesMultiplier)
{
    yell::ParameterBlock block;
    block.add("Scale", 1.0);
    yell::ExprPtr mult = block.add("mult", 0.6);

    auto rm = cubic_cell.reciprocal_metrical_matrix();
    auto uiso = yell::lit(0.01);
    yell::ExprPtr U[6];
    for (int i = 0; i < 6; ++i) U[i] = uiso * yell::lit(rm[i]);
    Atom a("C", XRay, mult,
           yell::lit(0.0), yell::lit(0.0), yell::lit(0.0),
           U[0], U[1], U[2], U[3], U[4], U[5]);

    a.update_caches(block.values());
    EXPECT_NEAR(0.6, a.occ_cache, 1e-12);

    block.set("mult", 0.4);
    a.update_caches(block.values());
    EXPECT_NEAR(0.4, a.occ_cache, 1e-12);
}

TEST_F(ParameterizedAtomFixture, ComponentProbAndMultAreIndependent)
{
    // Verify the split: component_prob is the variant probability (set via
    // set_occupancy), mult_expr is the per-atom refinable multiplier.
    // get_occupancy() returns component_prob only; occ_cache = both combined.
    yell::ParameterBlock block;
    block.add("Scale", 1.0);
    yell::ExprPtr mult = block.add("pCu", 0.7);

    Atom a("C", XRay, mult,
           yell::lit(0.0), yell::lit(0.0), yell::lit(0.0),
           yell::lit(0.0), yell::lit(0.0), yell::lit(0.0),
           yell::lit(0.0), yell::lit(0.0), yell::lit(0.0));

    // Before set_occupancy: component_prob = lit(1.0)
    Eigen::VectorXd p = block.values();
    EXPECT_NEAR(1.0, a.get_occupancy()->eval(p), 1e-12);
    EXPECT_NEAR(0.7, a.mult_expr->eval(p), 1e-12);

    a.update_caches(p);
    EXPECT_NEAR(0.7, a.occ_cache, 1e-12);  // 1.0 * 0.7

    // After set_occupancy (simulating Variant assignment of prob=0.5):
    a.set_occupancy(yell::lit(0.5));
    EXPECT_NEAR(0.5, a.get_occupancy()->eval(p), 1e-12);  // component_prob only
    EXPECT_NEAR(0.7, a.mult_expr->eval(p), 1e-12);        // mult unchanged

    // full_occupancy() = component_prob * mult_expr
    EXPECT_NEAR(0.35, a.full_occupancy()->eval(p), 1e-12);

    // occ_cache reflects the combined value after update_caches
    a.update_caches(p);
    EXPECT_NEAR(0.35, a.occ_cache, 1e-12);
}

// ─────────────────────────────────────────────────────────────────────────────
// 6. Model: parse-once integration tests
// ─────────────────────────────────────────────────────────────────────────────
//
// Uses a 2-component variant (C + Void) with SubstitutionalCorrelation so
// pairs are actually generated and the intensity map is non-trivially zero.
// Refinable variables are stored at indices 1..N (index 0 = Scale).
// So params passed to calculate() must be {Scale, x, Uiso}.
// Correlation parameter (0.5) is a literal in the input string.

static std::string simple_model_str(double x_val, double Uiso_val)
{
    std::ostringstream oss;
    oss << "Cell 4 4 4  90 90 90\n"
        << "DiffuseScatteringGrid -1 -1 -1  1 1 1  3 3 3\n"
        << "CalculationMethod direct\n"
        << "LaueSymmetry -1\n"
        << "RefinableVariables [ x = " << x_val << " Uiso = " << Uiso_val << " ]\n"
        << "UnitCell [\n"
        << "  V = Variant [ (p=0.5) C 1 x 0 0 Uiso (p=0.5) Void ]\n"
        << "]\n"
        << "Correlations [\n"
        << "  [ (1,0,0) SubstitutionalCorrelation(V,V,0.5) ]\n"
        << "]\n";
    return oss.str();
}

TEST(ModelParseOnce, CalculatesWithoutCrash)
{
    Model m(simple_model_str(0.25, 0.01));
    EXPECT_NO_THROW(m.calculate({1.0, 0.25, 0.01}));
}

TEST(ModelParseOnce, ParsesOnlyOnce)
{
    Model m(simple_model_str(0.25, 0.01));
    m.calculate({1.0, 0.25, 0.01});
    EXPECT_TRUE(m.model_parsed_);
    m.calculate({1.0, 0.30, 0.02});
    EXPECT_TRUE(m.model_parsed_);
}

TEST(ModelParseOnce, DifferentParamsDifferentResult)
{
    // Two models with different initial parameters give different intensity maps.
    Model m1(simple_model_str(0.25, 0.01));
    Model m2(simple_model_str(0.45, 0.05));

    m1.calculate({1.0, 0.25, 0.01});
    m2.calculate({1.0, 0.45, 0.05});

    bool differs = false;
    for (int i = 0; i < m1.intensity_map.size_1d(); i++)
        if (std::abs(m1.intensity_map.at(i) - m2.intensity_map.at(i)) > 1e-10)
            { differs = true; break; }
    EXPECT_TRUE(differs);
}

TEST(ModelParseOnce, SecondCallWithDifferentParamsChangesResult)
{
    // Single model called twice with different params → different intensity maps.
    Model m(simple_model_str(0.25, 0.01));
    m.calculate({1.0, 0.25, 0.01});
    IntensityMap first = m.intensity_map;

    m.calculate({1.0, 0.45, 0.05});
    IntensityMap second = m.intensity_map;

    bool differs = false;
    for (int i = 0; i < first.size_1d(); i++)
        if (std::abs(first.at(i) - second.at(i)) > 1e-10)
            { differs = true; break; }
    EXPECT_TRUE(differs);
}

TEST(ModelParseOnce, MixedDerivativesKeyword)
{
    std::string ms = "Cell 4 4 4 90 90 90\nDiffuseScatteringGrid 0 0 0 1 1 1 1 1 1\nLaueSymmetry -1\nDerivatives mixed\nUnitCell []\nCorrelations []\n";
    Model m(ms);
    EXPECT_EQ(m.derivatives_mode, MIXED);
}

// ─────────────────────────────────────────────────────────────────────────────
// PattersonPeak tests
// ─────────────────────────────────────────────────────────────────────────────

TEST(PattersonPeakTests, ScattererListNonEmpty)
{
    Model m(simple_model_str(0.25, 0.01));
    m.calculate({1.0, 0.25, 0.01});

    ScattererList sl;
    EXPECT_GT(sl.size(), 0);
}

TEST(PattersonPeakTests, PeaksFromPairsCountMatchesPairs)
{
    Model m(simple_model_str(0.25, 0.01));
    m.calculate({1.0, 0.25, 0.01});
    Eigen::VectorXd q = Eigen::VectorXd::Map(m.refinement_parameters.data(), m.refinement_parameters.size());

    ScattererList sl;
    std::vector<PattersonPeak> full_peaks, avg_peaks;
    peaks_from_pairs(m.atomic_pairs, q, sl, full_peaks, avg_peaks);

    EXPECT_EQ((int)full_peaks.size(), (int)m.atomic_pairs.size());
    EXPECT_EQ((int)avg_peaks.size(),  (int)m.atomic_pairs.size());
}

TEST(PattersonPeakTests, FullPeakCoefficientMatchesPairP)
{
    Model m(simple_model_str(0.25, 0.01));
    m.calculate({1.0, 0.25, 0.01});
    Eigen::VectorXd q = Eigen::VectorXd::Map(m.refinement_parameters.data(), m.refinement_parameters.size());

    ScattererList sl;
    std::vector<PattersonPeak> full_peaks, avg_peaks;
    peaks_from_pairs(m.atomic_pairs, q, sl, full_peaks, avg_peaks);

    for (int k = 0; k < (int)m.atomic_pairs.size(); ++k) {
        double expected = m.atomic_pairs[k].p(false)->eval(q) * m.atomic_pairs[k].multiplier;
        EXPECT_NEAR(full_peaks[k].coefficient, expected, 1e-12);
    }
}

TEST(PattersonPeakTests, ValidScattererIndices)
{
    Model m(simple_model_str(0.25, 0.01));
    m.calculate({1.0, 0.25, 0.01});
    Eigen::VectorXd q = Eigen::VectorXd::Map(m.refinement_parameters.data(), m.refinement_parameters.size());

    ScattererList sl;
    std::vector<PattersonPeak> full_peaks, avg_peaks;
    peaks_from_pairs(m.atomic_pairs, q, sl, full_peaks, avg_peaks);

    for (auto& pk : full_peaks) {
        EXPECT_GE(pk.type1_idx, 0);
        EXPECT_GE(pk.type2_idx, 0);
        EXPECT_LT(pk.type1_idx, sl.size());
        EXPECT_LT(pk.type2_idx, sl.size());
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Pair probability regression tests
//
// Regression: multiplier corruption bug.
// MultiplicityCorrelation used to multiply both pair.p(false) and pair.p(true)
// by its scalar, AND multiply pair.multiplier. This caused pair probabilities
// to grow with symmetry multiplicity (e.g. 48 × p for m-3m), making scale
// refine to 0 via VarProj.
//
// Invariant: pair.p(false)->eval * pair.multiplier is the physical coefficient.
//   - pair.p(false)->eval must equal the joint probability set by SubstitutionalCorrelation
//     (or product of marginals for pair.p(true)), NOT that value times any multiplier.
//   - pair.multiplier absorbs all counting factors (LaueSymmetry, MultiplicityCorrelation).
// ─────────────────────────────────────────────────────────────────────────────

// Helper: build a model string with m-3m and a SubstitutionalCorrelation at (1,0,0).
// joint_prob = occupancy^2 + delta, so it differs visibly from the average (occupancy^2).
static std::string cubic_variant_model(double occupancy, double joint_prob_100)
{
    std::ostringstream oss;
    oss << "Cell 5 5 5  90 90 90\n"
        << "DiffuseScatteringGrid -1 -1 -1  1 1 1  3 3 3\n"
        << "CalculationMethod direct\n"
        << "LaueSymmetry m-3m\n"
        << "UnitCell [\n"
        << "  V = Variant [ (p=" << occupancy << ") Au 1 0 0 0 0.01 (p=" << (1.0 - occupancy) << ") Void ]\n"
        << "]\n"
        << "Correlations [\n"
        << "  [ (1,0,0) SubstitutionalCorrelation(V,V," << joint_prob_100 << ") ]\n"
        << "]\n";
    return oss.str();
}

// After calculate(), pair.p(false)->eval and pair.p(true)->eval must NOT be
// multiplied by the LaueSymmetry group multiplicity (48 for m-3m).
// Before the fix, pair.p(false)->eval was 48*(4/9) ≈ 21.3 instead of 4/9.
TEST(PairProbabilityRegression, PairProbabilityNotMultipliedByLaueMultiplicity)
{
    const double p_occ = 0.67;                  // marginal occupancy
    const double c_100 = 0.40;                  // joint probability at (1,0,0)
    const double avg   = p_occ * p_occ;         // product of marginals ≈ 0.449

    Model m(cubic_variant_model(p_occ, c_100));
    m.calculate({1.0});

    Eigen::VectorXd q = Eigen::VectorXd::Map(m.refinement_parameters.data(),
                                             m.refinement_parameters.size());

    // All pair probabilities must be physical (≤ 1.0).
    // With the bug they would be ≈ 48 × probability (e.g. 21.3), which is >> 1.
    for (AtomicPair& pair : m.atomic_pairs) {
        double p_full = pair.p(false)->eval(q);
        double p_avg  = pair.p(true)->eval(q);
        EXPECT_LE(std::abs(p_full), 1.0 + 1e-9)
            << "pair.p(false) = " << p_full << " — looks like it was multiplied by symmetry factor";
        EXPECT_LE(std::abs(p_avg),  1.0 + 1e-9)
            << "pair.p(true) = "  << p_avg  << " — looks like it was multiplied by symmetry factor";
    }

    // We verify the coefficient p * multiplier is reasonable (within a small factor of c_100).
    for (AtomicPair& pair : m.atomic_pairs) {
        double coeff_full = pair.p(false)->eval(q) * pair.multiplier;
        double coeff_avg  = pair.p(true)->eval(q)  * pair.multiplier;
        // Coefficients may legitimately be larger than 1 (multiplier accounts for copies).
        // But they must not be 48× the actual probability — catch the old bug.
        EXPECT_LT(coeff_full, 48.0 * std::max(c_100, avg) + 1e-6)
            << "coefficient = " << coeff_full << " exceeds 48×probability — multiplicity bug?";
        EXPECT_LT(coeff_avg,  48.0 * std::max(c_100, avg) + 1e-6)
            << "avg coefficient = " << coeff_avg;
    }
}

// Focused unit test: MultiplicityCorrelation must not touch pair.p(), only pair.multiplier.
// This is the atomic regression for the bug fixed in ec114f14.
TEST(PairProbabilityRegression, MultiplicityCorrelationDoesNotChangeProbability)
{
    // Set up an atom with occupancy 0.4 and a self-pair.
    Atom a("Au", 0.4, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    AtomicPairPool pool;
    AtomicPair& pair = pool.get_pair(&a, &a);

    // Set a custom joint probability via SubstitutionalCorrelation.
    pair.p(false) = yell::lit(0.3);   // full: correlated value
    pair.p(true)  = yell::lit(0.16);  // avg:  p^2 = 0.4^2

    const double mc_factor = 6.0;
    MultiplicityCorrelation mcor(mc_factor);
    mcor.modify_pairs(&pool);

    Eigen::VectorXd zero_p;
    // Probabilities must be UNCHANGED — only multiplier grows.
    EXPECT_NEAR(0.3,  pair.p(false)->eval(zero_p), 1e-12)
        << "p(false) was multiplied by MultiplicityCorrelation factor — bug!";
    EXPECT_NEAR(0.16, pair.p(true)->eval(zero_p),  1e-12)
        << "p(true) was multiplied by MultiplicityCorrelation factor — bug!";
    EXPECT_NEAR(mc_factor, pair.multiplier, 1e-12)
        << "multiplier should equal the MultiplicityCorrelation factor";
}

// ─────────────────────────────────────────────────────────────────────────────
// Derivative tests
// ─────────────────────────────────────────────────────────────────────────────

TEST(DerivativeTests, CalculateDerivativeConsistentWithJacobian)
{
    Model m(simple_model_str(0.25, 0.01));
    m.refine_in_asu_val = false; // Compare full maps
    std::vector<double> params = {1.5, 0.25, 0.01};
    m.calculate(params);

    // Compute full Jacobian (n_obs x n_params)
    OptionalIntensityMap wts; // empty weights = 1.0
    IntensityMap exp_map = m.intensity_map; // dummy experiment
    Eigen::MatrixXd J = m.compute_analytical_jacobian_direct(params, exp_map, wts);

    for (int j = 0; j < (int)params.size(); ++j) {
        IntensityMap dI = m.calculate_derivative(params, j);
        
        // J(ii, j) = -d(Model_i)/dp_j * w_i
        // calculate_derivative returns d(Model_i)/dp_j
        for (int i = 0; i < dI.size_1d(); ++i) {
            EXPECT_NEAR(J(i, j), -dI.at(i), 1e-10) << "Mismatch in param " << j << " pixel " << i;
        }
    }
}

TEST(DerivativeTests, MixedJacobianConsistentWithDirect)
{
    Model m(simple_model_str(0.25, 0.01));
    m.refine_in_asu_val = false; 
    std::vector<double> params = {1.5, 0.25, 0.01};
    
    OptionalIntensityMap wts;
    IntensityMap exp_map = m.intensity_map; 

    Eigen::MatrixXd J_direct = m.compute_analytical_jacobian_direct(params, exp_map, wts);
    Eigen::MatrixXd J_mixed  = m.compute_jacobian_mixed(params, exp_map, wts);

    ASSERT_EQ(J_direct.rows(), J_mixed.rows());
    ASSERT_EQ(J_direct.cols(), J_mixed.cols());

    for (int j = 0; j < J_direct.cols(); ++j) {
        for (int i = 0; i < J_direct.rows(); ++i) {
            EXPECT_NEAR(J_direct(i, j), J_mixed(i, j), 1e-3) << "Mismatch in param " << j << " pixel " << i;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Gram–Charlier atom-line parsing (Phase 2): GramCharlier3/4 tail, a* baking,
// optionality, and harmonic-atom invariance.
// ─────────────────────────────────────────────────────────────────────────────

static std::string gc_model_str(const std::string& gc_tail)
{
    // Non-cubic cell so the per-component a* products differ: a*=0.25, b*=0.2, c*=0.125
    std::ostringstream oss;
    oss << "Cell 4 5 8  90 90 90\n"
        << "DiffuseScatteringGrid -1 -1 -1  1 1 1  3 3 3\n"
        << "CalculationMethod direct\n"
        << "LaueSymmetry -1\n"
        << "RefinableVariables [ ]\n"
        << "UnitCell [\n"
        << "  V = Variant [ (p=0.5) C 1 0 0 0  0.01 0.01 0.01 0 0 0 " << gc_tail
        << " (p=0.5) Void ]\n"
        << "]\n"
        << "Correlations [\n"
        << "  [ (1,0,0) SubstitutionalCorrelation(V,V,0.5) ]\n"
        << "]\n";
    return oss.str();
}

static Atom* only_atom(Model& m)
{
    for (Atom* a : m.all_atoms_) if (a->label == "C") return a;
    return nullptr;
}

TEST(GramCharlierParse, AnisoAtomWithoutTailIsHarmonic)
{
    Model m(gc_model_str(""));   // plain anisotropic atom, no GC tail
    Atom* a = only_atom(m);
    ASSERT_NE(a, nullptr);
    EXPECT_FALSE(a->anharmonic);
}

TEST(GramCharlierParse, Order3TailBakesAstarProducts)
{
    // c111 = 2 (n=0), c123 = 3 (n=4); rest 0
    Model m(gc_model_str("GramCharlier3[ 2 0 0 0 3 0 0 0 0 0 ]"));
    Atom* a = only_atom(m);
    ASSERT_NE(a, nullptr);
    EXPECT_TRUE(a->anharmonic);

    const double as = 0.25, bs = 0.2, cs = 0.125;
    EXPECT_NEAR(a->C_cache.c[0], 2.0 * as * as * as, 1e-12);   // 111
    EXPECT_NEAR(a->C_cache.c[4], 3.0 * as * bs * cs, 1e-12);   // 123
    EXPECT_NEAR(a->C_cache.c[1], 0.0, 1e-15);                  // 112 untouched
    // no GramCharlier4 → D stays zero
    for (int n = 0; n < yell::GC4_N; ++n) EXPECT_NEAR(a->D_cache.d[n], 0.0, 1e-15);
}

TEST(GramCharlierParse, Order4TailBakesAstarProducts)
{
    // d1111 = 5 (n=0), d1122 = 7 (n=3); rest 0   (no GramCharlier3)
    Model m(gc_model_str("GramCharlier4[ 5 0 0 7 0 0 0 0 0 0 0 0 0 0 0 ]"));
    Atom* a = only_atom(m);
    ASSERT_NE(a, nullptr);
    EXPECT_TRUE(a->anharmonic);

    const double as = 0.25, bs = 0.2;
    EXPECT_NEAR(a->D_cache.d[0], 5.0 * as * as * as * as, 1e-12);  // 1111
    EXPECT_NEAR(a->D_cache.d[3], 7.0 * as * as * bs * bs, 1e-12);  // 1122
    // no GramCharlier3 → C stays zero
    for (int n = 0; n < yell::GC3_N; ++n) EXPECT_NEAR(a->C_cache.c[n], 0.0, 1e-15);
}

TEST(GramCharlierParse, BothTailsParseAndCalculate)
{
    Model m(gc_model_str(
        "GramCharlier3[ 1 0 0 0 0 0 0 0 0 0 ] "
        "GramCharlier4[ 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ]"));
    Atom* a = only_atom(m);
    ASSERT_NE(a, nullptr);
    EXPECT_TRUE(a->anharmonic);
    // forward calc still runs (G(s) not yet wired, but parsing must not break it)
    EXPECT_NO_THROW(m.calculate({1.0}));
}

// ─────────────────────────────────────────────────────────────────────────────
// Gram–Charlier pair combination (Phase 3a): C = C2 − C1, D = D1 + D2.
// ─────────────────────────────────────────────────────────────────────────────

static Atom make_anh_atom(const std::string& label, double c111, double d1111, bool anh)
{
    Atom a(label, 1.0, 0, 0, 0, 0.01, 0.01, 0.01, 0, 0, 0);  // aniso literal ctor
    if (anh) {
        a.C.c[0] = yell::lit(c111);     // 111
        a.D.d[0] = yell::lit(d1111);    // 1111
        a.anharmonic = true;
    }
    return a;
}

TEST(GramCharlierPair, CombinesOddSubtractEvenAdd)
{
    Atom a1 = make_anh_atom("C", /*c111*/ 2.0, /*d1111*/ 3.0, true);
    Atom a2 = make_anh_atom("C", /*c111*/ 5.0, /*d1111*/ 7.0, true);
    AtomicPair pair(a1, a2);
    EXPECT_TRUE(pair.anharmonic_);

    Eigen::VectorXd p = Eigen::VectorXd::Zero(1);
    yell::tensor3 C = pair.C(false).eval(p);
    yell::tensor4 D = pair.D(false).eval(p);
    EXPECT_NEAR(C.c[0], 5.0 - 2.0, 1e-12);   // C2 − C1
    EXPECT_NEAR(D.d[0], 3.0 + 7.0, 1e-12);   // D1 + D2

    // average track seeded identically to full (atomic GC enters both)
    EXPECT_NEAR(pair.C(true).eval(p).c[0], 3.0, 1e-12);
    EXPECT_NEAR(pair.D(true).eval(p).d[0], 10.0, 1e-12);
}

TEST(GramCharlierPair, HarmonicAtomsGiveHarmonicPair)
{
    Atom a1 = make_anh_atom("C", 0, 0, false);
    Atom a2 = make_anh_atom("C", 0, 0, false);
    AtomicPair pair(a1, a2);
    EXPECT_FALSE(pair.anharmonic_);

    Eigen::VectorXd p = Eigen::VectorXd::Zero(1);
    yell::tensor3 C = pair.C(false).eval(p);
    for (int n = 0; n < yell::GC3_N; ++n) EXPECT_NEAR(C.c[n], 0.0, 1e-15);
}

TEST(GramCharlierPair, OneAnharmonicAtomMarksPair)
{
    Atom a1 = make_anh_atom("C", 2.0, 3.0, true);
    Atom a2 = make_anh_atom("C", 0, 0, false);  // harmonic
    AtomicPair pair(a1, a2);
    EXPECT_TRUE(pair.anharmonic_);

    Eigen::VectorXd p = Eigen::VectorXd::Zero(1);
    EXPECT_NEAR(pair.C(false).eval(p).c[0], 0.0 - 2.0, 1e-12);  // C2(0) − C1(2)
    EXPECT_NEAR(pair.D(false).eval(p).d[0], 3.0 + 0.0, 1e-12);  // D1 + D2(0)
}

// ─────────────────────────────────────────────────────────────────────────────
// AnharmonicCorrelation modifier (Phase 3b): left/right modes + combined-basis
// coefficients assembled into the pair's spatial κ(Δu) tensor.
// ─────────────────────────────────────────────────────────────────────────────

TEST(AnharmonicCorrelationModifier, SelfPairRank3HandCheck)
{
    // K=2 self-pair: left=[modeL on Cu, +x], right=[modeR on Cu, +x].
    // e_left = -x, e_right = +x, so Δu_x = q_right - q_left.
    // coeffs over combined modes (n=3, K=2): [k000,k001,k011,k111] = [2,3,5,7].
    // spatial xxx = κ3(Δu_x) = k111 - 3k011 + 3k001 - k000 = 7 -15 +9 -2 = -1.
    Atom cu("Cu", 1.0, 0, 0, 0, 0.01, 0.01, 0.01, 0, 0, 0);
    ADPMode modeL; modeL.add_atom(&cu, vec3<double>(1, 0, 0));
    ADPMode modeR; modeR.add_atom(&cu, vec3<double>(1, 0, 0));

    std::vector<ADPMode*> left{&modeL}, right{&modeR};
    std::vector<yell::ExprPtr> coeffs{yell::lit(2), yell::lit(3), yell::lit(5), yell::lit(7)};
    AnharmonicCorrelation corr(left, right, 3, coeffs);

    AtomicPairPool pool;
    corr.modify_pairs(&pool);
    AtomicPair& pair = pool.get_pair(&cu, &cu);
    EXPECT_TRUE(pair.anharmonic_);

    Eigen::VectorXd p = Eigen::VectorXd::Zero(1);
    yell::tensor3 C = pair.C(false).eval(p);
    EXPECT_NEAR(C.c[0], -1.0, 1e-12);                 // xxx
    for (int n = 1; n < yell::GC3_N; ++n)             // only x-direction populated
        EXPECT_NEAR(C.c[n], 0.0, 1e-12) << "comp " << n;
}

TEST(AnharmonicCorrelationModifier, WrongCoeffCountThrows)
{
    Atom cu("Cu", 1.0, 0, 0, 0, 0.01, 0.01, 0.01, 0, 0, 0);
    ADPMode m; m.add_atom(&cu, vec3<double>(1, 0, 0));
    std::vector<ADPMode*> left{&m}, right{&m};   // K=2 → expects 4 coeffs for n=3
    std::vector<yell::ExprPtr> coeffs{yell::lit(1), yell::lit(2)};  // only 2
    EXPECT_THROW(AnharmonicCorrelation(left, right, 3, coeffs), std::string);
}

TEST(AnharmonicCorrelationModifier, ParsesAndCalculates)
{
    std::ostringstream oss;
    oss << "Cell 4 4 4  90 90 90\n"
        << "DiffuseScatteringGrid -1 -1 -1  1 1 1  3 3 3\n"
        << "CalculationMethod direct\n"
        << "LaueSymmetry -1\n"
        << "RefinableVariables [ ]\n"
        << "UnitCell [\n"
        << "  V = Variant [ (p=0.5) Cu = Cu 1 0 0 0  0.01 0.01 0.01 0 0 0 (p=0.5) Void ]\n"
        << "]\n"
        << "Modes[\n"
        << "  Cu_x = TranslationalMode(Cu, x)\n"
        << "]\n"
        << "Correlations [\n"
        << "  [ (1,0,0)\n"
        << "    SubstitutionalCorrelation(V,V,0.5)\n"
        << "    AnharmonicCorrelation3([Cu_x],[Cu_x], [ 0 0 0 1 ])\n"  // k111 = 1
        << "  ]\n"
        << "]\n";
    Model m(oss.str());
    EXPECT_NO_THROW(m.calculate({1.0}));
}
