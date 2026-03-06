#pragma once
// Forward-mode automatic differentiation via expression trees.
//
// All differentiable model parameters are registered in a ParameterBlock and
// referenced through expression nodes (ExprPtr).  Each node evaluates to either
//   • a plain double  via  expr->eval(params.values())
//   • a Dual number   via  expr->eval_d(params.values())
//     where Dual carries (value, gradient w.r.t. every parameter).
//
// Example:
//   ParameterBlock params;
//   ExprPtr x = params.add("x", 0.5);
//   ExprPtr y = params.add("y", 0.1);
//   ExprPtr e = x * x + sin(y);
//
//   double v   = e->eval(params.values());      // plain evaluation
//   Dual   d   = e->eval_d(params.values());    // value + gradient
//   // d.value()           == v
//   // d.derivatives()[0]  == ∂e/∂x  at current values
//   // d.derivatives()[1]  == ∂e/∂y

#include <Eigen/Core>
#include <unsupported/Eigen/AutoDiff>
#include <string>
#include <memory>
#include <vector>
#include <unordered_map>
#include <array>
#include <cmath>
#include <stdexcept>

namespace yell {

// ── Dual number ───────────────────────────────────────────────────────────────
// Dynamic-size forward-mode dual: value + gradient w.r.t. all parameters.
// Backed by Eigen::AutoDiffScalar so matrix/vector operations with
// Eigen::Matrix<Dual, N, M> propagate derivatives automatically.
using Dual = Eigen::AutoDiffScalar<Eigen::VectorXd>;

// ── Expr base ─────────────────────────────────────────────────────────────────
class Expr;
using ExprPtr = std::shared_ptr<const Expr>;

class Expr {
public:
    virtual ~Expr() = default;

    // Fast evaluation — no gradient computation.
    virtual double eval(const Eigen::VectorXd& p) const = 0;

    // Evaluation with gradient.  p.size() sets the number of derivatives.
    virtual Dual   eval_d(const Eigen::VectorXd& p) const = 0;
};

// ── Literal (constant) ────────────────────────────────────────────────────────
class Literal : public Expr {
    double value_;
public:
    explicit Literal(double v) : value_(v) {}
    double eval(const Eigen::VectorXd&) const override { return value_; }
    Dual   eval_d(const Eigen::VectorXd& p) const override {
        return Dual(value_, Eigen::VectorXd::Zero(p.size()));
    }
};

inline ExprPtr lit(double v) { return std::make_shared<Literal>(v); }

// ── ParamRef (leaf — refers to one named parameter) ───────────────────────────
class ParamRef : public Expr {
    int idx_;
public:
    explicit ParamRef(int idx) : idx_(idx) {}
    double eval(const Eigen::VectorXd& p) const override { return p[idx_]; }
    Dual   eval_d(const Eigen::VectorXd& p) const override {
        // Constructor (value, nbDerivatives, activeIndex) sets gradient = eᵢ.
        return Dual(p[idx_], p.size(), idx_);
    }
};

// ── BinaryExpr ────────────────────────────────────────────────────────────────
enum class BinOp { ADD, SUB, MUL, DIV };

class BinaryExpr : public Expr {
    ExprPtr lhs_, rhs_;
    BinOp   op_;
public:
    BinaryExpr(ExprPtr l, ExprPtr r, BinOp op) : lhs_(l), rhs_(r), op_(op) {}

    double eval(const Eigen::VectorXd& p) const override {
        double l = lhs_->eval(p), r = rhs_->eval(p);
        switch (op_) {
            case BinOp::ADD: return l + r;
            case BinOp::SUB: return l - r;
            case BinOp::MUL: return l * r;
            case BinOp::DIV: return l / r;
        }
        return 0.0;
    }

    Dual eval_d(const Eigen::VectorXd& p) const override {
        // Eigen::AutoDiffScalar defines +, -, *, / with correct chain rules.
        Dual l = lhs_->eval_d(p), r = rhs_->eval_d(p);
        switch (op_) {
            case BinOp::ADD: return l + r;
            case BinOp::SUB: return l - r;
            case BinOp::MUL: return l * r;
            case BinOp::DIV: return l / r;
        }
        return Dual(0.0, Eigen::VectorXd::Zero(p.size()));
    }
};

// ── NegateExpr ────────────────────────────────────────────────────────────────
class NegateExpr : public Expr {
    ExprPtr arg_;
public:
    explicit NegateExpr(ExprPtr a) : arg_(a) {}
    double eval(const Eigen::VectorXd& p) const override { return -arg_->eval(p); }
    Dual   eval_d(const Eigen::VectorXd& p) const override { return -arg_->eval_d(p); }
};

// ── FuncExpr (unary math functions) ──────────────────────────────────────────
// All derivative rules are implemented manually using the
//   Dual(value, coeff * gradient)
// constructor pattern.  This avoids relying on which std:: or Eigen::
// overloads happen to be available for AutoDiffScalar in a given Eigen version.
//
// Corner cases (same as Ceres Jet):
//   sqrt(0)  → NaN gradient (derivative is genuinely infinite)
//   log(≤0)  → NaN gradient
//   abs(0)   → right-hand sub-gradient (+1 · g) — kink is not differentiable;
//               use sqrt(x*x + ε) if smoothness is needed
//   asin/acos(±1) → NaN gradient (singularity)

enum class MathFunc {
    SIN, COS, TAN,
    EXP, LOG, SQRT, ABS,
    ASIN, ACOS, ATAN,
    SINH, COSH, TANH
};

class FuncExpr : public Expr {
    ExprPtr  arg_;
    MathFunc func_;

    static Dual apply(MathFunc f, const Dual& x) {
        const double           v = x.value();
        const Eigen::VectorXd& g = x.derivatives();  // const ref — no copy

        switch (f) {
            case MathFunc::SIN:  return Dual(std::sin(v),   std::cos(v) * g);
            case MathFunc::COS:  return Dual(std::cos(v),  -std::sin(v) * g);
            case MathFunc::TAN: {
                double c = std::cos(v);
                return Dual(std::tan(v), g / (c * c));
            }
            case MathFunc::EXP: {
                double ev = std::exp(v);
                return Dual(ev, ev * g);
            }
            case MathFunc::LOG:  return Dual(std::log(v),  g / v);
            case MathFunc::SQRT: return Dual(std::sqrt(v), g / (2.0 * std::sqrt(v)));
            case MathFunc::ABS:  return Dual(std::abs(v),  (v >= 0.0 ? 1.0 : -1.0) * g);
            case MathFunc::ASIN: return Dual(std::asin(v),  g / std::sqrt(1.0 - v * v));
            case MathFunc::ACOS: return Dual(std::acos(v), -g / std::sqrt(1.0 - v * v));
            case MathFunc::ATAN: return Dual(std::atan(v),  g / (1.0 + v * v));
            case MathFunc::SINH: return Dual(std::sinh(v),  std::cosh(v) * g);
            case MathFunc::COSH: return Dual(std::cosh(v),  std::sinh(v) * g);
            case MathFunc::TANH: {
                double t = std::tanh(v);
                return Dual(t, (1.0 - t * t) * g);
            }
        }
        return Dual(0.0, Eigen::VectorXd::Zero(g.size()));
    }

    static double apply_plain(MathFunc f, double x) {
        switch (f) {
            case MathFunc::SIN:  return std::sin(x);
            case MathFunc::COS:  return std::cos(x);
            case MathFunc::TAN:  return std::tan(x);
            case MathFunc::EXP:  return std::exp(x);
            case MathFunc::LOG:  return std::log(x);
            case MathFunc::SQRT: return std::sqrt(x);
            case MathFunc::ABS:  return std::abs(x);
            case MathFunc::ASIN: return std::asin(x);
            case MathFunc::ACOS: return std::acos(x);
            case MathFunc::ATAN: return std::atan(x);
            case MathFunc::SINH: return std::sinh(x);
            case MathFunc::COSH: return std::cosh(x);
            case MathFunc::TANH: return std::tanh(x);
        }
        return 0.0;
    }

public:
    FuncExpr(ExprPtr a, MathFunc f) : arg_(a), func_(f) {}
    double eval(const Eigen::VectorXd& p)   const override { return apply_plain(func_, arg_->eval(p)); }
    Dual   eval_d(const Eigen::VectorXd& p) const override { return apply(func_, arg_->eval_d(p)); }
};

// ── PowExpr (base^constant) ───────────────────────────────────────────────────
// d/dx  x^n  =  n · x^(n−1)
class PowExpr : public Expr {
    ExprPtr base_;
    double  exp_;
public:
    PowExpr(ExprPtr b, double e) : base_(b), exp_(e) {}
    double eval(const Eigen::VectorXd& p) const override {
        return std::pow(base_->eval(p), exp_);
    }
    Dual eval_d(const Eigen::VectorXd& p) const override {
        Dual b = base_->eval_d(p);
        double v = b.value();
        return Dual(std::pow(v, exp_),
                    exp_ * std::pow(v, exp_ - 1.0) * b.derivatives());
    }
};

// ── GeneralPowExpr (base^exponent, both expressions) ─────────────────────────
// d(u^v) = u^v · (v/u · du + ln(u) · dv)   requires u > 0
class GeneralPowExpr : public Expr {
    ExprPtr base_, exp_;
public:
    GeneralPowExpr(ExprPtr b, ExprPtr e) : base_(b), exp_(e) {}
    double eval(const Eigen::VectorXd& p) const override {
        return std::pow(base_->eval(p), exp_->eval(p));
    }
    Dual eval_d(const Eigen::VectorXd& p) const override {
        Dual b = base_->eval_d(p);
        Dual e = exp_->eval_d(p);
        double u = b.value(), v = e.value();
        double uv = std::pow(u, v);
        Eigen::VectorXd deriv =
            uv * (v / u * b.derivatives() + std::log(u) * e.derivatives());
        return Dual(uv, deriv);
    }
};

// ── Operator overloads on ExprPtr ─────────────────────────────────────────────
// ExprPtr × ExprPtr
inline ExprPtr operator+(const ExprPtr& l, const ExprPtr& r) { return std::make_shared<BinaryExpr>(l, r, BinOp::ADD); }
inline ExprPtr operator-(const ExprPtr& l, const ExprPtr& r) { return std::make_shared<BinaryExpr>(l, r, BinOp::SUB); }
inline ExprPtr operator*(const ExprPtr& l, const ExprPtr& r) { return std::make_shared<BinaryExpr>(l, r, BinOp::MUL); }
inline ExprPtr operator/(const ExprPtr& l, const ExprPtr& r) { return std::make_shared<BinaryExpr>(l, r, BinOp::DIV); }
inline ExprPtr operator-(const ExprPtr& a)                   { return std::make_shared<NegateExpr>(a); }

// ExprPtr × double  (auto-promote constant to Literal)
inline ExprPtr operator+(const ExprPtr& l, double r) { return l + lit(r); }
inline ExprPtr operator+(double l, const ExprPtr& r) { return lit(l) + r; }
inline ExprPtr operator-(const ExprPtr& l, double r) { return l - lit(r); }
inline ExprPtr operator-(double l, const ExprPtr& r) { return lit(l) - r; }
inline ExprPtr operator*(const ExprPtr& l, double r) { return l * lit(r); }
inline ExprPtr operator*(double l, const ExprPtr& r) { return lit(l) * r; }
inline ExprPtr operator/(const ExprPtr& l, double r) { return l / lit(r); }
inline ExprPtr operator/(double l, const ExprPtr& r) { return lit(l) / r; }

// ── Math functions returning ExprPtr ─────────────────────────────────────────
inline ExprPtr sin (const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::SIN);  }
inline ExprPtr cos (const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::COS);  }
inline ExprPtr tan (const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::TAN);  }
inline ExprPtr exp (const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::EXP);  }
inline ExprPtr log (const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::LOG);  }
inline ExprPtr sqrt(const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::SQRT); }
inline ExprPtr abs (const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::ABS);  }
inline ExprPtr asin(const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::ASIN); }
inline ExprPtr acos(const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::ACOS); }
inline ExprPtr atan(const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::ATAN); }
inline ExprPtr sinh(const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::SINH); }
inline ExprPtr cosh(const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::COSH); }
inline ExprPtr tanh(const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::TANH); }

inline ExprPtr pow(const ExprPtr& b, double e)       { return std::make_shared<PowExpr>(b, e);        }
inline ExprPtr pow(const ExprPtr& b, const ExprPtr& e) { return std::make_shared<GeneralPowExpr>(b, e); }

// ── ParameterBlock ────────────────────────────────────────────────────────────
// Owns all named refineable parameters and their current values.
// Pass params.values() into eval / eval_d calls.
class ParameterBlock {
public:
    // Register a new named parameter; returns its expression node.
    // Throws std::invalid_argument if the name is already registered.
    ExprPtr add(const std::string& name, double initial_value) {
        if (name_to_idx_.count(name))
            throw std::invalid_argument("Parameter already registered: " + name);
        int idx = static_cast<int>(values_.size());
        name_to_idx_[name] = idx;
        names_.push_back(name);
        values_.conservativeResize(idx + 1);
        values_[idx] = initial_value;
        return std::make_shared<ParamRef>(idx);
    }

    // Retrieve the expression node for an already-registered parameter.
    ExprPtr operator[](const std::string& name) const {
        return std::make_shared<ParamRef>(index_of(name));
    }

    int index_of(const std::string& name) const {
        auto it = name_to_idx_.find(name);
        if (it == name_to_idx_.end())
            throw std::out_of_range("Unknown parameter: " + name);
        return it->second;
    }

    void set(const std::string& name, double v) { values_[index_of(name)] = v; }

    void set_values(const Eigen::VectorXd& v) {
        if (v.size() != values_.size())
            throw std::invalid_argument("Parameter vector size mismatch");
        values_ = v;
    }

    // Current parameter values — pass to eval / eval_d.
    const Eigen::VectorXd& values() const { return values_; }

    int size() const { return static_cast<int>(values_.size()); }

    const std::vector<std::string>& names() const { return names_; }

private:
    std::vector<std::string>             names_;
    std::unordered_map<std::string, int> name_to_idx_;
    Eigen::VectorXd                      values_;
};

}  // namespace yell
