#pragma once

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

class Expr;
using ExprPtr = std::shared_ptr<const Expr>;
using Dual = Eigen::AutoDiffScalar<Eigen::VectorXd>;

struct EvaluationCache {
    std::unordered_map<const Expr*, double> values;
    std::unordered_map<const Expr*, Dual>   duals;
    void clear() { values.clear(); duals.clear(); }
};

class Expr {
public:
    virtual ~Expr() = default;
    virtual double eval(const Eigen::VectorXd& p, EvaluationCache* cache = nullptr) const = 0;
    virtual Dual   eval_d(const Eigen::VectorXd& p, EvaluationCache* cache = nullptr) const = 0;
};

class Literal : public Expr {
    double value_;
public:
    explicit Literal(double v) : value_(v) {}
    double eval(const Eigen::VectorXd&, EvaluationCache* = nullptr) const override { return value_; }
    Dual   eval_d(const Eigen::VectorXd& p, EvaluationCache* = nullptr) const override {
        return Dual(value_, Eigen::VectorXd::Zero(p.size()));
    }
};

inline ExprPtr lit(double v) { return std::make_shared<Literal>(v); }

class ParamRef : public Expr {
    int idx_;
public:
    explicit ParamRef(int idx) : idx_(idx) {}
    double eval(const Eigen::VectorXd& p, EvaluationCache* = nullptr) const override { return p[idx_]; }
    Dual   eval_d(const Eigen::VectorXd& p, EvaluationCache* = nullptr) const override {
        return Dual(p[idx_], p.size(), idx_);
    }
};

enum class BinOp { ADD, SUB, MUL, DIV };

class BinaryExpr : public Expr {
    ExprPtr lhs_, rhs_;
    BinOp   op_;
public:
    BinaryExpr(ExprPtr l, ExprPtr r, BinOp op) : lhs_(l), rhs_(r), op_(op) {}

    double eval(const Eigen::VectorXd& p, EvaluationCache* cache = nullptr) const override {
        const Expr* self = this;
        if (cache && cache->values.count(self)) return cache->values.at(self);
        double l = lhs_->eval(p, cache), r = rhs_->eval(p, cache);
        double res = 0;
        switch (op_) {
            case BinOp::ADD: res = l + r; break;
            case BinOp::SUB: res = l - r; break;
            case BinOp::MUL: res = l * r; break;
            case BinOp::DIV: res = l / r; break;
        }
        if (cache) cache->values[self] = res;
        return res;
    }

    Dual eval_d(const Eigen::VectorXd& p, EvaluationCache* cache = nullptr) const override {
        const Expr* self = this;
        if (cache && cache->duals.count(self)) return cache->duals.at(self);
        Dual l = lhs_->eval_d(p, cache), r = rhs_->eval_d(p, cache);
        Dual res(0.0, Eigen::VectorXd::Zero(p.size()));
        switch (op_) {
            case BinOp::ADD: res = l + r; break;
            case BinOp::SUB: res = l - r; break;
            case BinOp::MUL: res = l * r; break;
            case BinOp::DIV: res = l / r; break;
        }
        if (cache) cache->duals[self] = res;
        return res;
    }
};

class NegateExpr : public Expr {
    ExprPtr arg_;
public:
    explicit NegateExpr(ExprPtr a) : arg_(a) {}
    double eval(const Eigen::VectorXd& p, EvaluationCache* cache = nullptr) const override {
        const Expr* self = this;
        if (cache && cache->values.count(self)) return cache->values.at(self);
        double res = -arg_->eval(p, cache);
        if (cache) cache->values[self] = res;
        return res;
    }
    Dual eval_d(const Eigen::VectorXd& p, EvaluationCache* cache = nullptr) const override {
        const Expr* self = this;
        if (cache && cache->duals.count(self)) return cache->duals.at(self);
        Dual res = -arg_->eval_d(p, cache);
        if (cache) cache->duals[self] = res;
        return res;
    }
};

enum class MathFunc {
    FSIN, FCOS, FTAN,
    FEXP, FLOG, FSQRT, FABS,
    FASIN, FACOS, FATAN,
    FSINH, FCOSH, FTANH
};

class FuncExpr : public Expr {
    ExprPtr  arg_;
    MathFunc func_;

    static Dual apply(MathFunc f, const Dual& x) {
        const double           v = x.value();
        const Eigen::VectorXd& g = x.derivatives();
        switch (f) {
            case MathFunc::FSIN:  return Dual(std::sin(v),   std::cos(v) * g);
            case MathFunc::FCOS:  return Dual(std::cos(v),  -std::sin(v) * g);
            case MathFunc::FTAN: { double c = std::cos(v); return Dual(std::tan(v), g / (c * c)); }
            case MathFunc::FEXP: { double ev = std::exp(v); return Dual(ev, ev * g); }
            case MathFunc::FLOG:  return Dual(std::log(v),  g / v);
            case MathFunc::FSQRT: return Dual(std::sqrt(v), g / (2.0 * std::sqrt(v)));
            case MathFunc::FABS:  return Dual(std::abs(v),  (v >= 0.0 ? 1.0 : -1.0) * g);
            case MathFunc::FASIN: return Dual(std::asin(v),  g / std::sqrt(1.0 - v * v));
            case MathFunc::FACOS: return Dual(std::acos(v), -g / std::sqrt(1.0 - v * v));
            case MathFunc::FATAN: return Dual(std::atan(v),  g / (1.0 + v * v));
            case MathFunc::FSINH: return Dual(std::sinh(v),  std::cosh(v) * g);
            case MathFunc::FCOSH: return Dual(std::cosh(v),  std::sinh(v) * g);
            case MathFunc::FTANH: { double t = std::tanh(v); return Dual(t, (1.0 - t * t) * g); }
        }
        return Dual(0.0, Eigen::VectorXd::Zero(g.size()));
    }

public:
    FuncExpr(ExprPtr a, MathFunc f) : arg_(a), func_(f) {}
    double eval(const Eigen::VectorXd& p, EvaluationCache* cache = nullptr) const override {
        const Expr* self = this;
        if (cache && cache->values.count(self)) return cache->values.at(self);
        double val = arg_->eval(p, cache);
        double res = 0;
        switch (func_) {
            case MathFunc::FSIN: res = std::sin(val); break; case MathFunc::FCOS: res = std::cos(val); break;
            case MathFunc::FTAN: res = std::tan(val); break; case MathFunc::FEXP: res = std::exp(val); break;
            case MathFunc::FLOG: res = std::log(val); break; case MathFunc::FSQRT: res = std::sqrt(val); break;
            case MathFunc::FABS: res = std::abs(val); break; case MathFunc::FASIN: res = std::asin(val); break;
            case MathFunc::FACOS: res = std::acos(val); break; case MathFunc::FATAN: res = std::atan(val); break;
            case MathFunc::FSINH: res = std::sinh(val); break; case MathFunc::FCOSH: res = std::cosh(val); break;
            case MathFunc::FTANH: res = std::tanh(val); break;
        }
        if (cache) cache->values[self] = res;
        return res;
    }
    Dual eval_d(const Eigen::VectorXd& p, EvaluationCache* cache = nullptr) const override {
        const Expr* self = this;
        if (cache && cache->duals.count(self)) return cache->duals.at(self);
        Dual res = apply(func_, arg_->eval_d(p, cache));
        if (cache) cache->duals[self] = res;
        return res;
    }
};

class PowExpr : public Expr {
    ExprPtr base_;
    double  exp_;
public:
    PowExpr(ExprPtr b, double e) : base_(b), exp_(e) {}
    double eval(const Eigen::VectorXd& p, EvaluationCache* cache = nullptr) const override {
        const Expr* self = this;
        if (cache && cache->values.count(self)) return cache->values.at(self);
        double res = std::pow(base_->eval(p, cache), exp_);
        if (cache) cache->values[self] = res;
        return res;
    }
    Dual eval_d(const Eigen::VectorXd& p, EvaluationCache* cache = nullptr) const override {
        const Expr* self = this;
        if (cache && cache->duals.count(self)) return cache->duals.at(self);
        Dual b = base_->eval_d(p, cache);
        double v = b.value();
        Dual res = Dual(std::pow(v, exp_), exp_ * std::pow(v, exp_ - 1.0) * b.derivatives());
        if (cache) cache->duals[self] = res;
        return res;
    }
};

class GeneralPowExpr : public Expr {
    ExprPtr base_, exp_;
public:
    GeneralPowExpr(ExprPtr b, ExprPtr e) : base_(b), exp_(e) {}
    double eval(const Eigen::VectorXd& p, EvaluationCache* cache = nullptr) const override {
        const Expr* self = this;
        if (cache && cache->values.count(self)) return cache->values.at(self);
        double res = std::pow(base_->eval(p, cache), exp_->eval(p, cache));
        if (cache) cache->values[self] = res;
        return res;
    }
    Dual eval_d(const Eigen::VectorXd& p, EvaluationCache* cache = nullptr) const override {
        const Expr* self = this;
        if (cache && cache->duals.count(self)) return cache->duals.at(self);
        Dual b = base_->eval_d(p, cache);
        Dual e = exp_->eval_d(p, cache);
        double u = b.value(), v = e.value();
        double uv = std::pow(u, v);
        Eigen::VectorXd deriv = uv * (v / u * b.derivatives() + std::log(u) * e.derivatives());
        Dual res(uv, deriv);
        if (cache) cache->duals[self] = res;
        return res;
    }
};

inline ExprPtr operator+(const ExprPtr& l, const ExprPtr& r) { return std::make_shared<BinaryExpr>(l, r, BinOp::ADD); }
inline ExprPtr operator-(const ExprPtr& l, const ExprPtr& r) { return std::make_shared<BinaryExpr>(l, r, BinOp::SUB); }
inline ExprPtr operator*(const ExprPtr& l, const ExprPtr& r) { return std::make_shared<BinaryExpr>(l, r, BinOp::MUL); }
inline ExprPtr operator/(const ExprPtr& l, const ExprPtr& r) { return std::make_shared<BinaryExpr>(l, r, BinOp::DIV); }
inline ExprPtr operator-(const ExprPtr& a)                   { return std::make_shared<NegateExpr>(a); }

inline ExprPtr operator+(const ExprPtr& l, double r) { return l + lit(r); }
inline ExprPtr operator+(double l, const ExprPtr& r) { return lit(l) + r; }
inline ExprPtr operator-(const ExprPtr& l, double r) { return l - lit(r); }
inline ExprPtr operator-(double l, const ExprPtr& r) { return lit(l) - r; }
inline ExprPtr operator*(const ExprPtr& l, double r) { return l * lit(r); }
inline ExprPtr operator*(double l, const ExprPtr& r) { return lit(l) * r; }
inline ExprPtr operator/(const ExprPtr& l, double r) { return l / lit(r); }
inline ExprPtr operator/(double l, const ExprPtr& r) { return lit(l) / r; }

inline ExprPtr sin (const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::FSIN);  }
inline ExprPtr cos (const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::FCOS);  }
inline ExprPtr tan (const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::FTAN);  }
inline ExprPtr exp (const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::FEXP);  }
inline ExprPtr log (const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::FLOG);  }
inline ExprPtr sqrt(const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::FSQRT); }
inline ExprPtr abs (const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::FABS);  }
inline ExprPtr asin(const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::FASIN); }
inline ExprPtr acos(const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::FACOS); }
inline ExprPtr atan(const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::FATAN); }
inline ExprPtr sinh(const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::FSINH); }
inline ExprPtr cosh(const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::FCOSH); }
inline ExprPtr tanh(const ExprPtr& a) { return std::make_shared<FuncExpr>(a, MathFunc::FTANH); }

inline ExprPtr pow(const ExprPtr& b, double e)       { return std::make_shared<PowExpr>(b, e);        }
inline ExprPtr pow(const ExprPtr& b, const ExprPtr& e) { return std::make_shared<GeneralPowExpr>(b, e); }

class ParameterBlock {
public:
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

    const Eigen::VectorXd& values() const { return values_; }
    int size() const { return static_cast<int>(values_.size()); }
    const std::vector<std::string>& names() const { return names_; }

private:
    std::vector<std::string>             names_;
    std::unordered_map<std::string, int> name_to_idx_;
    Eigen::VectorXd                      values_;
};

}  // namespace yell
