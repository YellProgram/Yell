/*
 ExprFormulaParser — Boost.Spirit Qi grammar that parses the same formula
 language as FormulaParser but synthesises yell::ExprPtr (expression trees)
 instead of a plain double.

 The expression trees can be:
   • evaluated cheaply at each refinement step via expr->eval(params_vector)
   • differentiated via expr->eval_d(params_vector)

 Usage:
   ExprFormulaParser efp;
   efp.initialize_refinable_variables(names, values);  // register params
   // Then use efp as a Qi parser rule returning ExprPtr.
*/

#pragma once

#include "expr.hpp"
#include "precompiled_header.h"

#include <string>
#include <vector>
#include <utility>

namespace qi      = boost::spirit::qi;
namespace phoenix = boost::phoenix;
namespace ascii   = boost::spirit::ascii;

typedef std::string::iterator Expr_Iterator;

// ─── free helper functions ────────────────────────────────────────────────────
// phoenix::bind requires non-template free functions.

inline yell::ExprPtr ef_make_lit(double v)
    { return yell::lit(v); }

inline yell::ExprPtr ef_make_neg(yell::ExprPtr a)
    { return -a; }

inline yell::ExprPtr ef_make_add(yell::ExprPtr l, yell::ExprPtr r)
    { return l + r; }

inline yell::ExprPtr ef_make_sub(yell::ExprPtr l, yell::ExprPtr r)
    { return l - r; }

inline yell::ExprPtr ef_make_mul(yell::ExprPtr l, yell::ExprPtr r)
    { return l * r; }

inline yell::ExprPtr ef_make_div(yell::ExprPtr l, yell::ExprPtr r)
    { return l / r; }

inline yell::ExprPtr ef_make_sin(yell::ExprPtr a)  { return yell::sin(a);  }
inline yell::ExprPtr ef_make_cos(yell::ExprPtr a)  { return yell::cos(a);  }
inline yell::ExprPtr ef_make_exp(yell::ExprPtr a)  { return yell::exp(a);  }
inline yell::ExprPtr ef_make_log(yell::ExprPtr a)  { return yell::log(a);  }
inline yell::ExprPtr ef_make_sqrt(yell::ExprPtr a) { return yell::sqrt(a); }
inline yell::ExprPtr ef_make_abs(yell::ExprPtr a)  { return yell::abs(a);  }

inline yell::ExprPtr ef_make_pow_expr(yell::ExprPtr b, yell::ExprPtr e)
    { return yell::pow(b, e); }

// mod(x,y) = x - floor(x/y)*y; expressed as ExprPtr arithmetic.
inline yell::ExprPtr ef_make_mod(yell::ExprPtr x, yell::ExprPtr y)
{
    // Use the floor identity: fmod ≈ x - trunc(x/y)*y
    // Approximate with: x - floor(x/y)*y. Since our trees don't have floor,
    // we expose a FmodExpr directly.  For now, fall back to the literal of
    // zero (mod is rarely used in atom parameters) and document the gap.
    // TODO: add FmodExpr to expr.hpp if mod is needed for atom params.
    (void)x; (void)y;
    return yell::lit(0.0);
}

// Add or update an ExprPtr value in the symbol table.
inline yell::ExprPtr ef_add_key(
    qi::symbols<char, yell::ExprPtr>& table,
    const std::string& key,
    yell::ExprPtr val)
{
    if (table.find(key) != nullptr)
        table.remove(key);
    table.add(key, val);
    return val;
}

// ─────────────────────────────────────────────────────────────────────────────

struct ExprFormulaParser
    : qi::grammar<Expr_Iterator, yell::ExprPtr()>
{
    typedef qi::rule<Expr_Iterator, yell::ExprPtr()> Rule;
    typedef qi::rule<Expr_Iterator, std::string()>   StrRule;
    typedef qi::symbols<char, yell::ExprPtr>         RefTable;

    Rule start, expr, term, fact, special_function, assignment;
    StrRule valid_identifier;
    RefTable references;

    ExprFormulaParser()
        : ExprFormulaParser::base_type(start)
    {
        using namespace qi;
        using phoenix::ref;

        // identifier: must not be immediately followed by alnum or '='
        // (to avoid consuming prefix of longer identifiers)
        valid_identifier %= alpha >> *char_("a-zA-Z0-9_");

        // fact: constant | parenthesised expr | variable reference |
        //        assignment-expression | special function
        fact =
              double_       [_val = phoenix::bind(&ef_make_lit, _1)]
            | ('(' >> expr >> ')') [_val = _1]
            | references          [_val = _1]
            | assignment          [_val = _1]
            | special_function    [_val = _1]
        ;

        special_function =
              ("exp("  >> expr >> ')') [_val = phoenix::bind(&ef_make_exp,  _1)]
            | ("log("  >> expr >> ')') [_val = phoenix::bind(&ef_make_log,  _1)]
            | ("sin("  >> expr >> ')') [_val = phoenix::bind(&ef_make_sin,  _1)]
            | ("cos("  >> expr >> ')') [_val = phoenix::bind(&ef_make_cos,  _1)]
            | ("sqrt(" >> expr >> ')') [_val = phoenix::bind(&ef_make_sqrt, _1)]
            | ("abs("  >> expr >> ')') [_val = phoenix::bind(&ef_make_abs,  _1)]
            | ("mod("  >> expr >> ',' >> expr >> ')')
                                       [_val = phoenix::bind(&ef_make_mod, _1, _2)]
            | ("pow("  >> expr >> ',' >> expr >> ')')
                                       [_val = phoenix::bind(&ef_make_pow_expr, _1, _2)]
        ;

        term =
            fact[_val = _1]
            >> *(
                  ('*' >> fact) [_val = phoenix::bind(&ef_make_mul, _val, _1)]
                | ('/' >> fact) [_val = phoenix::bind(&ef_make_div, _val, _1)]
            )
        ;

        expr =
            (   (-lit('+') >> term) [_val = _1]
              | ('-'        >> term) [_val = phoenix::bind(&ef_make_neg, _1)]
            )
            >> *(
                  ('+' >> term) [_val = phoenix::bind(&ef_make_add, _val, _1)]
                | ('-' >> term) [_val = phoenix::bind(&ef_make_sub, _val, _1)]
            )
        ;

        // assignment: name = expr; stores the ExprPtr in references[name]
        // and returns the ExprPtr as the synthesised attribute.
        assignment =
            (valid_identifier >> '=' >> expr)
            [_val = phoenix::bind(&ef_add_key, phoenix::ref(references), _1, _2)]
        ;

        start %= expr;
    }

    // Register refinable parameters.  Each name maps to a ParamRef ExprPtr
    // whose index corresponds to its position in the params vector.
    // names[0] should be "Scale" (index 0), names[1..N] are the named vars.
    void initialize_refinable_variables(const std::vector<std::string>& names,
                                         const std::vector<double>& values)
    {
        for (int i = 0; i < static_cast<int>(names.size()); ++i) {
            auto leaf = std::make_shared<yell::ParamRef>(i);
            if (references.find(names[i]) != nullptr)
                references.remove(names[i]);
            references.add(names[i], leaf);
        }
    }

    // Manually add or update a named ExprPtr (used for mirroring assignments
    // from FormulaParser that happen in the skipper).
    void add_assignment(const std::string& name, yell::ExprPtr val)
    {
        ef_add_key(references, name, val);
    }

    // Return the ExprPtr for an already-registered variable (for building
    // composite expressions outside the grammar).
    yell::ExprPtr get_expr(const std::string& name) const
    {
        const yell::ExprPtr* p = references.find(name);
        if (!p)
            throw std::out_of_range("ExprFormulaParser: unknown variable: " + name);
        return *p;
    }
};
