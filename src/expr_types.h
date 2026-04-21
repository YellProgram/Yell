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

// Expression-based vec3 / sym_mat3 types — used by AtomicPair parameters and
// Atom ADP/position expressions.  Extracted from AtomicPairs.h so that
// ChemicalStructure.h can include this without creating a circular dependency.

#pragma once

#include "expr.hpp"
#include "utils.h"

using namespace scitbx;

struct vec3_expr {
    yell::ExprPtr x, y, z;
    vec3_expr(yell::ExprPtr _x = yell::lit(0), yell::ExprPtr _y = yell::lit(0), yell::ExprPtr _z = yell::lit(0))
        : x(_x), y(_y), z(_z) {}
    vec3_expr(vec3<double> v) : x(yell::lit(v[0])), y(yell::lit(v[1])), z(yell::lit(v[2])) {}

    yell::ExprPtr& operator[](int i) {
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }
    const yell::ExprPtr& operator[](int i) const {
        if (i == 0) return x;
        if (i == 1) return y;
        return z;
    }

    void operator+=(const vec3_expr& other) {
        x = x + other.x; y = y + other.y; z = z + other.z;
    }
    void operator+=(const vec3<double>& other) {
        x = x + other[0]; y = y + other[1]; z = z + other[2];
    }
    void operator-=(const vec3_expr& other) {
        x = x - other.x; y = y - other.y; z = z - other.z;
    }
    void operator-=(const vec3<double>& other) {
        x = x - other[0]; y = y - other[1]; z = z - other[2];
    }
    vec3_expr operator-() const {
        return vec3_expr(-x, -y, -z);
    }
};

inline yell::ExprPtr& operator/=(yell::ExprPtr& lhs, double rhs) {
    lhs = lhs / rhs;
    return lhs;
}

struct sym_mat3_expr {
    yell::ExprPtr u11, u22, u33, u12, u13, u23;
    sym_mat3_expr(yell::ExprPtr _11 = yell::lit(0), yell::ExprPtr _22 = yell::lit(0), yell::ExprPtr _33 = yell::lit(0),
                  yell::ExprPtr _12 = yell::lit(0), yell::ExprPtr _13 = yell::lit(0), yell::ExprPtr _23 = yell::lit(0))
        : u11(_11), u22(_22), u33(_33), u12(_12), u13(_13), u23(_23) {}
    sym_mat3_expr(sym_mat3<double> v)
        : u11(yell::lit(v[0])), u22(yell::lit(v[1])), u33(yell::lit(v[2])),
          u12(yell::lit(v[3])), u13(yell::lit(v[4])), u23(yell::lit(v[5])) {}

    yell::ExprPtr& operator[](int i) {
        if (i == 0) return u11; if (i == 1) return u22; if (i == 2) return u33;
        if (i == 3) return u12; if (i == 4) return u13; return u23;
    }
    const yell::ExprPtr& operator[](int i) const {
        if (i == 0) return u11; if (i == 1) return u22; if (i == 2) return u33;
        if (i == 3) return u12; if (i == 4) return u13; return u23;
    }

    void operator+=(const sym_mat3_expr& other) {
        u11 = u11 + other.u11; u22 = u22 + other.u22; u33 = u33 + other.u33;
        u12 = u12 + other.u12; u13 = u13 + other.u13; u23 = u23 + other.u23;
    }
    void operator+=(const sym_mat3<double>& other) {
        u11 = u11 + other[0]; u22 = u22 + other[1]; u33 = u33 + other[2];
        u12 = u12 + other[3]; u13 = u13 + other[4]; u23 = u23 + other[5];
    }
};

inline vec3_expr operator*(const mat3<double>& m, const vec3_expr& v) {
    return vec3_expr(
        m[0]*v.x + m[1]*v.y + m[2]*v.z,
        m[3]*v.x + m[4]*v.y + m[5]*v.z,
        m[6]*v.x + m[7]*v.y + m[8]*v.z
    );
}

inline sym_mat3_expr operator*(const mat3<double>& m, const sym_mat3_expr& U) {
    auto get_U = [&](int i, int j) {
        if (i == j) return U[i];
        if (i == 0 && j == 1) return U[3]; if (i == 1 && j == 0) return U[3];
        if (i == 0 && j == 2) return U[4]; if (i == 2 && j == 0) return U[4];
        if (i == 1 && j == 2) return U[5]; if (i == 2 && j == 1) return U[5];
        return U[0];
    };

    auto calc_res = [&](int r, int c) {
        yell::ExprPtr res = yell::lit(0);
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                res = res + m[r*3 + i] * get_U(i, j) * m[c*3 + j];
        return res;
    };

    return sym_mat3_expr(
        calc_res(0,0), calc_res(1,1), calc_res(2,2),
        calc_res(0,1), calc_res(0,2), calc_res(1,2)
    );
}

inline sym_mat3_expr trusted_mat_to_sym_mat(const sym_mat3_expr& inp) {
    return inp;
}

inline sym_mat3_expr operator+(const sym_mat3_expr& a, const sym_mat3_expr& b) {
    return sym_mat3_expr(a.u11+b.u11, a.u22+b.u22, a.u33+b.u33,
                         a.u12+b.u12, a.u13+b.u13, a.u23+b.u23);
}

inline sym_mat3_expr operator*(const yell::ExprPtr& s, const sym_mat3<double>& m) {
    return sym_mat3_expr(s*m[0], s*m[1], s*m[2], s*m[3], s*m[4], s*m[5]);
}

inline vec3_expr operator*(const yell::ExprPtr& s, const vec3<double>& v) {
    return vec3_expr(s*v[0], s*v[1], s*v[2]);
}

inline bool operator==(const vec3_expr& lhs, const vec3<double>& rhs) {
    Eigen::VectorXd zero_p;
    return almost_equal(vec3<double>(lhs.x->eval(zero_p), lhs.y->eval(zero_p), lhs.z->eval(zero_p)), rhs);
}
inline bool operator==(const vec3<double>& lhs, const vec3_expr& rhs) { return rhs == lhs; }
inline bool operator!=(const vec3_expr& lhs, const vec3<double>& rhs) { return !(lhs == rhs); }
inline bool operator!=(const vec3<double>& lhs, const vec3_expr& rhs) { return !(lhs == rhs); }

inline bool operator==(const sym_mat3_expr& lhs, const sym_mat3<double>& rhs) {
    Eigen::VectorXd zero_p;
    return almost_equal(sym_mat3<double>(lhs.u11->eval(zero_p), lhs.u22->eval(zero_p), lhs.u33->eval(zero_p),
                                         lhs.u12->eval(zero_p), lhs.u13->eval(zero_p), lhs.u23->eval(zero_p)), rhs);
}
inline bool operator==(const sym_mat3<double>& lhs, const sym_mat3_expr& rhs) { return rhs == lhs; }
inline bool operator!=(const sym_mat3_expr& lhs, const sym_mat3<double>& rhs) { return !(lhs == rhs); }
inline bool operator!=(const sym_mat3<double>& lhs, const sym_mat3_expr& rhs) { return !(lhs == rhs); }
