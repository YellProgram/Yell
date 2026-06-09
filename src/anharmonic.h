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

// ─────────────────────────────────────────────────────────────────────────────
// anharmonic.h — Gram–Charlier tensor toolkit (Phase 1).
//
// Fully-symmetric rank-3 / rank-4 tensors over the 3 reciprocal directions, their
// contraction with the reciprocal vector s=(h,k,l), the Gram–Charlier multiplier
// G(s), and the outer-product assembly used to build a pair's spatial tensor from
// mode amplitudes (see GRAM_CHARLIER_PLAN.md §1.1, §8.1, §8.5).
//
// Convention (ITA §6.1.1.6, with Yell's explicit-2π / s=(h,k,l) basis):
//     T(s) = T_harmonic(s) · G(s)
//     G(s) = 1 − (4π³ i/3) C(s,s,s) + (2π⁴/3) D(s,s,s,s)
//     C(s,s,s)   = Σ_{jkl}  C_{jkl}   s_j s_k s_l        (10 unique comps)
//     D(s,s,s,s) = Σ_{jklm} D_{jklm}  s_j s_k s_l s_m    (15 unique comps)
//
// Component order is CIF lexicographic over the sorted index tuple:
//     rank3:  111 112 113 122 123 133 222 223 233 333
//     rank4:  1111 1112 1113 1122 1123 1133 1222 1223 1233 1333
//             2222 2223 2233 2333 3333
// (indices stored 0-based: 1→0, 2→1, 3→2).
// ─────────────────────────────────────────────────────────────────────────────

#ifndef YELL_ANHARMONIC_H
#define YELL_ANHARMONIC_H

#include <complex>
#include <cmath>
#include <vector>

#include <scitbx/vec3.h>

#include "expr.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace yell {

// ── Canonical sorted index tuples (0-based) and permutation multiplicities ──────
// MULT is the number of distinct index permutations of the tuple, so that
//   Σ_{all 3^n index tuples} A·s…  ==  Σ_{unique n} MULT[n] · A[n] · s…(sorted)
// for a fully-symmetric tensor A. (Σ MULT3 = 27, Σ MULT4 = 81.)

inline constexpr int GC3_N = 10;
inline constexpr int GC4_N = 15;

inline constexpr int GC3_IDX[GC3_N][3] = {
    {0,0,0},{0,0,1},{0,0,2},{0,1,1},{0,1,2},{0,2,2},{1,1,1},{1,1,2},{1,2,2},{2,2,2}
};
inline constexpr int GC3_MULT[GC3_N] = { 1,3,3,3,6,3,1,3,3,1 };

inline constexpr int GC4_IDX[GC4_N][4] = {
    {0,0,0,0},{0,0,0,1},{0,0,0,2},{0,0,1,1},{0,0,1,2},{0,0,2,2},
    {0,1,1,1},{0,1,1,2},{0,1,2,2},{0,2,2,2},
    {1,1,1,1},{1,1,1,2},{1,1,2,2},{1,2,2,2},{2,2,2,2}
};
inline constexpr int GC4_MULT[GC4_N] = { 1,4,4,6,12,6,4,12,12,4,1,4,6,4,1 };

// ── Symmetric multi-index enumeration over K modes ──────────────────────────────
// Non-decreasing n-tuples over {0..K-1} in lexicographic (CIF) order — the layout
// of an AnharmonicCorrelation coefficient list (combinations with replacement,
// count = C(K+n-1, n)). Used to map a coefficient position to its mode tuple.

inline std::vector<std::vector<int>> enumerate_sym(int K, int n) {
    std::vector<std::vector<int>> out;
    std::vector<int> t(n, 0);
    while (true) {
        out.push_back(t);
        int i = n - 1;
        while (i >= 0 && t[i] == K - 1) --i;
        if (i < 0) break;
        int v = t[i] + 1;
        for (int j = i; j < n; ++j) t[j] = v;   // keep non-decreasing
    }
    return out;
}

// ── Value tensors ───────────────────────────────────────────────────────────────
// Aggregates: `tensor3 t{};` zero-initialises all components.

struct tensor3 { double c[GC3_N]; };   ///< symmetric rank-3, spatial (reciprocal) basis
struct tensor4 { double d[GC4_N]; };   ///< symmetric rank-4, spatial (reciprocal) basis

/// C(s,s,s) = Σ_{jkl} C_{jkl} s_j s_k s_l.
inline double contract3(const tensor3& C, const scitbx::vec3<double>& s) {
    double r = 0.0;
    for (int n = 0; n < GC3_N; ++n) {
        const int* I = GC3_IDX[n];
        r += GC3_MULT[n] * C.c[n] * s[I[0]] * s[I[1]] * s[I[2]];
    }
    return r;
}

/// D(s,s,s,s) = Σ_{jklm} D_{jklm} s_j s_k s_l s_m.
inline double contract4(const tensor4& D, const scitbx::vec3<double>& s) {
    double r = 0.0;
    for (int n = 0; n < GC4_N; ++n) {
        const int* I = GC4_IDX[n];
        r += GC4_MULT[n] * D.d[n] * s[I[0]] * s[I[1]] * s[I[2]] * s[I[3]];
    }
    return r;
}

/// Gram–Charlier multiplier G(s) (real from rank-4 kurtosis, imaginary from rank-3 skewness).
inline std::complex<double> gram_charlier_factor(const scitbx::vec3<double>& s,
                                                 const tensor3& C, const tensor4& D) {
    const double p3 = 4.0 * M_PI * M_PI * M_PI / 3.0;        // 4π³/3
    const double p4 = 2.0 * M_PI * M_PI * M_PI * M_PI / 3.0; // 2π⁴/3
    return std::complex<double>(1.0 + p4 * contract4(D, s), -p3 * contract3(C, s));
}

// ── Outer-product assembly (mode cumulants → spatial tensor) ─────────────────────
// AnharmonicCorrelation builds  C_spatial = Σ_{modes a,b,c} κ_{abc} d_a⊗d_b⊗d_c
// by summing one ordered mode-tuple at a time. Because the full mode sum covers
// every ordering, the accumulated tensor is symmetric and each component may be
// stored at its sorted representative. `w` is κ for that mode-tuple.

inline void accumulate_outer3(tensor3& acc, double w,
                              const scitbx::vec3<double>& a,
                              const scitbx::vec3<double>& b,
                              const scitbx::vec3<double>& c) {
    for (int n = 0; n < GC3_N; ++n) {
        const int* I = GC3_IDX[n];
        acc.c[n] += w * a[I[0]] * b[I[1]] * c[I[2]];
    }
}

inline void accumulate_outer4(tensor4& acc, double w,
                              const scitbx::vec3<double>& a,
                              const scitbx::vec3<double>& b,
                              const scitbx::vec3<double>& c,
                              const scitbx::vec3<double>& d) {
    for (int n = 0; n < GC4_N; ++n) {
        const int* I = GC4_IDX[n];
        acc.d[n] += w * a[I[0]] * b[I[1]] * c[I[2]] * d[I[3]];
    }
}

// ── Expression-tree mirrors (built at parse, evaluated at bake) ──────────────────

struct tensor3_expr {
    ExprPtr c[GC3_N];
    tensor3_expr() { for (int n = 0; n < GC3_N; ++n) c[n] = lit(0); }
    tensor3 eval(const Eigen::VectorXd& p, EvaluationCache* cache = nullptr) const {
        tensor3 t{};
        for (int n = 0; n < GC3_N; ++n) t.c[n] = c[n]->eval(p, cache);
        return t;
    }
};

struct tensor4_expr {
    ExprPtr d[GC4_N];
    tensor4_expr() { for (int n = 0; n < GC4_N; ++n) d[n] = lit(0); }
    tensor4 eval(const Eigen::VectorXd& p, EvaluationCache* cache = nullptr) const {
        tensor4 t{};
        for (int n = 0; n < GC4_N; ++n) t.d[n] = d[n]->eval(p, cache);
        return t;
    }
};

inline tensor3_expr operator+(const tensor3_expr& a, const tensor3_expr& b) {
    tensor3_expr r; for (int n = 0; n < GC3_N; ++n) r.c[n] = a.c[n] + b.c[n]; return r;
}
inline tensor3_expr operator-(const tensor3_expr& a, const tensor3_expr& b) {
    tensor3_expr r; for (int n = 0; n < GC3_N; ++n) r.c[n] = a.c[n] - b.c[n]; return r;
}
inline tensor3_expr operator-(const tensor3_expr& a) {
    tensor3_expr r; for (int n = 0; n < GC3_N; ++n) r.c[n] = -a.c[n]; return r;
}

inline tensor4_expr operator+(const tensor4_expr& a, const tensor4_expr& b) {
    tensor4_expr r; for (int n = 0; n < GC4_N; ++n) r.d[n] = a.d[n] + b.d[n]; return r;
}
inline tensor4_expr operator-(const tensor4_expr& a, const tensor4_expr& b) {
    tensor4_expr r; for (int n = 0; n < GC4_N; ++n) r.d[n] = a.d[n] - b.d[n]; return r;
}
inline tensor4_expr operator-(const tensor4_expr& a) {
    tensor4_expr r; for (int n = 0; n < GC4_N; ++n) r.d[n] = -a.d[n]; return r;
}

}  // namespace yell

#endif // YELL_ANHARMONIC_H
