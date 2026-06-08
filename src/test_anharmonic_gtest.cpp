/*
 Tests for anharmonic.h — Gram–Charlier tensor toolkit (Phase 1).

 Strategy: every contraction / assembly is checked against an independent
 brute-force loop over the full 3^n index space, so a wrong multiplicity,
 index table, or prefactor fails loudly.
*/

#include <gtest/gtest.h>

#include "anharmonic.h"

#include <array>
#include <complex>
#include <cmath>

using yell::tensor3;
using yell::tensor4;
using scitbx::vec3;

// ── helpers: expand a unique-component tensor into a full symmetric array ────────

static int find3(int i, int j, int k) {
    int t[3] = {i, j, k};
    std::sort(t, t + 3);
    for (int n = 0; n < yell::GC3_N; ++n)
        if (yell::GC3_IDX[n][0] == t[0] && yell::GC3_IDX[n][1] == t[1] && yell::GC3_IDX[n][2] == t[2])
            return n;
    return -1;
}

static int find4(int i, int j, int k, int l) {
    int t[4] = {i, j, k, l};
    std::sort(t, t + 4);
    for (int n = 0; n < yell::GC4_N; ++n)
        if (yell::GC4_IDX[n][0] == t[0] && yell::GC4_IDX[n][1] == t[1] &&
            yell::GC4_IDX[n][2] == t[2] && yell::GC4_IDX[n][3] == t[3])
            return n;
    return -1;
}

static double brute3(const tensor3& C, const vec3<double>& s) {
    double r = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                r += C.c[find3(i, j, k)] * s[i] * s[j] * s[k];
    return r;
}

static double brute4(const tensor4& D, const vec3<double>& s) {
    double r = 0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    r += D.d[find4(i, j, k, l)] * s[i] * s[j] * s[k] * s[l];
    return r;
}

// distinct component values so a swapped index table is caught
static tensor3 sample3() {
    tensor3 C{};
    for (int n = 0; n < yell::GC3_N; ++n) C.c[n] = 0.1 + 0.07 * n;
    return C;
}
static tensor4 sample4() {
    tensor4 D{};
    for (int n = 0; n < yell::GC4_N; ++n) D.d[n] = -0.05 + 0.03 * n;
    return D;
}

// ── index-table integrity ───────────────────────────────────────────────────────

TEST(Anharmonic, MultiplicitiesSumToFullSpace) {
    int s3 = 0; for (int n = 0; n < yell::GC3_N; ++n) s3 += yell::GC3_MULT[n];
    int s4 = 0; for (int n = 0; n < yell::GC4_N; ++n) s4 += yell::GC4_MULT[n];
    EXPECT_EQ(s3, 27);
    EXPECT_EQ(s4, 81);
}

TEST(Anharmonic, MultiplicityEqualsPermutationCount) {
    // independent count: number of (i,j,k) that sort to this representative
    for (int n = 0; n < yell::GC3_N; ++n) {
        int cnt = 0;
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k)
            if (find3(i, j, k) == n) ++cnt;
        EXPECT_EQ(cnt, yell::GC3_MULT[n]) << "rank3 comp " << n;
    }
    for (int n = 0; n < yell::GC4_N; ++n) {
        int cnt = 0;
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l)
                if (find4(i, j, k, l) == n) ++cnt;
        EXPECT_EQ(cnt, yell::GC4_MULT[n]) << "rank4 comp " << n;
    }
}

// ── contractions vs brute force ──────────────────────────────────────────────────

TEST(Anharmonic, Contract3MatchesBruteForce) {
    tensor3 C = sample3();
    for (vec3<double> s : { vec3<double>(1, 0, 0), vec3<double>(0, 2, 0),
                            vec3<double>(1, 2, 3), vec3<double>(-2, 0.5, 1.7) })
        EXPECT_NEAR(yell::contract3(C, s), brute3(C, s), 1e-12);
}

TEST(Anharmonic, Contract4MatchesBruteForce) {
    tensor4 D = sample4();
    for (vec3<double> s : { vec3<double>(1, 0, 0), vec3<double>(0, 2, 0),
                            vec3<double>(1, 2, 3), vec3<double>(-2, 0.5, 1.7) })
        EXPECT_NEAR(yell::contract4(D, s), brute4(D, s), 1e-12);
}

TEST(Anharmonic, SingleComponentContraction) {
    // C with only c123 (n=4, mult 6) = 1 → contract = 6·s_x s_y s_z
    tensor3 C{}; C.c[4] = 1.0;
    vec3<double> s(1, 2, 3);
    EXPECT_NEAR(yell::contract3(C, s), 6.0 * 1 * 2 * 3, 1e-12);
    // C with only c111 (n=0, mult 1) = 1 → contract = s_x^3
    tensor3 C2{}; C2.c[0] = 1.0;
    EXPECT_NEAR(yell::contract3(C2, s), 1.0, 1e-12);  // s_x = 1
}

// ── Gram–Charlier factor ─────────────────────────────────────────────────────────

TEST(Anharmonic, GramCharlierFactorFormula) {
    tensor3 C = sample3();
    tensor4 D = sample4();
    vec3<double> s(0.3, -0.7, 1.1);
    const double p3 = 4.0 * M_PI * M_PI * M_PI / 3.0;
    const double p4 = 2.0 * M_PI * M_PI * M_PI * M_PI / 3.0;
    std::complex<double> expected(1.0 + p4 * brute4(D, s), -p3 * brute3(C, s));
    std::complex<double> got = yell::gram_charlier_factor(s, C, D);
    EXPECT_NEAR(got.real(), expected.real(), 1e-12);
    EXPECT_NEAR(got.imag(), expected.imag(), 1e-12);
}

TEST(Anharmonic, GramCharlierIsUnitWhenZero) {
    tensor3 C{}; tensor4 D{};
    auto g = yell::gram_charlier_factor(vec3<double>(1, 2, 3), C, D);
    EXPECT_NEAR(g.real(), 1.0, 1e-15);
    EXPECT_NEAR(g.imag(), 0.0, 1e-15);
}

// ── outer-product assembly (mode cumulants → spatial tensor) ─────────────────────

TEST(Anharmonic, OuterProduct3AssemblyContracts) {
    // modes: two arbitrary vectors with a symmetric 3rd-cumulant κ over them.
    // Build C_spatial = Σ_{a,b,c∈{0,1}} κ[a][b][c] d_a⊗d_b⊗d_c, then check
    // contract3(C_spatial, s) == Σ κ[a][b][c] (d_a·s)(d_b·s)(d_c·s).
    vec3<double> d0(1, 0.5, -0.2), d1(0.3, 1.0, 0.7);
    vec3<double> d[2] = { d0, d1 };
    // symmetric 3rd cumulant over the two modes (fill every permutation equal)
    double kappa[2][2][2] = {};
    kappa[0][0][0] = 0.4; kappa[1][1][1] = 0.9;
    int perms001[3][3] = {{0,0,1},{0,1,0},{1,0,0}};
    for (auto& pr : perms001) kappa[pr[0]][pr[1]][pr[2]] = 0.15;
    int perms011[3][3] = {{0,1,1},{1,0,1},{1,1,0}};
    for (auto& pr : perms011) kappa[pr[0]][pr[1]][pr[2]] = -0.25;

    tensor3 C{};
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
            for (int c = 0; c < 2; ++c)
                yell::accumulate_outer3(C, kappa[a][b][c], d[a], d[b], d[c]);

    vec3<double> s(0.6, -0.4, 1.3);
    double expected = 0;
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
            for (int c = 0; c < 2; ++c)
                expected += kappa[a][b][c] * (d[a] * s) * (d[b] * s) * (d[c] * s);

    EXPECT_NEAR(yell::contract3(C, s), expected, 1e-12);
}

TEST(Anharmonic, OuterProduct4AssemblyContracts) {
    vec3<double> d0(0.8, -0.3, 0.4), d1(0.2, 0.9, -0.5);
    vec3<double> d[2] = { d0, d1 };
    // simple diagonal 4th cumulant over the two modes
    double kappa[2] = { 0.5, 1.2 };  // κ_aaaa only
    tensor4 D{};
    for (int a = 0; a < 2; ++a)
        yell::accumulate_outer4(D, kappa[a], d[a], d[a], d[a], d[a]);

    vec3<double> s(1.1, 0.7, -0.9);
    double expected = 0;
    for (int a = 0; a < 2; ++a) {
        double ds = d[a] * s;
        expected += kappa[a] * ds * ds * ds * ds;
    }
    EXPECT_NEAR(yell::contract4(D, s), expected, 1e-12);
}

// ── expr mirrors ─────────────────────────────────────────────────────────────────

TEST(Anharmonic, Tensor3ExprEvalAndArithmetic) {
    Eigen::VectorXd p(2); p << 3.0, 5.0;
    yell::tensor3_expr A;  // all zero by default
    A.c[0] = std::make_shared<yell::ParamRef>(0);   // = p[0]
    A.c[4] = std::make_shared<yell::ParamRef>(1);   // = p[1]

    tensor3 a = A.eval(p);
    EXPECT_NEAR(a.c[0], 3.0, 1e-15);
    EXPECT_NEAR(a.c[4], 5.0, 1e-15);
    EXPECT_NEAR(a.c[1], 0.0, 1e-15);

    yell::tensor3_expr S = A + A;
    EXPECT_NEAR(S.eval(p).c[0], 6.0, 1e-15);
    yell::tensor3_expr Neg = -A;
    EXPECT_NEAR(Neg.eval(p).c[4], -5.0, 1e-15);
    yell::tensor3_expr Diff = A - A;
    EXPECT_NEAR(Diff.eval(p).c[0], 0.0, 1e-15);
}

TEST(Anharmonic, Tensor4ExprDefaultsZero) {
    Eigen::VectorXd p(1); p << 1.0;
    yell::tensor4_expr A;
    tensor4 a = A.eval(p);
    for (int n = 0; n < yell::GC4_N; ++n) EXPECT_NEAR(a.d[n], 0.0, 1e-15);
}
