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

#ifndef YELL_ATOMIC_PAIRS_H
#define YELL_ATOMIC_PAIRS_H

#include "ChemicalStructure.h"
#include "IntensityMap.h"
#include "utils.h"
#include "expr.hpp"

#include <Eigen/Core>
#include <vector>
#include <string>

using namespace std;
using namespace scitbx;

// Forward declarations needed for cross-references within this header
class ADPMode;
class SubstitutionalCorrelation;
class AtomicPairPool;

// ─────────────────────────────────────────────────────────────────────────────
// Free function declarations
// ─────────────────────────────────────────────────────────────────────────────

void add_pair_to_appropriate_place(IntensityMap& small_piece, IntensityMap& accumulator,
                                   vec3<int> r, vector<bool> periodic);

ADPMode* translational_mode(ChemicalUnit* unit, int direction, sym_mat3<double>);
ADPMode  z_rot_mode(ChemicalUnit* unit);
ADPMode* rot_mode(ChemicalUnit* unit, vec3<double> axis, vec3<double> point_on_axis,
                  sym_mat3<double> metrical_matrix);
ADPMode  combine_modes(vector<ADPMode> modes);
vector<SubstitutionalCorrelation*> correlators_from_cuns(ChemicalUnitNode* node1,
                                                          ChemicalUnitNode* node2,
                                                          vector<yell::ExprPtr> corr);

// ─────────────────────────────────────────────────────────────────────────────
// PairModifier — abstract base
// ─────────────────────────────────────────────────────────────────────────────

class PairModifier {
public:
    virtual ~PairModifier() {}
    virtual bool generates_pairs() = 0;
    virtual void modify_pairs(AtomicPairPool* const) = 0;
    virtual void update(const Eigen::VectorXd&) {}
    virtual PairModifier* clone() const = 0;
};

// ─────────────────────────────────────────────────────────────────────────────
// Expression-based parameter structures
// ─────────────────────────────────────────────────────────────────────────────

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

struct ParameterizedParams {
    yell::ExprPtr occupancy;
    vec3_expr r;
    sym_mat3_expr U;

    ParameterizedParams() : occupancy(yell::lit(0)) {}
    ParameterizedParams(double occ, vec3<double> _r, sym_mat3<double> _U)
        : occupancy(yell::lit(occ)), r(_r), U(_U) {}
};

// ─────────────────────────────────────────────────────────────────────────────
// AtomicPair
// ─────────────────────────────────────────────────────────────────────────────
class AtomicPair {
public:
    AtomicPair() {}

    AtomicPair(Atom& _atom1, Atom& _atom2)
        : real(_atom1.occupancy * _atom2.occupancy, _atom2.r - _atom1.r, _atom1.U + _atom2.U),
          average(_atom1.occupancy * _atom2.occupancy, _atom2.r - _atom1.r, _atom1.U + _atom2.U),
          atomic_type1(_atom1.atomic_type),
          atomic_type2(_atom2.atomic_type),
          multiplier(_atom1.multiplier * _atom2.multiplier),
          atom1(&_atom1),
          atom2(&_atom2)
    {}

    yell::ExprPtr&   p(bool average_flag = false) { return params(average_flag).occupancy; }
    vec3_expr&       r(bool average_flag = false) { return params(average_flag).r; }
    sym_mat3_expr&   U(bool average_flag = false) { return params(average_flag).U; }

    yell::ExprPtr&   average_p() { return p(true); }
    yell::ExprPtr&   real_p()    { return p(false); }
    vec3_expr&       average_r() { return r(true); }
    vec3_expr&       real_r()    { return r(false); }
    sym_mat3_expr&   average_U() { return U(true); }
    sym_mat3_expr&   real_U()    { return U(false); }

    Scatterer* atomic_type1;
    Scatterer* atomic_type2;
    double multiplier; ///< 1/symmetry_multiplicity
    Atom* atom1;
    Atom* atom2;

    bool pair_is_withing(Grid grid, const Eigen::VectorXd& p_vals) {
        vec3<double> r_val(average_r().x->eval(p_vals), average_r().y->eval(p_vals), average_r().z->eval(p_vals));
        for (int i = 0; i < 3; ++i)
            if (abs(r_val[i]) > abs(grid.lower_limits[i]))
                return false;
        return true;
    }

    string to_string(const Eigen::VectorXd& p_vals) {
        std::ostringstream oss;
        oss << atom1->label << ' ' << atom2->label << ' ' << multiplier
            << ' ' << real_p()->eval(p_vals)
            << ' ' << real_r().x->eval(p_vals) << ' ' << real_r().y->eval(p_vals) << ' ' << real_r().z->eval(p_vals)
            << ' ' << real_U().u11->eval(p_vals) << ' ' << real_U().u22->eval(p_vals) << ' ' << real_U().u33->eval(p_vals)
            << ' ' << real_U().u12->eval(p_vals) << ' ' << real_U().u13->eval(p_vals) << ' ' << real_U().u23->eval(p_vals)
            << ' ' << average_p()->eval(p_vals)
            << ' ' << average_r().x->eval(p_vals) << ' ' << average_r().y->eval(p_vals) << ' ' << average_r().z->eval(p_vals)
            << ' ' << average_U().u11->eval(p_vals) << ' ' << average_U().u22->eval(p_vals) << ' ' << average_U().u33->eval(p_vals)
            << ' ' << average_U().u12->eval(p_vals) << ' ' << average_U().u13->eval(p_vals) << ' ' << average_U().u23->eval(p_vals);
        return oss.str();
    }

    bool operator==(const AtomicPair& inp) const {
        return atom1 == inp.atom1 && atom2 == inp.atom2
            && almost_equal(multiplier, inp.multiplier);
    }
    bool operator!=(const AtomicPair& inp) const { return !(*this == inp); }

private:
    ParameterizedParams& params(bool average_flag) {
        return average_flag ? average : real;
    }
    ParameterizedParams real;
    ParameterizedParams average;
};

// ─────────────────────────────────────────────────────────────────────────────
// AtomicPairPool
// ─────────────────────────────────────────────────────────────────────────────

class AtomicPairPool {
public:
    AtomicPairPool() {}
    AtomicPairPool(const AtomicPairPool& other) {
        pairs = other.pairs;
        for (int i = 0; i < other.modifiers.size(); i++)
            modifiers.push_back(other.modifiers[i].clone());
    }
    AtomicPairPool& operator=(const AtomicPairPool& other) {
        if (this == &other) return *this;
        pairs = other.pairs;
        modifiers.clear();
        for (int i = 0; i < other.modifiers.size(); i++)
            modifiers.push_back(other.modifiers[i].clone());
        return *this;
    }

    AtomicPair& get_pair(Atom* atom1, Atom* atom2) {
        for (vector<AtomicPair>::iterator pair = pairs.begin(); pair != pairs.end(); pair++)
            if (pair->atom1 == atom1 && pair->atom2 == atom2)
                return *pair;
        pairs.push_back(AtomicPair(*atom1, *atom2));
        return pairs.back();
    }

    void add_modifiers(vector<SubstitutionalCorrelation*> _modifiers) {
        for (vector<SubstitutionalCorrelation*>::iterator it = _modifiers.begin(); it != _modifiers.end(); it++)
            modifiers.push_back((PairModifier*)(*it));
    }

    void add_modifier(PairModifier* modifier) { modifiers.push_back(modifier); }

    void invoke_correlators(const Eigen::VectorXd& p) {
        for (int i = 0; i < modifiers.size(); i++)
            modifiers[i].update(p);

        vector<PairModifier*> non_generating;
        for (int i = 0; i < modifiers.size(); i++) {
            if (modifiers[i].generates_pairs())
                modifiers[i].modify_pairs(this);
            else
                non_generating.push_back(&modifiers[i]);
        }
        for (size_t i = 0; i < non_generating.size(); i++)
            non_generating[i]->modify_pairs(this);
    }

    vector<AtomicPair> pairs;
    p_vector<PairModifier> modifiers;
};

// ─────────────────────────────────────────────────────────────────────────────
// AtomicDisplacement and ADPMode
// ─────────────────────────────────────────────────────────────────────────────

class AtomicDisplacement {
public:
    Atom* atom;
    vec3<double> displacement_vector;

    bool operator==(const AtomicDisplacement& inp) const {
        return atom == inp.atom
            && almost_equal(displacement_vector, inp.displacement_vector);
    }
    bool operator!=(const AtomicDisplacement& inp) const { return !(*this == inp); }
};

class ADPMode {
public:
    ADPMode() {}

    void add_atom(Atom* _atom, vec3<double> r) {
        AtomicDisplacement t;
        t.atom = _atom;
        t.displacement_vector = r;
        atomic_displacements.push_back(t);
    }

    vector<AtomicDisplacement> atomic_displacements;

    bool operator==(const ADPMode& inp) const {
        if (inp.atomic_displacements.size() != atomic_displacements.size())
            return false;
        for (int i = 0; i < atomic_displacements.size(); i++)
            if (!(atomic_displacements[i] == inp.atomic_displacements[i]))
                return false;
        return true;
    }

    bool operator!=(const ADPMode& inp) const { return !(*this == inp); }

    ADPMode* clone() const { return new ADPMode(*this); }
};

// ─────────────────────────────────────────────────────────────────────────────
// PairModifier subclasses
// ─────────────────────────────────────────────────────────────────────────────

class CellShifter : public PairModifier {
public:
    CellShifter(double r1, double r2, double r3) : shift(r1, r2, r3) {}

    bool generates_pairs() { return false; }

    void modify_pairs(AtomicPairPool* const pool) {
        for (vector<AtomicPair>::iterator pair = pool->pairs.begin(); pair != pool->pairs.end(); pair++) {
            pair->r()        += shift;
            pair->average_r() += shift;
        }
    }

    bool operator==(const CellShifter& inp) const { return almost_equal(shift, inp.shift); }
    bool operator!=(const CellShifter& inp) const { return !(*this == inp); }

    PairModifier* clone() const override { return new CellShifter(*this); }

    vec3<double> shift;
};

class SubstitutionalCorrelation : public PairModifier {
public:
    SubstitutionalCorrelation(ChemicalUnit* unit1, ChemicalUnit* unit2, yell::ExprPtr expr)
        : joint_probability_expr(expr)
    {
        chemical_units[0] = unit1;
        chemical_units[1] = unit2;
    }

    bool generates_pairs() { return true; }

    void modify_pairs(AtomicPairPool* const pool) {
        vector<Atom*> atoms1 = chemical_units[0]->get_atoms();
        vector<Atom*> atoms2 = chemical_units[1]->get_atoms();
        for (vector<Atom*>::iterator atom1 = atoms1.begin(); atom1 != atoms1.end(); atom1++)
            for (vector<Atom*>::iterator atom2 = atoms2.begin(); atom2 != atoms2.end(); atom2++) {
                AtomicPair& pair = pool->get_pair(*atom1, *atom2);
                pair.p() = joint_probability_expr;
            }
    }

    bool operator==(const SubstitutionalCorrelation& inp) const {
        return chemical_units[0] == inp.chemical_units[0]
            && chemical_units[1] == inp.chemical_units[1];
    }
    bool operator!=(const SubstitutionalCorrelation& inp) const { return !(*this == inp); }

    PairModifier* clone() const override { return new SubstitutionalCorrelation(*this); }

    ChemicalUnit* chemical_units[2];
    yell::ExprPtr joint_probability_expr;
};

class DoubleADPMode : public PairModifier {
public:
    DoubleADPMode(ADPMode* mode1, ADPMode* mode2, double _amplitude)
        : amplitude(_amplitude)
    {
        modes[0] = mode1;
        modes[1] = mode2;
    }

    bool generates_pairs() { return true; }

    void modify_pairs(AtomicPairPool* const pool) {
        for (vector<AtomicDisplacement>::iterator disp1 = modes[0]->atomic_displacements.begin();
             disp1 != modes[0]->atomic_displacements.end(); disp1++)
            for (vector<AtomicDisplacement>::iterator disp2 = modes[1]->atomic_displacements.begin();
                 disp2 != modes[1]->atomic_displacements.end(); disp2++)
                pool->get_pair(disp1->atom, disp2->atom).U() +=
                    -amplitude * (outer_product(disp1->displacement_vector, disp2->displacement_vector)
                                + outer_product(disp2->displacement_vector, disp1->displacement_vector));
    }

    bool operator==(const DoubleADPMode& inp) const {
        return inp.modes[0] == modes[0] && inp.modes[1] == modes[1]
            && almost_equal(amplitude, inp.amplitude);
    }
    bool operator!=(const DoubleADPMode& inp) const { return !(*this == inp); }

    PairModifier* clone() const override { return new DoubleADPMode(*this); }

private:
    ADPMode* modes[2];
    double amplitude;
};

class StaticShift : public PairModifier {
public:
    StaticShift() {}

    bool generates_pairs() { return true; }

    struct ModeAndAmplitude {
        ADPMode* mode;
        double amplitude;
    };

    void add_displacement(const int& StartOrEnd, ADPMode* mode, double amplitude) {
        if (StartOrEnd != 0 && StartOrEnd != 1)
            throw std::string("StartOrEnd flag should be 0 for start or 1 for end");
        ModeAndAmplitude ma = {mode, amplitude};
        modes_and_amplitudes[StartOrEnd].push_back(ma);
    }

    void modify_pairs(AtomicPairPool* const pool) {
        vector<Atom*> atoms[2];
        vector<vec3<double> > displacements[2];

        for (int st_end = 0; st_end < 2; st_end++)
            for (typename vector<ModeAndAmplitude>::iterator ma = modes_and_amplitudes[st_end].begin();
                 ma != modes_and_amplitudes[st_end].end(); ma++)
                for (vector<AtomicDisplacement>::iterator at_disp = ma->mode->atomic_displacements.begin();
                     at_disp != ma->mode->atomic_displacements.end(); at_disp++) {
                    atoms[st_end].push_back(at_disp->atom);
                    displacements[st_end].push_back(at_disp->displacement_vector * ma->amplitude);
                }

        vector<Atom*> unique0 = unique_elements(atoms[0]);
        for (vector<Atom*>::iterator atom1 = unique0.begin(); atom1 != unique0.end(); atom1++) {
            vector<Atom*>::iterator atom2 = atoms[1].begin();
            vector<vec3<double> >::iterator d2 = displacements[1].begin();
            for (; atom2 != atoms[1].end(); atom2++, d2++)
                pool->get_pair(*atom1, *atom2).r() += *d2;
        }

        vector<Atom*> unique1 = unique_elements(atoms[1]);
        vector<Atom*>::iterator atom1 = atoms[0].begin();
        vector<vec3<double> >::iterator d1 = displacements[0].begin();
        for (; atom1 != atoms[0].end(); atom1++, d1++)
            for (vector<Atom*>::iterator atom2 = unique1.begin(); atom2 != unique1.end(); atom2++)
                pool->get_pair(*atom1, *atom2).r() -= *d1;
    }

    PairModifier* clone() const override { return new StaticShift(*this); }

    vector<ModeAndAmplitude> modes_and_amplitudes[2]; ///< [0]=start [1]=end of PDF vector
};

class SizeEffect : public PairModifier {
public:
    bool generates_pairs() { return true; }

    SizeEffect(ADPMode* _mode, ChemicalUnit* _cu, double _amplitude)
        : cu(_cu), mode(_mode), amplitude(_amplitude), cu_to_adp_mode(false) {}
    SizeEffect(ChemicalUnit* _cu, ADPMode* _mode, double _amplitude)
        : cu(_cu), mode(_mode), amplitude(_amplitude), cu_to_adp_mode(true) {}

    void modify_pairs(AtomicPairPool* const pool) {
        vector<Atom*> atoms = cu->get_atoms();
        for (vector<Atom*>::iterator atom = atoms.begin(); atom != atoms.end(); atom++)
            for (vector<AtomicDisplacement>::iterator ad = mode->atomic_displacements.begin();
                 ad != mode->atomic_displacements.end(); ++ad)
                if (cu_to_adp_mode)
                    pool->get_pair(*atom, ad->atom).r() += ad->displacement_vector * amplitude;
                else
                    pool->get_pair(ad->atom, *atom).r() -= ad->displacement_vector * amplitude;
    }

    bool operator==(const SizeEffect& inp) const {
        return cu == inp.cu && mode == inp.mode
            && almost_equal(amplitude, inp.amplitude)
            && cu_to_adp_mode == inp.cu_to_adp_mode;
    }
    bool operator!=(const SizeEffect& inp) const { return !(*this == inp); }

    PairModifier* clone() const override { return new SizeEffect(*this); }

    bool cu_to_adp_mode; ///< true if CU is at start and ADP mode at end of vector
    ChemicalUnit* cu;
    ADPMode* mode;
    double amplitude;
};

class ZeroVectorCorrelation : public PairModifier {
public:
    ZeroVectorCorrelation() {}
    bool generates_pairs() { return false; }

    void modify_pairs(AtomicPairPool* const pool) {
        Eigen::VectorXd zero_p; // only for evaluating literals during init
        for (vector<AtomicPair>::iterator pair = pool->pairs.begin(); pair != pool->pairs.end(); pair++) {
            vec3<double> r_val(pair->r().x->eval(zero_p), pair->r().y->eval(zero_p), pair->r().z->eval(zero_p));
            if (r_val.length() < 0.0001 && pair->atom1 == pair->atom2) {
                pair->U() = sym_mat3<double>(0,0,0,0,0,0);
                pair->p() = yell::lit(pair->atom1->occupancy);
            }
            if (r_val.length() < 0.0001 && pair->atom1 != pair->atom2)
                pair->p() = yell::lit(0);
        }
    }

    PairModifier* clone() const override { return new ZeroVectorCorrelation(*this); }
};

class MultiplicityCorrelation : public PairModifier {
public:
    MultiplicityCorrelation(double _multiplier) : multiplier(_multiplier) {}

    bool generates_pairs() { return false; }

    void modify_pairs(AtomicPairPool* const pool) {
        for (vector<AtomicPair>::iterator pair = pool->pairs.begin(); pair != pool->pairs.end(); pair++) {
            pair->p(false) = pair->p(false) * multiplier;
            pair->p(true)  = pair->p(true)  * multiplier;
            pair->multiplier *= multiplier;
        }
    }

    bool operator==(const MultiplicityCorrelation& inp) const {
        return almost_equal(multiplier, inp.multiplier);
    }
    bool operator!=(const MultiplicityCorrelation& inp) const { return !(*this == inp); }

    PairModifier* clone() const override { return new MultiplicityCorrelation(*this); }

    double multiplier;
};

// ─────────────────────────────────────────────────────────────────────────────
// PattersonPeak PODs and Baking functions
// ─────────────────────────────────────────────────────────────────────────────

struct PattersonPeak {
    int              type1_idx;
    int              type2_idx;
    double           coefficient;
    vec3<double>     r;
    sym_mat3<double> U;
};

struct PeakSusceptibility {
    double           d_coefficient;
    vec3<double>     d_r;
    sym_mat3<double> d_U;
};

inline void peaks_from_pairs(
    vector<AtomicPair>&        pairs,
    const Eigen::VectorXd&     params,
    const ScattererList&       scatterers,
    vector<PattersonPeak>&     full_peaks,
    vector<PattersonPeak>&     avg_peaks)
{
    full_peaks.clear();
    avg_peaks.clear();
    full_peaks.reserve(pairs.size());
    avg_peaks.reserve(pairs.size());

    yell::EvaluationCache cache;

    for (AtomicPair& pair : pairs) {
        int idx1 = scatterers.index_of(pair.atomic_type1);
        int idx2 = scatterers.index_of(pair.atomic_type2);

        PattersonPeak pk;
        pk.type1_idx  = idx1;
        pk.type2_idx  = idx2;

        pk.coefficient = pair.p(false)->eval(params, &cache) * pair.multiplier;
        pk.r = vec3<double>(pair.r(false).x->eval(params, &cache), pair.r(false).y->eval(params, &cache), pair.r(false).z->eval(params, &cache));
        pk.U = sym_mat3<double>(pair.U(false).u11->eval(params, &cache), pair.U(false).u22->eval(params, &cache), pair.U(false).u33->eval(params, &cache),
                                pair.U(false).u12->eval(params, &cache), pair.U(false).u13->eval(params, &cache), pair.U(false).u23->eval(params, &cache));
        full_peaks.push_back(pk);

        pk.coefficient = pair.p(true)->eval(params, &cache) * pair.multiplier;
        pk.r = vec3<double>(pair.r(true).x->eval(params, &cache), pair.r(true).y->eval(params, &cache), pair.r(true).z->eval(params, &cache));
        pk.U = sym_mat3<double>(pair.U(true).u11->eval(params, &cache), pair.U(true).u22->eval(params, &cache), pair.U(true).u33->eval(params, &cache),
                                pair.U(true).u12->eval(params, &cache), pair.U(true).u13->eval(params, &cache), pair.U(true).u23->eval(params, &cache));
        avg_peaks.push_back(pk);
    }
}

inline void susceptibilities_from_pairs(
    vector<AtomicPair>&        pairs,
    const Eigen::VectorXd&     params,
    int                        param_idx,
    vector<PeakSusceptibility>& full_susc,
    vector<PeakSusceptibility>& avg_susc)
{
    full_susc.clear();
    avg_susc.clear();
    full_susc.reserve(pairs.size());
    avg_susc.reserve(pairs.size());

    yell::EvaluationCache cache;

    for (AtomicPair& pair : pairs) {
        auto bake_susc = [&](bool avg) {
            PeakSusceptibility s;
            yell::Dual p_d   = pair.p(avg)->eval_d(params, &cache);
            yell::Dual rx_d  = pair.r(avg).x->eval_d(params, &cache);
            yell::Dual ry_d  = pair.r(avg).y->eval_d(params, &cache);
            yell::Dual rz_d  = pair.r(avg).z->eval_d(params, &cache);
            yell::Dual u11_d = pair.U(avg).u11->eval_d(params, &cache);
            yell::Dual u22_d = pair.U(avg).u22->eval_d(params, &cache);
            yell::Dual u33_d = pair.U(avg).u33->eval_d(params, &cache);
            yell::Dual u12_d = pair.U(avg).u12->eval_d(params, &cache);
            yell::Dual u13_d = pair.U(avg).u13->eval_d(params, &cache);
            yell::Dual u23_d = pair.U(avg).u23->eval_d(params, &cache);

            s.d_coefficient = p_d.derivatives()[param_idx] * pair.multiplier;
            s.d_r = vec3<double>(rx_d.derivatives()[param_idx], ry_d.derivatives()[param_idx], rz_d.derivatives()[param_idx]);
            s.d_U = sym_mat3<double>(u11_d.derivatives()[param_idx], u22_d.derivatives()[param_idx], u33_d.derivatives()[param_idx],
                                     u12_d.derivatives()[param_idx], u13_d.derivatives()[param_idx], u23_d.derivatives()[param_idx]);
            return s;
        };

        full_susc.push_back(bake_susc(false));
        avg_susc.push_back(bake_susc(true));
    }
}

#endif // YELL_ATOMIC_PAIRS_H
