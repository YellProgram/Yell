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
                                   vec3<int> r, vector<bool> periodic,
                                   vec3<int> border_pixels = vec3<int>(0, 0, 0));

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

// vec3_expr, sym_mat3_expr and their operators are defined in expr_types.h,
// which is included transitively via ChemicalStructure.h.

// ─────────────────────────────────────────────────────────────────────────────
// ParameterizedParams
// ─────────────────────────────────────────────────────────────────────────────

struct ParameterizedParams {
    yell::ExprPtr occupancy;
    vec3_expr r;
    sym_mat3_expr U;
    yell::tensor3_expr C;   ///< anharmonic 3rd order (spatial/reciprocal basis); default lit(0)
    yell::tensor4_expr D;   ///< anharmonic 4th order; default lit(0)

    ParameterizedParams() : occupancy(yell::lit(0)) {}
    ParameterizedParams(double occ, vec3<double> _r, sym_mat3<double> _U)
        : occupancy(yell::lit(occ)), r(_r), U(_U) {}
    ParameterizedParams(double occ, vec3<double> _r, sym_mat3_expr _U)
        : occupancy(yell::lit(occ)), r(_r), U(std::move(_U)) {}
    ParameterizedParams(yell::ExprPtr occ, vec3<double> _r, sym_mat3_expr _U)
        : occupancy(occ), r(_r), U(std::move(_U)) {}
    ParameterizedParams(yell::ExprPtr occ, vec3_expr _r, sym_mat3_expr _U)
        : occupancy(occ), r(std::move(_r)), U(std::move(_U)) {}
};

// ─────────────────────────────────────────────────────────────────────────────
// AtomicPair
// ─────────────────────────────────────────────────────────────────────────────
class AtomicPair {
public:
    AtomicPair() {}

    AtomicPair(Atom& _atom1, Atom& _atom2)
        : real(_atom1.full_occupancy() * _atom2.full_occupancy(),
               vec3_expr(_atom2.r[0]-_atom1.r[0],
                         _atom2.r[1]-_atom1.r[1],
                         _atom2.r[2]-_atom1.r[2]),
               sym_mat3_expr(_atom1.U[0]+_atom2.U[0],
                             _atom1.U[1]+_atom2.U[1],
                             _atom1.U[2]+_atom2.U[2],
                             _atom1.U[3]+_atom2.U[3],
                             _atom1.U[4]+_atom2.U[4],
                             _atom1.U[5]+_atom2.U[5])),
          average(_atom1.full_occupancy() * _atom2.full_occupancy(),
               vec3_expr(_atom2.r[0]-_atom1.r[0],
                         _atom2.r[1]-_atom1.r[1],
                         _atom2.r[2]-_atom1.r[2]),
               sym_mat3_expr(_atom1.U[0]+_atom2.U[0],
                             _atom1.U[1]+_atom2.U[1],
                             _atom1.U[2]+_atom2.U[2],
                             _atom1.U[3]+_atom2.U[3],
                             _atom1.U[4]+_atom2.U[4],
                             _atom1.U[5]+_atom2.U[5])),
          atomic_type1(_atom1.atomic_type),
          atomic_type2(_atom2.atomic_type),
          multiplier(1.0),  // atom occupancy lives in atom.occupancy (ExprPtr); this holds LaueSymmetry factor only
          atom1(&_atom1),
          atom2(&_atom2)
    {
        // Anharmonic cumulants combine like r and U: odd order subtracts, even adds
        // (interatomic vector u = u2 − u1). Atomic GC enters both the full and average
        // tracks; correlations (AnharmonicCorrelation) later add to the full track only.
        if (_atom1.anharmonic || _atom2.anharmonic) {
            real.C = _atom2.C - _atom1.C;   // 3rd order: C2 − C1
            real.D = _atom1.D + _atom2.D;   // 4th order: D1 + D2
            average.C = real.C;
            average.D = real.D;
            anharmonic_ = true;
        }
    }

    yell::ExprPtr&     p(bool average_flag = false) { return params(average_flag).occupancy; }
    vec3_expr&         r(bool average_flag = false) { return params(average_flag).r; }
    sym_mat3_expr&     U(bool average_flag = false) { return params(average_flag).U; }
    yell::tensor3_expr& C(bool average_flag = false) { return params(average_flag).C; }
    yell::tensor4_expr& D(bool average_flag = false) { return params(average_flag).D; }

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
    bool anharmonic_ = false; ///< true if either constituent atom (or a correlation) is anharmonic

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
        // joint_probability_expr lives at the variant-probability level.
        // Scale it by each atom's refinable multiplier so that per-atom
        // occupancy parameters (e.g. pCu) flow through to the pair probability
        // and appear in analytical derivatives.
        vector<Atom*> atoms1 = chemical_units[0]->get_atoms();
        vector<Atom*> atoms2 = chemical_units[1]->get_atoms();
        for (vector<Atom*>::iterator atom1 = atoms1.begin(); atom1 != atoms1.end(); atom1++)
            for (vector<Atom*>::iterator atom2 = atoms2.begin(); atom2 != atoms2.end(); atom2++) {
                AtomicPair& pair = pool->get_pair(*atom1, *atom2);
                pair.real_p() = joint_probability_expr
                                * (*atom1)->mult_expr * (*atom2)->mult_expr;
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
    DoubleADPMode(ADPMode* mode1, ADPMode* mode2, yell::ExprPtr _amplitude)
        : amplitude(_amplitude)
    {
        modes[0] = mode1;
        modes[1] = mode2;
    }

    bool generates_pairs() { return true; }

    void modify_pairs(AtomicPairPool* const pool) {
        for (auto disp1 = modes[0]->atomic_displacements.begin();
             disp1 != modes[0]->atomic_displacements.end(); ++disp1)
            for (auto disp2 = modes[1]->atomic_displacements.begin();
                 disp2 != modes[1]->atomic_displacements.end(); ++disp2) {
                sym_mat3<double> outer =
                    -(outer_product(disp1->displacement_vector, disp2->displacement_vector)
                    + outer_product(disp2->displacement_vector, disp1->displacement_vector));
                pool->get_pair(disp1->atom, disp2->atom).U() += amplitude * outer;
            }
    }

    bool operator==(const DoubleADPMode& inp) const {
        return inp.modes[0] == modes[0] && inp.modes[1] == modes[1];
    }
    bool operator!=(const DoubleADPMode& inp) const { return !(*this == inp); }

    PairModifier* clone() const override { return new DoubleADPMode(*this); }

private:
    ADPMode* modes[2];
    yell::ExprPtr amplitude;
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

    SizeEffect(ADPMode* _mode, ChemicalUnit* _cu, yell::ExprPtr _amplitude)
        : cu(_cu), mode(_mode), amplitude(_amplitude), cu_to_adp_mode(false) {}
    SizeEffect(ChemicalUnit* _cu, ADPMode* _mode, yell::ExprPtr _amplitude)
        : cu(_cu), mode(_mode), amplitude(_amplitude), cu_to_adp_mode(true) {}

    void modify_pairs(AtomicPairPool* const pool) {
        vector<Atom*> atoms = cu->get_atoms();
        for (auto atom = atoms.begin(); atom != atoms.end(); ++atom)
            for (auto ad = mode->atomic_displacements.begin();
                 ad != mode->atomic_displacements.end(); ++ad)
                if (cu_to_adp_mode)
                    pool->get_pair(*atom, ad->atom).r() += amplitude * ad->displacement_vector;
                else
                    pool->get_pair(ad->atom, *atom).r() -= amplitude * ad->displacement_vector;
    }

    bool operator==(const SizeEffect& inp) const {
        return cu == inp.cu && mode == inp.mode && cu_to_adp_mode == inp.cu_to_adp_mode;
    }
    bool operator!=(const SizeEffect& inp) const { return !(*this == inp); }

    PairModifier* clone() const override { return new SizeEffect(*this); }

    bool cu_to_adp_mode; ///< true if CU is at start and ADP mode at end of vector
    ChemicalUnit* cu;
    ADPMode* mode;
    yell::ExprPtr amplitude;
};


class MultiplicityCorrelation : public PairModifier {
public:
    MultiplicityCorrelation(double _multiplier) : multiplier(_multiplier) {}

    bool generates_pairs() { return false; }

    void modify_pairs(AtomicPairPool* const pool) {
        for (vector<AtomicPair>::iterator pair = pool->pairs.begin(); pair != pool->pairs.end(); pair++) {
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
// AnharmonicCorrelation — free Gram–Charlier pair coefficients from mode cumulants.
//
// Left modes act on atom 1, right modes on atom 2 (like the bundled ADPCorrelation).
// coeffs is the fully-symmetric rank-`order` tensor over the COMBINED basis
// [left ++ right] in CIF order (combinations-with-replacement), i.e. the mode-amplitude
// cumulants ⟨q_a q_b …⟩. The pair's spatial tensor (κ of the relative displacement
// Δu = u2 − u1) is assembled as Σ coeff · e_a⊗e_b⊗…, with e_a = −d_a(atom1) for a left
// mode and +d_a(atom2) for a right mode. Only κ(Δu) is observable, so the combined
// coefficient set is over-complete — the user ties redundant components symbolically
// (see GRAM_CHARLIER_INDEPENDENT_COMPONENTS.md for auto-detection).
// ─────────────────────────────────────────────────────────────────────────────

class AnharmonicCorrelation : public PairModifier {
public:
    AnharmonicCorrelation(vector<ADPMode*> left, vector<ADPMode*> right,
                          int order, vector<yell::ExprPtr> coeffs)
        : left_(std::move(left)), right_(std::move(right)),
          order_(order), coeffs_(std::move(coeffs))
    {
        int K = (int)(left_.size() + right_.size());
        size_t expected = yell::enumerate_sym(K, order_).size();
        if (coeffs_.size() != expected)
            throw string("AnharmonicCorrelation" + std::to_string(order_) + " expects " +
                         std::to_string(expected) + " coefficients for " + std::to_string(K) +
                         " combined modes, got " + std::to_string(coeffs_.size()));
    }

    bool generates_pairs() override { return false; }

    void modify_pairs(AtomicPairPool* const pool) override {
        vector<ADPMode*> comb = left_;
        comb.insert(comb.end(), right_.begin(), right_.end());
        const int K  = (int)comb.size();
        const int Lc = (int)left_.size();

        vector<vector<int>> tuples = yell::enumerate_sym(K, order_);

        for (Atom* i : side_atoms(left_))
            for (Atom* j : side_atoms(right_)) {
                AtomicPair& pair = pool->get_pair(i, j);
                // mode a's contribution to Δu = u_j − u_i (left enters with −, right with +)
                vector<vec3<double>> e(K);
                for (int a = 0; a < K; ++a)
                    e[a] = (a < Lc) ? -disp(comb[a], i) : disp(comb[a], j);

                if (order_ == 3) assemble3(pair.C(false), tuples, e);
                else             assemble4(pair.D(false), tuples, e);
                pair.anharmonic_ = true;
            }
    }

    PairModifier* clone() const override { return new AnharmonicCorrelation(*this); }

private:
    static vec3<double> disp(ADPMode* m, Atom* at) {
        for (auto& ad : m->atomic_displacements)
            if (ad.atom == at) return ad.displacement_vector;
        return vec3<double>(0, 0, 0);
    }
    static vector<Atom*> side_atoms(const vector<ADPMode*>& modes) {
        vector<Atom*> out;
        for (ADPMode* m : modes)
            for (auto& ad : m->atomic_displacements)
                if (std::find(out.begin(), out.end(), ad.atom) == out.end())
                    out.push_back(ad.atom);
        return out;
    }

    // acc.c[sc] += Σ_tuples coeff · Σ_{perms of tuple} ∏ e over the spatial axes of sc
    void assemble3(yell::tensor3_expr& acc, const vector<vector<int>>& tuples,
                   const vector<vec3<double>>& e) {
        for (size_t ti = 0; ti < tuples.size(); ++ti) {
            double G[yell::GC3_N] = {0};
            vector<int> perm = tuples[ti];
            do {
                for (int sc = 0; sc < yell::GC3_N; ++sc) {
                    const int* I = yell::GC3_IDX[sc];
                    G[sc] += e[perm[0]][I[0]] * e[perm[1]][I[1]] * e[perm[2]][I[2]];
                }
            } while (std::next_permutation(perm.begin(), perm.end()));
            for (int sc = 0; sc < yell::GC3_N; ++sc)
                if (G[sc] != 0.0) acc.c[sc] = acc.c[sc] + coeffs_[ti] * G[sc];
        }
    }
    void assemble4(yell::tensor4_expr& acc, const vector<vector<int>>& tuples,
                   const vector<vec3<double>>& e) {
        for (size_t ti = 0; ti < tuples.size(); ++ti) {
            double G[yell::GC4_N] = {0};
            vector<int> perm = tuples[ti];
            do {
                for (int sc = 0; sc < yell::GC4_N; ++sc) {
                    const int* I = yell::GC4_IDX[sc];
                    G[sc] += e[perm[0]][I[0]] * e[perm[1]][I[1]] * e[perm[2]][I[2]] * e[perm[3]][I[3]];
                }
            } while (std::next_permutation(perm.begin(), perm.end()));
            for (int sc = 0; sc < yell::GC4_N; ++sc)
                if (G[sc] != 0.0) acc.d[sc] = acc.d[sc] + coeffs_[ti] * G[sc];
        }
    }

    vector<ADPMode*>      left_, right_;
    int                   order_;
    vector<yell::ExprPtr> coeffs_;
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
    bool             anharmonic = false;  ///< true → apply the Gram–Charlier factor G(s)
    yell::tensor3    C{};                  ///< baked spatial 3rd-order tensor (valid if anharmonic)
    yell::tensor4    D{};                  ///< baked spatial 4th-order tensor
};

struct PeakSusceptibility {
    double           d_coefficient;
    vec3<double>     d_r;
    sym_mat3<double> d_U;
    yell::tensor3    d_C{};   ///< ∂C/∂p (zero unless the peak is anharmonic)
    yell::tensor4    d_D{};   ///< ∂D/∂p
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

        pk.anharmonic = pair.anharmonic_;

        pk.coefficient = pair.p(false)->eval(params, &cache) * pair.multiplier;
        pk.r = vec3<double>(pair.r(false).x->eval(params, &cache), pair.r(false).y->eval(params, &cache), pair.r(false).z->eval(params, &cache));
        pk.U = sym_mat3<double>(pair.U(false).u11->eval(params, &cache), pair.U(false).u22->eval(params, &cache), pair.U(false).u33->eval(params, &cache),
                                pair.U(false).u12->eval(params, &cache), pair.U(false).u13->eval(params, &cache), pair.U(false).u23->eval(params, &cache));
        if (pair.anharmonic_) { pk.C = pair.C(false).eval(params, &cache); pk.D = pair.D(false).eval(params, &cache); }
        full_peaks.push_back(pk);

        pk.coefficient = pair.p(true)->eval(params, &cache) * pair.multiplier;
        pk.r = vec3<double>(pair.r(true).x->eval(params, &cache), pair.r(true).y->eval(params, &cache), pair.r(true).z->eval(params, &cache));
        pk.U = sym_mat3<double>(pair.U(true).u11->eval(params, &cache), pair.U(true).u22->eval(params, &cache), pair.U(true).u33->eval(params, &cache),
                                pair.U(true).u12->eval(params, &cache), pair.U(true).u13->eval(params, &cache), pair.U(true).u23->eval(params, &cache));
        if (pair.anharmonic_) { pk.C = pair.C(true).eval(params, &cache); pk.D = pair.D(true).eval(params, &cache); }
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
            if (pair.anharmonic_) {       // ∂C/∂p, ∂D/∂p (no multiplier — matches pk.C/pk.D)
                for (int n = 0; n < yell::GC3_N; ++n)
                    s.d_C.c[n] = pair.C(avg).c[n]->eval_d(params, &cache).derivatives()[param_idx];
                for (int n = 0; n < yell::GC4_N; ++n)
                    s.d_D.d[n] = pair.D(avg).d[n]->eval_d(params, &cache).derivatives()[param_idx];
            }
            return s;
        };

        full_susc.push_back(bake_susc(false));
        avg_susc.push_back(bake_susc(true));
    }
}

#endif // YELL_ATOMIC_PAIRS_H
