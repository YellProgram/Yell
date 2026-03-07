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

    // Accessors — average_flag selects between real and average parameters
    double&          p(bool average_flag = false) { return params(average_flag).occupancy; }
    vec3<double>&    r(bool average_flag = false) { return params(average_flag).r; }
    sym_mat3<double>& U(bool average_flag = false) { return params(average_flag).U; }

    double&          average_p() { return p(true); }
    double&          real_p()    { return p(false); }
    vec3<double>&    average_r() { return r(true); }
    vec3<double>&    real_r()    { return r(false); }
    sym_mat3<double>& average_U() { return U(true); }
    sym_mat3<double>& real_U()   { return U(false); }

    Scatterer* atomic_type1;
    Scatterer* atomic_type2;
    double multiplier; ///< 1/symmetry_multiplicity
    Atom* atom1;
    Atom* atom2;
    yell::ExprPtr p_real_expr; ///< ExprPtr for real occupancy; null if not set

    bool operator==(const AtomicPair& inp) {
        return average == inp.average && real == inp.real
            && atom1 == inp.atom1 && atom2 == inp.atom2
            && almost_equal(multiplier, inp.multiplier);
    }

    bool pair_is_withing(Grid grid) {
        for (int i = 0; i < 3; ++i)
            if (abs(average_r()[i]) > abs(grid.lower_limits[i]))
                return false;
        return true;
    }

    string to_string() {
        std::ostringstream oss;
        oss << atom1->label << ' ' << atom2->label << ' ' << multiplier
            << ' ' << p() << ' ' << r()[0] << ' ' << r()[1] << ' ' << r()[2];
        for (int j = 0; j < 6; ++j) oss << ' ' << U()[j];
        oss << ' ' << average_p()
            << ' ' << average_r()[0] << ' ' << average_r()[1] << ' ' << average_r()[2];
        for (int j = 0; j < 6; ++j) oss << ' ' << average_U()[j];
        return oss.str();
    }

private:
    AtomicParams& params(bool average_flag) {
        return average_flag ? average : real;
    }
    AtomicParams real;
    AtomicParams average;
};

// ─────────────────────────────────────────────────────────────────────────────
// AtomicPairPool
// ─────────────────────────────────────────────────────────────────────────────

class AtomicPairPool {
public:
    AtomicPairPool() {}

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
        // Update all modifiers with current parameters before generating pairs.
        for (int i = 0; i < modifiers.size(); i++)
            modifiers[i].update(p);

        vector<PairModifier*> non_generating;
        for (int i = 0; i < modifiers.size(); i++) {
            if (modifiers[i].generates_pairs())
                modifiers[i].modify_pairs(this);
            else
                non_generating.push_back(&modifiers[i]);
        }
        for (vector<PairModifier*>::iterator it = non_generating.begin(); it != non_generating.end(); it++)
            (*it)->modify_pairs(this);
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

    bool operator==(const AtomicDisplacement& inp) {
        return atom == inp.atom
            && almost_equal(displacement_vector, inp.displacement_vector);
    }
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

    bool operator==(const ADPMode& inp) {
        if (inp.atomic_displacements.size() != atomic_displacements.size())
            return false;
        for (int i = 0; i < atomic_displacements.size(); i++)
            if (!(atomic_displacements[i] == inp.atomic_displacements[i]))
                return false;
        return true;
    }

    bool operator!=(const ADPMode& inp) { return !operator==(inp); }
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

    bool operator==(const CellShifter& inp) { return almost_equal(shift, inp.shift); }

    vec3<double> shift;
};

class SubstitutionalCorrelation : public PairModifier {
public:
    SubstitutionalCorrelation(ChemicalUnit* unit1, ChemicalUnit* unit2, yell::ExprPtr expr)
        : joint_probability_expr(expr), joint_probability(0.0)
    {
        chemical_units[0] = unit1;
        chemical_units[1] = unit2;
    }

    void update(const Eigen::VectorXd& p) override {
        joint_probability = joint_probability_expr->eval(p);
    }

    bool generates_pairs() { return true; }

    void modify_pairs(AtomicPairPool* const pool) {
        vector<Atom*> atoms1 = chemical_units[0]->get_atoms();
        vector<Atom*> atoms2 = chemical_units[1]->get_atoms();
        for (vector<Atom*>::iterator atom1 = atoms1.begin(); atom1 != atoms1.end(); atom1++)
            for (vector<Atom*>::iterator atom2 = atoms2.begin(); atom2 != atoms2.end(); atom2++) {
                AtomicPair& pair = pool->get_pair(*atom1, *atom2);
                pair.p() = joint_probability;
                pair.p_real_expr = joint_probability_expr;
            }
    }

    bool operator==(const SubstitutionalCorrelation& inp) {
        return almost_equal(joint_probability, inp.joint_probability)
            && chemical_units[0] == inp.chemical_units[0]
            && chemical_units[1] == inp.chemical_units[1];
    }

    ChemicalUnit* chemical_units[2];
    double joint_probability;
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

    bool operator==(const DoubleADPMode& inp) {
        return inp.modes[0] == modes[0] && inp.modes[1] == modes[1]
            && almost_equal(amplitude, inp.amplitude);
    }

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
            throw string("StartOrEnd flag should be 0 for start or 1 for end");
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

    bool operator==(const SizeEffect& inp) {
        return cu == inp.cu && mode == inp.mode
            && almost_equal(amplitude, inp.amplitude)
            && cu_to_adp_mode == inp.cu_to_adp_mode;
    }

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
        for (vector<AtomicPair>::iterator pair = pool->pairs.begin(); pair != pool->pairs.end(); pair++) {
            if (pair->r().length() < 0.0001 && pair->atom1 == pair->atom2) {
                pair->U() = sym_mat3<double>(0,0,0,0,0,0);
                pair->p() = pair->atom1->occupancy;
            }
            if (pair->r().length() < 0.0001 && pair->atom1 != pair->atom2)
                pair->p() = 0;
        }
    }
};

class MultiplicityCorrelation : public PairModifier {
public:
    MultiplicityCorrelation(double _multiplier) : multiplier(_multiplier) {}

    bool generates_pairs() { return false; }

    void modify_pairs(AtomicPairPool* const pool) {
        for (vector<AtomicPair>::iterator pair = pool->pairs.begin(); pair != pool->pairs.end(); pair++)
            pair->multiplier *= multiplier;
    }

    bool operator==(const MultiplicityCorrelation& inp) {
        return almost_equal(multiplier, inp.multiplier);
    }

    double multiplier;
};

// ─────────────────────────────────────────────────────────────────────────────
// PattersonPeak — a single, flat Patterson-space contribution.
//
// Unlike AtomicPair (which carries both real and average parameters),
// a PattersonPeak represents one term of either the full or the average
// intensity sum.  Scatterer types are stored as integer indices into a
// ScattererList, making the inner per-pixel loop cache-friendly.
//
// Fields:
//   type1_idx / type2_idx  — indices into ScattererList::f()
//   coefficient            — p * N  (occupancy × symmetry multiplicity)
//   multiplier             — N alone (for scaling p_expr derivatives)
//   r                      — displacement vector (fractional coords)
//   U                      — combined ADP tensor (fractional)
//   p_expr                 — ExprPtr for p; null if not a refined parameter
// ─────────────────────────────────────────────────────────────────────────────

struct PattersonPeak {
    int              type1_idx;
    int              type2_idx;
    double           coefficient;
    double           multiplier;
    vec3<double>     r;
    sym_mat3<double> U;
    yell::ExprPtr    p_expr; ///< null for average peaks or non-parameterized real peaks
};

/// Convert a list of AtomicPairs into two PattersonPeak lists using the
/// given ScattererList for scatterer→index mapping.
/// full_peaks: uses real (non-average) pair parameters, carries p_expr.
/// avg_peaks:  uses average pair parameters, p_expr is always null.
inline void peaks_from_pairs(
    vector<AtomicPair>&        pairs,
    const ScattererList&       scatterers,
    vector<PattersonPeak>&     full_peaks,
    vector<PattersonPeak>&     avg_peaks)
{
    full_peaks.clear();
    avg_peaks.clear();
    full_peaks.reserve(pairs.size());
    avg_peaks.reserve(pairs.size());

    for (AtomicPair& pair : pairs) {
        int idx1 = scatterers.index_of(pair.atomic_type1);
        int idx2 = scatterers.index_of(pair.atomic_type2);

        PattersonPeak pk;
        pk.type1_idx  = idx1;
        pk.type2_idx  = idx2;
        pk.multiplier = pair.multiplier;

        pk.coefficient = pair.p(false) * pair.multiplier;
        pk.r           = pair.r(false);
        pk.U           = pair.U(false);
        pk.p_expr      = pair.p_real_expr; // nullable
        full_peaks.push_back(pk);

        pk.coefficient = pair.p(true) * pair.multiplier;
        pk.r           = pair.r(true);
        pk.U           = pair.U(true);
        pk.p_expr      = nullptr; // average occupancy is fixed
        avg_peaks.push_back(pk);
    }
}

#endif // YELL_ATOMIC_PAIRS_H
