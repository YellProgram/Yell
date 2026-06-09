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

#ifndef YELL_CHEMICAL_STRUCTURE_H
#define YELL_CHEMICAL_STRUCTURE_H

#include "utils.h"
#include "Scatterers.h"
#include "LaueSymmetry.h"
#include "expr.hpp"
#include "expr_types.h"
#include "anharmonic.h"

#include <cctbx/uctbx.h>
#include <scitbx/array_family/tiny.h>
#include <vector>
#include <string>

using namespace std;
using namespace scitbx;

class Atom;
class AtomicPair;

// ─────────────────────────────────────────────────────────────────────────────
// ChemicalUnit hierarchy
// ─────────────────────────────────────────────────────────────────────────────

class ChemicalUnit {
public:
    virtual ~ChemicalUnit() {}
    virtual yell::ExprPtr  get_occupancy() = 0;
    virtual void           set_occupancy(yell::ExprPtr) = 0;
    // Non-virtual double overload: converts to lit(d) and delegates.
    // Called by the Variant parser and backward-compat code.
    void set_occupancy(double d) { set_occupancy(yell::lit(d)); }
    virtual vector<Atom*>  get_atoms() = 0;
    virtual ChemicalUnit*  create_symmetric(mat3<double>, vec3<double>) = 0;
    virtual ChemicalUnit*  operator[](int) = 0;
    virtual ChemicalUnit*  clone() const = 0;
};

class AtomicAssembly : public ChemicalUnit {
public:
    AtomicAssembly() : occupancy_expr(yell::lit(1.0)) {}

    AtomicAssembly(vector<ChemicalUnit*> units) : occupancy_expr(yell::lit(1.0)) {
        for (int i = 0; i < (int)units.size(); i++)
            chemical_units.push_back(units[i]);
    }

    void add_chemical_unit(ChemicalUnit* unit) { chemical_units.push_back(unit); }

    yell::ExprPtr get_occupancy() override { return occupancy_expr; }

    void set_occupancy(yell::ExprPtr e) override {
        occupancy_expr = e;
        for (int i = 0; i < (int)chemical_units.size(); ++i)
            chemical_units[i].set_occupancy(e);
    }

    vector<Atom*> get_atoms() {
        vector<Atom*> atoms, t;
        for (int i = 0; i < (int)chemical_units.size(); i++) {
            t = chemical_units[i].get_atoms();
            atoms.insert(atoms.end(), t.begin(), t.end());
        }
        return atoms;
    }

    AtomicAssembly* create_symmetric(mat3<double> sym_matrix, vec3<double> translation) {
        AtomicAssembly* result = new AtomicAssembly();
        result->occupancy_expr = occupancy_expr;
        for (int i = 0; i < (int)chemical_units.size(); i++)
            result->add_chemical_unit(chemical_units[i].create_symmetric(sym_matrix, translation));
        return result;
    }

    ChemicalUnit* operator[](int n) { return &chemical_units[n]; }

    ChemicalUnit* clone() const override {
        vector<ChemicalUnit*> units;
        for (int i = 0; i < (int)chemical_units.size(); i++)
            units.push_back(chemical_units[i].clone());
        AtomicAssembly* a = new AtomicAssembly(units);
        a->occupancy_expr = occupancy_expr;
        return a;
    }

    yell::ExprPtr occupancy_expr;
    p_vector<ChemicalUnit> chemical_units;
};

class ChemicalUnitNode {
public:
    ChemicalUnitNode() {}

    void add_chemical_unit(ChemicalUnit* unit) { chemical_units.push_back(unit); }

    p_vector<ChemicalUnit> chemical_units;

    bool complain_if_sum_of_occupancies_is_not_one() {
        Eigen::VectorXd zero_p;
        double sum = 0;
        for (int i = 0; i < (int)chemical_units.size(); ++i)
            sum += chemical_units[i].get_occupancy()->eval(zero_p);

        if (almost_equal(1, sum))
            return true;
        REPORT(ERROR) << "The sum of occupancies should be 1. Here it is " << sum << "\n";
        return false;
    }

    virtual ChemicalUnitNode* clone() const {
        ChemicalUnitNode* res = new ChemicalUnitNode();
        res->chemical_units = chemical_units; // triggers p_vector deep copy
        return res;
    }
    virtual ~ChemicalUnitNode() {}
};

class UnitCell {
public:
    cctbx::uctbx::unit_cell cell;
    LaueSymmetry laue_symmetry;

    UnitCell() {}

    UnitCell(double a, double b, double c, double alpha, double beta, double gamma)
        : cell(scitbx::af::tiny<double,6>(a, b, c, alpha, beta, gamma)) {}

    UnitCell(const cctbx::uctbx::unit_cell& cell) : cell(cell) {}

    void set_laue_symmetry(string sym) { laue_symmetry = LaueSymmetry(sym); }

    void add_node(ChemicalUnitNode* node) { chemical_unit_nodes.push_back(node); }

    bool operator==(const UnitCell& inp) const {
        return almost_equal(cell.parameters(), inp.cell.parameters());
    }
    bool operator!=(const UnitCell& inp) const { return !(*this == inp); }

    string to_string() {
        std::ostringstream res;
        for (int i = 0; i < 6; ++i)
            res << cell.parameters()[i] << " ";
        return res.str();
    }

    p_vector<ChemicalUnitNode> chemical_unit_nodes;
};

// ─────────────────────────────────────────────────────────────────────────────
// Atom
// ─────────────────────────────────────────────────────────────────────────────
/// TODO: rename U to beta
class Atom : public ChemicalUnit {
public:
    string label;
    Scatterer* atomic_type;

    /// Variant component probability — set by set_occupancy() when the atom is
    /// placed inside a Variant node.  Pure probability, independent of per-atom
    /// multiplier.  get_occupancy() returns this so consistency checks and
    /// correlators_from_cuns see the marginal probability, not the full occupancy.
    yell::ExprPtr component_prob;

    /// Refinable per-atom occupancy multiplier — from param_exprs[0] in
    /// construct_atom*.  Represents disorder that is not resolved via the Variant
    /// mechanism (e.g. a refined site occupancy pCu).  Default: lit(1.0).
    yell::ExprPtr mult_expr;

    /// Position (fractional coords) as expression trees.
    /// Set by construct_atom*; transformed by create_symmetric().
    yell::ExprPtr r[3];

    /// ADP tensor in fractional coordinates as expression trees.
    /// Populated by construct_atom*; ADP conversion (Å² → frac) baked at parse time.
    yell::ExprPtr U[6];

    /// Anharmonic Gram–Charlier coefficients (reciprocal basis, a* baked at parse
    /// time like U). Default lit(0). Only meaningful when `anharmonic` is true.
    yell::tensor3_expr C;   ///< 3rd-order (skewness), 10 components
    yell::tensor4_expr D;   ///< 4th-order (kurtosis), 15 components
    bool anharmonic = false;

    // ── Double caches for MolecularScatterer inner loop and ADPMode ─────────
    double           occ_cache;   // = component_prob * mult_expr evaluated
    vec3<double>     r_cache;
    sym_mat3<double> U_cache;
    yell::tensor3    C_cache{};    // valid after update_caches() when anharmonic
    yell::tensor4    D_cache{};

    /// Full occupancy ExprPtr used by AtomicPair and SubstitutionalCorrelation.
    yell::ExprPtr full_occupancy() const { return component_prob * mult_expr; }

    /// Re-evaluate all ExprPtr trees and write to double caches.
    void update_caches(const Eigen::VectorXd& p, yell::EvaluationCache* cache = nullptr) {
        occ_cache = component_prob->eval(p, cache) * mult_expr->eval(p, cache);
        for (int i = 0; i < 3; ++i) r_cache[i] = r[i]->eval(p, cache);
        for (int i = 0; i < 6; ++i) U_cache[i] = U[i]->eval(p, cache);
        if (anharmonic) {          // skip the 25 evals for the harmonic majority
            C_cache = C.eval(p, cache);
            D_cache = D.eval(p, cache);
        }
    }

    Atom() : occ_cache(1.0) {
        component_prob = yell::lit(1.0);
        mult_expr      = yell::lit(1.0);
        for (int i = 0; i < 6; ++i) U[i] = yell::lit(0);
        for (int i = 0; i < 3; ++i) r[i] = yell::lit(0);
        atomic_type = nullptr;
    }

    /// ExprPtr constructor — used by construct_atom* in model.h.
    /// component_prob starts as lit(1.0); set_occupancy() updates it when the
    /// atom is placed in a Variant.  Call update_caches(params) after construction.
    Atom(string const& _label, ScatteringType st,
         yell::ExprPtr mult,
         yell::ExprPtr rx, yell::ExprPtr ry, yell::ExprPtr rz,
         yell::ExprPtr U0, yell::ExprPtr U1, yell::ExprPtr U2,
         yell::ExprPtr U3, yell::ExprPtr U4, yell::ExprPtr U5)
        : label(_label), component_prob(yell::lit(1.0)), mult_expr(mult), occ_cache(0)
    {
        atomic_type = AtomicTypeCollection::get(_label, st);
        r[0] = rx; r[1] = ry; r[2] = rz;
        U[0] = U0; U[1] = U1; U[2] = U2;
        U[3] = U3; U[4] = U4; U[5] = U5;
    }

    /// Literal constructor with isotropic ADP (used in tests).
    /// _multiplier is the per-atom mult (default 1.0); _occupancy is NOT used —
    /// component_prob stays lit(1.0) until set_occupancy() is called.
    Atom(string const& _label, double _multiplier,
         double r1, double r2, double r3,
         double Uiso,
         sym_mat3<double> reciprocal_metric_tensor,
         ScatteringType scattering_type = XRay)
        : label(_label),
          component_prob(yell::lit(1.0)),
          mult_expr(yell::lit(_multiplier)),
          occ_cache(_multiplier),
          r_cache(r1, r2, r3),
          U_cache(Uiso * reciprocal_metric_tensor)
    {
        atomic_type = AtomicTypeCollection::get(_label, scattering_type);
        r[0] = yell::lit(r1); r[1] = yell::lit(r2); r[2] = yell::lit(r3);
        for (int i = 0; i < 6; ++i) U[i] = yell::lit(U_cache[i]);
    }

    /// Literal constructor with anisotropic ADP (used in tests and direct construction).
    Atom(string const& _label, double _multiplier,
         double r1, double r2, double r3,
         double U11, double U22, double U33, double U12, double U13, double U23,
         ScatteringType scattering_type = XRay)
        : label(_label),
          component_prob(yell::lit(1.0)),
          mult_expr(yell::lit(_multiplier)),
          occ_cache(_multiplier),
          r_cache(r1, r2, r3),
          U_cache(U11, U22, U33, U12, U13, U23)
    {
        atomic_type = AtomicTypeCollection::get(_label, scattering_type);
        r[0] = yell::lit(r1); r[1] = yell::lit(r2); r[2] = yell::lit(r3);
        U[0] = yell::lit(U11); U[1] = yell::lit(U22); U[2] = yell::lit(U33);
        U[3] = yell::lit(U12); U[4] = yell::lit(U13); U[5] = yell::lit(U23);
    }

    /// Returns the pure variant component probability (not multiplied by mult_expr).
    /// Used by consistency checks and correlators_from_cuns.
    yell::ExprPtr get_occupancy() override { return component_prob; }

    /// Called by the Variant parser to stamp the component probability onto the atom.
    /// Does NOT update occ_cache — that happens in update_caches() with real params.
    void set_occupancy(yell::ExprPtr e) override {
        component_prob = e;
    }

    vector<Atom*> get_atoms() {
        vector<Atom*> v;
        v.push_back(this);
        return v;
    }

    Atom* create_symmetric(mat3<double> M, vec3<double> t) {
        Atom* result = new Atom(*this);
        // Transform double caches
        result->r_cache = M * r_cache + t;
        result->U_cache = trusted_mat_to_sym_mat(M * U_cache * M.transpose());
        // Transform ExprPtr r
        yell::ExprPtr new_r[3];
        for (int i = 0; i < 3; ++i)
            new_r[i] = yell::lit(M[i*3+0])*r[0] + yell::lit(M[i*3+1])*r[1]
                     + yell::lit(M[i*3+2])*r[2] + yell::lit(t[i]);
        for (int i = 0; i < 3; ++i)
            result->r[i] = new_r[i];
        // Transform ExprPtr U: compute M * U_expr * M^T element-wise
        auto get_u = [&](int k, int l) -> yell::ExprPtr {
            if (k == 0 && l == 0) return U[0];
            if (k == 1 && l == 1) return U[1];
            if (k == 2 && l == 2) return U[2];
            if ((k == 0 && l == 1) || (k == 1 && l == 0)) return U[3];
            if ((k == 0 && l == 2) || (k == 2 && l == 0)) return U[4];
            return U[5]; // (1,2) and (2,1)
        };
        auto calc_ij = [&](int ii, int jj) -> yell::ExprPtr {
            yell::ExprPtr res = yell::lit(0);
            for (int k = 0; k < 3; ++k)
                for (int l = 0; l < 3; ++l)
                    res = res + M[ii*3+k] * get_u(k,l) * M[jj*3+l];
            return res;
        };
        result->U[0] = calc_ij(0,0); result->U[1] = calc_ij(1,1); result->U[2] = calc_ij(2,2);
        result->U[3] = calc_ij(0,1); result->U[4] = calc_ij(0,2); result->U[5] = calc_ij(1,2);
        return result;
    }

    ChemicalUnit* operator[](int) { throw "atom is a leaf node"; }

    ChemicalUnit* clone() const override { return new Atom(*this); }

    bool operator==(const Atom& inp) const {
        return almost_equal(r_cache, inp.r_cache)
            && almost_equal(U_cache, inp.U_cache)
            && atomic_type == inp.atomic_type;
    }
    bool operator!=(const Atom& inp) const { return !(*this == inp); }
};

// ─────────────────────────────────────────────────────────────────────────────
// MolecularScatterer — defined here because it needs Atom
// ─────────────────────────────────────────────────────────────────────────────

class MolecularScatterer : public Scatterer {
public:
    MolecularScatterer(ChemicalUnit* molecule) {
        constituent_atoms = molecule->get_atoms();
    }

    double form_factor_at(double) { return 0; }

    complex<double> form_factor_at_c(vec3<double> s, double d_star_sq);

    vector<Atom*> constituent_atoms;

    inline static complex<double> form_factor_in_a_point(complex<double> f,
                                                          double occ,
                                                          vec3<double> r,
                                                          sym_mat3<double> U,
                                                          vec3<double> s);
};

#endif // YELL_CHEMICAL_STRUCTURE_H
