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
    virtual double        get_occupancy() = 0;
    virtual void          set_occupancy(double) = 0;
    virtual vector<Atom*> get_atoms() = 0;
    virtual ChemicalUnit* create_symmetric(mat3<double>, vec3<double>) = 0;
    virtual ChemicalUnit* operator[](int) = 0;
    virtual ChemicalUnit* clone() const = 0;
};

class AtomicAssembly : public ChemicalUnit {
public:
    AtomicAssembly() {}

    AtomicAssembly(vector<ChemicalUnit*> units) {
        for (int i = 0; i < units.size(); i++)
            chemical_units.push_back(units[i]);
    }

    void add_chemical_unit(ChemicalUnit* unit) { chemical_units.push_back(unit); }

    double get_occupancy()           { return occupancy;   }
    void   set_occupancy(double _occ) {
        occupancy = _occ;
        for (int i = 0; i < (int)chemical_units.size(); ++i)
            chemical_units[i].set_occupancy(_occ);
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
        for (int i = 0; i < (int)chemical_units.size(); i++)
            result->add_chemical_unit(chemical_units[i].create_symmetric(sym_matrix, translation));
        return result;
    }

    ChemicalUnit* operator[](int n) { return &chemical_units[n]; }

    ChemicalUnit* clone() const override {
        vector<ChemicalUnit*> units;
        for (int i = 0; i < chemical_units.size(); i++)
            units.push_back(chemical_units[i].clone());
        return new AtomicAssembly(units);
    }

    double occupancy;
    p_vector<ChemicalUnit> chemical_units;
};

class ChemicalUnitNode {
public:
    ChemicalUnitNode() {}

    void add_chemical_unit(ChemicalUnit* unit) { chemical_units.push_back(unit); }

    p_vector<ChemicalUnit> chemical_units;

    bool complain_if_sum_of_occupancies_is_not_one() {
        double sum = 0;
        for (int i = 0; i < chemical_units.size(); ++i)
            sum += chemical_units[i].get_occupancy();

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
// AtomicParams and Atom
// ─────────────────────────────────────────────────────────────────────────────
//TODO: check if this structure is ever used, if not - delete.
class AtomicParams {
public:
    double occupancy;
    vec3<double> r;
    sym_mat3<double> U;

    AtomicParams(double _occupancy, vec3<double> _r, sym_mat3<double> _Uinp)
        : occupancy(_occupancy), r(_r), U(_Uinp) {}
    AtomicParams() {}

    bool operator==(const AtomicParams& inp) const {
        return almost_equal(occupancy, inp.occupancy)
            && almost_equal(r, inp.r)
            && almost_equal(U, inp.U);
    }
};

/// TODO: make sure that precursor constructors are not accessible
/// TODO: rename U to beta
class Atom : public ChemicalUnit {
public:
    string label;
    Scatterer* atomic_type;
    double occupancy;
    double multiplier; ///< multiplier from symmetry; unlike occupancy, unaffected by SubstitutionalCorrelation
    vec3<double> r;
    sym_mat3<double> U;
    /// ExprPtr representations of U[0..5] (fractional coords).
    /// Populated by ParameterizedAtomData::update() so analytical derivatives
    /// flow through pair.U() ExprPtr trees.  Initialised to lit(U[i]) here.
    yell::ExprPtr U_expr[6];
    /// Live expression for the atom's total neutral occupancy = comp_prob_expr * mult_expr.
    /// Initialized to lit(occupancy*multiplier) in constructors; in the parse-once path
    /// construct_atom sets it to param_exprs[0] (mult ExprPtr), then set_occupancy(p)
    /// multiplies by lit(p).  Evaluated by AtomicPair and ZeroVectorCorrelation.
    yell::ExprPtr occupancy_expr;

    Atom() {
        for (int i = 0; i < 6; ++i) U_expr[i] = yell::lit(0);
        occupancy_expr = yell::lit(1.0);
    }

    /// Atom with Uiso and metric tensor.
    Atom(string const& _label, double _multiplier, double _occupancy,
         double r1, double r2, double r3,
         double Uiso,
         sym_mat3<double> reciprocal_metric_tensor,
         ScatteringType scattering_type = XRay)
        : multiplier(_multiplier), occupancy(_occupancy), label(_label),
          r(r1, r2, r3), U(Uiso * reciprocal_metric_tensor)
    {
        atomic_type = AtomicTypeCollection::get(_label, scattering_type);
        for (int i = 0; i < 6; ++i) U_expr[i] = yell::lit(U[i]);
        occupancy_expr = yell::lit(occupancy * multiplier);
    }

    /// Constructor for tests only.
    Atom(string const& _label, double _occupancy,
         double r1, double r2, double r3,
         double U11, double U22, double U33, double U12, double U13, double U23)
        : occupancy(_occupancy), r(r1, r2, r3), U(U11, U22, U33, U12, U13, U23),
          label(_label), multiplier(1)
    {
        atomic_type = AtomicTypeCollection::get(_label, XRay);
        for (int i = 0; i < 6; ++i) U_expr[i] = yell::lit(U[i]);
        occupancy_expr = yell::lit(occupancy * multiplier);
    }

    /// General constructor. Uij parameters are Uij/ai*aj.
    Atom(string const& _label, double _multiplier, double _occupancy,
         double r1, double r2, double r3,
         double U11, double U22, double U33, double U12, double U13, double U23,
         ScatteringType scattering_type = XRay)
        : multiplier(_multiplier), occupancy(_occupancy),
          r(r1, r2, r3), U(U11, U22, U33, U12, U13, U23), label(_label)
    {
        atomic_type = AtomicTypeCollection::get(_label, scattering_type);
        for (int i = 0; i < 6; ++i) U_expr[i] = yell::lit(U[i]);
        occupancy_expr = yell::lit(occupancy * multiplier);
    }

    double get_occupancy()           { return occupancy; }
    void   set_occupancy(double _o)  {
        occupancy = _o;
        // Multiply the existing expression (mult_expr) by the component probability.
        // After construct_atom sets occupancy_expr = param_exprs[0] (mult ExprPtr),
        // this produces comp_prob * mult_expr so derivatives flow through both.
        occupancy_expr = yell::lit(_o) * occupancy_expr;
    }

    vector<Atom*> get_atoms() {
        vector<Atom*> v;
        v.push_back(this);
        return v;
    }

    Atom* create_symmetric(mat3<double> transformation_matrix, vec3<double> translation) {
        Atom* result = new Atom(*this);
        result->r = transformation_matrix * r + translation;
        result->U = trusted_mat_to_sym_mat(transformation_matrix * U * transformation_matrix.transpose());
        return result;
    }

    ChemicalUnit* operator[](int) { throw "atom is a leaf node"; }

    ChemicalUnit* clone() const override { return new Atom(*this); }

    bool operator==(const Atom& inp) const {
        return almost_equal(occupancy, inp.occupancy)
            && almost_equal(r, inp.r)
            && almost_equal(U, inp.U)
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
                                                          double p,
                                                          double N,
                                                          vec3<double> r,
                                                          sym_mat3<double> U,
                                                          vec3<double> s);
};

#endif // YELL_CHEMICAL_STRUCTURE_H
