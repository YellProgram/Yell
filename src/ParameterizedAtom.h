/*
 ParameterizedAtomData — stores an Atom's parameters as yell::ExprPtr trees
 so that the atom values can be re-evaluated at each refinement step without
 re-parsing the input file.

 Ownership: atom_ptr is a NON-OWNING raw pointer into the ChemicalUnit tree
 (owned by UnitCell::chemical_unit_nodes).  The ExprPtr trees in param_exprs
 are reference-counted and outlive the parser.

 Layout of param_exprs:
   isotropic == true  (5 exprs)  : [mult, x, y, z, Uiso_ang]
   isotropic == false (10 exprs) : [mult, x, y, z, U11, U22, U33, U12, U13, U23]
                                    ADPs in Å²; converted to fractional on update.
*/

#pragma once

#include "expr.hpp"
#include "ChemicalStructure.h"  // Atom, UnitCell

#include <vector>
#include <Eigen/Core>

struct ParameterizedAtomData
{
    std::vector<yell::ExprPtr> param_exprs;  // 5 or 10 expressions
    bool                       isotropic;    // true → 5-param Uiso form
    Atom*                      atom_ptr;     // non-owning
    cctbx::uctbx::unit_cell    unit_cell;

    // Re-evaluate all expression trees and write results into *atom_ptr.
    // p must be the current refinement parameter vector (same indexing as
    // when initialize_refinable_variables was called).
    //
    // TODO: RACE CONDITION SUSPECT — PARALLEL JACOBIAN SAFETY
    //
    // Model::clone() copy-constructs parameterized_atoms_, which copies the
    // Atom* pointers by value.  All 16 Model clones therefore share the SAME
    // underlying Atom objects.  When CeresMinimizer spawns 16 threads that each
    // call clone->calculate_derivative() → pad.update(), every thread writes
    // through *atom_ptr to the same heap-allocated Atom struct simultaneously.
    //
    // Affected fields written per call: multiplier, r[0..2], U[0..5], U_expr[0..5].
    //
    // This data race is currently benign only because the tricarboxamide test
    // model does NOT use parameterized atoms.  Any model that combines
    //   • ParameterizedAtom (e.g. "Atom C1 x+dp1 ...")
    //   • MaxProcessors > 1
    //   • Derivatives analytical
    // will produce silent corruption or a crash.
    //
    // Fix: give each Model clone its own Atom copies instead of sharing pointers.
    // The fix should happen in Model::clone(): deep-copy the atom tree
    // (cell.chemical_unit_nodes already deep-copies Atom objects via p_vector,
    // but parameterized_atoms_[i].atom_ptr still points into the ORIGINAL tree).
    // After cloning, remap each atom_ptr to the corresponding Atom in the clone's
    // cell.chemical_unit_nodes by matching atom label or pointer identity.
    void update(const Eigen::VectorXd& p, yell::EvaluationCache* cache = nullptr) const
    {
        atom_ptr->multiplier = param_exprs[0]->eval(p, cache);
        atom_ptr->r[0]       = param_exprs[1]->eval(p, cache);
        atom_ptr->r[1]       = param_exprs[2]->eval(p, cache);
        atom_ptr->r[2]       = param_exprs[3]->eval(p, cache);

        if (isotropic) {
            double Uiso = param_exprs[4]->eval(p, cache);
            auto rm = unit_cell.reciprocal_metrical_matrix();
            atom_ptr->U = Uiso * rm;
            for (int i = 0; i < 6; ++i)
                atom_ptr->U_expr[i] = param_exprs[4] * rm[i];
        } else {
            // ADP conversion: U_frac[i][j] = U_ang[i][j] / (a_i * a_j)
            const double a = unit_cell.parameters()[0];
            const double b = unit_cell.parameters()[1];
            const double c = unit_cell.parameters()[2];
            atom_ptr->U[0] = param_exprs[4]->eval(p, cache) / (a * a);
            atom_ptr->U[1] = param_exprs[5]->eval(p, cache) / (b * b);
            atom_ptr->U[2] = param_exprs[6]->eval(p, cache) / (c * c);
            atom_ptr->U[3] = param_exprs[7]->eval(p, cache) / (a * b);
            atom_ptr->U[4] = param_exprs[8]->eval(p, cache) / (a * c);
            atom_ptr->U[5] = param_exprs[9]->eval(p, cache) / (b * c);
            atom_ptr->U_expr[0] = param_exprs[4] / (a * a);
            atom_ptr->U_expr[1] = param_exprs[5] / (b * b);
            atom_ptr->U_expr[2] = param_exprs[6] / (c * c);
            atom_ptr->U_expr[3] = param_exprs[7] / (a * b);
            atom_ptr->U_expr[4] = param_exprs[8] / (a * c);
            atom_ptr->U_expr[5] = param_exprs[9] / (b * c);
        }
    }
};
