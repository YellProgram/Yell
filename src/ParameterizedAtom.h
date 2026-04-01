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
    // NOTE: atom_ptr race condition — FIXED in Model::clone() (see model.h).
    //
    // Previously, Model::clone() copy-constructed parameterized_atoms_ but left
    // atom_ptr pointing into the original model's atom tree, so all N parallel
    // Jacobian threads wrote through the same Atom objects simultaneously.
    //
    // The fix in Model::clone() builds a map (original Atom* → clone Atom*) from
    // the p_vector deep-copy of cell.chemical_unit_nodes and remaps every atom_ptr
    // before returning the clone.  Each clone thread now writes to its own Atom.
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
            // ADP conversion: U_stored[ij] = U_ang[ij] * a*_i * a*_j
            // Uses scalar reciprocal lattice lengths a* = sqrt(G*[i,i])
            // so the DWF matches the SHELX/CIF convention:
            //   T = exp(-2pi^2 (U11 h^2 a*^2 + U22 k^2 b*^2 + U33 l^2 c*^2
            //                 + 2U12 hk a*b* + 2U13 hl a*c* + 2U23 kl b*c*))
            // For orthogonal cells a*=1/a so this reduces to the old formula.
            auto rm   = unit_cell.reciprocal_metrical_matrix();
            const double astar = std::sqrt(rm[0]);
            const double bstar = std::sqrt(rm[1]);
            const double cstar = std::sqrt(rm[2]);
            atom_ptr->U[0] = param_exprs[4]->eval(p, cache) * (astar * astar);
            atom_ptr->U[1] = param_exprs[5]->eval(p, cache) * (bstar * bstar);
            atom_ptr->U[2] = param_exprs[6]->eval(p, cache) * (cstar * cstar);
            atom_ptr->U[3] = param_exprs[7]->eval(p, cache) * (astar * bstar);
            atom_ptr->U[4] = param_exprs[8]->eval(p, cache) * (astar * cstar);
            atom_ptr->U[5] = param_exprs[9]->eval(p, cache) * (bstar * cstar);
            atom_ptr->U_expr[0] = param_exprs[4] * (astar * astar);
            atom_ptr->U_expr[1] = param_exprs[5] * (bstar * bstar);
            atom_ptr->U_expr[2] = param_exprs[6] * (cstar * cstar);
            atom_ptr->U_expr[3] = param_exprs[7] * (astar * bstar);
            atom_ptr->U_expr[4] = param_exprs[8] * (astar * cstar);
            atom_ptr->U_expr[5] = param_exprs[9] * (bstar * cstar);
        }
    }
};
