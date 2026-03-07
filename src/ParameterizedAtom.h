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
    void update(const Eigen::VectorXd& p, yell::EvaluationCache* cache = nullptr) const
    {
        atom_ptr->multiplier = param_exprs[0]->eval(p, cache);
        atom_ptr->r[0]       = param_exprs[1]->eval(p, cache);
        atom_ptr->r[1]       = param_exprs[2]->eval(p, cache);
        atom_ptr->r[2]       = param_exprs[3]->eval(p, cache);

        if (isotropic) {
            double Uiso = param_exprs[4]->eval(p, cache);
            atom_ptr->U = Uiso * unit_cell.reciprocal_metrical_matrix();
        } else {
            // ADP conversion: U_frac[i][j] = U_ang[i][j] / (a_i * a_j)
            const double a = unit_cell.parameters()[0];
            const double b = unit_cell.parameters()[1];
            const double c = unit_cell.parameters()[2];
            atom_ptr->U[0] = param_exprs[4]->eval(p, cache) / (a * a);  // U11
            atom_ptr->U[1] = param_exprs[5]->eval(p, cache) / (b * b);  // U22
            atom_ptr->U[2] = param_exprs[6]->eval(p, cache) / (c * c);  // U33
            atom_ptr->U[3] = param_exprs[7]->eval(p, cache) / (a * b);  // U12
            atom_ptr->U[4] = param_exprs[8]->eval(p, cache) / (a * c);  // U13
            atom_ptr->U[5] = param_exprs[9]->eval(p, cache) / (b * c);  // U23
        }
    }
};
