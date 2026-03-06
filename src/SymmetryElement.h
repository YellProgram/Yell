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

#ifndef YELL_SYMMETRY_ELEMENT_H
#define YELL_SYMMETRY_ELEMENT_H

#include "diffuser_core.h"

class SymmetryElement {
public:
    SymmetryElement(mat3<double> _permutation_matrix, vec3<double> _displacement)
        : permutation_matrix(_permutation_matrix), displacement(_displacement) {}

    SymmetryElement()
        : permutation_matrix(mat3<double>(1,0,0,0,1,0,0,0,1)),
          displacement(vec3<double>(0,0,0)) {}

    bool operator==(SymmetryElement inp)
    {
        return almost_equal(displacement, inp.displacement)
            && almost_equal(permutation_matrix, inp.permutation_matrix);
    }

    vec3<double> displacement;
    mat3<double> permutation_matrix;
};

#endif // YELL_SYMMETRY_ELEMENT_H
