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

// Legacy umbrella header — prefer including specific headers directly.

#ifndef basic_classes_H
#define basic_classes_H

#include <cctbx/sgtbx/rt_mx.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <scitbx/fftpack/complex_to_complex_3d.h>

#include "utils.h"
#include "Scatterers.h"
#include "ChemicalStructure.h"
#include "AtomicPairs.h"
#include "Calculator.h"
#include "SymmetryElement.h"
#include "CeresMinimizer.h"

using namespace cctbx::sgtbx;

extern OutputHandler report;

class Error {
public:
    Error(string inp) : message(inp) {}
    string message;
};

#endif // basic_classes_H
