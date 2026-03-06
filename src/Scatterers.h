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

#ifndef YELL_SCATTERERS_H
#define YELL_SCATTERERS_H

#include "diffuser_core.h"
#include "IntensityMap.h"
#include "Grid.h"
#include "OutputHandler.h"

#include <cctbx/eltbx/neutron.h>
#include <cctbx/eltbx/electron_scattering.h>
#include <cctbx/eltbx/xray_scattering.h>
#include <cctbx/eltbx/xray_scattering/gaussian.h>

#include <complex>
#include <string>
#include <map>

using namespace std;
using namespace scitbx;

enum ScatteringType { XRay, Neutron, Electrons }; // TODO: change to enum class

class Scatterer;

class AtomicTypeCollection {
public:
    AtomicTypeCollection() {}
    ~AtomicTypeCollection();
    static AtomicTypeCollection& get();
    map<string, Scatterer*> types;

    static void calculate_form_factors_of_all_atoms_on_grid(vec3<int> map_size, Grid grid);
    static void update_current_form_factors(vec3<double> s, double d_star_square);

    /** Creates and registers a Scatterer for label/type. Returns existing one if already registered. */
    static Scatterer* get(string const& label, ScatteringType t);
    static string strip_label(string const& label);
    static void add(string const& label, Scatterer* s);
};

class Scatterer {
public:
    virtual complex<double> form_factor_at_c(vec3<double> s, double d_star_sq) = 0;
    virtual double form_factor_at(double d_star_sq) = 0;
    virtual ~Scatterer() {}

    void update_current_form_factor(vec3<double> s, double d_star_sq) {
        current_form_factor = form_factor_at_c(s, d_star_sq);
    }

    void calculate_form_factors_on_grid(vec3<int> map_size, Grid grid) {
        gridded_form_factors = IntensityMap(map_size);
        gridded_form_factors.set_grid(grid);

        gridded_form_factors.init_iterator();
        while (gridded_form_factors.next()) {
            AtomicTypeCollection::update_current_form_factors(
                gridded_form_factors.current_s(),
                gridded_form_factors.current_d_star_square());
            gridded_form_factors.current_array_value_c() = current_form_factor;
        }
    }

    IntensityMap gridded_form_factors;
    complex<double> current_form_factor; ///< cached form factor for direct calculation
};

/// X-ray atomic form factor
class AtomicType : public Scatterer {
public:
    AtomicType() {}

    AtomicType(string const& _label) {
        cctbx::eltbx::xray_scattering::wk1995 wk(_label);
        gauss = wk.fetch();
    }

    double form_factor_at(double d_star_sq) {
        return gauss.at_d_star_sq(d_star_sq);
    }

    complex<double> form_factor_at_c(vec3<double> s, double d_star_sq) {
        return form_factor_at(d_star_sq);
    }

private:
    cctbx::eltbx::xray_scattering::gaussian gauss;
};

// TODO: add tests
// TODO: remove copy-pastes in this code
/// Neutron atomic form factor
class NeutronScattererAtom : public Scatterer {
public:
    NeutronScattererAtom() {}

    NeutronScattererAtom(string const& _label)
        : scat_table_entry(_label) {}

    double form_factor_at(double) {
        return 0; // TODO: check why this is needed
    }

    complex<double> form_factor_at_c(vec3<double> s, double d_star_sq) {
        return scat_table_entry.bound_coh_scatt_length();
    }

private:
    cctbx::eltbx::neutron::neutron_news_1992_table scat_table_entry;
};

/// Electron atomic form factor
class ElectronScattererAtom : public Scatterer {
public:
    ElectronScattererAtom() {}

    ElectronScattererAtom(string const& _label) {
        cctbx::eltbx::electron_scattering::peng1996 wk(_label);
        gauss = wk.fetch();
    }

    double form_factor_at(double d_star_sq) {
        return gauss.at_d_star_sq(d_star_sq);
    }

    complex<double> form_factor_at_c(vec3<double> s, double d_star_sq) {
        return form_factor_at(d_star_sq);
    }

private:
    cctbx::eltbx::xray_scattering::gaussian gauss;
};

#endif // YELL_SCATTERERS_H
