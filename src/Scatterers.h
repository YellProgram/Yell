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
#include <unordered_map>

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

    /// Legacy: sets current_form_factor on this shared object.
    /// Used only by the FFT path (single-threaded).  Do not call from
    /// multi-model code; use ScattererList::update() instead.
    void update_current_form_factor(vec3<double> s, double d_star_sq) {
        current_form_factor = form_factor_at_c(s, d_star_sq);
    }

    /// Legacy: pre-computes gridded form factors onto this shared object.
    /// Used only when the old AtomicPair-based FFT path is active.
    /// Prefer ScattererList::compute_form_factors_on_grid() for new code.
    void calculate_form_factors_on_grid(vec3<int> map_size, Grid grid) {
        gridded_form_factors = IntensityMap(map_size);
        gridded_form_factors.set_grid(grid);
        gridded_form_factors.init_iterator();
        while (gridded_form_factors.next()) {
            gridded_form_factors.current_array_value_c() = form_factor_at_c(
                gridded_form_factors.current_s(),
                gridded_form_factors.current_d_star_square());
        }
    }

    /// Legacy: mutable cache fields on the shared Scatterer object.
    /// New code must use ScattererList instead.
    IntensityMap gridded_form_factors;
    complex<double> current_form_factor;
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

// ─────────────────────────────────────────────────────────────────────────────
// ScattererList — flat, index-addressable list of scatterers.
//
// Built once from the global AtomicTypeCollection after all atoms have been
// registered.  Stores pre-computed per-pixel form factors in a contiguous
// vector so the inner pair loop can use integer indices instead of pointer
// dereferences.
// ─────────────────────────────────────────────────────────────────────────────

class ScattererList {
public:
    /// Snapshot the current AtomicTypeCollection into a flat vector.
    ScattererList() {
        for (auto& kv : AtomicTypeCollection::get().types)
            scatterers_.push_back(kv.second);
        form_factors_.resize(scatterers_.size());
    }

    int size() const { return (int)scatterers_.size(); }

    /// Return the index of scatterer s, or -1 if not found.
    int index_of(Scatterer* s) const {
        for (int i = 0; i < (int)scatterers_.size(); ++i)
            if (scatterers_[i] == s) return i;
        return -1;
    }

    /// Recompute per-pixel form factors for the given reciprocal-space point.
    /// Stores results in form_factors_[i] — no writes to any shared Scatterer state.
    void update(vec3<double> s_vec, double d_star_sq) {
        for (int i = 0; i < (int)scatterers_.size(); ++i)
            form_factors_[i] = scatterers_[i]->form_factor_at_c(s_vec, d_star_sq);
    }

    /// Per-pixel form factor at index i (valid after update()).
    complex<double> f(int i) const { return form_factors_[i]; }

    /// Register a per-clone override: when computing form factors for scatterer
    /// `orig`, use `replacement->form_factor_at_c()` instead.  The original pointer
    /// is still used for index_of() so PattersonPeak type indices remain valid.
    void add_override(Scatterer* orig, Scatterer* replacement) {
        overrides_[orig] = replacement;
    }

    /// Pre-compute gridded form factors for the FFT path.
    /// Replaces AtomicTypeCollection::calculate_form_factors_of_all_atoms_on_grid():
    /// results are stored per-model in gridded_form_factors_[i] instead of on the
    /// shared Scatterer objects, making the FFT path safe for multi-model use.
    /// Per-clone overrides registered via add_override() are used instead of the
    /// global scatterer — this ensures MolecularScatterer form factors are computed
    /// from per-clone private atoms, eliminating the data race.
    void compute_form_factors_on_grid(vec3<int> map_size, Grid grid) {
        gridded_form_factors_.resize(scatterers_.size());
        for (int i = 0; i < (int)scatterers_.size(); ++i) {
            Scatterer* s = scatterers_[i];
            auto it = overrides_.find(s);
            if (it != overrides_.end()) s = it->second;

            IntensityMap ff_map(map_size);
            ff_map.set_grid(grid);
            ff_map.init_iterator();
            while (ff_map.next())
                ff_map.current_array_value_c() = s->form_factor_at_c(
                    ff_map.current_s(), ff_map.current_d_star_square());
            gridded_form_factors_[i] = std::move(ff_map);
        }
    }

    /// Gridded form factor for scatterer idx at grid index (valid after
    /// compute_form_factors_on_grid()).
    complex<double> f_gridded(int idx, af::c_grid<3,int>::index_type index) {
        return gridded_form_factors_[idx].at_c(index);
    }

    // Public to allow Model::clone() to remap MolecularScatterer constituent atoms.
    vector<Scatterer*>        scatterers_;

private:
    vector<complex<double>>   form_factors_;
    vector<IntensityMap>      gridded_form_factors_;
    // Per-clone form-factor overrides: orig ptr → replacement ptr.
    unordered_map<Scatterer*, Scatterer*> overrides_;
};

#endif // YELL_SCATTERERS_H
