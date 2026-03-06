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

#ifndef YELL_UTILS_H
#define YELL_UTILS_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <assert.h>

#include <scitbx/vec3.h>
#include <scitbx/sym_mat3.h>

using namespace std;
using namespace scitbx;

const int VectorStart = 0;
const int VectorEnd   = 1;

/// Skips check for matrix being symmetric. Just takes upper-triangular part.
sym_mat3<double> trusted_mat_to_sym_mat(mat3<double> inp);

sym_mat3<double> outer_product(vec3<double> v1, vec3<double> v2);

/**
 * \brief Unsafe owning pointer vector.
 *
 * Holds a vector of pointers, gives access via operator[], and destroys
 * the pointed-to objects on destruction.
 * TODO: replace with std::vector<std::unique_ptr<Obj>>
 */
template<class Obj>
class p_vector {
public:
    p_vector(p_vector const &) { assert(false); } // copy not implemented
    p_vector() {}

    ~p_vector() {
        for (int i = 0; i < pointer_vector.size(); i++)
            delete pointer_vector[i];
    }

    void push_back(Obj* a)             { pointer_vector.push_back(a); }
    void concat(vector<Obj*> a)        { pointer_vector.insert(pointer_vector.end(), a.begin(), a.end()); }
    int  size()                        { return pointer_vector.size(); }
    Obj& operator[](int i)             { return (*pointer_vector[i]); }

private:
    vector<Obj*> pointer_vector;
};

/// Prints a vec3 to stdout (for debugging).
template<class NumType>
void print_vector(vec3<NumType> v)
{
    cout << v[0] << " " << v[1] << " " << v[2] << endl;
}

/// Returns a sorted, deduplicated copy of inp.
template<class T>
vector<T> unique_elements(vector<T> inp)
{
    vector<T> res;
    sort(inp.begin(), inp.end());
    unique_copy(inp.begin(), inp.end(), back_inserter(res));
    return res;
}

#endif // YELL_UTILS_H
