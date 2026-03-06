# Yell — Claude Code Guide

## Project Overview

Yell is a C++ diffuse scattering analysis and refinement tool. It parses a `model.txt`
input file (Boost.Spirit Qi grammar), builds a crystal structure model, computes
diffuse scattering intensity maps, and refines structural parameters against
experimental data using Ceres Solver.

- **Source root**: `/Users/asimonov/ag/C++/Yell/Yell`
- **Build dir**: `cmake-build-debug` (CLion default; CMake out-of-source)
- **Active branch**: `parameterized_model` (parse-once + ExprPtr; branched from `replacing_minimizer`)
- **C++ standard**: C++20
- **Current version**: 1.2.9c (in `src/main.cpp`)

---

## Build

```bash
cd /Users/asimonov/ag/C++/Yell/Yell

# Configure (first time or after CMakeLists changes)
cmake -B cmake-build-debug -DCMAKE_BUILD_TYPE=Debug

# Build the main binary
cmake --build cmake-build-debug --target yell -j16

# Build CxxTest suite (test_diffuser.h)
cmake --build cmake-build-debug --target test-yell -j16

# Build Google Test suite (test_expr_gtest.cpp)
cmake --build cmake-build-debug --target test-yell-gtest -j16
```

First build downloads external dependencies (Eigen, glog, gflags, Ceres, googletest)
via ExternalProject/FetchContent — takes ~10 minutes. Subsequent builds are fast.

### Key CMake variables

| Variable | Value |
|---|---|
| `BOOST_INCLUDEDIR` | `../libs/boost_1_88_0` (forced, prevents CLion cache clobber) |
| `CCTBX_PATH` | `lib/cctbx_stubs` |
| `LEVMAR_PATH` | `lib/levmar` |
| Eigen | Downloaded by ExternalProject to `cmake-build-debug/eigen-install` |
| Ceres | Downloaded to `cmake-build-debug/ceres-solver` |

### macOS LAPACK
On Apple: `-framework accelerate` (set automatically in CMakeLists).

---

## Run Tests

```bash
# CxxTest — test_diffuser.h (110 tests, existing diffuse-scattering logic)
cmake-build-debug/test-yell

# Google Test — test_expr_gtest.cpp (expr/ExprFormulaParser/ParameterizedAtom/Model)
cmake-build-debug/test-yell-gtest

# Run a specific gtest filter
cmake-build-debug/test-yell-gtest --gtest_filter="ModelParseOnce.*"
```

---

## Run Yell

Yell reads `model.txt` from the **current working directory** and writes output files
(intensity maps as HDF5, `refined_parameters.txt`, `refinement_trajectory.json`) there.

```bash
# Example: single_pair refinement test
cd "/Users/asimonov/ag/C++/Yell/Yell-jax-stuff/yell-test-example/1. single_pair"
/Users/asimonov/ag/C++/Yell/Yell/cmake-build-debug/yell
```

Expected output for `single_pair`:
- `Scale` refined to ~1.0
- `p_AuAu_100` refined to ~0.25
- Small standard errors on both parameters
- Output files: `cli_out_diffuse.h5`, `refined_parameters.txt`, `refinement_trajectory.json`

---

## Key Source Files

| File | Purpose |
|---|---|
| `src/main.cpp` | Entry point, version string, top-level orchestration |
| `src/model.h` / `model.cpp` | `Model` class: owns parsed structure, runs `calculate()` |
| `src/InputFileParser.h/.cpp` | Boost.Spirit Qi grammar for `model.txt` |
| `src/FormulaParser.h/.cpp` | Formula sub-grammar; `named_assignment` returns `pair<string,double>` |
| `src/ExprFormulaParser.h` | Expression-tree formula parser producing `yell::ExprPtr` |
| `src/expr.hpp` | `ExprPtr = shared_ptr<const Expr>`, `ParameterBlock`, `ParamRef`, `Literal`, autodiff via Eigen |
| `src/ParameterizedAtom.h` | `ParameterizedAtomData`: atom pointer + ExprPtr trees; `update(VectorXd)` |
| `src/basic_classes.h` | Umbrella header; includes `AtomicPairs.h`, `ChemicalStructure.h`, etc. |
| `src/AtomicPairs.h` | `AtomicPairPool`, `SubstitutionalCorrelation`, `CellShifter`, `ADPMode`, etc. |
| `src/ChemicalStructure.h` | `Atom`, `ChemicalUnit`, `ChemicalUnitNode` (Variant), `UnitCell` |
| `src/CeresMinimizer.h/.cpp` | Ceres-based nonlinear least-squares minimizer |
| `src/diffuser_core.cpp` | FFT-based and direct diffuse scattering calculation |
| `src/test_diffuser.h` | CxxTest suite (existing tests) |
| `src/test_expr_gtest.cpp` | Google Test suite (expr, ExprFormulaParser, parse-once model) |

---

## Architecture

### Parse-once model (replacing_minimizer branch)

`Model(string)` constructor calls `parse_model_()`, which runs the Boost.Spirit parser
once — building the full structure (atoms, pools, correlators, ExprPtr trees).

`Model::calculate(params)` is the hot path called by the minimizer at every step:

- Updates atoms via `ParameterizedAtomData::update(p)` — re-evaluates ExprPtr trees.
- Calls `pool->invoke_correlators(p)` which calls `modifier.update(p)` on each
  `PairModifier` before generating pairs. `SubstitutionalCorrelation` stores an
  `ExprPtr` for its joint probability; `update(p)` evaluates it.
- Clears and rebuilds pair caches each call. No re-parsing.

### ExprPtr expression trees

`yell::ExprPtr = shared_ptr<const Expr>`. Built at parse time, evaluated at refinement
time:

```cpp
ExprPtr e = b.add("x", 0.5);          // ParamRef node, index into VectorXd
ExprPtr f = 2.0 * e + yell::lit(1.0); // arithmetic tree
double v = f->eval(p);                 // fast eval
auto d = f->eval_d(p);                 // Eigen autodiff: d.value(), d.derivatives()
```

`ExprFormulaParser` is a Spirit Qi grammar that parses formula strings into ExprPtr.
**Important**: it has no internal whitespace handling — use `"x+1.0"` not `"x + 1.0"`.
It is used inside `lexeme[]` context in the atom rule.

### Pair generation

Only these modifiers generate pairs (i.e. `generates_pairs() == true`):
- `SubstitutionalCorrelation`
- `DoubleADPMode`
- `StaticShift`
- `SizeEffect`

`CellShifter` and `MultiplicityCorrelation` do **not** generate pairs — they only
modify existing ones. A model that has only `Multiplicity` correlators will produce
zero diffuse scattering.

### Two-component Variant syntax

To get non-zero diffuse scattering you need a 2-component Variant with
`SubstitutionalCorrelation`:

```
V = Variant [ (p=0.5) C 1 x 0 0 Uiso (p=0.5) Void ]
Correlations [
  [ (1,0,0) SubstitutionalCorrelation(V,V,0.5) ]
]
```

`Void` is the keyword for an empty chemical unit (not `[]`).
For `s1=s2=2`, `SubstitutionalCorrelation` takes `(s1-1)*(s2-1) = 1` parameter.

---

## Boost.Spirit Qi — Common Pitfalls

### phoenix::bind vs std::bind (C++20)
`using namespace std;` (from AtomicPairs.h) + `using phoenix::bind;` → std::bind wins.
**Fix**: always write `phoenix::bind(&func, ...)` explicitly. Never import with `using`.

### phoenix::bind with free functions
Must use `&`: `phoenix::bind(&free_function, ...)` — not `phoenix::bind(free_function, ...)`.

### phoenix::try_/catch_ broken in Boost ≥ 1.88
Replace with a plain C++ wrapper function called via `phoenix::bind`:
```cpp
static void safe_foo(Out& out, ..., bool& _pass) {
    try { out = foo(...); } catch (...) { _pass = false; }
}
// in rule:
rule[phoenix::bind(&safe_foo, _val, ..., _pass)]
```

### Grammar copy assignment
Boost.Spirit grammars are non-copyable (rule members deleted copy assignment).
Never write `efp = ExprFormulaParser();` — keep one instance and call
`initialize_refinable_variables()` on it to reset state.

### CLion CMake cache
CLion re-runs cmake on edit, which can clobber `CACHE PATH` variables.
The `BOOST_INCLUDEDIR` uses `FORCE` to prevent this.

---

## Refinement Pipeline

```
model.txt
    ↓ (InputParser + FormulaParser + ExprFormulaParser)
Model (pools, atoms, correlators, parameterized_atoms_)
    ↓ calculate(params) → intensity_map
CeresMinimizer (DynamicNumericDiffCostFunction)
    ↓ iterates params
refined_parameters.txt  +  refinement_trajectory.json
```

The minimizer calls `Model::calculate(params)` as a black box at every step.
Currently uses Ceres `DynamicNumericDiffCostFunction` (numerical Jacobian via Ceres).

---

## Roadmap (replacing_minimizer branch)

See task list for full detail. High-level sequence:

1. Move model parsing to constructor (remove lazy parse in `calculate`)
2. Audit: ensure pools, correlators, MolecularScatterers, Variants all owned by Model
3. Check if ADPModes can be parameterized with ExprPtr (if not: roadmap item)
4. End-to-end test: `single_pair` refines to Scale=1, p_AuAu_100=0.25
5. Add `Derivatives finite_difference|analytical` input option
6. Add `JacobianMultiplier` input option (scale Jacobian columns written to disk)
7. Finite-difference derivatives computed up to pairs; skip pairs with no derivatives
8. Introduce `PattersonPeak` class (integer scatterer indices, no real/average flag)
9. Reimplement intensity calculation over list of `PattersonPeak`
10. Derivatives through FFT path; test atom scatterer power derivatives
11. Full refinement loop with analytical derivatives
