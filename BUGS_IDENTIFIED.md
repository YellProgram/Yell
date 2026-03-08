# Bugs Identified — 2026-03-07

The following issues were identified during the refinement engine investigation and need to be addressed starting from "ground zero" (commit `853e2870`).

### 1. `MultiplicityCorrelation` Double Scaling
- **Location**: `src/AtomicPairs.h`
- **Issue**: `modify_pairs` was incorrectly multiplying the pair probabilities (`p` and `average_p`) by the multiplier.
- **Result**: Double-scaling of intensities because the `multiplier` field is already used in calculations. It should only modify the `multiplier` field.

### 2. Obsolete `filter_pairs` code in `LaueSymmetry`
- **Location**: `src/LaueSymmetry.cpp` / `apply_patterson_symmetry`
- **Issue**: The code was calling `filter_pairs_from_asymmetric_unit`, which is part of an older, dead path.
- **Result**: This function used `decrease_multiplicity` macros that divided probabilities at special positions (like the origin), causing "flattened" origin intensities (making them the same as general peaks).

### 3. Incomplete `4/mmm` Symmetry Generators
- **Location**: `src/LaueSymmetry.h`
- **Issue**: The generator list for `4/mmm` incorrectly consisted of only four order-2 mirrors/inversions.
- **Result**: Only 8 equivalent positions were generated instead of the required 16. The four-fold rotation `y,-x,z` was missing.

### 4. Unconditional Averaging in `apply_generator`
- **Location**: `src/LaueSymmetry.cpp` / `apply_generator(IntensityMap&)`
- **Issue**: A final 4-fold rotation averaging block was executed for *every* generator symbol because it was missing an `else` statement.
- **Result**: Incorrect over-averaging of intensities for all symmetry operations.

### 5. Partial Map Iteration in `apply_generator`
- **Location**: `src/LaueSymmetry.cpp`
- **Issue**: Several symmetry operations (like `mx` and `my`) used loops starting at `size/2`.
- **Result**: Assumes a specific map centering and can leave parts of the map un-averaged or zeroed out if the grid indexing doesn't match that assumption perfectly.

### 6. ADP Sandwich Transformation Semantics
- **Location**: `src/AtomicPairs.h`
- **Refactor**: The `operator*` for ADP expression transformation ($R \cdot U \cdot R^T$) is semantically confusing.
- **Recommendation**: Refactor to a named function `sandwich_product` to maintain clarity while keeping the high-level expression tree generation.
