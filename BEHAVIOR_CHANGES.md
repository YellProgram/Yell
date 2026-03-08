# Yell Behavioral Changes and Bug Fixes (Investigation since b1c6cf21)

This document tracks changes in Yell's calculation logic, unintended bug fixes, and behavioral shifts observed during the development from commit `b1c6cf21`.

## 1. MultiplicityCorrelation logic fix (Commit ec114f14)
*   **Previous Behavior:** `MultiplicityCorrelation` was multiplying both the occupancy `p` and the `multiplier` field of an `AtomicPair` by the user-provided multiplier.
*   **Current Behavior:** Only the `multiplier` field is modified.
*   **Impact:** In the old version, if both `p` and `multiplier` were used in intensity calculation (which is standard), the multiplicity was effectively being applied twice. This fix corrected a long-standing double-counting bug.

## 2. SubstitutionalCorrelation joint probability handling
*   **Previous Behavior:** In the early `parameterized_model` branch (near `0b5dddf8`), the migration to `ExprPtr` might have caused issues where the joint probability was not correctly updated for all pairs in a pool if they shared the same chemical units.
*   **Current Behavior:** Properly uses `joint_probability_expr` to set `pair.p()` during `modify_pairs`.
*   **Impact:** Correct propagation of refined variables to structural correlations.

## 3. Scale keyword persistence (Commit 853e2870)
*   **Previous Behavior:** If a `Scale` keyword appeared in `model.txt` *before* the `RefinableVariables` block, it was silently reset to 1.0 when the refinable variables were parsed.
*   **Current Behavior:** The `Scale` value is preserved and used as the initial value even after `RefinableVariables` are processed.
*   **Impact:** Fixed an unintended reset of the user-specified scale.

## 4. Multiplicity in `correlators_from_cuns` (Commit ec114f14)
*   **Change:** Identified that the code was silently ignoring the last column and row of joint probability matrices if they were provided by the user (as they are redundant for 2x2 or 3x3 matrices where they must sum to 1).
*   **Impact:** Added a TODO to detect and warn/error on inconsistent input rather than silently fixing it.

## 5. Implementation of Parse-Once Model (Commit 0b5dddf8)
*   **Change:** Shifted from re-parsing `model.txt` on every refinement step to a single parse that builds an expression tree (`ExprPtr`).
*   **Impact:** Fixed potential performance bottlenecks and ensured consistency. The old version's re-parsing might have been susceptible to state changes if the parser was not perfectly stateless.

## 6. PattersonPeak and Cache-Friendly Intensity Loops (Commit c49958a1 and others)
*   **Change:** Replaced the direct use of `AtomicPair` in calculation loops with `PattersonPeak`, which separates "real" and "average" contributions into flat lists.
*   **Impact:** This cleanup fixed inconsistencies in how symmetry and multipliers were applied to different types of peaks (real vs average).

## 7. R-factor calculation behavior
*   **Observation:** While `Rw` was present in `b1c6cf21`, the refinement logic has been significantly updated (Ceres solver, pseudo-inverse for covariance).
*   **Impact:** The final reported `Rw` might differ slightly due to more robust convergence in the new refinement engine.
