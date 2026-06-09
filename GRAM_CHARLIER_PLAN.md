# Anharmonic Gram–Charlier PDF — Implementation Plan

Status: on `parallel_refinement`. **Phases 1–8 done** (toolkit, atom coefficients,
pair combination + AnharmonicCorrelation, forward G(s) in direct & FFT, symmetry
guard, refinement recovery, analytical derivatives, docs). Only the deferred niceties
remain: bundled ADPCorrelation matrix sugar; negative-PDF monitoring; auto-detection
of independent components (`GRAM_CHARLIER_INDEPENDENT_COMPONENTS.md`). Refinement tests
(both derivative modes) live in test_parallel_gtest.cpp (links the Ceres minimizer).

---

## 0. The one-sentence reason this is feasible

Yell's FFT engine computes the scattering of **each pair on its own small reciprocal-space
grid and inverse-FFTs it into PDF space** (`IntnsityCalculator::calculate_patterson_map_from_pairs_f`,
`Calculator.h:188`). In reciprocal space the Gram–Charlier anharmonic temperature factor is
just the harmonic Gaussian **multiplied by a low-order polynomial in the reciprocal vector
`s`**. So adding anharmonicity = multiplying the per-pixel value by one extra complex factor
`G(s)`, for the handful of pairs that need it. No new transform, no new data flow. That is
exactly what "the algorithm is sort of ready for Gram–Charlier" means.

---

## 1. Physics & conventions (pin these first)

### 1.1 The temperature factor

Yell's per-pair scattering at reciprocal point `s = (h,k,l)` (continuous fractional reciprocal
coords — confirmed: `Grid::s_at` feeds `s` straight into `cell.d_star_sq(miller::index)`):

```
T_pair(s) = p · conj(f1) f2 · exp( +2πi (s·r)  − 2π² (s·U·s) )      // current harmonic form
```

The displacement PDF's characteristic function is `E[exp(2πi s·u)]`, whose **cumulant
(log) expansion is exact and additive under convolution**:

```
log T(s) = 2πi (s·r)                          // n=1, mean      → r
           − 2π²   (s·U·s)                     // n=2, U         (κ2)
           − (4π³ i / 3)  C(s,s,s)             // n=3, skewness  (κ3)
           + (2π⁴ / 3)    D(s,s,s,s)           // n=4, kurtosis  (κ4)
           + ...
```

(the n-th term is `(2πi)^n/n! · κ_n` contracted with `s`; check: n=2 gives `−2π² κ2`,
matching the existing `M2PISQ = −2π²` exponent and `M_2PI` phase exactly).

The **Gram–Charlier series** is the expansion of the exponential of the higher terms:

```
T(s) = T_harmonic(s) · G(s),    G(s) = 1 − (4π³ i/3) C(s,s,s) + (2π⁴/3) D(s,s,s,s) + ...
```

- `C(s,s,s) = Σ_{jkl} C_{jkl} s_j s_k s_l`  — fully symmetric rank-3 tensor, **10** unique comps.
- `D(s,s,s,s) = Σ_{jklm} D_{jklm} s_j s_k s_l s_m` — fully symmetric rank-4, **15** unique comps.

**Grounded in ITA §6.1.1.6** (`references/ITC_sec6o1o1_incl_gram_Charlie.pdf`):
the GC temperature factor (6.1.1.47) and the cumulant form (6.1.1.63) are

```
GC:        T(S) = T₀(S)·[ 1 + (i³/3!) c^{jkl} SⱼSₖSₗ + (i⁴/4!) c^{jklm} SⱼSₖSₗSₘ + … ]
cumulant:  T(S) = exp[ iκʲSⱼ − ½κ^{jk}SⱼSₖ − (i/6)κ^{jkl}SⱼSₖSₗ + (1/24)κ^{jklm}SⱼSₖSₗSₘ + … ]
```

ITA's `S` folds in 2π (6.1.1.92, `F=∫ρ e^{iS·r}`); Yell keeps 2π explicit with `s=(h,k,l)`.
Substituting `S = 2π s` reproduces the prefactors above exactly (`−4π³i/3`, `+2π⁴/3`).
**Key truncation result:** equating the two forms and dropping terms below `S⁶` (the `(κ³)²`
cross term is 6th order) gives `c^{jkl} = κ^{jkl}` and `c^{jklm} = κ^{jklm}` — i.e. **to 4th
order the GC coefficients *are* the cumulants.** This is what makes the pair combination (§1.2)
linear and exact.

`G(s)` is **complex**: the rank-3 term is imaginary (skewness ⇒ asymmetric PDF), the rank-4
term is real.

> **Validation anchor:** still confirm the overall sign against an independent Python/numpy
> reference for one atom (§7). The cctbx in `lib/cctbx_stubs` is a *stub* (only uctbx/sgtbx/eltbx)
> — it ships no anharmonic ADP code, so we own the convention.

### 1.2 Atom → pair combination rule (uncorrelated baseline)

The interatomic vector is `u = u2 − u1`. For **independent** displacements the pair
characteristic function factorizes, `T_pair(s) = T₂(s)·T₁(−s)`; cumulants are additive under
convolution and the `−s` flips **odd** terms (ITA 6.1.1.63). Adding the logs:

| order | quantity | pair value (already / proposed)        | sign |
|------:|----------|----------------------------------------|------|
| 1     | `r`      | `r2 − r1`   (exists, `AtomicPairs.h:103`) | subtract (odd) |
| 2     | `U`      | `U1 + U2`   (exists, `AtomicPairs.h:105`) | add (even) |
| 3     | `C`      | `C2 − C1`   (new)                       | subtract (odd) |
| 4     | `D`      | `D1 + D2`   (new)                       | add (even) |

Because `c=κ` at the truncation (§1.1), the GC coefficients combine by the same add/subtract
rule — exact to 4th order (the product `G₂(s)·G₁(−s)` only generates cross terms at `s⁶`). So
any pair touching an anharmonic atom inherits `±` that atom's tensor (sign by which end it sits),
rotated into the pair frame. This is the natural continuation of the existing `r2−r1`, `U1+U2`,
and feeds both the full and average peak tracks.

### 1.3 Free anharmonic correlations (harmonic atoms ≠ zero pair GC)

The pair tensors are cumulants of the **relative** displacement `Δu = u2 − u1`. Marginal
(single-atom) Gaussianity does **not** force the *joint* two-atom distribution to be Gaussian:
an anharmonic interatomic potential / nonlinear neighbour coupling gives `κ3(Δu), κ4(Δu) ≠ 0`
even when every atom is on-average harmonic. So **pair GC coefficients are legitimate free
refinable parameters**, constrained only by the pair-vector symmetry — they vanish only if the
joint is genuinely Gaussian (independence is the special case where they reduce to §1.2).

This maps onto Yell's existing two-layer design for `U`:

```
atomic U  →  baseline pair U = U1+U2   ,  ADPCorrelation/DoubleADPMode adds correlated U on top
atomic c  →  baseline pair c = c2∓c1   ,  AnharmonicCorrelation (new) adds FREE pair c on top
```

Both sources write into the same `tensor3_expr/tensor4_expr` on the `AtomicPair`, so the forward
`G(s)` and derivative code are agnostic to which filled them.

### 1.4 Symmetry constraints (user-supplied)

Pair tensors are contravariant rank-3/4 tensors at the pair-vector site in the **Patterson
(Laue) group**: invariant under that site's point group (`R⊗R⊗R · C = C`). ITA Table 6.1.1.10
lists allowed indices per site symmetry; the Patterson is always centrosymmetric, so a
self-centrosymmetric pair kills all odd (3rd-order) components while 4th survives. **We do not
derive these** — the user supplies the allowed independent components (dependents tied via
expressions). The engine only needs to *rotate* `C`,`D` when the Laue expansion replicates a
pair (§4).

### 1.5 Units / basis — Ångströms everywhere

**Decision: Å for everything, a\* baked internally — uniform with U.** The user enters `U` in Å²,
`C` in Å³, `D` in Å⁴; Yell bakes the matching `a*` products at parse time so all tensors contract
with `s=(h,k,l)`. `construct_atom` already does this for `U` (`U_ij → U_ij·a*_i·a*_j`,
`model.h:249`); extend identically: `C_{jkl} → C·a*_j a*_k a*_l`, `D_{jklm} → D·a*_j a*_k a*_l a*_m`.

> Rejected alternative: copy CIF/JANA `C_ijk/D_ijkl` verbatim (already a\*-folded, no conversion).
> That would make U Å² but C/D a different basis — inconsistent. We keep one rule: Å in, a\* baked.
> Component **order** still follows CIF lexicographic (`111 112 113 122 123 133 222 223 233 333`)
> for paste-friendliness; that is independent of units. Phase 0 confirms the per-component a\*
> products against one known JANA atom.

---

## 2. Where each change goes (code map)

| Concern | File / symbol | Change |
|---|---|---|
| tensor types + contraction + symmetry transform | new `src/anharmonic.h` | `tensor3`,`tensor4`, `tensor3_expr`,`tensor4_expr`; `contract3(C,s)`,`contract4(D,s)`; rank-3/4 `operator*(mat3, …)` |
| atom storage | `ChemicalStructure.h` `Atom` | add `yell::ExprPtr C[10], D[15]`; default `lit(0)`; extend `update_caches`, `reset` |
| atom construction + a\* bake | `model.h` `construct_atom*` | accept optional anharmonic exprs, bake conversion |
| input syntax | `InputFileParser.cpp` atom rule (`:440`) | optional `GramCharlier3[…] GramCharlier4[…]` suffix |
| pair storage | `AtomicPairs.h` `ParameterizedParams` / `AtomicPair` | add `tensor3_expr C; tensor4_expr D;` + accessors; ctor sets `C2−C1`, `D1+D2` |
| bundled ADPCorrelation | `InputFileParser.cpp` grammar | sugar: `[L…]×[R…] + matrix` → `Σ DoubleADPMode`; matrix avoids U ordering |
| free pair correlation | `AtomicPairs.h` new `AnharmonicCorrelation : PairModifier` | non-pair-generating modifier (like `DoubleADPMode`); mode list + symmetric mode-cumulant tensor, assembled via outer products onto a pair's spatial `C`,`D` |
| pair symmetry | `LaueSymmetry` setup / `apply_patterson_symmetry` | **guard only** (Option B): error if `generators_on_vectors` non-empty AND any pair anharmonic. No tensor rotation — map track handles `G(s)` for free. See §4 / §8.6. |
| baked peak | `AtomicPairs.h` `PattersonPeak` + new side-table | add lean `int anh_idx` (−1 = harmonic); `vector<AnharmonicData>` holds `C[10],D[15]` |
| bake peaks | `AtomicPairs.h` `peaks_from_pairs` | fill side-table only for anharmonic pairs |
| forward FFT | `Calculator.h` new `…_from_pairs_anharmonic_f` | clone of the hot loop × `G(s)`; hot loop untouched |
| forward direct | `Calculator.h` `calculate_scattering_from_patterson_peaks` | `× G(s)` when `anh_idx≥0` (ground-truth path) |
| derivatives (phase 7) | `PeakSusceptibility`, both derivative loops | add `d_C`, `d_D`; product rule on `G` |
| output | `model.cpp` reporting | write refined `C`,`D` + esds |

---

## 3. Keeping the harmonic path fast (non-negotiable constraints)

1. **Partition peaks into two lists**: harmonic (the 1000s) and anharmonic (the dozens).
   The existing hot routine `calculate_patterson_map_from_pairs_f` stays **byte-for-byte
   unchanged** and runs on the harmonic list. A new sibling routine handles the anharmonic
   list and accumulates into the **same** map under the same `accum_mutex`. ⇒ zero perf and
   zero behavioural regression for purely-harmonic models (the overwhelming majority).
2. **Lean `PattersonPeak`**: do **not** inline 25 doubles into the peak struct (would ~triple
   it and thrash cache for the harmonic 1000s). Keep one `int anh_idx`; C/D live in a side
   `vector<AnharmonicData>` touched only for the dozens.
3. **Symmetry transform**: `multiply_pairs_by_matrix` skips the rank-3/4 tensor rotation unless
   the pair carries anharmonic terms (cheap per-pair flag).
4. **Lazy allocation**: anharmonic side-structures only allocated when ≥1 anharmonic atom is
   parsed. Models with none pay nothing, allocate nothing.
5. `G(s)` for an anharmonic peak is ~10–15 FMAs/pixel on dozens of peaks — negligible.

---

## 4. Implementation order (each phase independently testable, builds green)

**Phase 0 — Prototype the math, no Yell.** numpy/scratch script: build `G(s)` and the
atom→pair combination, FFT one pair both "analytic T(s)" and "Gaussian × G(s)", confirm they
agree. Freezes the sign/prefactor convention before any C++.

**Phase 1 — Tensor toolkit (`anharmonic.h`), pure, unit-tested in isolation.**
`tensor3/tensor4` (double + expr), `contract3/contract4` with correct **multinomial
multiplicities** (1·xxx, 3·xxy, 6·xyz for rank-3; 1/4/6/12 for rank-4), rank-3/4 symmetry
transforms. gtest vs brute-force 27-/81-term sums. No wiring yet.

**Phase 2 — Atom plumbing + parser.** Fields on `Atom`, `construct_atom*` accepts optional
anharmonic exprs with a\* bake, new grammar suffix. Test: parse a model, round-trip the stored
exprs. Harmonic models still parse identically.

**Phase 3 — Pair combination + free-correlation modifier.** `AtomicPair` carries `C = C2−C1`,
`D = D1+D2`; `peaks_from_pairs` fills the anharmonic side-table; `anh_idx` set. Add the
`AnharmonicCorrelation` modifier (§1.3) that adds free pair-level `C`,`D` on top, plus its
grammar (a `Correlations [...]` entry, analogous to `ADPCorrelation`). Test: combination signs;
harmonic pairs ⇒ `anh_idx == −1`; modifier adds onto the baseline.

**Phase 4 — Symmetry guard (Option B, decided).** Laue symmetry is applied either on the map
(fast, default) or — when the grid is incompatible with a generator — on the pairs by physical
rotation (`generators_on_vectors`, `multiply_pairs_by_matrix`). The **map track needs no tensor
rotation**: a mate pair's contribution equals the ASU pair evaluated at `R⁻¹s`, so map averaging
absorbs the rotation for `U` *and* for the anharmonic `G(s)` automatically — zero new code. The
pair track *would* need rank-3/4 rotation, but **we do not implement it now.** Instead, at setup,
if `generators_on_vectors` is non-empty **and** any pair is anharmonic, throw a clear error:
"anharmonic terms require a symmetry-compatible grid (Laue symmetry must be applied on the map);
adjust your DiffuseScatteringGrid." The rank-3/4 `operator*` stays specced in §8.1 but
**unimplemented** — promoting to full rotation later is additive and isolated. The user supplies
symmetry-allowed components directly (no constraint derivation, ever).

**Phase 5 — Forward intensity (the milestone).** Add `G(s)` to the **direct** path (ground
truth) and the new **anharmonic FFT** routine. Test: **direct vs FFT must match** on a small
anharmonic model — this is the primary correctness harness and catches every sign/multiplicity
bug. Also compare a single anharmonic pair against the Phase-0 numpy reference.

**Phase 6 — Refinement, finite-difference.** With `Derivatives finite_difference` the Ceres
numeric-diff path calls `calculate()` as a black box, so anharmonic parameters refine with **no
extra derivative code**. Test: simulate data with known `C`,`D`, refine from zero, recover them.

**Phase 7 — Analytical derivatives (optional perf).** Extend `PeakSusceptibility` with
`d_C`,`d_D`; apply product rule `∂(base·G) = ∂base·G + base·∂G` in both derivative loops
(`Calculator.h:299`, `:404`). Test: analytical Jacobian columns vs finite-difference. Run under
TSan on disorder-s01 (the new FFT routine uses the same clone machinery; see CLAUDE.md).

**Phase 8 — Output & docs.** Report refined `C`,`D` with esds; document syntax; warn about
non-positive PDF (Gram–Charlier can go negative — monitor and optionally flag).

Forward correctness (Phases 0–6) is the bulk of the value; analytical derivatives (7) are a
speed optimisation that can land later.

---

## 5. Debugging strategy

- **Direct vs FFT is the workhorse.** Yell keeps both a direct (`calculate_scattering_from_patterson_peaks`)
  and an FFT engine. Implement `G(s)` in both; any divergence = a bug in the new FFT routine.
  (No existing direct-vs-FFT gtest was found — add one as part of Phase 5; it pays off forever.)
- **Numeric vs analytical derivatives** (Phase 7) — column-by-column compare.
- **`GramCharlierOrder 0|3|4` switch** to disable anharmonicity at runtime for bisecting a
  suspected regression without recompiling.
- **Single-pair analytic check** against the Phase-0 reference.
- **TSan** on the new parallel routine (clone path), per the existing CLAUDE.md workflow.

---

## 6. Decisions for you (defaults chosen, easy to change)

Symmetry constraints are **out of scope** — the user provides symmetry-allowed components
(independent comps as parameters, dependents via expressions). Remaining choices:

1. **Where to declare anharmonicity** — *Recommended:* support **both** — per-**atom**
   coefficients that auto-seed the pair baseline (§1.2), **and** a per-pair `AnharmonicCorrelation`
   modifier for free correlated anharmonicity (§1.3). Both write the same pair tensors.
2. **Order to support** — *Recommended:* both 3rd **and** 4th from the start (15+10 comps, same
   machinery). Could ship 3rd-only first if you want a smaller first cut.
3. **Syntax (settled — see `gram_charlier_dummy_model.txt`):** atom-line tail
   `GramCharlier3[10] GramCharlier4[15]` (CIF order, Å); bundled
   `ADPCorrelation([L…],[R…],[matrix])`; free `AnharmonicCorrelation3/4(unit,unit,[tensor])`.
   Units Å everywhere, a\* baked internally.
4. **Laue expansion of anharmonic pairs (settled — Option B):** no rank-3/4 rotation; the map
   symmetry track handles `G(s)` for free, and a grid-incompatible generator + an anharmonic pair
   throws a clear error. Full rotation (Option A) is deferred and additive.

---

## 7. Reference for the convention (use to freeze signs in Phase 0)

Gram–Charlier anharmonic ADP background and the `C^{jkl}/D^{jklm}` formalism:

- Kuhs (1992/2023) review — *Gram–Charlier approach for anharmonic atomic displacements*,
  Crystallography Reviews 29(3). https://www.tandfonline.com/doi/abs/10.1080/0889311X.2023.2266400
- Johnson & Levy treatment / IT Vol B §1.2 (cumulant ↔ Gram–Charlier).
- Worked refinement example (component counts, special positions):
  https://pmc.ncbi.nlm.nih.gov/articles/PMC8196607/
- Negative-density caveat in GC refinements:
  https://www.researchgate.net/publication/327425948

---

## 8. Class-level implementation plan

This maps the phases onto concrete classes, mirroring existing patterns (`sym_mat3_expr`,
`PairModifier`, `StructurePartRef` factories). New code is **additive** — every struct gains
optional fields that default to the harmonic case, so existing models compile and run unchanged.

### 8.1 New file `src/anharmonic.h` — tensor toolkit (Phase 1)

The value types and the contraction/rotation algebra. Pure, no Yell deps beyond `expr.hpp`.

```cpp
// Canonical unique-component order + multinomial multiplicities (the only place these live):
//   rank3 (10): 111 222 333 | 112 122 113 133 223 233 | 123
//               mult:  1  1  1 |  3   3   3   3   3   3 |  6
//   rank4 (15): 1111 2222 3333 | 1112 1113 1222 1333 2223 2333 | 1122 1133 2233 | 1123 1223 1233
//               mult:   1   1   1 |   4    4    4    4    4    4 |    6    6    6 |  12   12   12
struct tensor3 { double c[10]; /* idx→(i,j,k) table + mult[] as constexpr */ };
struct tensor4 { double d[15]; };

double contract3(const tensor3& C, const vec3<double>& s);  // Σ mult·C·sᵢsⱼsₖ
double contract4(const tensor4& D, const vec3<double>& s);

// The Gram–Charlier multiplier G(s) for a baked peak (§1.1):
std::complex<double> gram_charlier_factor(const vec3<double>& s,
                                          const tensor3& C, const tensor4& D);
//   = 1 + complex(0, -4*M_PI*M_PI*M_PI/3 * contract3) + (2*pow(M_PI,4)/3 * contract4)

// Expr-tree mirrors of sym_mat3_expr (built at parse, evaluated at bake):
struct tensor3_expr { yell::ExprPtr c[10]; /* default lit(0) */ };
struct tensor4_expr { yell::ExprPtr d[15]; };
tensor3_expr operator+(…), operator-(…), operator-(const tensor3_expr&); // and tensor4_expr
// DEFERRED (future Option A, not built in this pass — see §8.6):
//   tensor3_expr operator*(const mat3<double>& R, const tensor3_expr&);   // rank-3 rotation
//   tensor4_expr operator*(const mat3<double>& R, const tensor4_expr&);   // rank-4 rotation

// Symmetrized outer products for AnharmonicCorrelation assembly (§8.5):
tensor3 outer3(const vec3<double>& a, const vec3<double>& b, const vec3<double>& c);
tensor4 outer4(const vec3<double>&, const vec3<double>&, const vec3<double>&, const vec3<double>&);
// (+ _expr-weighted variants accumulating coeff·d⊗d⊗… into a tensor*_expr)
```

The two `operator*(mat3,…)` are the rank-3/4 analogues of `expr_types.h:106`: expand the unique
component to the full 3³/3⁴ array, apply `R` on each index, refold — emitting `ExprPtr` sums
with `double` coefficients (EvaluationCache dedups). gtest each contraction and rotation against
a brute-force 27-/81-term loop; round-trip a full-order rotation to identity.

### 8.2 `Atom` (ChemicalStructure.h) — atomic coefficients (Phase 2)

```cpp
yell::ExprPtr C[10], D[15];      // default lit(0); fractional-reciprocal basis (a* baked)
bool          anharmonic = false;
double        C_cache[10], D_cache[15];   // filled by update_caches alongside U_cache
```

`update_caches()` evaluates `C/D` into the caches; `reset()` sets them to `lit(0)` and
`anharmonic=false`. A setter `set_anharmonic(C_exprs, D_exprs)` flips the flag. Untouched atoms
stay byte-identical to today.

### 8.3 `model.h` construct/factory + parser (Phase 2)

- `construct_atom` gains an optional anharmonic tail; bake `a*` products (`C_jkl·a*_j a*_k a*_l`,
  rank-4 analogously) exactly like the existing `U` bake (`model.h:249`), then `set_anharmonic`.
- Grammar (`InputFileParser.cpp` atom rule `:440`): optional suffix
  `GramCharlier3[ 10×expr ] GramCharlier4[ 15×expr ]` (both omittable). Harmonic atom lines
  parse identically.

### 8.4 `AtomicPair` / `ParameterizedParams` (AtomicPairs.h) — pair baseline (Phase 3)

```cpp
struct ParameterizedParams { … existing …; tensor3_expr C; tensor4_expr D; };  // default lit(0)
// AtomicPair ctor, both real & average tracks, when either atom is anharmonic:
//   C = a2.C(rank3_expr) - a1.C ;   D = a1.D + a2.D            // §1.2 odd-sub / even-add
//   anharmonic_ = a1.anharmonic || a2.anharmonic
bool anharmonic_ = false;
tensor3_expr& C(bool avg=false);  tensor4_expr& D(bool avg=false);   // accessors like U()
```

`get_pair` is unchanged; only the ctor seeds tensors. The `anharmonic_` flag gates all later
per-pair anharmonic work so the 1000s of harmonic pairs pay nothing.

### 8.5 Bundled correlations (Phase 3) — `ADPCorrelation` matrix + `AnharmonicCorrelation3/4`

Two related changes. **(a) Bundle the existing `ADPCorrelation`** into mode-list × mode-list +
amplitude matrix (sugar over the current per-component form), and **(b) add free pair
`AnharmonicCorrelation3/4`** as direct crystal-axis tensors.

**(a) Bundled ADPCorrelation** — *deferred* (independent sugar, not yet built). Pure parse-time
desugaring, no new runtime class: `ADPCorrelation([L…],[R…],[matrix])` → `Σ matrix[a][b]·
DoubleADPMode(L_a,R_b)`. Matrix is direction-indexed, so no U 6-vector ordering. Legacy scalar
form stays valid.

**(b) AnharmonicCorrelation3/4** — the higher-cumulant analogue of `ADPCorrelation`, in the
**same mode basis**. The unifying picture: mode amplitudes `q_m` have cumulants, and each
correlation supplies one cumulant order, assembled into a spatial tensor via the mode vectors
`d_m`:

```
ADPCorrelation        : U = Σ_ab   ⟨q_a q_b⟩      d_a⊗d_b           (rank-2, exists as DoubleADPMode)
AnharmonicCorrelation3: C = Σ_abc  ⟨q_a q_b q_c⟩  d_a⊗d_b⊗d_c       (rank-3, NEW)
AnharmonicCorrelation4: D = Σ_abcd ⟨q…⟩           d_a⊗d_b⊗d_c⊗d_d   (rank-4, NEW)
```

Input is a **single mode list + a fully-symmetric rank-n tensor over those M modes** (the
mode-space cumulants, *not* the raw crystal tensor — they coincide only for unit `[x,y,z]`
modes). Coefficient count `C(M+n−1, n)` (M=3 → 10/15), CIF order over 1-based mode positions.
Parameter count therefore follows the modes — one mode ⇒ one number (1-D anharmonicity / a
libration). This dissolves the "rank-3 vs 2-atom" worry: modes carry the atoms (like
`DoubleADPMode`), and the assembled pair tensor is just a symmetric spatial rank-n object.

```cpp
class AnharmonicCorrelation : public PairModifier {   // order = 3 or 4
  vector<ADPMode*> modes;                 // the direction basis (M modes)
  vector<yell::ExprPtr> coeffs;           // C(M+n-1,n) symmetric cumulants, CIF order
  int order;                              // 3 or 4
  bool generates_pairs() override { return false; }
  void modify_pairs(AtomicPairPool* pool) override {
    // for each ordered n-tuple of modes, accumulate coeff · (d⊗d⊗…) onto the
    // spatial tensor of pool->get_pair(atom_a, atom_b); set p.anharmonic_ = true.
    // Needs symmetrized outer-product helpers vec3⊗vec3⊗… → tensor3_expr/tensor4_expr.
  }
  PairModifier* clone() const override { return new AnharmonicCorrelation(*this); }
};
```

- Outer-product assembly helpers (`anharmonic.h`): `outer3(d_a,d_b,d_c) → tensor3`, `outer4 → tensor4`
  (and their `_expr` forms), summed over all ordered n-tuples ⇒ symmetric spatial tensor.
- Factory `Model::create_anharmonic_correlation(vector<StructurePartRef> modes, int order, vector<ExprPtr> coeffs)`
  (sibling of `create_double_adp_mode`, `model.h:491`).
- Grammar: bundled `adp_correlation` (mode-list × mode-list + matrix) **and**
  `anharmonic_correlation` (`[modes] , [coeffs]`) added to the `atomic_pair_pool` alternatives
  (`InputFileParser.cpp:242`), via `do_*` safe-wrappers.
- The pair still stores the **spatial** `tensor3_expr/tensor4_expr` (§8.4); §8.6 symmetry rotation
  and §8.7 baking are unchanged — only the *source* of the spatial tensor differs.

### 8.6 Symmetry guard (Option B) — Phase 4

No tensor rotation is implemented. In `apply_patterson_symmetry` (or at model setup), if
`generators_on_vectors` is non-empty **and** any pair has `anharmonic_ == true`, throw with the
"anharmonic needs a symmetry-compatible grid" message. The map-symmetry track (the default)
handles anharmonic `G(s)` correctly with no code — see §4. `multiply_pairs_by_matrix` is left
untouched; the rank-3/4 `operator*(mat3, tensor*_expr)` of §8.1 is specced but **not built** yet
(future Option A, additive).

### 8.7 Baking: `PattersonPeak`, side-table, susceptibilities (Phases 3 & 7)

```cpp
struct PattersonPeak { … existing …; int anh_idx = -1; };      // lean: one int, no inline tensors
struct AnharmonicData { tensor3 C; tensor4 D; };               // side array, dozens of entries

// peaks_from_pairs(...) extra out-params:
//   vector<AnharmonicData>& full_anh, avg_anh
//   if (pair.anharmonic_) { eval 10+15 exprs → AnharmonicData; anh_idx = (int)anh.size()-1; }
//   else anh_idx = -1

struct PeakSusceptibility { … existing …; tensor3 d_C; tensor4 d_D; };  // Phase 7 only
```

### 8.8 Forward calc — partition + new routine (Phase 5)

In `Model::calculate_from_peaks` (`model.cpp:174`): split the paired `full`/`avg` arrays by
`anh_idx<0` into a **harmonic** sublist and an **anharmonic** sublist (O(n), once per call).

- Harmonic sublist → existing `calculate_patterson_map_from_pairs_f` **unchanged**.
- Anharmonic sublist → new `calculate_patterson_map_from_pairs_anharmonic_f`: a copy of the hot
  loop whose only added line multiplies the per-pixel value by
  `gram_charlier_factor(ppm.current_s(), anh.C, anh.D)`. Accumulates into the **same** map under
  the same `accum_mutex`.
- Direct path `calculate_scattering_from_patterson_peaks` (`Calculator.h:373`) multiplies by the
  same factor when `pk.anh_idx>=0` — the ground-truth for the **direct-vs-FFT** test (Phase 5).

### 8.9 Derivatives (Phase 7)

`susceptibilities_from_pairs` also bakes `d_C, d_D` (autodiff column `param_idx`, same `eval_d`
pattern). The two derivative loops (`Calculator.h:299`, `:404`) apply the product rule
`∂(base·G) = ∂base·G + base·∂G`, with `∂G = -(4π³i/3)·contract3(d_C,s) + (2π⁴/3)·contract4(d_D,s)`.
Mirror the §8.8 harmonic/anharmonic partition for the derivative maps.

### 8.10 `Model::clone()` (thread safety)

`clone()` already deep-copies atoms and remaps pointers; the new `Atom` `C/D` ExprPtr and the
pair `tensor*_expr` are value/shared-ptr members that copy with the existing machinery — **no new
clone logic**, provided no clone writes shared tensor state (it doesn't; tensors are baked into
per-call peak side-tables). Re-run TSan on an anharmonic model (Phase 7) to confirm.

### Build/test order summary

`anharmonic.h` + gtest (1) → `Atom` + parser round-trip (2) → `AtomicPair`/`AnharmonicCorrelation`
+ combination test (3) → symmetry rotation test (4) → direct-vs-FFT forward test (5) →
finite-diff refinement recovery (6) → analytical derivatives vs finite-diff + TSan (7) → output
& docs (8). Each step builds green and is covered by `test-yell-gtest`.
