# Finding the independent components of a mode-based correlation

Status: **design note for the future.** Not implemented. The current
`ADPCorrelation` / `AnharmonicCorrelation3/4` take the full, explicit list of
combined-basis coefficients and the user ties redundant ones by hand (symbolic
expressions). This note describes how that hand-tying could be *derived
automatically* for arbitrary modes (atomic / collective / rotational / anything).

---

## 1. The problem

A correlation is specified by **left modes** `{L_a}` (act on atom 1 of a pair) and
**right modes** `{R_b}` (act on atom 2), plus a fully-symmetric rank-`n`
**mode-amplitude cumulant** tensor `κ` over the combined basis `M = L ∪ R`
(`K = |M|` modes). The user types `C(K+n−1, n)` numbers.

But the **observable** is the pair PDF, which depends only on the cumulants of the
**relative displacement** `Δu = u₂ − u₁` — a symmetric rank-`n` tensor in the **3
spatial directions** (10 for n=3, 15 for n=4). The assembly is linear:

```
T_spatial = A · κ ,     T_{i…} = Σ_{a…∈M} κ_{a…} (e_a)_i (e_b)_j …      (n factors)
```

where `e_a` is mode `a`'s contribution to `Δu`:

```
e_a = −d_a(atom1)         for a left mode      (enters Δu = u₂−u₁ with −)
e_a = +d_a(atom2)         for a right mode
e_a = d_a(atom2) − d_a(atom1)   general / collective (differential across the pair)
```

`d_a(atom)` is the mode's displacement of that atom — works for atomic translations,
collective molecular translations, **rotations** (`d_a = ω_a × (r_atom − r_centre)`),
or anything that yields a per-atom 3-vector. The map `A` does not care how `e_a` was
produced; it only needs the `3×K` matrix `E = [e_1 … e_K]`.

Because `K` modes generally exceed the 3 spatial directions, `A` is **rank-deficient**:
many `κ` map to the same `T_spatial`. The redundant directions (`ker A`) are
**unobservable** and must be fixed/tied, or the refinement is singular. This is the
n=2 `ADPCorrelation(Cu_x,Au_y)` vs `(Cu_y,Au_x)` degeneracy, one order up.

---

## 2. The structure: `A = Sym^n(E)`

`T_{ijk} = Σ_{abc} κ_{abc} E_{ia} E_{jb} E_{kc}` is the symmetric `n`-th power of the
linear map `E : ℝ^K → ℝ³`. Two consequences:

- **Image** (= the genuinely refinable observables) = `Sym^n(col E)`. Its dimension is
  `C(d+n−1, n)` where `d = rank E ≤ 3`. So **at most 10 (n=3) / 15 (n=4)** independent
  parameters, fewer when the modes span fewer than 3 directions (e.g. planar or 1-D
  motion), and fewer still under site symmetry.
- **Kernel** (= the redundant directions to tie/fix) = `ker(Sym^n(E))`, of dimension
  `C(K+n−1,n) − dim(image)`.

So "how many independent components" is exactly `dim Sym^n(col E)`, and "which ties" is
a basis of `ker`.

### Collective modes touch several atom pairs

A collective mode modifies many atom pairs `(i,j)` at once, each with its **own** `E^{(ij)}`
(different `e_a` per pair). The true observable is the **stack** over all generated pairs:

```
A = [ Sym^n(E^{(ij)}) ]  stacked over all (i,j) in the correlation
```

`ker A` is the intersection of the per-pair kernels — collective modes are *more*
constrained (more pairs ⇒ fewer redundancies) than a single atomic pair. Symmetry-
equivalent pairs add no new information (their rows are linearly dependent), so one
representative per orbit suffices, but including all is harmless.

---

## 3. Algorithm

```
INPUT:  left modes, right modes, order n, the generated pairs {(i,j)}
OUTPUT: r = #independent params; a pivot set (free comps); ties for the rest

1. Build E^{(ij)} (3×K) for every generated pair from the e_a rule (§1).
2. Form A by stacking Sym^n(E^{(ij)}):
     - rows  = (pair, spatial component)         [#pairs × {10 or 15}]
     - cols  = combined-basis κ components        [C(K+n−1, n)]
     - entry = Σ over distinct permutations of the column's mode-tuple of
               ∏ e over the row's spatial axes   (the geometric factor; same
               multiplicity bookkeeping as the assembler already does)
3. r = rank(A).                       (SVD; threshold singular values)
4. Independent set: column-pivoted QR / RREF of A → r pivot columns = the comps
   to leave free.
5. Ties: for each non-pivot column c, RREF gives c = Σ α_k · pivot_k. Emit the
   symbolic relation  κ_c = Σ α_k κ_{pivot_k}.  (Most α are 0/±1 for axis-aligned
   atomic modes — reproducing today's hand-tying — but rational in general.)
6. (Optional) site symmetry of the pair vector can be folded in as extra rows
   (R⊗…⊗R − I) so the pivots/ties already respect it. We currently leave this to
   the user.
```

### Exact vs numeric

- **Numeric (SVD/QR):** robust, easy, gives `r` and a pivot set immediately. Ties come
  out as floating-point `α`; round to rationals for display. Good enough to *report*
  "you have `r` independent params; here is a maximal free set."
- **Exact (rational RREF):** needed if we want to *emit* clean symbolic ties as
  parameter expressions. The `E` entries are algebraic (mode vectors in the lattice
  basis); rational arithmetic on `Sym^n(E)` is feasible for the small sizes here
  (≤126 columns, ≤15·#pairs rows).

### Cost

Tiny: one `(rows)×(≤126)` SVD/RREF per correlation at parse time. Negligible next to a
single intensity evaluation.

---

## 4. What this would buy

- Auto-report the true DOF count and flag over-parameterised input.
- Auto-generate the redundancy ties the user types by hand today (the `U_xy`
  symmetrisation and its rank-3/4 analogues), uniformly for atomic, collective, and
  rotational modes.
- A single code path for all orders (`n=2` reproduces `ADPCorrelation`).

It does **not** replace crystallographic *site-symmetry* constraints unless step 6 is
enabled; those remain the user's responsibility per the main plan.

---

## 5. Reuse

Step 2's geometric factor is exactly the assembler in `AnharmonicCorrelation::modify_pairs`
(the per-`(pair, spatial-comp, mode-tuple)` product). Factor it into a routine that can
emit either an `ExprPtr` (assembly) or a `double` matrix column (this analysis), so the
two never drift. See `GRAM_CHARLIER_PLAN.md` §8.5.
