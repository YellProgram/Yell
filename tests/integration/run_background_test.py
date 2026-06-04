#!/usr/bin/env python3
"""
Integration test for isotropic background refinement (RefineBackground).

Strategy (self-contained — the generic refinement harness cannot inject a
background because the generated experiment.h5 = model.h5 = Scale*(Ifull-Iavg)
is structural only):

  1. Generate a structural "experiment" E0 from a target model (Refine false).
  2. Inject a known isotropic background B(|q|) into E0  ->  E1 = E0 + B.
  3. Refine E1 from a perturbed starting model with RefineBackground enabled.
  4. Assert that
       - the structural parameters (Scale, p1, p2) are recovered to the target,
         i.e. the injected background did NOT bias the structural fit, and
       - the refined background coefficients match the injection.

Two cases:
  * constant  (degree 0): B = C.  Exact, needs no |q|; checks Bkg_0 == C.
  * linear    (degree 1): B = A + G*|q|.  |q| from the cell metric.  Chebyshev
    degree-1 spans {1, |q|} under any affine normalisation, so this must be
    absorbed exactly -> verifies structural recovery despite a |q|-shaped
    background, and exercises the |q| basis.

Run:  python3 run_background_test.py [--bin /path/to/yell]
"""
import os, sys, re, shutil, subprocess, argparse
import numpy as np
import h5py

BASE = os.path.dirname(os.path.abspath(__file__))

def default_bin():
    for rel in ("../../cmake-build-release/yell", "../../cmake-build-debug/yell"):
        p = os.path.abspath(os.path.join(BASE, rel))
        if os.path.exists(p):
            return p
    return os.path.abspath(os.path.join(BASE, "../../cmake-build-debug/yell"))

# Structural model shared by all cases: Scale 132, p1=0.113, p2=0.298.
def model_text(refine, scale, p1, p2, extra=""):
    return f"""Cell 3.101 3.101 1  90 90 90
LaueSymmetry 4/mmm
DiffuseScatteringGrid -30 -30 0 0.1 0.1 1 600 600 1

CalculationMethod exact
Refine {refine}
{extra}
Scale {scale}
RefinableVariables[
p1={p1};
p2={p2};
]

UnitCell
[
  AuAg = Variant[
    (p=1/2)
    Au 1  0 0 0  0.002
    (p=1/2)
    Ag 1  0 0 0  0.002
  ]
]

Modes[
]

Correlations
[
  [(0,0,0)
   Multiplicity 1
   SubstitutionalCorrelation(AuAg,AuAg,1/2)
  ]
  [(1,0,0)
   Multiplicity 4
   SubstitutionalCorrelation(AuAg,AuAg, p1)
  ]
  [(1,1,0)
   Multiplicity 4
   SubstitutionalCorrelation(AuAg,AuAg, p2)
  ]
]
"""

TARGET = dict(scale=132, p1=0.113, p2=0.298)

def run_yell(yell, work):
    r = subprocess.run([yell], cwd=work, capture_output=True, text=True)
    if r.returncode != 0:
        print("  yell failed:\n" + r.stdout[-2000:] + "\n" + r.stderr[-2000:])
        raise RuntimeError("yell failed")
    return r.stdout

def q_length(h5path):
    """|q| per voxel via the reciprocal metric tensor, matching Grid::d_star_square_at
    for this orthorhombic (all-90deg) cell: d*^2 = (h/a)^2 + (k/b)^2 + (l/c)^2,
    (h,k,l) = lower_limits + step_sizes * index, in the same C-order as data."""
    with h5py.File(h5path, "r") as f:
        a, b, c = f["unit_cell"][0:3]
        low = np.array(f["lower_limits"][...], float)
        step = np.array(f["step_sizes"][...], float)
        shape = f["data"].shape
    ii, jj, kk = np.meshgrid(np.arange(shape[0]), np.arange(shape[1]),
                             np.arange(shape[2]), indexing="ij")
    h = low[0] + step[0] * ii
    k = low[1] + step[1] * jj
    l = low[2] + step[2] * kk
    dstar2 = (h / a) ** 2 + (k / b) ** 2 + (l / c) ** 2
    return np.sqrt(np.maximum(dstar2, 0.0))

def parse_params(path):
    with open(path) as f:
        txt = f.read()
    params = {}
    m = re.search(r'Scale\s+([0-9.eE+-]+)(?:\([0-9]+\))?', txt)
    if m: params['Scale'] = float(m.group(1))
    blk = re.search(r'RefinableVariables\s*\[(.*?)\]', txt, re.DOTALL)
    if blk:
        for name, val in re.findall(r'(\w+)\s*=\s*([0-9.eE+-]+)(?:\([0-9]+\))?\s*;', blk.group(1)):
            params[name] = float(val)
    for name, val in re.findall(r'#\s*(Bkg_\d+)\s*=\s*([0-9.eE+-]+)', txt):
        params[name] = float(val)
    return params

def run_case(name, yell, degree, inject):
    """inject: callable(E0_array, qlen_array) -> (E1_array, expected_bkg_dict)."""
    work = os.path.join(BASE, f"work_background_{name}")
    if os.path.exists(work): shutil.rmtree(work)
    os.makedirs(work)
    print(f"Testing background: {name}")

    # 1. structural experiment
    with open(os.path.join(work, "model.txt"), "w") as f:
        f.write(model_text("false", **TARGET))
    run_yell(yell, work)
    shutil.move(os.path.join(work, "model.h5"), os.path.join(work, "experiment.h5"))

    # 2. inject background
    qlen = q_length(os.path.join(work, "experiment.h5"))
    with h5py.File(os.path.join(work, "experiment.h5"), "r+") as f:
        E0 = f["data"][...]
        E1, expected = inject(E0, qlen)
        f["data"][...] = E1

    # 3. refine with background, structural params perturbed to 0 / Scale 1
    extra = f"RefineBackground true\nBackgroundDegree {degree}"
    with open(os.path.join(work, "model.txt"), "w") as f:
        f.write(model_text("true", scale=1, p1=0, p2=0, extra=extra))
    run_yell(yell, work)

    # 4. compare
    refined = parse_params(os.path.join(work, "refined_parameters.txt"))
    ok = True
    checks = {"Scale": TARGET["scale"], "p1": TARGET["p1"], "p2": TARGET["p2"], **expected}
    # structural params must be tight; background coeff tolerance a touch looser
    for key, want in checks.items():
        got = refined.get(key)
        tol = 5e-3 if key.startswith("Bkg") else 1e-3
        if got is None:
            print(f"    {key}: MISSING (expected {want})"); ok = False; continue
        diff = abs(got - want)
        status = "OK" if diff <= tol else "FAIL"
        if diff > tol: ok = False
        print(f"    {key}: target={want:g}, refined={got:g}, diff={diff:g} [{status}]")
    print("  SUCCESS" if ok else "  FAILURE")
    return ok

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin", default=default_bin())
    args = ap.parse_args()
    if not os.path.exists(args.bin):
        print(f"ERROR: yell binary not found at {args.bin}"); sys.exit(1)

    def const_inject(E0, q):
        C = 50.0
        E1 = E0 + C
        return E1, {"Bkg_0": C}

    def linear_inject(E0, q):
        # B = A + G*|q|; degree-1 Chebyshev spans {1,|q|} so it is absorbed exactly.
        # Coefficients depend on the internal normalisation, so we only assert
        # structural recovery here (no Bkg_* expectation).
        A, G = 200.0, -8.0
        E1 = E0 + A + G * q
        return E1, {}

    results = []
    results.append(run_case("constant", args.bin, 0, const_inject))
    results.append(run_case("linear",   args.bin, 1, linear_inject))

    if all(results):
        print("\nAll background integration tests PASSED!")
    else:
        print("\nBackground integration tests FAILED!"); sys.exit(1)

if __name__ == "__main__":
    main()
