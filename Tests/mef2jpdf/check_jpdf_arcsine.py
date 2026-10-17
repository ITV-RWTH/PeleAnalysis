#!/usr/bin/env python3
"""
check_jpdf_arcsine.py <matrix_dat> <centers_dat> <center> <radius> [tol_frac]

Checks that the X-marginal (sum of each column over all rows) of a mef2jpdf
output matrix follows the arcsine distribution for a cylinder of given <center>
and <radius>:

    p(x) = 1 / (pi * sqrt(R^2 - (x - c)^2))   for |x - c| < R

<centers_dat>  path to the x bin-centers file (Pdf_<xvar_out>_x.dat)
<center>       cylinder axis center in the x direction (float)
<radius>       cylinder radius (float)
<tol_frac>     max allowed relative error per bin (default 0.20 = 20%)

Also verifies left-right symmetry: marginal[i] ~ marginal[N-1-i].

Exits 0 on success, 1 on failure.
"""
import sys
import numpy as np


def arcsine_weight(x_lo, x_hi, c, R):
    """Integral of arcsine density over [x_lo, x_hi] for cylinder c, R."""
    u_lo = np.clip((x_lo - c) / R, -1.0, 1.0)
    u_hi = np.clip((x_hi - c) / R, -1.0, 1.0)
    return (np.arcsin(u_hi) - np.arcsin(u_lo)) / np.pi


def main():
    if len(sys.argv) < 5:
        sys.exit(
            "Usage: check_jpdf_arcsine.py <matrix_dat> <centers_dat> "
            "<center> <radius> [tol_frac]"
        )

    matrix_path = sys.argv[1]
    centers_path = sys.argv[2]
    c = float(sys.argv[3])
    R = float(sys.argv[4])
    tol_frac = float(sys.argv[5]) if len(sys.argv) > 5 else 0.20

    try:
        data = np.loadtxt(matrix_path)
        centers = np.loadtxt(centers_path)
    except Exception as e:
        print(f"Cannot read input: {e}")
        sys.exit(1)

    nBins = len(centers)
    half_w = (centers[-1] - centers[0]) / (2 * (nBins - 1))
    edges = np.concatenate([[centers[0] - half_w],
                             centers[:-1] + np.diff(centers) / 2,
                             [centers[-1] + half_w]])

    # xvar-marginal: sum over yvar columns (axis=1) → one value per xvar bin (row)
    marginal = data.sum(axis=1)

    # Expected arcsine weights
    expected = np.array([arcsine_weight(edges[i], edges[i + 1], c, R)
                         for i in range(nBins)])
    total_exp = expected.sum()
    if total_exp < 0.5:
        print(f"Expected arcsine weights sum to {total_exp:.3f} -- bins likely outside cylinder range")
        sys.exit(1)
    expected /= total_exp  # normalise in case boundary clips

    # Relative error check (skip bins with expected weight < 1e-6)
    active = expected > 1e-6
    rel_err = np.where(active, np.abs(marginal - expected) / expected, 0.0)
    max_rel = float(rel_err.max())

    # Symmetry check: |marginal[i] - marginal[N-1-i]| relative to mean
    sym_err = np.abs(marginal - marginal[::-1]) / (marginal + marginal[::-1] + 1e-15) * 2
    max_sym = float(sym_err[:nBins // 2].max())

    ok = True
    if max_rel > tol_frac:
        worst = int(rel_err.argmax())
        print(
            f"Arcsine check FAIL: bin {worst} marginal={marginal[worst]:.5f} "
            f"expected={expected[worst]:.5f} rel_err={max_rel:.3f} > {tol_frac}"
        )
        ok = False
    else:
        print(f"Arcsine shape OK: max rel error = {max_rel:.3f} <= {tol_frac}")

    if max_sym > tol_frac:
        print(f"Symmetry check FAIL: max asymmetry = {max_sym:.3f} > {tol_frac}")
        ok = False
    else:
        print(f"Symmetry OK: max asymmetry = {max_sym:.3f}")

    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
