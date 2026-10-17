#!/usr/bin/env python3
"""
check_jpdf_marginal_uniform.py <matrix_dat> <axis> [tol_frac]

Reads a mef2jpdf output matrix and asserts that the marginal distribution along
<axis> is uniform.

  axis=0  sum over rows    -> checks that each COLUMN sum = 1/N
  axis=1  sum over columns -> checks that each ROW sum = 1/N

<tol_frac> is the allowed relative deviation per bin (default 0.10 = 10%).

Exits 0 on success, 1 on failure.
"""
import sys
import numpy as np


def main():
    if len(sys.argv) < 3:
        sys.exit("Usage: check_jpdf_marginal_uniform.py <matrix_dat> <axis> [tol_frac]")

    path = sys.argv[1]
    axis = int(sys.argv[2])
    tol_frac = float(sys.argv[3]) if len(sys.argv) > 3 else 0.10

    if axis not in (0, 1):
        sys.exit("axis must be 0 or 1")

    try:
        data = np.loadtxt(path)
    except Exception as e:
        print(f"Cannot read {path}: {e}")
        sys.exit(1)

    nBins = data.shape[0]
    marginal = data.sum(axis=axis)
    expected = 1.0 / nBins
    max_rel = float(np.max(np.abs(marginal - expected) / expected))

    if max_rel <= tol_frac:
        print(f"marginal (axis={axis}) uniform: max rel deviation = {max_rel:.3f} <= {tol_frac} -- OK")
        sys.exit(0)
    else:
        worst = int(np.argmax(np.abs(marginal - expected)))
        print(
            f"marginal (axis={axis}) NOT uniform: bin {worst} = {marginal[worst]:.6f}, "
            f"expected {expected:.6f}, rel error {max_rel:.3f} > {tol_frac} -- FAIL"
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
