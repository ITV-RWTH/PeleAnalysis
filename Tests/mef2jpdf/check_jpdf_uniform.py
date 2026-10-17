#!/usr/bin/env python3
"""
check_jpdf_uniform.py <matrix_dat> [tol_abs]

Reads a mef2jpdf output matrix and asserts that every bin value equals
1/N^2 (where N = number of bins per axis) within absolute tolerance <tol_abs>
(default 1e-10).

Exits 0 on success, 1 on failure.
"""
import sys
import numpy as np


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: check_jpdf_uniform.py <matrix_dat> [tol_abs]")

    path = sys.argv[1]
    tol = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-10

    try:
        data = np.loadtxt(path)
    except Exception as e:
        print(f"Cannot read {path}: {e}")
        sys.exit(1)

    nBins = data.shape[0]
    expected = 1.0 / (nBins * nBins)
    max_err = float(np.max(np.abs(data - expected)))

    if max_err <= tol:
        print(f"all {nBins}x{nBins} bins = {expected:.8f} (max deviation {max_err:.2e}) -- OK")
        sys.exit(0)
    else:
        print(f"max |bin - {expected:.8f}| = {max_err:.2e} > tol {tol:.2e} -- FAIL")
        sys.exit(1)


if __name__ == "__main__":
    main()
