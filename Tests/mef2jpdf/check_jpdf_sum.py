#!/usr/bin/env python3
"""
check_jpdf_sum.py <matrix_dat> <expected> [tol]

Reads a mef2jpdf output matrix (space-separated floats, one row per line) and
asserts that the sum of all values equals <expected> within <tol> (default 1e-6).

Exits 0 on success, 1 on failure.
"""
import sys
import numpy as np


def main():
    if len(sys.argv) < 3:
        sys.exit("Usage: check_jpdf_sum.py <matrix_dat> <expected> [tol]")

    path = sys.argv[1]
    expected = float(sys.argv[2])
    tol = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-6

    try:
        data = np.loadtxt(path)
    except Exception as e:
        print(f"Cannot read {path}: {e}")
        sys.exit(1)

    total = float(data.sum())
    if abs(total - expected) <= tol:
        print(f"sum = {total:.10f} (expected {expected}, tol {tol}) -- OK")
        sys.exit(0)
    else:
        print(f"sum = {total:.10f}, expected {expected} +/- {tol} -- FAIL")
        sys.exit(1)


if __name__ == "__main__":
    main()
