#!/usr/bin/env python3
"""
check_plt_bounds.py <plt_dir> <var_name> <lower> <upper>

Reads every cell of <var_name> across all AMR levels of an AMReX plotfile and
asserts that each value is finite and lies in [lower, upper].  Unlike
check_plt_range.py in Tests/generateTestPlt, NaN or Inf is always a failure.
Exits 0 on success, 1 on failure.
"""
import math
import os
import re
import struct
import sys


def parse_fab_hdr(line):
    little = "8 7 6 5 4 3 2 1" in line
    triples = re.findall(r"\((-?\d+),(-?\d+),(-?\d+)\)", line)
    if len(triples) < 2:
        raise ValueError(f"Cannot parse box from FAB header: {line!r}")
    lo = tuple(int(x) for x in triples[0])
    hi = tuple(int(x) for x in triples[1])
    return lo, hi, little


def main():
    if len(sys.argv) != 5:
        sys.exit("Usage: check_plt_bounds.py <plt_dir> <var_name> <lower> <upper>")

    plt_dir, var_name = sys.argv[1], sys.argv[2]
    lower, upper = float(sys.argv[3]), float(sys.argv[4])

    with open(os.path.join(plt_dir, "Header")) as f:
        hdr_lines = f.readlines()
    ncomp = int(hdr_lines[1])
    var_names = [hdr_lines[2 + i].strip() for i in range(ncomp)]
    if var_name not in var_names:
        sys.exit(f"Variable '{var_name}' not in plotfile.  Available: {var_names}")
    comp_idx = var_names.index(var_name)

    total = n_bad = 0
    obs_min, obs_max = math.inf, -math.inf
    lev = 0
    while os.path.isdir(os.path.join(plt_dir, f"Level_{lev}")):
        lev_dir = os.path.join(plt_dir, f"Level_{lev}")
        with open(os.path.join(lev_dir, "Cell_H")) as f:
            cell_h = f.read()
        for m in re.finditer(r"FabOnDisk:\s+(\S+)\s+(\d+)", cell_h):
            with open(os.path.join(lev_dir, m.group(1)), "rb") as f:
                f.seek(int(m.group(2)))
                lo, hi, little = parse_fab_hdr(f.readline().decode("latin-1"))
                n = math.prod(b - a + 1 for a, b in zip(lo, hi))
                f.seek(comp_idx * n * 8, 1)
                raw = f.read(n * 8)
                if len(raw) < n * 8:
                    sys.exit(f"Truncated FAB data in {lev_dir}/{m.group(1)}")
                vals = struct.unpack(("<" if little else ">") + "d" * n, raw)
            for v in vals:
                if not (math.isfinite(v) and lower <= v <= upper):
                    n_bad += 1
                if math.isfinite(v):
                    obs_min, obs_max = min(obs_min, v), max(obs_max, v)
            total += n
        lev += 1

    if total == 0:
        sys.exit(f"No cells found in {plt_dir}")
    print(
        f"{total} cells, min={obs_min:.6g}, max={obs_max:.6g}, "
        f"{n_bad} outside [{lower}, {upper}] or non-finite"
    )
    sys.exit(0 if n_bad == 0 else 1)


if __name__ == "__main__":
    main()
