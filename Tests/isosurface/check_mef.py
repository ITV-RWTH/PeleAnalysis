#!/usr/bin/env python3
"""Measure an isosurface MEF file and check it for periodic defects.

Prints the total length (2D) or area (3D) of the surface, with every element
measured in the minimum image of the periodic directions. Fails (exit 1) if

- a node lies outside the domain,
- an element edge is longer than 1.5 finest cell diagonals (marching squares
  and cubes cut one cell, so no edge can be longer than a cell diagonal), or
- with --expect, the measure differs from the expected value by more than
  --rtol.

Usage:
  check_mef.py FILE --lo 0 0 --hi 1 1 --periodic 1 1 --dx 0.015625 \
               [--expect 1.2345 --rtol 1e-9]
"""

import argparse
import math
import re
import struct
import sys


def read_mef(path):
    """Return the node coordinates and the 0-based element connectivity."""
    with open(path, "rb") as f:
        f.readline()  # time
        names = f.readline().decode().split()
        n_elts, per_elt = map(int, f.readline().split())
        header = f.readline().decode()
        box = re.search(r"\(\(([-\d,]+)\) \(([-\d,]+)\) \([-\d,]+\)\) (\d+)", header)
        n_nodes = int(box.group(2).split(",")[0]) - int(box.group(1).split(",")[0]) + 1
        n_comp = int(box.group(3))
        data = struct.unpack(f"<{n_nodes * n_comp}d", f.read(8 * n_nodes * n_comp))
        conn = struct.unpack(f"<{n_elts * per_elt}i", f.read(4 * n_elts * per_elt))
    dim = sum(1 for n in names if n in ("X", "Y", "Z"))
    nodes = [data[i * n_comp : i * n_comp + dim] for i in range(n_nodes)]
    elts = [
        [c - 1 for c in conn[i * per_elt : (i + 1) * per_elt]] for i in range(n_elts)
    ]
    return nodes, elts


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("mef")
    parser.add_argument("--lo", type=float, nargs="+", required=True)
    parser.add_argument("--hi", type=float, nargs="+", required=True)
    parser.add_argument("--periodic", type=int, nargs="+", required=True)
    parser.add_argument("--dx", type=float, required=True, help="finest cell size")
    parser.add_argument("--expect", type=float)
    parser.add_argument("--rtol", type=float, default=1e-9)
    args = parser.parse_args()

    nodes, elts = read_mef(args.mef)
    dim = len(args.lo)
    length = [h - l for l, h in zip(args.lo, args.hi)]

    def diff(a, b):
        d = [b[k] - a[k] for k in range(dim)]
        for k in range(dim):
            if args.periodic[k]:
                d[k] -= length[k] * round(d[k] / length[k])
        return d

    tol = 1e-12 * max(length)
    outside = sum(
        1
        for p in nodes
        if any(p[k] < args.lo[k] - tol or p[k] > args.hi[k] + tol for k in range(dim))
    )

    measure, max_edge = 0.0, 0.0
    for e in elts:
        p = [nodes[i] for i in e]
        edges = [diff(p[i], p[(i + 1) % len(p)]) for i in range(len(p))]
        max_edge = max(max_edge, *(math.hypot(*d) for d in edges))
        if dim == 2:
            measure += math.hypot(*edges[0])
        else:
            a, b = edges[0], diff(p[0], p[2])
            cross = (
                a[1] * b[2] - a[2] * b[1],
                a[2] * b[0] - a[0] * b[2],
                a[0] * b[1] - a[1] * b[0],
            )
            measure += 0.5 * math.hypot(*cross)

    diagonal = args.dx * math.sqrt(dim)
    print(f"measure = {measure:.12g}")
    print(f"elements = {len(elts)}, nodes = {len(nodes)}, outside the domain = {outside}")
    print(f"longest edge = {max_edge / diagonal:.3g} cell diagonals")

    ok = outside == 0 and max_edge <= 1.5 * diagonal
    if args.expect is not None:
        rel = abs(measure - args.expect) / abs(args.expect)
        print(f"expected = {args.expect:.12g}, relative difference = {rel:.3g}")
        ok = ok and rel <= args.rtol
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
