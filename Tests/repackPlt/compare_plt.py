#!/usr/bin/env python3
"""Compare the cell data of two AMReX plotfiles box by box (exact match).

Usage: compare_plt.py PLT_A PLT_B [--tol TOL] [--max-files N]

Reads every Level_*/<mf>_H MultiFab header and its Cell_D data files
directly (VisMF formats v1..v4), so it does not depend on any AMReX tool.
Exits 0 when all levels, MultiFabs, boxes and values match, 1 otherwise.
--max-files N additionally asserts that PLT_B has at most N data files per
MultiFab.
"""

import argparse
import glob
import os
import re
import sys

import numpy as np

BOX_RE = re.compile(r"\(\(([-\d,\s]+)\)\s*\(([-\d,\s]+)\)\s*\(([-\d,\s]+)\)\)")
RD_RE = re.compile(r"\(\((\d+),\s*\(([\d\s]+)\)\),\((\d+),\s*\(([\d\s]+)\)\)\)")


def ints(s):
    return [int(v) for v in s.replace(",", " ").split()]


def dtype_from_rd(nbytes, order):
    # order "1 2 ... n" = big endian, "n ... 1" = little endian
    endian = ">" if order[0] == 1 else "<"
    return np.dtype(f"{endian}f{nbytes}")


def read_mf(path):
    """Return (ncomp, list of (lo, hi, ndarray)) for MultiFab <path>."""
    with open(path + "_H") as f:
        text = f.read()
    lines = text.splitlines()
    vers = int(lines[0])
    ncomp = int(lines[2])
    ngrow = ints(lines[3].strip("()")) if "(" in lines[3] else [int(lines[3])]
    fods = re.findall(r"FabOnDisk:\s+(\S+)\s+(\d+)", text)
    rd = None
    if vers != 1:
        m = RD_RE.findall(text)
        nbytes, order = int(m[-1][2]), ints(m[-1][3])
        rd = dtype_from_rd(nbytes, order)
    # BoxArray boxes appear before the FabOnDisk list
    ba_text = text.split("FabOnDisk:")[0]
    boxes = [(ints(a), ints(b)) for a, b, _ in BOX_RE.findall(ba_text)]
    fabs = []
    base = os.path.dirname(path)
    for i, (fname, off) in enumerate(fods):
        with open(os.path.join(base, fname), "rb") as df:
            df.seek(int(off))
            if vers == 1:
                hdr = df.readline().decode()
                m = RD_RE.search(hdr)
                dt = dtype_from_rd(int(m.group(3)), ints(m.group(4)))
                b = BOX_RE.search(hdr)
                lo, hi = ints(b.group(1)), ints(b.group(2))
                nc = int(hdr.split(")")[-1])
            else:
                dt = rd
                lo = [l - g for l, g in zip(boxes[i][0], ngrow * 3)]
                hi = [h + g for h, g in zip(boxes[i][1], ngrow * 3)]
                nc = ncomp
            npts = int(np.prod([h - l + 1 for l, h in zip(lo, hi)]))
            data = np.fromfile(df, dtype=dt, count=npts * nc)
        fabs.append((tuple(lo), tuple(hi), data))
    return ncomp, fabs


def mf_names(plt):
    out = []
    for h in sorted(glob.glob(os.path.join(plt, "Level_*", "*_H"))):
        out.append(os.path.relpath(h, plt)[:-2])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--tol", type=float, default=0.0)
    ap.add_argument("--max-files", type=int, default=0)
    args = ap.parse_args()

    ok = True
    na, nb = mf_names(args.a), mf_names(args.b)
    if na != nb:
        print(f"FAIL: MultiFab lists differ: {na} vs {nb}")
        return 1
    for mf in na:
        ca, fa = read_mf(os.path.join(args.a, mf))
        cb, fb = read_mf(os.path.join(args.b, mf))
        if ca != cb or len(fa) != len(fb):
            print(
                f"FAIL: {mf}: ncomp/nboxes differ ({ca},{len(fa)}) vs ({cb},{len(fb)})"
            )
            ok = False
            continue
        maxerr = 0.0
        for (loa, hia, da), (lob, hib, db) in zip(fa, fb):
            if loa != lob or hia != hib or da.size != db.size:
                print(f"FAIL: {mf}: box mismatch {loa}-{hia} vs {lob}-{hib}")
                ok = False
                break
            if da.size:
                maxerr = max(
                    maxerr,
                    float(
                        np.max(np.abs(da.astype(np.float64) - db.astype(np.float64)))
                    ),
                )
        if maxerr > args.tol:
            print(f"FAIL: {mf}: max abs difference {maxerr:g} > {args.tol:g}")
            ok = False
        if args.max_files:
            nfiles = len(glob.glob(os.path.join(args.b, mf + "_D_*")))
            if nfiles > args.max_files:
                print(f"FAIL: {mf}: {nfiles} data files > {args.max_files}")
                ok = False
        print(f"{mf}: {len(fa)} boxes, {ca} comps, max diff {maxerr:g}")
    print("PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
