# mef2jpdf test suite

## Pipeline

```text
generateTestPlt3d.gnu.ex  →  AMReX plotfile
isosurface3d.gnu.ex       →  .mef surface file
mef2jpdf3d.gnu.ex         →  MEF_JPDFAverage_<X>_<Y>/ output directory
Python checkers           →  pass / fail
```

Run the full suite from this directory:

```bash
./run_tests.sh
# or, if the executables are already built:
./run_tests.sh --no-compile
```

All artefacts are written to `testrun/`.

---

## Test matrix

| ID | Geometry | xvar | yvar | nBins | Analytical result | Checker |
| -- | -------- | ---- | ---- | ----- | ----------------- | ------- |
| T1 | Sphere | X | Y | 20 | sum = 1.0 | `check_jpdf_sum.py` |
| T2 | Flat plane (z=0.5) | X | Y | 8 | both marginals uniform (15% tol) | `check_jpdf_sum.py`, `check_jpdf_marginal_uniform.py` |
| T3 | Sphere | X | Z | 20 | Z-marginal uniform (20% tol, belt theorem) | `check_jpdf_sum.py`, `check_jpdf_marginal_uniform.py` |
| T4 | Sphere + const\_field=3.14 | const\_field | Z | 10 | one column holds all weight; Z-marginal of that column uniform | inline Python |
| T5 | Sphere MEF fed twice | X | Y | 20 | sum = 1.0 (doubling cancels on normalisation) | `check_jpdf_sum.py` |
| T6 | Flat plane + sine\_x=sin(2πX) | sine\_x | Y | 20 | sine\_x-marginal follows arcsine (15% tol) | `check_jpdf_sum.py`, `check_jpdf_arcsine.py` |
| T7 | Sphere | X | Y | 20 | outer bin > center bin; X-marginal symmetric | `check_jpdf_sum.py`, inline Python |

---

## Analytical foundations

### T2 — Uniform marginals (flat plane)

The marching-cubes triangulation of an axis-aligned plane on a uniform non-periodic
Cartesian grid produces approximately equal-area triangles, but not exactly so:
boundary cells and the specific marching-cubes diagonal orientation introduce ~10%
variation in the marginals.  Both the X-marginal and Y-marginal are checked for
uniformity within 15% relative tolerance.

### T3 — Uniform Z-marginal (Archimedes belt theorem)

The lateral surface area of a spherical zone of height dz equals 2πR dz, independent of
the latitude.  Consequently, every Z bin receives exactly 1/nBins of the total area.
The discrete approximation on a 32³ grid introduces ~16% discretization error; the check
allows 20%.

### T4 — Constant mapped field

`const_field = 3.14` places every triangle in the same X bin (ix = 1 for the range
[3.0, 4.0) with nBins = 10).  All other columns must be zero.  Within column ix = 1,
the belt theorem ensures the Z distribution is uniform (20% tolerance).

Providing explicit `xmin = 3.0` / `xmax = 4.0` skips the first pass over the data that
`mef2jpdf` would otherwise use to compute the range — exercising the explicit-bounds
code path.

### T6 — Arcsine marginal (sine field on flat plane)

On the flat plane, `sine_x = sin(2πX)` where X is uniform on [0, 1].
The distribution of `u = sin(2πX)` is the arcsine distribution:

```text
p(u) = 1 / (π √(1 − u²)),   u ∈ (−1, 1)
```

The expected weight in bin i with edges [u\_lo, u\_hi] is:

```text
w_i = [arcsin(u_hi) − arcsin(u_lo)] / π
```

The formula `sin(arg_x) * cos(arg_y) * sin(arg_z)` gives zero when `fz=0` and
`phase_z=0` (since `sin(0)=0`). Setting `phase_z=π/2` makes the z-factor
`sin(π/2)=1`, so `sine_x = sin(2πX)` as intended.

The check allows 15% relative error per bin.  A 256³ grid is needed to keep the
discretization error below this threshold: with fewer cells, adjacent triangle
sine\_x values are too coarsely spaced to approximate the arcsine density well in
the middle bins.

### T7 — Non-uniform sphere jPDF

The sphere (X, Y) jPDF is not uniform: surface area concentrates near the equatorial
ring (large |Z|), so the outer X bins at mid-Y have higher density than the central
bins.  The test asserts:

- `P[0, mid] > P[mid, mid]` (outer X bin at center Y exceeds center bin)
- X-marginal is symmetric about X = 0.5 (sphere symmetry), max relative asymmetry < 5%

---

## Python helper scripts

| Script | Usage | Description |
| ------ | ----- | ----------- |
| `check_jpdf_sum.py` | `<matrix.dat> <expected> [tol]` | Asserts `sum(matrix) ≈ expected` within `tol` (default 1e-6) |
| `check_jpdf_uniform.py` | `<matrix.dat> [tol_abs]` | Asserts every bin = 1/N² within absolute tolerance (default 1e-10) |
| `check_jpdf_marginal_uniform.py` | `<matrix.dat> <axis> [tol_frac]` | Asserts marginal along `axis` (0=column sums, 1=row sums) is uniform within relative tolerance (default 10%) |
| `check_jpdf_arcsine.py` | `<matrix.dat> <centers.dat> <center> <radius> [tol_frac]` | Asserts xvar-marginal (row sums) matches arcsine distribution; checks left-right symmetry |

All scripts exit 0 on success and 1 on failure, printing a one-line summary.

---

## Input files

### Plotfile generation (`generateTestPlt`)

| File | Geometry | Fields | n\_cell |
| ---- | -------- | ------ | ------- |
| `gen_sphere.inp` | Sphere step at (0.5, 0.5, 0.5) R=0.4 | `sphere_step` | 32³ |
| `gen_plane_z.inp` | Plane step on z at 0.5 | `z_iso` | 32³ |
| `gen_sphere_const.inp` | Sphere step + constant field 3.14 | `sphere_step`, `const_field` | 32³ |
| `gen_cylinder.inp` | Flat plane at z=0.5 + sine\_x=sin(2πX) | `z_iso`, `sine_x` | 256³ |

### Isosurface extraction (`isosurface`)

| File | Plotfile | Iso component | Iso value | Extra comps |
| ---- | -------- | ------------- | --------- | ----------- |
| `iso_sphere.inp` | `plt_sphere` | `sphere_step` | 0.5 | — |
| `iso_plane_z.inp` | `plt_plane_z` | `z_iso` | 0.5 | — |
| `iso_sphere_const.inp` | `plt_sphere_const` | `sphere_step` | 0.5 | `const_field` |
| `iso_cylinder.inp` | `plt_cylinder` | `z_iso` | 0.5 | `sine_x` |

### jPDF computation (`mef2jpdf`)

| File | Input MEF | xvar | yvar | xvar\_out | yvar\_out |
| ---- | --------- | ---- | ---- | --------- | --------- |
| `jpdf_t1.inp` | `surf_sphere.mef` | X | Y | `X_t1` | `Y_t1` |
| `jpdf_t2.inp` | `surf_plane_z.mef` | X | Y | `X_t2` | `Y_t2` |
| `jpdf_t3.inp` | `surf_sphere.mef` | X | Z | `X_t3` | `Z_t3` |
| `jpdf_t4.inp` | `surf_sphere_const.mef` | `const_field` | Z | `cf_t4` | `Z_t4` |
| `jpdf_t5.inp` | `surf_sphere.mef` (×2) | X | Y | `X_t5` | `Y_t5` |
| `jpdf_t6.inp` | `surf_cylinder.mef` | `sine_x` | Y | `X_t6` | `Z_t6` |
| `jpdf_t7.inp` | `surf_sphere.mef` | X | Y | `X_t7` | `Y_t7` |

Unique `xvar_out`/`yvar_out` values prevent output-directory collisions when multiple
tests use the same variable pair.
