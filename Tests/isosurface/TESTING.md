# isosurface Test Suite

Tests for `Src/isosurface.cpp` on periodic domains: iso lines (2D) and iso
surfaces (3D) that cross periodic boundaries, on uniform and two-level grids.

## Running

```bash
cd PeleAnalysis/Tests/isosurface
./run_tests.sh              # build + run
./run_tests.sh --no-compile # skip the build, reuse the executables in Src/
```

Artefacts are written to `testrun/`. `python3` is required for `check_mef.py`.
`ISO2D` and `ISO3D` select other isosurface executables.

## Expected values

All tests use the iso line `c = 0.5` of `c = cos(2 pi x) cos(2 pi y)`, which has
two congruent closed lines per unit square.

- **Reference:** `gen_ref_2d.inp` holds one line alone in the non-periodic window
  `[0.25, 0.75]^2` (`dx = 1/64`). Its measured length `L_ref` is within 0.2% of
  the exact length 1.0810134.
- **Periodic tests:** the field is shifted by 10 cells. The cell centres
  coincide with those of the reference, and one line reaches less than a cell
  across both periodic boundaries. The measure is `2 L_ref` in 2D and
  `0.5 L_ref` in 3D (extruded over `z` in `[0, 0.25]`), to a relative `1e-9`.
- **Two-level variants:** `dx = 1/32` on level 0, and level 1 covers both lines.

`check_mef.py` measures the MEF in the minimum image of the periodic directions,
and also fails if a node lies outside the domain or an edge is longer than 1.5
finest cell diagonals.

## Test matrix

| # | Input | What is checked |
|---|---|---|
| Ref | `gen_ref_2d.inp` | Within 0.2% of the exact length; nodes inside; no long edge |
| per_2d | `gen_per_2d.inp` | `2 L_ref` (MEF and `computeArea`); nodes inside the domain; no long edge |
| per_amr_2d | `gen_per_amr_2d.inp` | As per_2d, two levels |
| per_3d | `gen_per_3d.inp` | `0.5 L_ref` (MEF and `computeArea`); nodes inside the domain; no long edge |
| per_amr_3d | `gen_per_amr_3d.inp` | As per_3d, two levels |
