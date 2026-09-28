# sequentialFilterPlt Test Suite

Tests for `Src/sequentialFilterPlt.cpp` — the box filter applied as three 1D passes — and for the fixes it shares with `Src/filterPlt.cpp`.

## Running

```bash
cd PeleAnalysis/Tests/sequentialFilterPlt
./run_tests.sh              # build + run
./run_tests.sh --no-compile # skip build, reuse existing executables
```

Output artefacts (plotfiles, build logs) are written to `testrun/` and are not tracked by git. `testrun/` is wiped at the start of each run so that stale output cannot make a failed run pass.

## Oracles

No part of the tool is reimplemented for testing. Results are checked against:

- **`filterPlt`**, an independent implementation (PelePhysics `Filter`, one 3D stencil) of the same box filter. A box filter is separable, so both tools must agree to round-off wherever no coarse/fine interface is involved.
- **AMReX `fcompare`** (built from `$AMREX_HOME/Tools/Plotfile`), which compares two plotfiles cell by cell and fails on any NaN. `$AMREX_HOME` defaults to the PelePhysics submodule.
- **Invariants of the filter**: the weights are positive and sum to one, so a constant field is preserved and a field bounded by [a, b] stays in [a, b]. The result must not depend on the box layout, the thread or rank count, or on how memory was initialised.
- **The input plotfile Header**, for the AMR layout of the output.

`check_plt_bounds.py` reads plotfiles directly; unlike `Tests/generateTestPlt/check_plt_range.py` it treats NaN and Inf as failures.

## Test plotfiles

Generated with `generateTestPlt`. All contain `sine_f` (range [-1, 1]), `sph` (smooth sphere, range [0, 1]) and `cst` (constant 3.5).

| Input | Layout |
|-------|--------|
| `gen_single.inp` | 64³, one level, non-periodic |
| `gen_single_periodic.inp` | 64³, one level, periodic |
| `gen_full_rr2.inp` / `gen_full_rr4.inp` | 16³ + a level 1 covering the whole domain, ref_ratio 2 / 4 |
| `gen_nested_rr4.inp` | 32³ + a central level 1, ref_ratio 4 |
| `gen_tight_3lev.inp` | 32³, three levels, level 2 only two cells inside level 1 |

## Test matrix

### Phase 3 — Agreement with filterPlt

| Test | What is checked |
|------|-----------------|
| T1–T4 | Single level, `base_fgr` 2, 8, 16 and periodic 8: matches filterPlt to 1e-12 everywhere, including at domain boundaries |
| T5–T7 | Level 1 covers the domain, so there is no coarse/fine interface: matches filterPlt on every level for ref_ratio 2 and 4 and for `same_fgr_all_levels = true` |
| T8 | The output differs from the input (the filter is not a no-op) |

### Phase 4 — Wide filters on tightly nested levels

For `base_fgr` 2 and 8 the level-2 filter half-width (4 and 16 cells) exceeds the two-cell nesting margin, so ghost cells must come from level 0. Runs use `amrex.init_snan=1`.

| Test | What is checked |
|------|-----------------|
| T9 | No NaN on any level |
| T10 | Output identical with and without `init_snan` |
| T11 | `cst` stays 3.5 |
| T12–T13 | `sph` stays in [0, 1], `sine_f` in [-1, 1] |
| T14–T15 | The same for filterPlt |

### Phase 5 — Output metadata and options

| Test | What is checked |
|------|-----------------|
| T16–T18 | Finest level, ref ratios and level domains in the output Header equal the input's (ref_ratio 4 for both tools, and the 3-level file) |
| T19 | Default output name is `<infile>_filtered` in the working directory |
| T20 | `max_filter_level=1` writes levels 0–1 only |
| T21 | `variables=sph` gives the same `sph` as a full run and writes one variable |
| T22–T23 | `max_grid_size=8` (smaller than the filter width) gives the same result as the input grids, single level and 3 levels |
| T24 | `interp_type=0` runs without NaN and stays bounded |

### Phase 6 — Error handling

| Test | What is checked |
|------|-----------------|
| T25 | `filter_type=2` aborts and writes no output |
| T26–T27 | `base_fgr` odd or zero aborts |
| T28 | Unknown variable aborts |
| T29–T30 | Missing or nonexistent `infile` aborts |

### Phase 7 — OpenMP

| Test | What is checked |
|------|-----------------|
| T31 | sequentialFilterPlt with 4 threads matches the serial build (3 levels) |
| T32 | filterPlt with 4 threads matches the serial build |

### Phase 8 — MPI (skipped if mpirun/mpiexec is not found)

| Test | What is checked |
|------|-----------------|
| T33 | 2 and 4 ranks match serial, single level |
| T34 | 2 and 4 ranks match serial, 3 levels with `max_grid_size=8` |
