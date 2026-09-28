# sequentialFilterPlt Test Suite

Tests for `Src/sequentialFilterPlt.cpp` — the box filter applied as three 1D passes — and for the fixes it shares with `Src/filterPlt.cpp`.

## Running

```bash
cd PeleAnalysis/Tests/sequentialFilterPlt
./run_tests.sh              # build + run
./run_tests.sh --no-compile # skip build, reuse existing executables
```

Build logs are written next to the script; plotfiles are written to `testrun/`. Neither is tracked by git. `testrun/` is wiped at the start of each run so that stale output cannot make a failed run pass.

## Oracles

No part of the tool is reimplemented for testing. Results are checked against:

- **`filterPlt`**, an independent implementation of the stencil (PelePhysics `Filter`, one 3D pass). A box filter is separable, so both tools must agree to round-off wherever no coarse/fine interface is involved. The two tools share the plotfile reading, the per-level filter width and the ghost-cell filling, so those are checked separately (Phases 4–6).
- **AMReX `fcompare` and `fnan`** (built from `$AMREX_HOME/Tools/Plotfile`; `$AMREX_HOME` defaults to the PelePhysics submodule). `fcompare` compares plotfiles cell by cell and fails on NaN; the suite also makes it fail when a variable is missing from either file.
- **Invariants**: the filter weights are positive and sum to one, and the limited conservative interpolation adds no new extrema, so a constant field is preserved and a field bounded by [a, b] stays in [a, b]. The result must not depend on the box layout, the thread count, the rank count or how memory was initialised, and filtering a level must not depend on finer levels.
- **The input plotfile Header**, for time, domain, level count, ref ratios and level domains of the output.

`check_plt_bounds.py` reads plotfiles directly; unlike `Tests/generateTestPlt/check_plt_range.py` it treats NaN and Inf as failures.

**Known limit:** levels 0–1 are checked exactly, but on the tightly nested level 2, where ghost cells come from level 0, values are only checked for NaN, bounds, the constant field and independence of memory, box layout and thread count. An exact check there would need a field that interpolation reproduces exactly (a linear one), which `generateTestPlt` does not provide.

## Test plotfiles

Generated with `generateTestPlt`. All contain `sine_f` (range [-1, 1]), `sph` (smooth sphere, range [0, 1]) and `cst` (constant 3.5).

| Input | Layout |
|-------|--------|
| `gen_single.inp` | 64³, one level, non-periodic |
| `gen_single_periodic.inp` | 64³, one level, periodic |
| `gen_full_rr2.inp` / `gen_full_rr4.inp` | 16³ + a level 1 covering the whole domain, ref_ratio 2 / 4 |
| `gen_nested_rr4.inp` | 32³ + a central level 1, ref_ratio 4 |
| `gen_tight_3lev.inp` | 32³, three levels; level 2 is two level-2 cells inside level 1 |

## Test matrix

Tolerances: 1e-12 between the two tools (their summation orders differ); bit-identical (0) where the tool is compared with itself.

### Phase 3 — Agreement with filterPlt

| Test | What is checked |
|------|-----------------|
| T1–T4 | Single level, `base_fgr` 2, 8, 16 and periodic 8: matches filterPlt everywhere, including at domain boundaries |
| T5–T7 | Level 1 covers the domain, so there is no coarse/fine interface: matches filterPlt on every level for ref_ratio 2 and 4 and for `same_fgr_all_levels = true` |
| T8 | The output differs from the input |

### Phase 4 — Filter width per level

| Test | What is checked |
|------|-----------------|
| T9 | Level 1 of `base_fgr=4` (scaled by ref_ratio 2) is bit-identical to level 1 of `base_fgr=8, same_fgr_all_levels=true` |
| T10 | `same_fgr_all_levels` changes level 1 |

### Phase 5 — Wide filters on tightly nested levels

For `base_fgr` 2 and 8 the level-2 filter half-width (4 and 16 cells) exceeds the two-cell nesting margin, so ghost cells must come from level 0. Runs use `amrex.init_snan=1`.

| Test | What is checked |
|------|-----------------|
| T11 | No NaN on any level |
| T12 | Bit-identical with and without `init_snan` |
| T13 | Levels 0–1 are bit-identical to a `max_filter_level=1` run |
| T14 | `cst` stays 3.5 |
| T15–T16 | `sph` stays in [0, 1], `sine_f` in [-1, 1] |
| T17–T18 | filterPlt: no NaN, `sph` stays in [0, 1] |

### Phase 6 — Output metadata and options

| Test | What is checked |
|------|-----------------|
| T19–T21 | Time, domain, finest level, ref ratios and level domains in the output Header equal the input's (ref_ratio 4 for both tools, and the 3-level file) |
| T22 | Default output name is `<infile>_filtered` in the working directory |
| T23 | `max_filter_level=1` writes levels 0–1 only |
| T24–T25 | `variables=sph` and `variables="cst sph"` write exactly those variables, bit-identical to a full run |
| T26 | `max_grid_size=8` with a 32-cell filter matches the input grids, single level |
| T27 | `max_grid_size=8` matches the input grids, 3 levels |
| T28 | `interp_type=0` runs without NaN and stays bounded |

### Phase 7 — Error handling

| Test | What is checked |
|------|-----------------|
| T29 | `filter_type=2` aborts and writes no output |
| T30–T31 | `base_fgr` odd or zero aborts |
| T32 | Unknown variable aborts |
| T33–T34 | Missing or nonexistent `infile` aborts |

### Phase 8 — OpenMP

| Test | What is checked |
|------|-----------------|
| T35 | sequentialFilterPlt with 4 threads is bit-identical to serial (3 levels) |
| T36 | filterPlt with 4 threads is bit-identical to serial |

### Phase 9 — MPI (skipped if mpirun/mpiexec is not found)

| Test | What is checked |
|------|-----------------|
| T37 | 2 and 4 ranks are bit-identical to serial, single level |
| T38 | 2 and 4 ranks are bit-identical to serial, 3 levels with `max_grid_size=8` |
