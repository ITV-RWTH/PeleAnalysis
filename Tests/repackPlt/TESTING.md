# repackPlt Test Suite

Tests for `repackPlt.cpp` using a synthetic 3-level 3D plotfile (918 boxes,
two analytic fields) generated with `generateTestPlt` on 4 MPI ranks, so the
input has 4 `Cell_D` files per level. Data equality is checked box by box with
the included `compare_plt.py`, which reads the VisMF headers and data files
directly (header versions 1–4, 32/64 bit) and does not depend on any AMReX
tool.

---

## Running the test suite

```bash
cd Tests/repackPlt
./run_tests.sh               # compile + run all tests
./run_tests.sh --no-compile  # skip build, use existing executables in Src/
```

The script compiles `repackPlt` (3D serial and MPI) and `generateTestPlt`,
generates three input plotfiles (default, VisMF header version 4, 32-bit
data) and runs all tests. All artefacts land in `testrun/`. Without
`mpirun`/`mpiexec` the input is generated serially (1 file per level) and the
MPI tests are skipped. Python 3 with numpy is required.

---

## Tests

| # | Test | Pass criterion |
|---|------|----------------|
| 1 | Serial, default `n_files=1`, default `outfile` | `plt_repack_nf1` has 1 `Cell_D` file per level; all boxes and values identical; `Header`, an extra top-level dir and an extra file in `Level_0` copied byte for byte |
| 2 | Input with VisMF header version 4 | Output `Cell_H` is version 4; data identical |
| 3 | 32-bit input | Output FAB data is 32-bit; data identical |
| 4 | `float_type=32` on 64-bit input | Output is 32-bit; values within 1e-5 of input |
| 5 | `outfile` = `infile` | Tool aborts; input unchanged |
| 6 | Missing `infile` | Tool aborts; no output written |
| 7 | `n_files=0` | Tool aborts |
| MPI | 4 ranks, `n_files` = 1, 2, 4 | Exactly `n_files` files per level; data identical |
| MPI | 2 ranks, `n_files=16` | Capped at 2 files per level; warning printed |
| MPI | 4 ranks, single-file input, `n_files=4` | 4 files per level; data identical to original |

---

## Manual comparison

```bash
python3 compare_plt.py plt_a plt_b [--tol 1e-5] [--max-files N]
```

Exits 0 when every level, MultiFab, box and value of the two plotfiles match
(within `--tol`), and, with `--max-files`, `plt_b` has at most `N` data files
per MultiFab.
