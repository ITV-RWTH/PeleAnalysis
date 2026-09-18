# stream Test Suite

Tests for the stream tool chain: `Src/partStream.cpp` (with `Src/StreamPC.cpp`)
and `Src/streamBinTubeStats.cpp`. Each case runs the whole chain

```
generateTestPlt  ->  isosurface  ->  partStream  ->  streamBinTubeStats
```

## Running

```bash
cd PeleAnalysis/Tests/stream
./run_tests.sh              # build + run
./run_tests.sh --no-compile # skip the build, reuse the executables in Src/
```

Artefacts (plotfiles, surfaces, streams, logs) are written to `testrun/` and can
be deleted freely. MPI tests are skipped when neither `mpirun` nor `mpiexec` is
found; one error-handling test is skipped without `python3`.

## Manual setup

```bash
cd PeleAnalysis/Src
make EBASE=generateTestPlt    DIM=3 DEBUG=FALSE USE_MPI=FALSE COMP=gnu
make EBASE=isosurface         DIM=3 DEBUG=FALSE USE_MPI=FALSE COMP=gnu
make EBASE=partStream         DIM=3 DEBUG=FALSE USE_MPI=FALSE COMP=gnu
make EBASE=partStream         DIM=3 DEBUG=FALSE USE_MPI=TRUE  COMP=gnu
make EBASE=streamBinTubeStats DIM=3 DEBUG=FALSE USE_MPI=FALSE COMP=gnu
```
plus the `DIM=2` variants of everything except the MPI build.

## Test plotfiles

| Plotfile | Domain | Periodicity | Vector field `ux uy uz` |
|---|---|---|---|
| `plt_stream_3d`      | `[0,1]^3`            | `1 1 1` | `(1,0,0)` |
| `plt_stream_3d_stop` | `[0,1]^3`            | `1 1 1` | `(1,0,0)` for `x < 0.9`, `(0,0,0)` above |
| `plt_stream_3d_crop` | `[0,1]x[0,0.9]x[0,1]`| `1 0 1` | `(1,0,0)` |
| `plt_stream_2d`      | `[0,1]^2`            | `1 1`   | `(1,0)` |

All of them carry `c` (a smoothed plane whose `0.5` contour is the seed
surface), `temp = 300` and `q = 1`. `temp` is not optional:
`streamBinTubeStats` looks up a variable of that name unconditionally.

## Where the expected numbers come from

With `n_cell = 32` on a unit domain, `dx = 1/32`. `partStream` is run with
`hRK = 0.5` and `cSpace = 0`, so every step advances by

```
step = hRK * dx = 1/64
```

`Nsteps = 20` gives `2*20 - 1 = 39` points per stream, i.e. 38 intervals:

```
stream length = 38/64 = 0.59375
```

`calcIntegral` divides the tube integral by the element area, so for the
constant field `q = 1` the reported `q_volInt` *is* the stream length, the same
for every element.

The `c` field rises from 0 to 1 around `x = 0.8`. In a periodic x-direction it
also drops back across the boundary, so the `c = 0.5` isosurface has two flat
sheets, each of unit extent:

```
total area   = 2
total volume = area * length = 2 * 0.59375 = 1.1875
```

Both totals are sums over elements. Each element is written once per corner
with identical values, so the runner sums every `AMREX_SPACEDIM`-th row.

Before the periodicity fix the streams were held at the periodic face and these
came out as 0.497–0.594 per tube and 1.090 in total.

## Test matrix

| # | Input | What is checked |
|---|---|---|
| T1 | `test1_periodic_3d.inp`, `stats1_3d.inp` | Streams cross a periodic boundary: `q_volInt = 0.59375` in every element, `q_avg = 1`, total area `2`, total volume `1.1875`, no NaN/Inf, mapped `q = 1` along every stream, no particles lost |
| T2 | `test2_stop_3d.inp`, `stats2_stop_3d.inp` | Where the vector field vanishes the streams stop: positions and statistics stay finite, no stream passes `x = 0.95`, tubes are shorter than in T1 (this used to divide by a zero magnitude and clamp the resulting NaN to the domain corner) |
| T4 | `test4_ghost_3d.inp` | Seeds between the last cell centre and a non-periodic wall map `temp = 300` and `q = 1` exactly (before the ghost cells outside that boundary were filled: `0.75` of the true value in serial, garbage under MPI) |
| T5 | `test5_periodic_2d.inp`, `stats3_2d.inp` | The 2D chain reproduces the same three numbers |
| T3 | `test3_outside_3d.inp` | Seeds outside the domain in a non-periodic direction abort, the message names the direction, and no output is written |
| T6 | `test6_nosurface.inp` | Asking for the stream binary without a seed surface aborts |
| T7 | `stats5_noheader.inp` | A stream directory without a readable `Header` aborts |
| T8 | `stats4_missing.inp` | A surface node whose stream is missing aborts (`break_stream_bin.py` removes one record from a good directory, which is the only way to reach this now that partStream refuses to produce such a set) |
| T-MPI2, T-MPI4 | `test1_periodic_3d.inp` | 2 and 4 ranks reproduce T1's numbers, and the surface file matches the serial one to 1e-12 |
| T-MPI4 wall | `test4_ghost_3d.inp` | The near-wall values are decomposition-independent |
| T-MPI4 abort | `test3_outside_3d.inp` | All ranks abort, under a timeout so a deadlock fails rather than hangs |

## Helper

| Script | Purpose |
|--------|---------|
| `break_stream_bin.py` | Removes the last stream record from a stream binary directory and decrements that file's count, leaving the `Header` untouched — what a truncated write looks like to the reader |
