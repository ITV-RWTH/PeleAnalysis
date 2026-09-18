#!/usr/bin/env bash
# run_tests.sh — compile and run all stream tests (2D/3D serial and 3D MPI)
#
# Usage:  ./run_tests.sh [--no-compile]
#
#   --no-compile   Skip the build step; assume executables already exist in Src/
#
# The suite runs the whole chain generateTestPlt -> isosurface -> partStream ->
# streamBinTubeStats.  All test artefacts are written to Tests/stream/testrun/
# so the source tree stays clean.  MPI tests require mpirun/mpiexec; they are
# skipped if neither is found.  Python 3 is required for one test that damages
# a stream file on purpose.

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$SCRIPT_DIR/../../Src"
RUN_DIR="$SCRIPT_DIR/testrun"
NPROC=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
COMPILE=true
for arg in "$@"; do
    case "$arg" in
        --no-compile) COMPILE=false ;;
        *) echo "Unknown argument: $arg"; exit 1 ;;
    esac
done

# ---------------------------------------------------------------------------
# Terminal colours
# ---------------------------------------------------------------------------
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
BOLD='\033[1m'; NC='\033[0m'

PASS=0; FAIL=0

section() { echo -e "\n${BOLD}=== $* ===${NC}"; }
pass()    { echo -e "  ${GREEN}PASS${NC} $*"; ((PASS++)) || true; }
fail()    { echo -e "  ${RED}FAIL${NC} $*"; ((FAIL++)) || true; }
skip()    { echo -e "  ${YELLOW}SKIP${NC} $*"; }

# The analytic answers, derived in TESTING.md:
#   step   = hRK * dx = 0.5 / 32                      = 1/64
#   length = (2*Nsteps - 2) * step = 38/64            = 0.59375
#   area   = two iso sheets of unit extent            = 2
#   volume = area * length                            = 1.1875
STREAM_LENGTH=0.59375
SURFACE_AREA=2.0
TUBE_VOLUME=1.1875

# ---------------------------------------------------------------------------
# Assertion helpers
# ---------------------------------------------------------------------------

check_file() {   # check_file <desc> <path>
    if [[ -f "$2" ]]; then pass "$1"; else fail "$1 — missing: $2"; fi
}

check_dir() {    # check_dir <desc> <path>
    if [[ -d "$2" ]]; then pass "$1"; else fail "$1 — missing dir: $2"; fi
}

# Tecplot files written by streamBinTubeStats start with
#   VARIABLES = X Y Z area volume q_avg q_volInt
#   ZONE ...
# so data starts on line 3 and column names can be looked up by name.
col_of() {       # col_of <file> <name>  -> column index, empty if absent
    awk -v n="$2" 'NR==1 { for (i = 3; i <= NF; i++) if ($i == n) { print i - 2; exit } }' "$1"
}

# Number of variables per row.  The connectivity written at the end of the file
# has only AMREX_SPACEDIM fields per line, so data rows are the ones carrying a
# full set of values.
nvars_of() {     # nvars_of <file>
    awk 'NR==1 { print NF - 2; exit }' "$1"
}

# col_all_close <desc> <file> <column name> <expected> <tol>
col_all_close() {
    local desc="$1" file="$2" name="$3" exp="$4" tol="$5" ci nv out
    if [[ ! -f "$file" ]]; then fail "$desc — missing: $file"; return; fi
    ci=$(col_of "$file" "$name"); nv=$(nvars_of "$file")
    if [[ -z "$ci" ]]; then fail "$desc — no column '$name' in $file"; return; fi
    if out=$(awk -v c="$ci" -v nv="$nv" -v e="$exp" -v t="$tol" '
        NR <= 2 { next }
        NF < nv { next }
        { n++; d = $c - e; if (d < 0) d = -d; if (d > m) m = d }
        END { if (n == 0) { printf "no data rows"; exit 1 }
              printf "%d rows, max |diff| = %.3g", n, m
              if (m > t) exit 1 }' "$file"); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out"
    fi
}

# col_sum_close <desc> <file> <column name> <values per element> <expected> <tol>
# Each element is written once per corner with identical values, so only every
# <values per element>-th data row is summed.
col_sum_close() {
    local desc="$1" file="$2" name="$3" per="$4" exp="$5" tol="$6" ci nv out
    if [[ ! -f "$file" ]]; then fail "$desc — missing: $file"; return; fi
    ci=$(col_of "$file" "$name"); nv=$(nvars_of "$file")
    if [[ -z "$ci" ]]; then fail "$desc — no column '$name' in $file"; return; fi
    if out=$(awk -v c="$ci" -v nv="$nv" -v p="$per" -v e="$exp" -v t="$tol" '
        NR <= 2 { next }
        NF < nv { next }
        { r++; if ((r - 1) % p == 0) { s += $c; n++ } }
        END { if (n == 0) { printf "no data rows"; exit 1 }
              d = s - e; if (d < 0) d = -d
              printf "%d elements, sum = %.10g", n, s
              if (d > t) exit 1 }' "$file"); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out"
    fi
}

# col_range_within <desc> <file> <column name> <min below> <max at most>
# Asserts that the smallest value is below <min below> (something did shorten)
# and no value exceeds <max at most>.
col_range_within() {
    local desc="$1" file="$2" name="$3" lo="$4" hi="$5" ci nv out
    if [[ ! -f "$file" ]]; then fail "$desc — missing: $file"; return; fi
    ci=$(col_of "$file" "$name"); nv=$(nvars_of "$file")
    if [[ -z "$ci" ]]; then fail "$desc — no column '$name' in $file"; return; fi
    if out=$(awk -v c="$ci" -v nv="$nv" -v lo="$lo" -v hi="$hi" '
        NR <= 2 { next }
        NF < nv { next }
        { n++; if (n == 1 || $c < mn) mn = $c; if (n == 1 || $c > mx) mx = $c }
        END { if (n == 0) { printf "no data rows"; exit 1 }
              printf "min = %.9g, max = %.9g", mn, mx
              if (mn >= lo || mx > hi) exit 1 }' "$file"); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out (wanted min < $lo and max <= $hi)"
    fi
}

# matrix_max_below <desc> <file> <limit>  — for the headerless matlab dumps
matrix_max_below() {
    local desc="$1" file="$2" lim="$3" out
    if [[ ! -f "$file" ]]; then fail "$desc — missing: $file"; return; fi
    if out=$(awk -v l="$lim" '
        { for (i = 1; i <= NF; i++) { n++; if (n == 1 || $i > m) m = $i } }
        END { if (n == 0) { printf "no values"; exit 1 }
              printf "max = %.9g", m
              if (m >= l) exit 1 }' "$file"); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out (limit $lim)"
    fi
}

# all_finite <desc> <file>
all_finite() {
    local desc="$1" file="$2" out
    if [[ ! -f "$file" ]]; then fail "$desc — missing: $file"; return; fi
    if out=$(awk '
        { for (i = 1; i <= NF; i++)
              if ($i ~ /[nN][aA][nN]/ || $i ~ /[iI][nN][fF]/) bad++ }
        END { printf "%d non-finite value(s)", bad + 0; if (bad > 0) exit 1 }' \
        "$file"); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out"
    fi
}

# matrix_all_close <desc> <file> <expected> <tol>
# For the whitespace matrices written by writeStreamsToMatlab and for the
# Tecplot stream dumps, where every number must equal the same value.
matrix_all_close() {
    local desc="$1" file="$2" exp="$3" tol="$4" out
    if [[ ! -f "$file" ]]; then fail "$desc — missing: $file"; return; fi
    if out=$(awk -v e="$exp" -v t="$tol" '
        /^VARIABLES/ || /^ZONE/ { next }
        { for (i = 1; i <= NF; i++) { n++; d = $i - e; if (d < 0) d = -d
                                      if (d > m) m = d } }
        END { if (n == 0) { printf "no values"; exit 1 }
              printf "%d values, max |diff| = %.3g", n, m
              if (m > t) exit 1 }' "$file"); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out"
    fi
}

# tec_col_all_close <desc> <file> <column> <expected> <tol>
# Tecplot stream dumps (writeStreams) carry quoted variable names and a ZONE
# line per stream, so the column is given by number: X Y Z <vars...>.
tec_col_all_close() {
    local desc="$1" file="$2" ci="$3" exp="$4" tol="$5" out
    if [[ ! -f "$file" ]]; then fail "$desc — missing: $file"; return; fi
    if out=$(awk -v c="$ci" -v e="$exp" -v t="$tol" '
        /^VARIABLES/ || /^ZONE/ { next }
        NF < c { next }
        { n++; d = $c - e; if (d < 0) d = -d; if (d > m) m = d }
        END { if (n == 0) { printf "no data rows"; exit 1 }
              printf "%d values, max |diff| = %.3g", n, m
              if (m > t) exit 1 }' "$file"); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out"
    fi
}

# files_close <desc> <file a> <file b> <tol>
files_close() {
    local desc="$1" a="$2" b="$3" tol="$4" out
    if [[ ! -f "$a" || ! -f "$b" ]]; then
        fail "$desc — missing $a or $b"; return
    fi
    if out=$(awk -v t="$tol" '
        FNR == NR { for (i = 1; i <= NF; i++) v[FNR, i] = $i; rows = FNR; next }
        { if (FNR > rows) { extra++; next }
          for (i = 1; i <= NF; i++) { d = $i - v[FNR, i]; if (d < 0) d = -d
                                      if (d > m) m = d } }
        END { printf "max |diff| = %.3g", m
              if (extra) { printf ", %d extra rows", extra; exit 1 }
              if (m > t) exit 1 }' "$a" "$b"); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out"
    fi
}

# expect_abort <desc> <command...>
# Passes if the command exits non-zero, as amrex::Abort does.
expect_abort() {
    local desc="$1"; shift
    if "$@" > expect_abort.log 2>&1; then
        fail "$desc — expected a non-zero exit, got 0"
    else
        pass "$desc"
    fi
}

# expect_abort_saying <desc> <pattern> <command...>
expect_abort_saying() {
    local desc="$1" pat="$2"; shift 2
    if "$@" > expect_abort.log 2>&1; then
        fail "$desc — expected a non-zero exit, got 0"
    elif grep -q "$pat" expect_abort.log; then
        pass "$desc"
    else
        fail "$desc — aborted, but the message does not mention '$pat'"
    fi
}

# run_with_timeout <seconds> <command...> — so a deadlock fails instead of hangs
run_with_timeout() {
    local secs="$1"; shift
    perl -e 'alarm shift; exec @ARGV' "$secs" "$@"
}

# ---------------------------------------------------------------------------
# Executable discovery
# ---------------------------------------------------------------------------
find_exe() {   # find_exe <glob>
    ls "$SRC_DIR"/$1 2>/dev/null | head -1 || true
}

# ---------------------------------------------------------------------------
# Phase 1 — Build
# ---------------------------------------------------------------------------
section "Phase 1: Build"

BUILD_OPTS="DEBUG=FALSE PRECISION=DOUBLE COMP=gnu -j$NPROC"

build_tool() {   # build_tool <label> <log> <make_args...>
    local label="$1" log="$2"; shift 2
    echo "  Building $label..."
    if make -C "$SRC_DIR" $BUILD_OPTS "$@" > "$SCRIPT_DIR/$log" 2>&1; then
        pass "$label built"
    else
        fail "$label build failed — see $log"
    fi
}

if $COMPILE; then
    build_tool "generateTestPlt 3D" build_genPlt3d.log \
        EBASE=generateTestPlt DIM=3 USE_MPI=FALSE
    build_tool "generateTestPlt 2D" build_genPlt2d.log \
        EBASE=generateTestPlt DIM=2 USE_MPI=FALSE
    build_tool "isosurface 3D" build_iso3d.log \
        EBASE=isosurface DIM=3 USE_MPI=FALSE
    build_tool "isosurface 2D" build_iso2d.log \
        EBASE=isosurface DIM=2 USE_MPI=FALSE
    build_tool "partStream 3D serial" build_ps3d.log \
        EBASE=partStream DIM=3 USE_MPI=FALSE
    build_tool "partStream 2D serial" build_ps2d.log \
        EBASE=partStream DIM=2 USE_MPI=FALSE
    build_tool "partStream 3D MPI" build_ps3dmpi.log \
        EBASE=partStream DIM=3 USE_MPI=TRUE
    build_tool "streamBinTubeStats 3D" build_sbts3d.log \
        EBASE=streamBinTubeStats DIM=3 USE_MPI=FALSE
    build_tool "streamBinTubeStats 2D" build_sbts2d.log \
        EBASE=streamBinTubeStats DIM=2 USE_MPI=FALSE
else
    skip "Build skipped (--no-compile)"
fi

GEN3D=$(find_exe "generateTestPlt3d.gnu.ex")
GEN2D=$(find_exe "generateTestPlt2d.gnu.ex")
ISO3D=$(find_exe "isosurface3d.gnu.ex")
ISO2D=$(find_exe "isosurface2d.gnu.ex")
PS3D=$(find_exe "partStream3d.gnu.ex")
PS2D=$(find_exe "partStream2d.gnu.ex")
PS3D_MPI=$(find_exe "partStream3d.gnu.MPI.ex")
SBTS3D=$(find_exe "streamBinTubeStats3d.gnu.ex")
SBTS2D=$(find_exe "streamBinTubeStats2d.gnu.ex")

check_file "generateTestPlt 3D"    "$GEN3D"
check_file "generateTestPlt 2D"    "$GEN2D"
check_file "isosurface 3D"         "$ISO3D"
check_file "isosurface 2D"         "$ISO2D"
check_file "partStream 3D serial"  "$PS3D"
check_file "partStream 2D serial"  "$PS2D"
check_file "partStream 3D MPI"     "$PS3D_MPI"
check_file "streamBinTubeStats 3D" "$SBTS3D"
check_file "streamBinTubeStats 2D" "$SBTS2D"

if command -v mpirun &>/dev/null; then
    MPI_AVAILABLE=true; MPI_CMD=mpirun
elif command -v mpiexec &>/dev/null; then
    MPI_AVAILABLE=true; MPI_CMD=mpiexec
else
    MPI_AVAILABLE=false
    skip "mpirun/mpiexec not found — MPI tests will be skipped"
fi

# ---------------------------------------------------------------------------
# Phase 2 — Plotfiles and seed surfaces
# ---------------------------------------------------------------------------
section "Phase 2: Generate plotfiles and seed surfaces"

rm -rf "$RUN_DIR"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"
T="$SCRIPT_DIR"

echo "  Generating plotfiles..."
"$GEN3D" "$T/gen_stream_3d.inp"      > gen_3d.log      2>&1
"$GEN3D" "$T/gen_stream_3d_stop.inp" > gen_3d_stop.log 2>&1
"$GEN3D" "$T/gen_stream_3d_crop.inp" > gen_3d_crop.log 2>&1
"$GEN2D" "$T/gen_stream_2d.inp"      > gen_2d.log      2>&1

check_dir "plt_stream_3d"      "plt_stream_3d"
check_dir "plt_stream_3d_stop" "plt_stream_3d_stop"
check_dir "plt_stream_3d_crop" "plt_stream_3d_crop"
check_dir "plt_stream_2d"      "plt_stream_2d"

echo "  Extracting seed surfaces..."
"$ISO3D" "$T/iso_3d.inp" > iso_3d.log 2>&1
"$ISO2D" "$T/iso_2d.inp" > iso_2d.log 2>&1

check_file "3D seed surface" "surf_3d.mef"
check_file "2D seed surface" "surf_2d.mef"

# ---------------------------------------------------------------------------
# Phase 3 — Streams across a periodic boundary (3D)
# ---------------------------------------------------------------------------
section "Phase 3: Periodic streams (3D)"

if "$PS3D" "$T/test1_periodic_3d.inp" > ps_3d.log 2>&1; then
    pass "T1 partStream ran"
else
    fail "T1 partStream failed — see testrun/ps_3d.log"
fi
check_file "T1 stream Header" "str_3d/Header"

if grep -q "LOST" ps_3d.log; then
    fail "T1 no particles lost"
else
    pass "T1 no particles lost"
fi

if "$SBTS3D" "$T/stats1_3d.inp" > sbts_3d.log 2>&1; then
    pass "T1 streamBinTubeStats ran"
else
    fail "T1 streamBinTubeStats failed — see testrun/sbts_3d.log"
fi

# A stream that is stopped at a periodic face instead of wrapping is shorter
# than 0.59375; this is the regression that motivated the periodicity fix.
col_all_close "T1 tube integral of q = stream length" \
    "str_3d_binVolInt.dat" "q_volInt" "$STREAM_LENGTH" 1e-9
col_all_close "T1 surface average of q" \
    "str_3d_binVolInt.dat" "q_avg" 1.0 1e-12
col_sum_close "T1 total surface area" \
    "str_3d_binVolInt.dat" "area" 3 "$SURFACE_AREA" 1e-9
col_sum_close "T1 total stream tube volume" \
    "str_3d_binVolInt.dat" "volume" 3 "$TUBE_VOLUME" 1e-9
all_finite "T1 no NaN or Inf in the surface file" "str_3d_binVolInt.dat"
matrix_all_close "T1 mapped q along every stream" "str_3d/q.dat" 1.0 1e-12
check_file "T1 tecplot surface from streams" "str_3d_surfTec.dat"

# ---------------------------------------------------------------------------
# Phase 4 — Vanishing vector field (3D)
# ---------------------------------------------------------------------------
section "Phase 4: Vanishing vector field (3D)"

if "$PS3D" "$T/test2_stop_3d.inp" > ps_3d_stop.log 2>&1; then
    pass "T2 partStream ran"
else
    fail "T2 partStream failed — see testrun/ps_3d_stop.log"
fi

if "$SBTS3D" "$T/stats2_stop_3d.inp" > sbts_3d_stop.log 2>&1; then
    pass "T2 streamBinTubeStats ran"
else
    fail "T2 streamBinTubeStats failed — see testrun/sbts_3d_stop.log"
fi

# Dividing by a zero magnitude used to put NaN into the position, which was
# then clamped to the corner of the domain.
all_finite "T2 no NaN or Inf in the stream positions" "str_3d_stop/X.dat"
all_finite "T2 no NaN or Inf in the surface file" "str_3d_stop_binVolInt.dat"
# ux is zero above x = 0.9, so no stream position may go far beyond it.  The
# old code divided by a zero magnitude and the NaN was clamped to the domain
# corner, which showed up as positions far away from where the field vanished.
matrix_max_below "T2 no stream runs past the vanishing field" \
    "str_3d_stop/X.dat" 0.95
# The seeds on the second (wrapped) sheet never reach the region where the
# field vanishes, so their tubes keep the full length; the ones approaching
# x = 0.9 from below must be shorter, and none may exceed the full length.
col_range_within "T2 tubes at the vanishing field are shorter" \
    "str_3d_stop_binVolInt.dat" "q_volInt" "$STREAM_LENGTH" 0.59375001

# ---------------------------------------------------------------------------
# Phase 5 — Interpolation next to a non-periodic boundary (3D)
# ---------------------------------------------------------------------------
section "Phase 5: Seeds next to a non-periodic wall (3D)"

if "$PS3D" "$T/test4_ghost_3d.inp" > ps_ghost.log 2>&1; then
    pass "T4 partStream ran"
else
    fail "T4 partStream failed — see testrun/ps_ghost.log"
fi

# Both mapped fields are constant.  Before the domain boundary was filled,
# these seeds interpolated against uninitialised ghost cells and returned
# 0.75 * the true value in serial, and garbage under MPI.
tec_col_all_close "T4 mapped temp next to the wall" \
    "str_ghost/str_00000.dat" 4 300.0 1e-9
tec_col_all_close "T4 mapped q next to the wall" \
    "str_ghost/str_00000.dat" 5 1.0 1e-12

# ---------------------------------------------------------------------------
# Phase 6 — 2D
# ---------------------------------------------------------------------------
section "Phase 6: Periodic streams (2D)"

if "$PS2D" "$T/test5_periodic_2d.inp" > ps_2d.log 2>&1; then
    pass "T5 partStream ran"
else
    fail "T5 partStream failed — see testrun/ps_2d.log"
fi

if "$SBTS2D" "$T/stats3_2d.inp" > sbts_2d.log 2>&1; then
    pass "T5 streamBinTubeStats ran"
else
    fail "T5 streamBinTubeStats failed — see testrun/sbts_2d.log"
fi

col_all_close "T5 tube integral of q = stream length" \
    "str_2d_binVolInt.dat" "q_volInt" "$STREAM_LENGTH" 1e-9
col_sum_close "T5 total surface length" \
    "str_2d_binVolInt.dat" "area" 2 "$SURFACE_AREA" 1e-9
col_sum_close "T5 total stream tube volume" \
    "str_2d_binVolInt.dat" "volume" 2 "$TUBE_VOLUME" 1e-9
all_finite "T5 no NaN or Inf in the surface file" "str_2d_binVolInt.dat"

# ---------------------------------------------------------------------------
# Phase 7 — Error handling
# ---------------------------------------------------------------------------
section "Phase 7: Error handling"

expect_abort_saying "T3 seeds outside a non-periodic domain abort" \
    "outside the domain" "$PS3D" "$T/test3_outside_3d.inp"
if [[ -d str_3d_outside ]]; then
    fail "T3 no stream output written on abort"
else
    pass "T3 no stream output written on abort"
fi

expect_abort "T6 stream binary without a seed surface aborts" \
    "$PS3D" "$T/test6_nosurface.inp"

expect_abort_saying "T7 missing stream Header aborts" \
    "Could not open" "$SBTS3D" "$T/stats5_noheader.inp"

if command -v python3 &>/dev/null; then
    rm -rf str_3d_broken
    cp -R str_3d str_3d_broken
    if python3 "$T/break_stream_bin.py" str_3d_broken > break.log 2>&1; then
        expect_abort_saying "T8 a surface node without a stream aborts" \
            "have no stream" "$SBTS3D" "$T/stats4_missing.inp"
    else
        fail "T8 could not remove a stream record — see testrun/break.log"
    fi
else
    skip "T8 a surface node without a stream (python3 not available)"
fi

# ---------------------------------------------------------------------------
# Phase 8 — MPI
# ---------------------------------------------------------------------------
section "Phase 8: MPI"

if ! $MPI_AVAILABLE; then
    skip "All MPI tests (mpirun/mpiexec not found)"
elif [[ -z "$PS3D_MPI" ]]; then
    skip "All MPI tests (partStream 3D MPI executable not found)"
else
    for np in 2 4; do
        echo "  partStream on $np ranks..."
        if run_with_timeout 600 $MPI_CMD -np $np "$PS3D_MPI" \
                "$T/test1_periodic_3d.inp" streamBinfile=str_3d_mpi$np \
                > ps_3d_mpi$np.log 2>&1; then
            pass "T-MPI$np partStream ran"
        else
            fail "T-MPI$np partStream failed — see testrun/ps_3d_mpi$np.log"
            continue
        fi

        if "$SBTS3D" "$T/stats1_3d.inp" infile=str_3d_mpi$np \
                > sbts_3d_mpi$np.log 2>&1; then
            pass "T-MPI$np streamBinTubeStats ran"
        else
            fail "T-MPI$np streamBinTubeStats failed"
            continue
        fi

        col_all_close "T-MPI$np tube integral of q" \
            "str_3d_mpi${np}_binVolInt.dat" "q_volInt" "$STREAM_LENGTH" 1e-9
        col_sum_close "T-MPI$np total stream tube volume" \
            "str_3d_mpi${np}_binVolInt.dat" "volume" 3 "$TUBE_VOLUME" 1e-9
        files_close "T-MPI$np surface file matches the serial one" \
            "str_3d_binVolInt.dat" "str_3d_mpi${np}_binVolInt.dat" 1e-12
    done

    # The ghost cells outside a non-periodic boundary are where a wrong answer
    # used to depend on the decomposition.
    if run_with_timeout 600 $MPI_CMD -np 4 "$PS3D_MPI" \
            "$T/test4_ghost_3d.inp" streamfile=str_ghost_mpi4 \
            > ps_ghost_mpi4.log 2>&1; then
        pass "T-MPI4 partStream next to a wall ran"
        cat str_ghost_mpi4/str_*.dat > str_ghost_mpi4_all.dat
        tec_col_all_close "T-MPI4 mapped temp next to the wall" \
            "str_ghost_mpi4_all.dat" 4 300.0 1e-9
        tec_col_all_close "T-MPI4 mapped q next to the wall" \
            "str_ghost_mpi4_all.dat" 5 1.0 1e-12
    else
        fail "T-MPI4 partStream next to a wall failed"
    fi

    # All ranks must abort, and none may hang waiting in a collective.
    expect_abort "T-MPI4 seeds outside a non-periodic domain abort" \
        run_with_timeout 300 $MPI_CMD -np 4 "$PS3D_MPI" \
        "$T/test3_outside_3d.inp"
fi

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
TOTAL=$((PASS + FAIL))
echo ""
echo -e "${BOLD}Results: ${GREEN}$PASS${NC}${BOLD}/$TOTAL passed${NC}"
if (( FAIL > 0 )); then
    echo -e "${RED}$FAIL test(s) failed.${NC}"
    exit 1
else
    echo -e "${GREEN}All tests passed.${NC}"
fi
