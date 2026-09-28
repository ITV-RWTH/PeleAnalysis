#!/usr/bin/env bash
# run_tests.sh — build and run all sequentialFilterPlt tests (3D serial, OpenMP and MPI)
#
# Usage:  ./run_tests.sh [--no-compile]
#
#   --no-compile   Skip the build step; assume executables already exist
#
# All test artefacts are written to Tests/sequentialFilterPlt/testrun/ so the
# source tree stays clean.  Plotfiles are compared with AMReX's fcompare, built
# from $AMREX_HOME/Tools/Plotfile.  MPI tests require mpirun/mpiexec and are
# skipped if neither is found.  Python 3 is required for the bounds checks.

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$SCRIPT_DIR/../../Src"
RUN_DIR="$SCRIPT_DIR/testrun"
AMREX_HOME="${AMREX_HOME:-$SCRIPT_DIR/../../Submodules/PelePhysics/Submodules/amrex}"
PLOTFILE_TOOLS="$AMREX_HOME/Tools/Plotfile"
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

# ---------------------------------------------------------------------------
# Assertion helpers
# ---------------------------------------------------------------------------

# Filtered values are O(1); differing summation orders give ~1e-15
EQ_TOL=1e-12

check_file() {   # check_file <desc> <path>
    if [[ -f "$2" ]]; then pass "$1"; else fail "$1 — missing: $2"; fi
}

check_dir() {    # check_dir <desc> <path>
    if [[ -d "$2" ]]; then pass "$1"; else fail "$1 — missing dir: $2"; fi
}

no_dir() {       # no_dir <desc> <path>  — assert directory does NOT exist
    if [[ ! -d "$2" ]]; then pass "$1"; else fail "$1 — unexpected dir: $2"; fi
}

# run_ok <desc> <cmd...>  — assert the command exits 0
run_ok() {
    local desc="$1"; shift
    if "$@" > /dev/null 2>&1; then pass "$desc"; else fail "$desc — exited non-zero"; fi
}

# run_aborts <desc> <cmd...>  — assert the command exits non-zero
run_aborts() {
    local desc="$1"; shift
    if "$@" > /dev/null 2>&1; then
        fail "$desc — expected abort, tool exited 0"
    else
        pass "$desc — aborted correctly"
    fi
}

# plt_agree <desc> <pltA> <pltB> [fcompare options...]
# fcompare fails on any NaN and on any difference above EQ_TOL
plt_agree() {
    local desc="$1" a="$2" b="$3"; shift 3
    if [[ ! -d "$a" || ! -d "$b" ]]; then
        fail "$desc — plotfile missing: $a or $b"; return
    fi
    local out
    if out=$("$FCOMPARE" --abs_tol "$EQ_TOL" "$@" "$a" "$b" 2>&1); then
        pass "$desc"
    else
        fail "$desc"; echo "$out" | sed 's/^/        /'
    fi
}

# plt_differ <desc> <pltA> <pltB>  — assert the plotfiles are NOT equal
plt_differ() {
    local desc="$1" a="$2" b="$3"
    if "$FCOMPARE" --abs_tol "$EQ_TOL" "$a" "$b" > /dev/null 2>&1; then
        fail "$desc — plotfiles unexpectedly agree"
    else
        pass "$desc"
    fi
}

# no_nan <desc> <plt>  — fcompare of a plotfile with itself fails iff it holds a NaN
no_nan() {
    local desc="$1" plt="$2"
    if [[ ! -d "$plt" ]]; then
        fail "$desc — plotfile missing: $plt"; return
    fi
    if "$FCOMPARE" "$plt" "$plt" > /dev/null 2>&1; then
        pass "$desc"
    else
        fail "$desc — NaN present"
    fi
}

# check_bounds <desc> <plt> <var> <lower> <upper>
check_bounds() {
    local desc="$1" plt="$2" var="$3" lo="$4" hi="$5"
    if [[ ! -d "$plt" ]]; then
        fail "$desc — plotfile missing: $plt"; return
    fi
    if ! command -v python3 &>/dev/null; then
        skip "$desc — python3 not available"; return
    fi
    local out
    if out=$(python3 "$SCRIPT_DIR/check_plt_bounds.py" "$plt" "$var" "$lo" "$hi" 2>&1); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out"
    fi
}

# header_line <plt> <offset>  — Header line <offset> after the variable names
header_line() {
    local ncomp
    ncomp=$(sed -n '2p' "$1/Header")
    sed -n "$((ncomp + 2 + $2))p" "$1/Header"
}

# same_amr_layout <desc> <input_plt> <output_plt>
# Compares finest level, ref ratios and level domains of the two Headers
same_amr_layout() {
    local desc="$1" in="$2" out="$3" off
    if [[ ! -d "$out" ]]; then
        fail "$desc — plotfile missing: $out"; return
    fi
    for off in 3 6 7; do
        if [[ "$(header_line "$in" $off)" != "$(header_line "$out" $off)" ]]; then
            fail "$desc — Header differs: '$(header_line "$in" $off)' vs '$(header_line "$out" $off)'"
            return
        fi
    done
    pass "$desc  (ref ratio: $(header_line "$out" 6))"
}

# ---------------------------------------------------------------------------
# Executable discovery
# ---------------------------------------------------------------------------
find_exe() {   # find_exe <dir> <name>
    ls "$1"/$2 2>/dev/null | head -1 || true
}

# ---------------------------------------------------------------------------
# Phase 1 — Build
# ---------------------------------------------------------------------------
section "Phase 1: Build"

BUILD_OPTS="DEBUG=FALSE PRECISION=DOUBLE COMP=gnu -j$NPROC"

build_tool() {   # build_tool <label> <log> <make_dir> <make_args...>
    local label="$1" log="$2" dir="$3"; shift 3
    echo "  Building $label..."
    if make -C "$dir" $BUILD_OPTS "$@" > "$SCRIPT_DIR/$log" 2>&1; then
        pass "$label built"
    else
        fail "$label build failed — see $log"
    fi
}

if $COMPILE; then
    build_tool "generateTestPlt 3D"         build_genPlt3d.log   "$SRC_DIR" EBASE=generateTestPlt     DIM=3 USE_MPI=FALSE USE_OMP=FALSE
    build_tool "sequentialFilterPlt 3D"     build_seq3d.log      "$SRC_DIR" EBASE=sequentialFilterPlt DIM=3 USE_MPI=FALSE USE_OMP=FALSE
    build_tool "sequentialFilterPlt 3D OMP" build_seq3domp.log   "$SRC_DIR" EBASE=sequentialFilterPlt DIM=3 USE_MPI=FALSE USE_OMP=TRUE
    build_tool "filterPlt 3D"               build_filter3d.log   "$SRC_DIR" EBASE=filterPlt           DIM=3 USE_MPI=FALSE USE_OMP=FALSE
    build_tool "filterPlt 3D OMP"           build_filter3domp.log "$SRC_DIR" EBASE=filterPlt          DIM=3 USE_MPI=FALSE USE_OMP=TRUE
    if command -v mpicxx &>/dev/null; then
        build_tool "sequentialFilterPlt 3D MPI" build_seq3dmpi.log "$SRC_DIR" EBASE=sequentialFilterPlt DIM=3 USE_MPI=TRUE USE_OMP=FALSE
    else
        skip "sequentialFilterPlt 3D MPI build (mpicxx not found)"
    fi
    build_tool "fcompare"                   build_fcompare.log   "$PLOTFILE_TOOLS" programs=fcompare  DIM=3 USE_MPI=FALSE USE_OMP=FALSE
else
    skip "Build skipped (--no-compile)"
fi

GEN3D=$(find_exe "$SRC_DIR" "generateTestPlt3d.gnu.ex")
SEQ3D=$(find_exe "$SRC_DIR" "sequentialFilterPlt3d.gnu.ex")
SEQ3D_OMP=$(find_exe "$SRC_DIR" "sequentialFilterPlt3d.gnu.OMP.ex")
SEQ3D_MPI=$(find_exe "$SRC_DIR" "sequentialFilterPlt3d.gnu.MPI.ex")
FILT3D=$(find_exe "$SRC_DIR" "filterPlt3d.gnu.ex")
FILT3D_OMP=$(find_exe "$SRC_DIR" "filterPlt3d.gnu.OMP.ex")
FCOMPARE=$(find_exe "$PLOTFILE_TOOLS" "fcompare*.ex")

check_file "generateTestPlt 3D"     "$GEN3D"
check_file "sequentialFilterPlt 3D" "$SEQ3D"
check_file "filterPlt 3D"           "$FILT3D"
check_file "fcompare"               "$FCOMPARE"

if [[ -z "$GEN3D" || -z "$SEQ3D" || -z "$FILT3D" || -z "$FCOMPARE" ]]; then
    echo -e "${RED}Cannot find required executables — aborting test run.${NC}"
    exit 1
fi

if command -v mpirun &>/dev/null; then
    MPI_AVAILABLE=true; MPI_CMD=mpirun
elif command -v mpiexec &>/dev/null; then
    MPI_AVAILABLE=true; MPI_CMD=mpiexec
else
    MPI_AVAILABLE=false
    skip "mpirun/mpiexec not found — MPI tests will be skipped"
fi

# ---------------------------------------------------------------------------
# Phase 2 — Generate test plotfiles
# ---------------------------------------------------------------------------
section "Phase 2: Generate test plotfiles"

rm -rf "$RUN_DIR"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

for name in single single_periodic full_rr2 full_rr4 nested_rr4 tight_3lev; do
    "$GEN3D" "$SCRIPT_DIR/gen_$name.inp" > /dev/null 2>&1
    check_file "plt_$name" "plt_$name/Header"
done

# filt <outfile> <filterPlt args...>  — filterPlt has no outfile option
filt() {
    local out="$1"; shift
    local in
    in=$(printf '%s\n' "$@" | sed -n 's/^infile=//p')
    "$FILT3D" "$@" > /dev/null 2>&1 && mv "${in##*/}_filtered" "$out"
}

# ---------------------------------------------------------------------------
# Phase 3 — Agreement with filterPlt (independent 3D stencil)
# ---------------------------------------------------------------------------
section "Phase 3: Agreement with filterPlt"

echo "  T1–T4 — single level, compared cell by cell (including domain boundaries)"
for spec in "single 2" "single 8" "single 16" "single_periodic 8"; do
    set -- $spec
    "$SEQ3D" infile=plt_$1 base_fgr=$2 outfile=seq_$1_$2 > /dev/null 2>&1 || true
    filt filt_$1_$2 infile=plt_$1 base_fgr=$2 max_grid_size=32 || true
    plt_agree "plt_$1 base_fgr=$2 matches filterPlt" seq_$1_$2 filt_$1_$2
done

echo "  T5–T7 — two levels, level 1 covers the domain (fine filter scaled by ref_ratio)"
for spec in "full_rr2 4 false" "full_rr4 2 false" "full_rr2 4 true"; do
    set -- $spec
    "$SEQ3D" infile=plt_$1 base_fgr=$2 same_fgr_all_levels=$3 outfile=seq_$1_$2_$3 > /dev/null 2>&1 || true
    filt filt_$1_$2_$3 infile=plt_$1 base_fgr=$2 same_fgr_all_levels=$3 max_grid_size=32 || true
    plt_agree "plt_$1 base_fgr=$2 same_fgr_all_levels=$3 matches filterPlt" seq_$1_$2_$3 filt_$1_$2_$3
done

echo "  T8 — filtering actually changes the data"
plt_differ "plt_single base_fgr=8 differs from its input" seq_single_8 plt_single

# ---------------------------------------------------------------------------
# Phase 4 — Coarse/fine ghost cells beyond the next coarser level
# ---------------------------------------------------------------------------
section "Phase 4: Wide filters on tightly nested levels"

# With init_snan, any ghost cell read before it is filled ends up as NaN
for fgr in 2 8; do
    "$SEQ3D" infile=plt_tight_3lev base_fgr=$fgr outfile=seq_tight_$fgr amrex.init_snan=1 > /dev/null 2>&1 || true
    "$SEQ3D" infile=plt_tight_3lev base_fgr=$fgr outfile=seq_tight_${fgr}_nosnan > /dev/null 2>&1 || true
    filt filt_tight_$fgr infile=plt_tight_3lev base_fgr=$fgr amrex.init_snan=1 || true
    no_nan       "T9  sequentialFilterPlt base_fgr=$fgr no NaN on any level" seq_tight_$fgr
    plt_agree    "T10 sequentialFilterPlt base_fgr=$fgr independent of memory init" seq_tight_$fgr seq_tight_${fgr}_nosnan
    check_bounds "T11 sequentialFilterPlt base_fgr=$fgr constant preserved" seq_tight_$fgr cst 3.4999999999 3.5000000001
    check_bounds "T12 sequentialFilterPlt base_fgr=$fgr sph stays in [0, 1]" seq_tight_$fgr sph 0 1
    check_bounds "T13 sequentialFilterPlt base_fgr=$fgr sine stays in [-1, 1]" seq_tight_$fgr sine_f -1 1
    no_nan       "T14 filterPlt base_fgr=$fgr no NaN on any level" filt_tight_$fgr
    check_bounds "T15 filterPlt base_fgr=$fgr sph stays in [0, 1]" filt_tight_$fgr sph 0 1
done

# ---------------------------------------------------------------------------
# Phase 5 — Output metadata and options
# ---------------------------------------------------------------------------
section "Phase 5: Output metadata and options"

"$SEQ3D" infile=plt_nested_rr4 outfile=seq_nested_rr4 > /dev/null 2>&1 || true
filt filt_nested_rr4 infile=plt_nested_rr4 || true
same_amr_layout "T16 sequentialFilterPlt keeps ref_ratio 4 and level domains" plt_nested_rr4 seq_nested_rr4
same_amr_layout "T17 filterPlt keeps ref_ratio 4 and level domains"           plt_nested_rr4 filt_nested_rr4
same_amr_layout "T18 sequentialFilterPlt keeps the 3-level layout"            plt_tight_3lev seq_tight_8

echo "  T19 — default output name <infile>_filtered in the working directory"
mkdir -p default_name
(cd default_name && "$SEQ3D" infile=../plt_single > /dev/null 2>&1) || true
check_dir "T19 default_name/plt_single_filtered" default_name/plt_single_filtered

echo "  T20 — max_filter_level=1 on a 3-level file"
"$SEQ3D" infile=plt_tight_3lev max_filter_level=1 outfile=seq_maxlev1 > /dev/null 2>&1 || true
check_dir "T20 Level_1 written"     seq_maxlev1/Level_1
no_dir    "T20 Level_2 not written" seq_maxlev1/Level_2

echo "  T21 — variables subset equals the same variable from a full run"
"$SEQ3D" infile=plt_single base_fgr=8 variables=sph outfile=seq_sph_only > /dev/null 2>&1 || true
plt_agree "T21 variables=sph matches full run" seq_sph_only seq_single_8
if [[ "$(sed -n '2p' seq_sph_only/Header 2>/dev/null)" == "1" ]]; then
    pass "T21 only one variable written"
else
    fail "T21 expected nComp=1 in seq_sph_only/Header"
fi

echo "  T22–T23 — max_grid_size smaller than the filter width does not change the result"
"$SEQ3D" infile=plt_single base_fgr=32 max_grid_size=8 outfile=seq_mgs8 > /dev/null 2>&1 || true
"$SEQ3D" infile=plt_single base_fgr=32 outfile=seq_mgs_default > /dev/null 2>&1 || true
plt_agree "T22 single level, max_grid_size=8, base_fgr=32" seq_mgs8 seq_mgs_default --allow_diff_grids
"$SEQ3D" infile=plt_tight_3lev base_fgr=8 max_grid_size=8 outfile=seq_tight_mgs8 > /dev/null 2>&1 || true
plt_agree "T23 3 levels, max_grid_size=8, base_fgr=8" seq_tight_mgs8 seq_tight_8_nosnan --allow_diff_grids

echo "  T24 — interp_type=0 (piecewise constant) runs without NaN"
"$SEQ3D" infile=plt_tight_3lev base_fgr=8 interp_type=0 outfile=seq_pc amrex.init_snan=1 > /dev/null 2>&1 || true
no_nan       "T24 interp_type=0 no NaN" seq_pc
check_bounds "T24 interp_type=0 sph stays in [0, 1]" seq_pc sph 0 1

# ---------------------------------------------------------------------------
# Phase 6 — Error handling
# ---------------------------------------------------------------------------
section "Phase 6: Error handling"

run_aborts "T25 filter_type=2"           "$SEQ3D" infile=plt_single filter_type=2 outfile=err1
no_dir     "T25 no output on abort"      err1
run_aborts "T26 odd base_fgr=3"          "$SEQ3D" infile=plt_single base_fgr=3 outfile=err2
run_aborts "T27 base_fgr=0"              "$SEQ3D" infile=plt_single base_fgr=0 outfile=err3
run_aborts "T28 unknown variable"        "$SEQ3D" infile=plt_single variables="sph nope" outfile=err4
run_aborts "T29 missing infile argument" "$SEQ3D" outfile=err5
run_aborts "T30 nonexistent infile"      "$SEQ3D" infile=no_such_plt outfile=err6

# ---------------------------------------------------------------------------
# Phase 7 — OpenMP
# ---------------------------------------------------------------------------
section "Phase 7: OpenMP"

if [[ -z "$SEQ3D_OMP" || -z "$FILT3D_OMP" ]]; then
    skip "OpenMP tests (OMP executables not found)"
else
    OMP_NUM_THREADS=4 "$SEQ3D_OMP" infile=plt_tight_3lev base_fgr=8 outfile=seq_omp > /dev/null 2>&1 || true
    plt_agree "T31 sequentialFilterPlt 4 threads matches serial" seq_omp seq_tight_8_nosnan
    OMP_NUM_THREADS=4 "$FILT3D_OMP" infile=plt_single base_fgr=8 max_grid_size=32 > /dev/null 2>&1 \
        && mv plt_single_filtered filt_omp || true
    plt_agree "T32 filterPlt 4 threads matches serial" filt_omp filt_single_8
fi

# ---------------------------------------------------------------------------
# Phase 8 — MPI
# ---------------------------------------------------------------------------
section "Phase 8: MPI"

if ! $MPI_AVAILABLE; then
    skip "All MPI tests (mpirun/mpiexec not found)"
elif [[ -z "$SEQ3D_MPI" ]]; then
    skip "All MPI tests (sequentialFilterPlt MPI executable not found)"
else
    for NP in 2 4; do
        $MPI_CMD -np $NP "$SEQ3D_MPI" infile=plt_single base_fgr=8 \
            outfile=seq_mpi${NP}_single > /dev/null 2>&1 || true
        plt_agree "T33 $NP ranks, single level matches serial" seq_mpi${NP}_single seq_single_8
        $MPI_CMD -np $NP "$SEQ3D_MPI" infile=plt_tight_3lev base_fgr=8 max_grid_size=8 \
            outfile=seq_mpi${NP}_tight amrex.init_snan=1 > /dev/null 2>&1 || true
        plt_agree "T34 $NP ranks, 3 levels matches serial" seq_mpi${NP}_tight seq_tight_8_nosnan --allow_diff_grids
    done
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
