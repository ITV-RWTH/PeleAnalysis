#!/usr/bin/env bash
# run_tests.sh — compile and run all repackPlt tests (3D serial and MPI)
#
# Usage:  ./run_tests.sh [--no-compile]
#
#   --no-compile   Skip the build step; assume executables already exist in Src/
#
# All test artefacts are written to Tests/repackPlt/testrun/ so the source
# tree stays clean.  The multi-file input plotfile and the MPI tests require
# mpirun/mpiexec; without it only the serial tests run.  Python 3 with numpy
# is required for the data comparison (compare_plt.py).

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$SCRIPT_DIR/../../Src"
RUN_DIR="$SCRIPT_DIR/testrun"
NPROC=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
CMP="$SCRIPT_DIR/compare_plt.py"

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
check_file() {   # check_file <desc> <path>
    if [[ -f "$2" ]]; then pass "$1"; else fail "$1 — missing: $2"; fi
}

no_dir() {       # no_dir <desc> <path>
    if [[ ! -d "$2" ]]; then pass "$1"; else fail "$1 — unexpected dir: $2"; fi
}

# check_same <desc> <plt_a> <plt_b> [compare_plt.py options]
check_same() {
    local desc="$1" a="$2" b="$3"; shift 3
    if [[ ! -d "$b" ]]; then fail "$desc — plotfile missing: $b"; return; fi
    local out
    if out=$(python3 "$CMP" "$a" "$b" "$@" 2>&1); then
        pass "$desc"
    else
        fail "$desc"; echo "$out" | sed 's/^/        /'
    fi
}

# check_nfiles <desc> <plt> <expected per level>
check_nfiles() {
    local desc="$1" plt="$2" exp="$3" ok=true lev n
    for lev in "$plt"/Level_*; do
        n=$(ls "$lev" | grep -c '^Cell_D_' || true)
        [[ "$n" == "$exp" ]] || { ok=false; echo "        $lev: $n files"; }
    done
    if $ok; then pass "$desc  ($exp Cell_D file(s) per level)"; else fail "$desc"; fi
}

# check_identical <desc> <file_a> <file_b>
check_identical() {
    if cmp -s "$2" "$3"; then pass "$1"; else fail "$1 — $2 and $3 differ"; fi
}

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

if command -v mpirun &>/dev/null; then
    MPI_AVAILABLE=true; MPI_CMD=mpirun
elif command -v mpiexec &>/dev/null; then
    MPI_AVAILABLE=true; MPI_CMD=mpiexec
else
    MPI_AVAILABLE=false
fi
# Open MPI refuses to run as root unless told otherwise (containers, CI)
export OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
MPI_FLAGS=""
if $MPI_AVAILABLE && $MPI_CMD --help 2>&1 | grep -q -- "--oversubscribe"; then
    MPI_FLAGS="--oversubscribe"
fi

if $COMPILE; then
    build_tool "repackPlt 3D serial"  build_repack3d.log    EBASE=repackPlt DIM=3 USE_MPI=FALSE
    if $MPI_AVAILABLE; then
        build_tool "generateTestPlt 3D MPI" build_gen3dmpi.log EBASE=generateTestPlt DIM=3 USE_MPI=TRUE
        build_tool "repackPlt 3D MPI"    build_repack3dmpi.log EBASE=repackPlt DIM=3 USE_MPI=TRUE
    else
        build_tool "generateTestPlt 3D"  build_gen3d.log    EBASE=generateTestPlt DIM=3 USE_MPI=FALSE
    fi
else
    skip "Build skipped (--no-compile)"
fi

REPACK=$(find_exe "repackPlt3d.gnu.ex")
REPACK_MPI=$(find_exe "repackPlt3d.gnu.MPI.ex")
GEN_MPI=$(find_exe "generateTestPlt3d.gnu.MPI.ex")
GEN=$(find_exe "generateTestPlt3d.gnu.ex")

check_file "repackPlt 3D serial" "$REPACK"
if ! $MPI_AVAILABLE; then
    skip "mpirun/mpiexec not found — MPI tests will be skipped"
fi
if ! python3 -c "import numpy" &>/dev/null; then
    echo "python3 with numpy is required"; exit 1
fi

# ---------------------------------------------------------------------------
# Phase 2 — Generate input plotfile (one Cell_D file per rank and level)
# ---------------------------------------------------------------------------
section "Phase 2: Generate test plotfile"

rm -rf "$RUN_DIR"; mkdir -p "$RUN_DIR"; cd "$RUN_DIR"

NIN=1
if $MPI_AVAILABLE && [[ -n "$GEN_MPI" ]]; then
    NIN=4
    $MPI_CMD $MPI_FLAGS -np 4 "$GEN_MPI" "$SCRIPT_DIR/generate_plt_repack.inp" > gen.log 2>&1
    $MPI_CMD $MPI_FLAGS -np 4 "$GEN_MPI" "$SCRIPT_DIR/generate_plt_repack.inp" \
        plotfile_name=plt_repack_v4 vismf.headerversion=4 > gen_v4.log 2>&1
    $MPI_CMD $MPI_FLAGS -np 4 "$GEN_MPI" "$SCRIPT_DIR/generate_plt_repack.inp" \
        plotfile_name=plt_repack_f32 fab.format=NATIVE_32 > gen_f32.log 2>&1
else
    "$GEN" "$SCRIPT_DIR/generate_plt_repack.inp" > gen.log 2>&1
    "$GEN" "$SCRIPT_DIR/generate_plt_repack.inp" \
        plotfile_name=plt_repack_v4 vismf.headerversion=4 > gen_v4.log 2>&1
    "$GEN" "$SCRIPT_DIR/generate_plt_repack.inp" \
        plotfile_name=plt_repack_f32 fab.format=NATIVE_32 > gen_f32.log 2>&1
fi
check_file "plt_repack/Header" "plt_repack/Header"
check_nfiles "input written with $NIN file(s) per level" plt_repack $NIN
# Extra content that must be carried over verbatim
mkdir -p plt_repack/extra_dir && echo "keep me" > plt_repack/extra_dir/note.txt
echo "level note" > plt_repack/Level_0/note.txt

# ---------------------------------------------------------------------------
# Phase 3 — Serial repack
# ---------------------------------------------------------------------------
section "Phase 3: Serial repack"

echo "  Test 1 — n_files=1 (default), default outdir repack_nf1"
"$REPACK" infiles=plt_repack > t1.log 2>&1 || fail "T1 tool exit status"
check_nfiles    "T1 one file per level"           repack_nf1/plt_repack 1
check_same      "T1 data identical"               plt_repack repack_nf1/plt_repack --max-files 1
check_identical "T1 Header copied verbatim"       plt_repack/Header repack_nf1/plt_repack/Header
check_identical "T1 extra dir copied"             plt_repack/extra_dir/note.txt repack_nf1/plt_repack/extra_dir/note.txt
check_identical "T1 extra level file copied"      plt_repack/Level_0/note.txt repack_nf1/plt_repack/Level_0/note.txt

echo "  Test 2 — header version 4 is preserved (user outdir)"
"$REPACK" infiles=plt_repack_v4 outdir=out_t2 > t2.log 2>&1 || fail "T2 tool exit status"
if [[ "$(head -1 out_t2/plt_repack_v4/Level_0/Cell_H 2>/dev/null)" == "4" ]]; then
    pass "T2 Cell_H version 4"; else fail "T2 Cell_H version"; fi
check_same      "T2 data identical"               plt_repack_v4 out_t2/plt_repack_v4 --max-files 1

echo "  Test 3 — 32-bit input stays 32-bit"
"$REPACK" infiles=plt_repack_f32 outdir=out_t3 > t3.log 2>&1 || fail "T3 tool exit status"
if head -c 60 out_t3/plt_repack_f32/Level_0/Cell_D_00000 2>/dev/null | grep -q "(4, (4 3 2 1))"; then
    pass "T3 output is 32-bit"; else fail "T3 output precision"; fi
check_same      "T3 data identical"               plt_repack_f32 out_t3/plt_repack_f32

echo "  Test 4 — float_type=32 converts 64-bit data"
"$REPACK" infiles=plt_repack outdir=out_t4 float_type=32 > t4.log 2>&1 || fail "T4 tool exit status"
if head -c 60 out_t4/plt_repack/Level_0/Cell_D_00000 2>/dev/null | grep -q "(4, (4 3 2 1))"; then
    pass "T4 output is 32-bit"; else fail "T4 output precision"; fi
check_same      "T4 data within float32 round-off" plt_repack out_t4/plt_repack --tol 1e-5

echo "  Test 5 — batch of three plotfiles into one outdir"
"$REPACK" infiles="plt_repack plt_repack_v4/ plt_repack_f32" outdir=out_t5/ > t5.log 2>&1 \
    || fail "T5 tool exit status"
for p in plt_repack plt_repack_v4 plt_repack_f32; do
    check_nfiles "T5 $p one file per level" out_t5/$p 1
    check_same   "T5 $p data identical"     $p out_t5/$p --max-files 1
done
n=$(ls out_t5 2>/dev/null | wc -l)
if [[ "$n" == "3" ]]; then pass "T5 outdir holds exactly 3 plotfiles"; else fail "T5 outdir holds $n entries"; fi

# ---------------------------------------------------------------------------
# Phase 4 — Error handling (validation happens before anything is written)
# ---------------------------------------------------------------------------
section "Phase 4: Error handling"

if "$REPACK" infiles=plt_repack outdir=. > t6.log 2>&1; then
    fail "T6 outdir = folder of the input — expected abort"; else pass "T6 outdir = input folder aborted"; fi
check_same "T6 input untouched" plt_repack repack_nf1/plt_repack
mkdir -p dup && cp -r plt_repack_v4 dup/plt_repack
if "$REPACK" infiles="plt_repack dup/plt_repack" outdir=out_t7 > t7.log 2>&1; then
    fail "T7 duplicate names — expected abort"; else pass "T7 duplicate names aborted"; fi
no_dir "T7 no output written" out_t7
if "$REPACK" infiles="plt_repack no_such_plt" outdir=out_t8 > t8.log 2>&1; then
    fail "T8 missing input in batch — expected abort"; else pass "T8 missing input in batch aborted"; fi
no_dir "T8 no output written for the valid input either" out_t8
if "$REPACK" n_files=2 > t9.log 2>&1; then
    fail "T9 no infiles — expected abort"; else pass "T9 no infiles aborted"; fi
if "$REPACK" infiles=plt_repack n_files=0 outdir=out_t10 > t10.log 2>&1; then
    fail "T10 n_files=0 — expected abort"; else pass "T10 n_files=0 aborted"; fi
no_dir "T10 no output written" out_t10

# ---------------------------------------------------------------------------
# Phase 5 — MPI repack
# ---------------------------------------------------------------------------
section "Phase 5: MPI repack"

if ! $MPI_AVAILABLE || [[ -z "$REPACK_MPI" ]]; then
    skip "All MPI tests (mpirun or MPI executable not found)"
else
    for NF in 1 2 4; do
        echo "  Test — 4 ranks, n_files=$NF"
        $MPI_CMD $MPI_FLAGS -np 4 "$REPACK_MPI" infiles=plt_repack n_files=$NF \
            > mpi_nf$NF.log 2>&1 || fail "T-MPI nf$NF exit status"
        check_nfiles "T-MPI nf$NF file count" repack_nf$NF/plt_repack $NF
        check_same   "T-MPI nf$NF data identical" plt_repack repack_nf$NF/plt_repack --max-files $NF
    done

    echo "  Test — n_files larger than the rank count is capped"
    $MPI_CMD $MPI_FLAGS -np 2 "$REPACK_MPI" infiles=plt_repack n_files=16 \
        outdir=mpi_cap > mpi_cap.log 2>&1 || fail "T-MPI cap exit status"
    check_nfiles "T-MPI cap: 2 ranks -> 2 files" mpi_cap/plt_repack 2
    if grep -q "WARNING: n_files=16" mpi_cap.log; then pass "T-MPI cap warning printed"; else fail "T-MPI cap warning"; fi

    echo "  Test — round trip 1 file -> 4 files"
    $MPI_CMD $MPI_FLAGS -np 4 "$REPACK_MPI" infiles=repack_nf1/plt_repack n_files=4 \
        outdir=mpi_rt > mpi_rt.log 2>&1 || fail "T-MPI round trip exit status"
    check_nfiles "T-MPI round trip file count" mpi_rt/plt_repack 4
    check_same   "T-MPI round trip data identical" plt_repack mpi_rt/plt_repack

    echo "  Test — batch of three plotfiles, 4 ranks, n_files=2"
    $MPI_CMD $MPI_FLAGS -np 4 "$REPACK_MPI" infiles="plt_repack plt_repack_v4 plt_repack_f32" \
        n_files=2 outdir=mpi_batch > mpi_batch.log 2>&1 || fail "T-MPI batch exit status"
    for p in plt_repack plt_repack_v4 plt_repack_f32; do
        check_nfiles "T-MPI batch $p file count"     mpi_batch/$p 2
        check_same   "T-MPI batch $p data identical" $p mpi_batch/$p --max-files 2
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
