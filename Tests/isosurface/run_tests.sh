#!/usr/bin/env bash
# run_tests.sh — compile and run the isosurface tests (2D and 3D, serial)
#
# Usage:  ./run_tests.sh [--no-compile]
#
#   --no-compile   Skip the build step; assume executables already exist in Src/
#
# Iso lines and surfaces that cross periodic boundaries are measured against a
# reference taken without periodicity, see TESTING.md.  All test artefacts are
# written to Tests/isosurface/testrun/ so the source tree stays clean.  Python 3
# is required (check_mef.py).  ISO2D and ISO3D may be set to test other
# isosurface executables.

set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SRC_DIR="$SCRIPT_DIR/../../Src"
RUN_DIR="$SCRIPT_DIR/testrun"
NPROC=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
CHECK="python3 $SCRIPT_DIR/check_mef.py"

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

# The answers, derived in TESTING.md:
#   L_exact = length of the closed line cos(2 pi x) cos(2 pi y) = 0.5
#   periodic 2D: two lines congruent to the reference line   -> 2 L_ref
#   periodic 3D: the 2D lines extruded over z in [0, 0.25]     -> 0.5 L_ref
L_EXACT=1.0810133978566598
DX=0.015625   # finest cell size in all tests

# ---------------------------------------------------------------------------
# Assertion helpers
# ---------------------------------------------------------------------------
check_file() {   # check_file <desc> <path>
    if [[ -f "$2" ]]; then pass "$1"; else fail "$1 — missing: $2"; fi
}

# measure_ok <desc> <mef> <expected> <rtol> <check_mef.py domain args...>
# Measure, node and edge checks of check_mef.py.
measure_ok() {
    local desc="$1" mef="$2" exp="$3" rtol="$4" out; shift 4
    if [[ ! -f "$mef" ]]; then fail "$desc — missing: $mef"; return; fi
    if out=$($CHECK "$mef" "$@" --dx "$DX" --expect "$exp" --rtol "$rtol"); then
        pass "$desc  ($(echo "$out" | tr '\n' ';' | sed 's/;$//; s/;/; /g'))"
    else
        fail "$desc — $(echo "$out" | tr '\n' ';' | sed 's/;$//; s/;/; /g')"
    fi
}

# area_ok <desc> <log> <expected> <rtol>: the area printed by computeArea
area_ok() {
    local desc="$1" log="$2" exp="$3" rtol="$4" a
    a=$(awk '/^Total area = / { print $4 }' "$log")
    if [[ -n "$a" ]] && awk -v a="$a" -v e="$exp" -v t="$rtol" \
        'BEGIN { d = (a - e) / e; if (d < 0) d = -d; exit !(d <= t) }'; then
        pass "$desc  (Total area = $a)"
    else
        fail "$desc — Total area = ${a:-missing}, expected $exp"
    fi
}

# iso <log/outfile base> <exe> <plotfile> <is_per...>
iso() {
    local base="$1" exe="$2" plt="$3"; shift 3
    "$exe" infile="$plt" isoCompName=c isoVal=0.5 computeArea=1 writeSurf=1 \
        surfFormat=MEF outfile_base="$base" is_per="$*" > "$base.log" 2>&1
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

if $COMPILE; then
    build_tool "generateTestPlt 2D" build_genPlt2d.log \
        EBASE=generateTestPlt DIM=2 USE_MPI=FALSE
    build_tool "generateTestPlt 3D" build_genPlt3d.log \
        EBASE=generateTestPlt DIM=3 USE_MPI=FALSE
    build_tool "isosurface 2D" build_iso2d.log \
        EBASE=isosurface DIM=2 USE_MPI=FALSE
    build_tool "isosurface 3D" build_iso3d.log \
        EBASE=isosurface DIM=3 USE_MPI=FALSE
else
    skip "Build skipped (--no-compile)"
fi

GEN2D=$(find_exe "generateTestPlt2d.gnu.ex")
GEN3D=$(find_exe "generateTestPlt3d.gnu.ex")
ISO2D=${ISO2D:-$(find_exe "isosurface2d.gnu.ex")}
ISO3D=${ISO3D:-$(find_exe "isosurface3d.gnu.ex")}

check_file "generateTestPlt 2D" "$GEN2D"
check_file "generateTestPlt 3D" "$GEN3D"
check_file "isosurface 2D"      "$ISO2D"
check_file "isosurface 3D"      "$ISO3D"

# ---------------------------------------------------------------------------
# Phase 2 — Test plotfiles
# ---------------------------------------------------------------------------
section "Phase 2: Test plotfiles"

rm -rf "$RUN_DIR"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

for name in ref_2d per_2d per_amr_2d; do
    "$GEN2D" "$SCRIPT_DIR/gen_$name.inp" > "gen_$name.log" 2>&1 || true
    check_file "plt_iso_$name" "plt_iso_$name/Header"
done
for name in per_3d per_amr_3d; do
    "$GEN3D" "$SCRIPT_DIR/gen_$name.inp" > "gen_$name.log" 2>&1 || true
    check_file "plt_iso_$name" "plt_iso_$name/Header"
done

# ---------------------------------------------------------------------------
# Phase 3 — Reference line without periodicity
# ---------------------------------------------------------------------------
section "Phase 3: Reference iso line (non-periodic window)"

iso ref_2d "$ISO2D" plt_iso_ref_2d 0 0 || fail "isosurface on plt_iso_ref_2d failed"
L_REF=$($CHECK ref_2d.mef --lo 0.25 0.25 --hi 0.75 0.75 --periodic 0 0 --dx "$DX" \
        | awk '/^measure/ { print $3 }')
echo "  L_ref = $L_REF"
# Marching squares cuts chords: L_ref is slightly shorter than the exact line
measure_ok "Reference line within 0.2% of the exact length" ref_2d.mef \
    "$L_EXACT" 2e-3 --lo 0.25 0.25 --hi 0.75 0.75 --periodic 0 0

TWO_L=$(awk -v l="$L_REF" 'BEGIN { printf "%.17g", 2 * l }')
EXTRUDED=$(awk -v l="$L_REF" 'BEGIN { printf "%.17g", 0.5 * l }')

# ---------------------------------------------------------------------------
# Phase 4 — Periodic iso lines (2D)
# ---------------------------------------------------------------------------
section "Phase 4: Periodic iso lines (2D)"

for name in per_2d per_amr_2d; do
    iso "$name" "$ISO2D" "plt_iso_$name" 1 1 || fail "isosurface on plt_iso_$name failed"
    measure_ok "$name: two lines of length L_ref" "$name.mef" "$TWO_L" 1e-9 \
        --lo 0 0 --hi 1 1 --periodic 1 1
    area_ok "$name: computeArea" "$name.log" "$TWO_L" 1e-8
done

# ---------------------------------------------------------------------------
# Phase 5 — Periodic iso surfaces (3D)
# ---------------------------------------------------------------------------
section "Phase 5: Periodic iso surfaces (3D)"

for name in per_3d per_amr_3d; do
    iso "$name" "$ISO3D" "plt_iso_$name" 1 1 1 || fail "isosurface on plt_iso_$name failed"
    measure_ok "$name: lines of length 2 L_ref times 0.25" "$name.mef" "$EXTRUDED" 1e-9 \
        --lo 0 0 0 --hi 1 1 0.25 --periodic 1 1 1
    area_ok "$name: computeArea" "$name.log" "$EXTRUDED" 1e-8
done

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
