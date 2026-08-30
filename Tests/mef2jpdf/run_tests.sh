#!/usr/bin/env bash
# run_tests.sh — compile and run all mef2jpdf tests
#
# Usage:  ./run_tests.sh [--no-compile]
#
#   --no-compile   Skip the build step; assume executables already exist in Src/
#
# Pipeline: generateTestPlt -> isosurface -> mef2jpdf
# All test artefacts are written to Tests/mef2jpdf/testrun/.
# Python 3 (with numpy) is required for jPDF assertions.

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

# ---------------------------------------------------------------------------
# Assertion helpers
# ---------------------------------------------------------------------------

check_file() {   # check_file <desc> <path>
    if [[ -f "$2" ]]; then pass "$1"; else fail "$1 — missing: $2"; fi
}

check_dir() {    # check_dir <desc> <path>
    if [[ -d "$2" ]]; then pass "$1"; else fail "$1 — missing dir: $2"; fi
}

# check_jpdf_sum <desc> <matrix_dat> <expected> [<tol>]
check_jpdf_sum() {
    local desc="$1" mat="$2" expected="$3" tol="${4:-1e-6}"
    if [[ ! -f "$mat" ]]; then
        fail "$desc — matrix file missing: $mat"; return
    fi
    if ! command -v python3 &>/dev/null; then
        skip "$desc — python3 not available"; return
    fi
    local out
    if out=$(python3 "$SCRIPT_DIR/check_jpdf_sum.py" "$mat" "$expected" "$tol" 2>&1); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out"
    fi
}

# check_jpdf_uniform <desc> <matrix_dat> [<tol_abs>]
check_jpdf_uniform() {
    local desc="$1" mat="$2" tol="${3:-1e-10}"
    if [[ ! -f "$mat" ]]; then
        fail "$desc — matrix file missing: $mat"; return
    fi
    if ! command -v python3 &>/dev/null; then
        skip "$desc — python3 not available"; return
    fi
    local out
    if out=$(python3 "$SCRIPT_DIR/check_jpdf_uniform.py" "$mat" "$tol" 2>&1); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out"
    fi
}

# check_jpdf_marginal_uniform <desc> <matrix_dat> <axis> [<tol_frac>]
check_jpdf_marginal_uniform() {
    local desc="$1" mat="$2" axis="$3" tol="${4:-0.10}"
    if [[ ! -f "$mat" ]]; then
        fail "$desc — matrix file missing: $mat"; return
    fi
    if ! command -v python3 &>/dev/null; then
        skip "$desc — python3 not available"; return
    fi
    local out
    if out=$(python3 "$SCRIPT_DIR/check_jpdf_marginal_uniform.py" \
                "$mat" "$axis" "$tol" 2>&1); then
        pass "$desc  ($out)"
    else
        fail "$desc — $out"
    fi
}

# check_jpdf_arcsine <desc> <matrix_dat> <centers_dat> <center> <radius> [<tol_frac>]
check_jpdf_arcsine() {
    local desc="$1" mat="$2" cen_dat="$3" center="$4" radius="$5" tol="${6:-0.20}"
    if [[ ! -f "$mat" ]] || [[ ! -f "$cen_dat" ]]; then
        fail "$desc — output files missing"; return
    fi
    if ! command -v python3 &>/dev/null; then
        skip "$desc — python3 not available"; return
    fi
    local out
    if out=$(python3 "$SCRIPT_DIR/check_jpdf_arcsine.py" \
                "$mat" "$cen_dat" "$center" "$radius" "$tol" 2>&1); then
        pass "$desc"
    else
        fail "$desc — $out"
    fi
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
        fail "$label build failed — see $SCRIPT_DIR/$log"
    fi
}

if $COMPILE; then
    build_tool "generateTestPlt 3D" build_genPlt3d.log  EBASE=generateTestPlt DIM=3 USE_MPI=FALSE
    build_tool "isosurface 3D"      build_iso3d.log     EBASE=isosurface      DIM=3 USE_MPI=FALSE
    build_tool "mef2jpdf 3D"        build_mef2jpdf3d.log EBASE=mef2jpdf       DIM=3 USE_MPI=FALSE
else
    skip "Build skipped (--no-compile)"
fi

GEN3D=$(find_exe "generateTestPlt3d.gnu.ex")
ISO3D=$(find_exe "isosurface3d.gnu.ex")
MEF2JPDF=$(find_exe "mef2jpdf3d.gnu.ex")

check_file "generateTestPlt 3D" "$GEN3D"
check_file "isosurface 3D"      "$ISO3D"
check_file "mef2jpdf 3D"        "$MEF2JPDF"

# ---------------------------------------------------------------------------
# Phase 2 — Generate test plotfiles
# ---------------------------------------------------------------------------
section "Phase 2: Generate test plotfiles"

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

for cfg in gen_sphere gen_plane_z gen_sphere_const gen_cylinder; do
    echo "  Generating $cfg..."
    "$GEN3D" "$SCRIPT_DIR/$cfg.inp" > /dev/null 2>&1
done

check_dir  "plt_sphere"             "plt_sphere"
check_file "plt_sphere/Header"      "plt_sphere/Header"
check_dir  "plt_plane_z"            "plt_plane_z"
check_file "plt_plane_z/Header"     "plt_plane_z/Header"
check_dir  "plt_sphere_const"       "plt_sphere_const"
check_file "plt_sphere_const/Header" "plt_sphere_const/Header"
check_dir  "plt_cylinder"           "plt_cylinder"
check_file "plt_cylinder/Header"    "plt_cylinder/Header"

# ---------------------------------------------------------------------------
# Phase 3 — Extract isosurfaces
# ---------------------------------------------------------------------------
section "Phase 3: Extract isosurfaces"

for cfg in iso_sphere iso_plane_z iso_sphere_const iso_cylinder; do
    echo "  Extracting $cfg..."
    "$ISO3D" "$SCRIPT_DIR/$cfg.inp" > /dev/null 2>&1
done

check_file "surf_sphere.mef"       "surf_sphere.mef"
check_file "surf_plane_z.mef"      "surf_plane_z.mef"
check_file "surf_sphere_const.mef" "surf_sphere_const.mef"
check_file "surf_cylinder.mef"     "surf_cylinder.mef"

# ---------------------------------------------------------------------------
# Phase 4 — mef2jpdf tests
# ---------------------------------------------------------------------------
section "Phase 4: T1 — Normalisation (sphere, X vs Y)"

"$MEF2JPDF" "$SCRIPT_DIR/jpdf_t1.inp" > /dev/null 2>&1
check_dir      "T1 output dir"   "MEF_JPDFAverage_X_t1_Y_t1"
check_jpdf_sum "T1 sum = 1.0"    "MEF_JPDFAverage_X_t1_Y_t1/Pdf_X_t1_Y_t1.dat" 1.0

# ---------------------------------------------------------------------------
section "Phase 4: T2 — Uniform marginals (flat plane, X vs Y)"

"$MEF2JPDF" "$SCRIPT_DIR/jpdf_t2.inp" > /dev/null 2>&1
check_dir                    "T2 output dir"         "MEF_JPDFAverage_X_t2_Y_t2"
check_jpdf_sum               "T2 sum = 1.0"          "MEF_JPDFAverage_X_t2_Y_t2/Pdf_X_t2_Y_t2.dat" 1.0
check_jpdf_marginal_uniform  "T2 X-marginal uniform" "MEF_JPDFAverage_X_t2_Y_t2/Pdf_X_t2_Y_t2.dat" 1 0.15
check_jpdf_marginal_uniform  "T2 Y-marginal uniform" "MEF_JPDFAverage_X_t2_Y_t2/Pdf_X_t2_Y_t2.dat" 0 0.15

# ---------------------------------------------------------------------------
section "Phase 4: T3 — Sphere belt theorem (X vs Z, Z-marginal uniform)"

"$MEF2JPDF" "$SCRIPT_DIR/jpdf_t3.inp" > /dev/null 2>&1
check_dir                    "T3 output dir"             "MEF_JPDFAverage_X_t3_Z_t3"
check_jpdf_sum               "T3 sum = 1.0"              "MEF_JPDFAverage_X_t3_Z_t3/Pdf_X_t3_Z_t3.dat" 1.0
check_jpdf_marginal_uniform  "T3 Z-marginal uniform"     "MEF_JPDFAverage_X_t3_Z_t3/Pdf_X_t3_Z_t3.dat" 0 0.20

# ---------------------------------------------------------------------------
section "Phase 4: T4 — Constant mapped field, explicit bounds"
# const_field = 3.14, xmin=3.0 xmax=4.0 nBins=10 -> ix=1 gets all weight

"$MEF2JPDF" "$SCRIPT_DIR/jpdf_t4.inp" > /dev/null 2>&1
check_dir      "T4 output dir"   "MEF_JPDFAverage_cf_t4_Z_t4"
check_jpdf_sum "T4 sum = 1.0"    "MEF_JPDFAverage_cf_t4_Z_t4/Pdf_cf_t4_Z_t4.dat" 1.0

# const_field=3.14 -> ix=1 (bin [3.1, 3.2)); all other rows must be zero;
# within row 1 the Z distribution must be uniform (belt theorem, 20% tol)
if command -v python3 &>/dev/null && [[ -f "MEF_JPDFAverage_cf_t4_Z_t4/Pdf_cf_t4_Z_t4.dat" ]]; then
    if out=$(python3 - <<'EOF'
import sys, numpy as np
m = np.loadtxt("MEF_JPDFAverage_cf_t4_Z_t4/Pdf_cf_t4_Z_t4.dat")
other = np.delete(m, 1, axis=0)   # rows != ix=1
leak = float(other.sum())
row1 = m[1, :]
expected = 1.0 / m.shape[1]
max_rel = float(np.max(np.abs(row1 - expected)) / expected)
if leak > 1e-6:
    print(f"FAIL leak={leak:.2e}")
    sys.exit(1)
if max_rel > 0.20:
    print(f"FAIL row1 non-uniform max_rel={max_rel:.3f}")
    sys.exit(1)
print(f"OK leak={leak:.2e} Z-uniform max_rel={max_rel:.3f}")
EOF
    ); then
        pass "T4 const_field concentrates in one bin, Z-marginal uniform  ($out)"
    else
        fail "T4 const_field check — $out"
    fi
fi

# ---------------------------------------------------------------------------
section "Phase 4: T5 — Multi-file (sphere MEF fed twice)"

"$MEF2JPDF" "$SCRIPT_DIR/jpdf_t5.inp" > /dev/null 2>&1
check_dir      "T5 output dir"   "MEF_JPDFAverage_X_t5_Y_t5"
check_jpdf_sum "T5 sum = 1.0"    "MEF_JPDFAverage_X_t5_Y_t5/Pdf_X_t5_Y_t5.dat" 1.0

# ---------------------------------------------------------------------------
section "Phase 4: T6 — Arcsine distribution (flat plane, sine_x = sin(2πX))"
# sine_x = sin(2pi*X) on a flat plane; X uniform -> sine_x follows arcsine distribution

"$MEF2JPDF" "$SCRIPT_DIR/jpdf_t6.inp" > /dev/null 2>&1
check_dir        "T6 output dir"   "MEF_JPDFAverage_X_t6_Z_t6"
check_jpdf_sum   "T6 sum = 1.0"    "MEF_JPDFAverage_X_t6_Z_t6/Pdf_X_t6_Z_t6.dat" 1.0
check_jpdf_arcsine \
    "T6 sine_x-marginal follows arcsine (c=0, R=1)" \
    "MEF_JPDFAverage_X_t6_Z_t6/Pdf_X_t6_Z_t6.dat" \
    "MEF_JPDFAverage_X_t6_Z_t6/Pdf_X_t6_x.dat" \
    0.0 1.0 0.15

# ---------------------------------------------------------------------------
section "Phase 4: T7 — Sphere (X,Y) joint distribution is non-uniform"
# The sphere (X,Y) jPDF is non-uniform: highest density near the equatorial
# ring (outer projection), lowest at the center (polar projection).
# Check: symmetry about both axes and outer bin > center bin.

"$MEF2JPDF" "$SCRIPT_DIR/jpdf_t7.inp" > /dev/null 2>&1
check_dir      "T7 output dir"   "MEF_JPDFAverage_X_t7_Y_t7"
check_jpdf_sum "T7 sum = 1.0"    "MEF_JPDFAverage_X_t7_Y_t7/Pdf_X_t7_Y_t7.dat" 1.0

if command -v python3 &>/dev/null && [[ -f "MEF_JPDFAverage_X_t7_Y_t7/Pdf_X_t7_Y_t7.dat" ]]; then
    if out=$(python3 - <<'EOF'
import sys, numpy as np
m = np.loadtxt("MEF_JPDFAverage_X_t7_Y_t7/Pdf_X_t7_Y_t7.dat")
N = m.shape[0]
mid = N // 2

# Non-uniformity: outer X-bin at center Y should exceed center X-bin
# (sphere density diverges at equatorial ring, minimum at poles)
outer = float(m[0, mid])
center = float(m[mid, mid])
if outer <= center:
    print(f"FAIL outer={outer:.5f} <= center={center:.5f}")
    sys.exit(1)

# Left-right symmetry of X-marginal (sphere is symmetric about X=0.5)
mx = m.sum(axis=1)   # X-marginal: sum over Y
sym_err = np.max(np.abs(mx - mx[::-1]) / (mx + mx[::-1] + 1e-15))
if sym_err > 0.05:
    print(f"FAIL X-marginal asymmetry={sym_err:.3f}")
    sys.exit(1)

print(f"OK outer={outer:.5f} > center={center:.5f}; X-sym={sym_err:.4f}")
EOF
    ); then
        pass "T7 sphere non-uniform: outer > center, symmetric  ($out)"
    else
        fail "T7 sphere non-uniform check — $out"
    fi
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
