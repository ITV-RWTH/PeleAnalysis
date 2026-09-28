#!/usr/bin/env bash
# Build every tool listed as "#EBASE = <tool>" in ./GNUmakefile (tools behind a
# double "##" are ignored) and fail on the first compilation error.
#
# Usage (from Src/ or Src/ModelSpecificAnalysis/):
#   build_tools.sh <extra make args...>
#
# Environment:
#   DIM      2 or 3 (required)
#   USE_EB   TRUE or FALSE (default FALSE)
#   BACKEND  CPU, CUDA, HIP or SYCL (default CPU)
#   NPROCS   parallel make jobs (default: nproc)
#
# Tools known not to compile for a given configuration are skipped. Each list
# below is a space-separated set of tool names; the skip list for a build is
# the union of every list that applies to it:
#   EXPECTED_<DIM>D_FAILURES           e.g. EXPECTED_2D_FAILURES
#   EXPECTED_EB_FAILURES               any DIM, when USE_EB=TRUE
#   EXPECTED_EB_<DIM>D_FAILURES        this DIM, when USE_EB=TRUE
#   EXPECTED_<BACKEND>_FAILURES        e.g. EXPECTED_CUDA_FAILURES
set +e

DIM=${DIM:?DIM must be set}
USE_EB=${USE_EB:-FALSE}
BACKEND=${BACKEND:-CPU}
NPROCS=${NPROCS:-$(nproc)}

list() { local name=$1; echo "${!name}" | tr -d "'\""; }

SKIP="$(list "EXPECTED_${DIM}D_FAILURES")"
if [ "${USE_EB}" = "TRUE" ]; then
  SKIP="${SKIP} $(list EXPECTED_EB_FAILURES) $(list "EXPECTED_EB_${DIM}D_FAILURES")"
fi
if [ "${BACKEND}" != "CPU" ]; then
  SKIP="${SKIP} $(list "EXPECTED_${BACKEND}_FAILURES")"
fi
SKIP=" $(echo ${SKIP}) "

EBASE_OPTIONS=$(grep -oP '^#(?!#)\s*EBASE\s*=\s*\K[^ ]+' GNUmakefile)

for TYPE in ${EBASE_OPTIONS}; do
  if [[ "${SKIP}" == *" ${TYPE} "* ]]; then
    echo "Skipping ${TYPE}: known unsupported for DIM=${DIM} USE_EB=${USE_EB} BACKEND=${BACKEND}."
    continue
  fi
  printf "\n-------- %s (DIM=%s, USE_EB=%s, BACKEND=%s) --------\n" "${TYPE}" "${DIM}" "${USE_EB}" "${BACKEND}"
  make -j "${NPROCS}" EBASE="${TYPE}" USE_CCACHE=TRUE DIM="${DIM}" USE_EB="${USE_EB}" "$@" 2>&1
  RESULT=$?
  make realclean > /dev/null 2>&1
  if [ "${RESULT}" -ne 0 ]; then
    echo "ERROR: Compilation failed for ${TYPE} (DIM=${DIM}, USE_EB=${USE_EB}, BACKEND=${BACKEND})."
    exit "${RESULT}"
  fi
done
