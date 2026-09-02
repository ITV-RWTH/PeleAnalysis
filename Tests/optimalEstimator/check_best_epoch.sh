#!/usr/bin/env bash
#
# Best-epoch check for optimalEstimatorTraining.
#
# The tool documents that the weights written to disk are those of the best
# epoch, not the last. This asserts the observable form of that promise: the R2
# reported on the closing "Restored weights from epoch N" line must equal the
# best R2 anywhere in the epoch log.
#
# It used to fail, and only on runs that converged well. The snapshot and the
# patience counter shared one condition,
#
#     if (val_loss < best_val - min_delta_abs)      min_delta_abs = min_delta*targVar
#
# so once best_val fell below min_delta_abs the right-hand side went negative
# and no non-negative loss could ever satisfy it again: the snapshot froze for
# the rest of the run. In R2 terms that is R2 > 1 - min_delta, i.e. 0.999 by
# default. The run then always stopped exactly `patience` epochs later, having
# discarded every improvement made in between - and since the irreducible error
# is read straight off that discarded-from estimator, it came out inflated, by
# a different amount for each feature set.
#
# Case A is therefore the regression: it must reach R2 > 0.999, so that it
# crosses the threshold that used to freeze the snapshot, and it must still
# restore its best epoch. Case B fits a target against a feature that carries
# no information about it, lands at R2 near 0, never goes near the threshold,
# and must satisfy the same property - the point being that the fix is not
# meant to change what a low-R2 run does.
#
# Usage:  ./check_best_epoch.sh <training-exe> <plotfile>
# e.g.    ./check_best_epoch.sh ../../Src/optimalEstimatorTraining2d.gnu.ex plt_oe_high_r2
#
# Build the plotfile from gen_oe_high_r2.inp in this directory:
#     ../../Src/generateTestPlt2d.gnu.ex gen_oe_high_r2.inp

set -u

EXE=${1:-}
PLT=${2:-}
NEURONS=${NEURONS:-"16 16"}
NEPOCHS=${NEPOCHS:-800}
TOL=${TOL:-1e-9}

if [ -z "$EXE" ] || [ -z "$PLT" ]; then
  sed -n '2,34p' "$0"
  exit 1
fi
if [ ! -x "$EXE" ]; then echo "not executable: $EXE"; exit 1; fi
if [ ! -d "$PLT" ]; then echo "no such plotfile: $PLT"; exit 1; fi

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

STATUS=0

# $1 label, $2 features, $3 target, $4 "high"|"low"
run_case() {
  local label=$1 feat=$2 targ=$3 regime=$4
  local log="$WORK/$label.log"

  echo "--- $label: features=$feat targets=$targ ---"
  "$EXE" infile="$PLT" features="$feat" targets="$targ" neurons="$NEURONS" \
    nEpochs="$NEPOCHS" use_double=1 print_every=1 \
    model_path="$WORK/m_$label" minmax_path="$WORK/mm_$label" > "$log" 2>&1
  if ! grep -q "^Restored weights from epoch" "$log"; then
    echo "FAIL: no \"Restored weights\" line (run it by hand to see why)"
    tail -5 "$log"
    STATUS=1
    return
  fi

  awk -v tol="$TOL" -v regime="$regime" -v label="$label" '
    # Epoch [N/M], Training Loss: t, Validation Loss: v, R2: r, lr: l
    /^Epoch \[/ {
      ep = $2; gsub(/[^0-9\/]/, "", ep); split(ep, e, "/")
      for (i = 1; i <= NF; i++)
        if ($i == "R2:") {
          r = $(i+1) + 0
          if (n++ == 0 || r > best) { best = r; best_ep = e[1] + 0 }
        }
      next
    }
    # Restored weights from epoch N (validation loss v, R2 r).
    /^Restored weights from epoch/ {
      got_ep = $5 + 0
      for (i = 1; i <= NF; i++)
        if ($i == "R2") { got = $(i+1) + 0 }
      seen = 1
    }
    END {
      if (!seen) { print "FAIL: could not parse the restored line"; exit 1 }
      d = got - best; if (d < 0) d = -d
      printf "  best R2 %.10f at epoch %d; restored epoch %d, R2 %.10f\n", best, best_ep, got_ep, got

      if (regime == "high" && best <= 0.999) {
        printf "FAIL: %s never passed R2 = 0.999 (best %.10f), so it does not\n", label, best
        print  "      exercise the threshold this check exists for. Train longer"
        print  "      (NEPOCHS) or widen the network (NEURONS)."
        exit 1
      }
      if (regime == "low" && best >= 0.9) {
        printf "FAIL: %s reached R2 %.10f; it is meant to be an uninformative fit\n", label, best
        exit 1
      }
      if (d > tol) {
        printf "FAIL: restored R2 %.10f is not the best R2 %.10f (differs by %g)\n", got, best, d
        exit 1
      }
      printf "PASS: restored the best epoch (R2 %.10f)\n", got
    }' "$log" || STATUS=1
}

# A: targ is an exact function of feat - converges past R2 = 0.999, which is
#    where the snapshot used to freeze. This is the regression case.
run_case high feat targ high

# B: noise varies only in y and targ only in x, so no estimator can do better
#    than the unconditional mean: R2 ~ 0, nowhere near the threshold.
run_case low noise targ low

exit $STATUS
