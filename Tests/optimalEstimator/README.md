# optimalEstimator checks

Not a full suite in the sense of `Tests/jpdf` — CI is build-only and has no
libtorch, so nothing here runs automatically. These are checks to run by hand
after touching the training tool.

## `check_best_epoch.sh`

Asserts that the network written to disk is the one from the best epoch: the R²
on the closing `Restored weights from epoch N` line must equal the best R²
anywhere in the epoch log.

```sh
# build both executables first
cd ../../Src
make EBASE=generateTestPlt DIM=2
make EBASE=optimalEstimatorTraining DIM=2

# generate the plotfile, then run the check
cd ../Tests/optimalEstimator
../../Src/generateTestPlt2d.gnu.ex gen_oe_high_r2.inp
./check_best_epoch.sh ../../Src/optimalEstimatorTraining2d.gnu.ex plt_oe_high_r2
```

It runs two fits against that plotfile:

- **high-R²** — `targ` is an exact deterministic function of `feat` (both are
  functions of `x`, and `feat` is strictly monotone in `x` over the whole
  domain), so the irreducible error is zero and the fit runs up past
  R² = 0.999. The check fails if it does not get there, because the point is to
  cross that threshold;
- **low-R²** — `noise` varies only in `y` and `targ` only in `x`, so no
  estimator beats the unconditional mean and R² stays near 0.

Both must restore their best epoch. This is a regression check: the snapshot and
the early-stopping counter used to share the condition

```
if (val_loss < best_val - min_delta_abs)      // min_delta_abs = min_delta * targVar
```

so once `best_val` fell below `min_delta_abs` the threshold went negative and no
non-negative loss could satisfy it again — the snapshot froze permanently. In R²
terms that is R² > 1 − `min_delta`, 0.999 by default. A converged run then
discarded every improvement made over its last `patience` epochs and wrote out a
worse estimator, which inflated the irreducible error it exists to measure — and
by a different amount for each feature set, so it corrupted the comparison
between them too. The high-R² case reproduces that; the low-R² case is there to
show the two regimes are not being conflated.

Environment overrides: `NEURONS`, `NEPOCHS`, `TOL`.

## `check_mpi_invariance.sh`

Asserts that the training trajectory does not depend on the MPI rank count.

```sh
# build an MPI training executable first
cd ../../Src
make EBASE=optimalEstimatorTraining DIM=2 USE_MPI=TRUE

# then, from a directory holding a MULTI-LEVEL plotfile
cd /path/to/run
.../Tests/optimalEstimator/check_mpi_invariance.sh \
    .../Src/optimalEstimatorTraining2d.gnu.MPI.ex plt_test
```

It runs 1, 2 and 4 ranks with `batch_size=-1` (so each rank's single mini-batch
is its whole local training set) and `use_double=1`, then requires every epoch's
Training Loss, Validation Loss and R² to agree to 1e-10 relative.

Two properties make that identity hold, and the check exists to defend both:

- each rank normalises its loss by the mini-batch weight summed over **all**
  ranks, so the gradients are shares that sum to the exact global gradient;
- the train/validation split is drawn over the global box numbering with a fixed
  seed, so the training set does not move with the rank count.

The check is only meaningful on a **multi-level** plotfile. On a single level
every cell carries weight 1, the ranks are automatically balanced in weight, and
the identity holds even under the older, incorrect scheme.

It must be run on a plotfile whose refinement is uneven enough that ranks end up
with different shares of coarse and fine cells — which is the normal case, since
the distribution mapping is built independently per level. Generate one with
`Src/generateTestPlt.cpp` using `amr.max_level = 2` and a refinement box covering
part of the domain.

Environment overrides: `FEATURES`, `TARGETS`, `NEURONS`, `RANKS`, `TOL`.

Note that convergence curves are *not* a substitute for this check: Adam's
per-parameter normalisation is invariant to a constant rescaling of the gradient,
so a mis-weighted reduction can be invisible in the loss history while still
converging on the wrong objective.
