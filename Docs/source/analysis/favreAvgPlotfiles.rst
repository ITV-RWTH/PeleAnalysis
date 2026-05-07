favreAvgPlotfiles
==================

Overview
########

Utility to average AMReX plotfiles on the same domain but with non-matching AMR refinement structures. This tool computes ensemble or time-averaged statistics from multiple plotfiles, including Favre-averaged means and turbulent variances.

Usage
#####

::

   ./favreAvgPlotfiles3d.gnu.ex infiles=<file1 file2 ...> [options]

Example
#######

::

   ./favreAvgPlotfiles3d.gnu.ex infiles="plt29000 plt30000 plt31000" outfile="avg_result"

   ./favreAvgPlotfiles3d.gnu.ex infiles="plt29000 plt30000" do_variance=0 do_divide=1

Tool Options
############

::

   #------------------- IO CONTROL -----------------------------------------------------------
   infiles = plt29000 plt30000 plt31000     # List of AMReX plotfiles to average [REQUIRED]
   outfile = plt_averaged                   # DEF: plt_averaged; Output averaged plotfile name

`infiles` specifies the list of input AMReX plotfiles to be averaged. All files must have the same domain geometry. Multiple files can be provided as a space-separated list.

`outfile` specifies the name of the output plotfile containing the averaged statistics. If not provided, defaults to ``plt_averaged``.

::

   #------------------- VARIABLE SELECTION --------------------------------------------------
   variables = x_velocity y_velocity z_velocity density    # Variables to extract (DEF: all)

`variables` allows selective extraction of specific variables from the input plotfiles. If not specified, all variables present in the input files are processed. Variable names must match exactly those in the input plotfiles.

::

   #------------------- AVERAGING OPTIONS ---------------------------------------------------
   do_average = 1                           # DEF: 1; Compute ensemble/time-averaged means
   do_variance = 1                          # DEF: 1; Compute turbulent variances
   do_divide = 1                            # DEF: 1; Divide by rho for Favre averaging

These flags control what statistics are computed in the output file.

`do_average` enables computation of the mean fields. When set to 1, the ensemble/time-averaged values are computed. Set to 0 to skip mean computation.

`do_variance` enables computation of turbulent variances. When set to 1, the variance around the mean is computed. The variance is calculated using:

.. math::

   \text{var}(\phi) = \langle \phi^2 \rangle - \langle \phi \rangle^2

Set to 0 to skip variance computation. Note: At least one of `do_average` or `do_variance` must be 1.

`do_divide` controls whether to apply Favre averaging (density-weighted). When set to 1, means are divided by the mean density to produce Favre-averaged quantities:

.. math::

   \langle \phi \rangle_f = \frac{\langle \rho \phi \rangle}{\langle \rho \rangle}

When set to 0, rho-weighted means are output instead:

.. math::

   \langle \rho \phi \rangle

::

   #------------------- AMR CONTROL ---------------------------------------------------------
   output_max_level = 1000                  # DEF: 1000; Maximum refinement level to keep
   output_max_grid_size = 32                # DEF: 32; Maximum grid size in output
   interp_type = 1                          # DEF: 1; 0=piecewise const, 1=linear interp

`output_max_level` specifies the finest AMR level to include in the output (zero-indexed). Levels are refined up to this value. Setting this lower reduces memory usage and I/O for highly refined datasets. The default value is large enough to include all available levels in most cases.

`output_max_grid_size` sets the maximum grid size for boxes in the output plotfile when the input files have different grid structures. This is ignored if all input files have identical box arrays on a given level. Smaller values increase the number of boxes but may improve load balancing.

`interp_type` controls the interpolation scheme when filling patches between different AMR levels. Use 0 for piecewise constant (no interpolation) or 1 for linear interpolation (default, recommended).

Output Variables
#################

The output variable names depend on the selected averaging options:

**With do_divide=1 (Favre-averaged, default):**

When Favre averaging is enabled, the output contains:

- ``rho_mean`` — mean density :math:`\langle \rho \rangle`
- ``<variable>_favre_mean`` — Favre-averaged mean (if do_average=1): :math:`\langle \phi \rangle_f = \frac{\langle \rho \phi \rangle}{\langle \rho \rangle}`
- ``<variable>_favre_variance`` — Favre-averaged variance (if do_variance=1): :math:`\text{var}_f(\phi) = \langle \phi^2 \rangle_f - \langle \phi \rangle_f^2`

**With do_divide=0 (rho-weighted):**

When Favre averaging is disabled, the output contains rho-weighted (conservative) quantities:

- ``rho_mean`` — mean density :math:`\langle \rho \rangle`
- ``<variable>_mean`` — rho-weighted mean (if do_average=1): :math:`\langle \rho \phi \rangle`
- ``<variable>_variance`` — variance (if do_variance=1): :math:`\text{var}(\phi) = \langle \phi^2 \rangle - \langle \phi \rangle^2`

**Example output with default settings (do_average=1, do_variance=1, do_divide=1):**

::

   0   rho_mean
   1   x_velocity_favre_mean
   2   y_velocity_favre_mean
   3   z_velocity_favre_mean
   4   density_favre_mean
   5   rhoh_favre_mean
   6   temp_favre_mean
   7   x_velocity_favre_variance
   8   y_velocity_favre_variance
   9   z_velocity_favre_variance
   ...

Constraints and Notes
#####################

- All input plotfiles must have the same physical domain geometry and extent.
- All files must contain density (named ``density`` in the variable list) for Favre-averaging operations.
- At least one of `do_average` or `do_variance` must be set to 1.
- Output file size depends on the number of variables and grid refinement. Consider using `output_max_level` to reduce file size for highly refined datasets.
- For compressible turbulent flows, Favre-averaged quantities (do_divide=1) are typically more physically meaningful than rho-weighted averages. Favre averaging properly weights the statistics by density, making it the standard approach for compressible flows where :math:`\rho` varies significantly.
