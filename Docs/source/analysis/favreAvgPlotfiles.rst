.. highlight:: bash

favreAvgPlotfiles
*****************

Average multiple AMReX plot files defined on the same physical domain and write the result as a new plot file. The input plot files may have non-matching AMR box layouts. In that case, the tool constructs a combined AMR grid, fill-patches each input plot file onto this common grid, and accumulates mean quantities.

By default, the tool writes Favre-averaged quantities. It can also be used to write conservative mean quantities based on ``rho * variable``.

Usage: ::

   ./favreAvgPlotfiles.gnu.MPI.ex infiles=<s1 s2 s3> [options]

Example: ::

   ./favreAvgPlotfiles.gnu.MPI.ex infiles="plt00010 plt00020 plt00030" outfile=plt_favre

Tool Options
#############

::

   #------------------- IO CONTROL -----------------------------------------------------------
   infiles = plt00010 plt00020 plt00030       # Input AMReX plot files
   outfile = plt_averaged                     # DEF: plt_averaged; Output plot file name

`infiles` specifies the input AMReX plot files to combine and average. This argument is required. All input plot files must represent the same physical domain and compatible AMR hierarchy.

`outfile` specifies the name of the output plot file. If `outfile` is not provided, the output is named ``plt_averaged`` by default.

::

   #------------------- Variable Selection ---------------------------------------------------
   variables = density temp x_velocity        # DEF: all variables; Variables to include in output

`variables` specifies the variables to read from the input plot files and write to the averaged output. If `variables` is not provided, all variables are used. If a variable list is provided, each listed variable must be present in every input plot file.

::

   #------------------- AMR Control ----------------------------------------------------------
   output_max_level = 1000                    # DEF: 1000; Maximum AMR level to include, zero-indexed
   output_max_grid_size = 32                  # DEF: 32; Maximum grid size for combined output grids

`output_max_level` sets the maximum refinement level to include in the averaged output. The value is zero-indexed. For example, ``output_max_level = 0`` writes only the base level, while ``output_max_level = 1`` writes levels 0 and 1.

`output_max_grid_size` controls the maximum grid size used when the tool has to construct a new combined BoxArray because the input plot files do not have identical AMR box layouts. If all input BoxArrays are identical on a level, this option is ignored for that level.

::

   #------------------- FillPatch / Interpolation --------------------------------------------
   interp_type = 1                            # DEF: 1; 0 = piecewise constant, 1 = cell conservative linear

`interp_type` controls the interpolation used when fill-patching data from each input plot file onto the combined output grid. Use ``interp_type = 0`` for piecewise constant interpolation and ``interp_type = 1`` for cell conservative linear interpolation.

::

   #------------------- Averaging Mode -------------------------------------------------------
   favre_average = 1                          # DEF: 1; 0 = conservative mean output, 1 = Favre-average output

`favre_average` controls whether the accumulated quantities are normalized by ``rho_mean``.

If ``favre_average = 1``, the tool outputs Favre-averaged quantities. The Favre mean is computed as:

::

   mean(rho * variable) / mean(rho)

If ``favre_average = 0``, the tool outputs conservative mean quantities:

::

   mean(rho * variable)
   mean(rho * variable^2)

Output Variables
################

The output plot file always contains:

::

   rho_mean

For ``favre_average = 1``:

::

   <variable>_favre_mean
   <variable>_favre_variance

For ``favre_average = 0``:

::

   rho_<variable>_mean
   rho_<variable>_2_mean
