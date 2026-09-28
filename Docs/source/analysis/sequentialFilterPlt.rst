.. highlight:: bash

sequentialFilterPlt
******************************************
Apply a box filter to an AMReX plot file as three successive one-dimensional
passes (first in *x*, then *y*, then *z*) rather than as a single
three-dimensional stencil. A box filter is separable, so the result is the same
filter as the one :doc:`filterPlt` applies with ``filter_type = 1``, but the
grow (ghost) region is only ever needed in one direction at a time.

This matters for memory. For a box of :math:`N^3` cells and a filter half-width
of :math:`g` cells, a single-pass 3D filter has to hold a grown box of
:math:`(N+2g)^3` cells, whereas each pass here only needs :math:`(N+2g)N^2`.
Two such arrays are cycled through as source and destination of the three
passes, so the peak footprint grows linearly rather than cubically in the filter
width. Use this tool when :doc:`filterPlt` runs out of memory at the filter
width you want; for narrow filters :doc:`filterPlt` is simpler and cheaper.

Like :doc:`filterPlt`, this tool uses the PelePhysics ``PltFileManager`` utility
to read plot files, so ``AMREX_HOME`` and ``PELE_PHYSICS_HOME`` must be defined
in the ``GNUmakefile``. For multilevel plot files the filter can be applied
either at a constant absolute filter width or at a constant filter-to-grid ratio
across levels. The filter-to-grid ratio on the base level must be even.

Only the box filter is implemented here. Requesting any other ``filter_type``
aborts rather than silently box-filtering; use :doc:`filterPlt` for the
other filter types offered by
`PelePhysics <https://amrex-combustion.github.io/PelePhysics/Utility.html#filter>`_.

This tool must be compiled with ``DIM = 3``.

Usage: ::

   ./sequentialFilterPlt3d.gnu.MPI.ex infile=FILE [OPTIONS]

Help: ::

   ./sequentialFilterPlt3d.gnu.MPI.ex -h

Example: ::

   ./sequentialFilterPlt3d.gnu.MPI.ex ./InputSamples/sequentialFilterPlt.inp

Tool Options
#############
::

   #------------------- IO CONTROL -----------------------------------------------------------
   infile  = plt00000                         # pltfile to filter
   outfile = plt00000_filtered                # DEF: "<infile>_filtered"; name of output file
   n_files = 64                               # DEF: AMReX default; cap on the number of plotfile data files (VisMF), clamped to the MPI rank count

`infile` is required. `outfile` defaults to the base name of `infile` with
``_filtered`` appended, written into the current working directory. The output is
a standard multi-level AMReX plot file carrying the geometry and simulation time
of the input.
::

   #------------------- Operation control ----------------------------------------------------
   variables = temp HeatRelease               # DEF: all variables in the file; names of the variables to filter
   max_filter_level = 4                       # DEF: 1000; max refinement level to filter, zero-indexed

`variables` selects which fields to filter and write; the names must match those
in the input plot file exactly, and an unknown name is a fatal error. Omit it to
filter every variable in the file. `max_filter_level` caps the number of AMR
levels read and written; it is clamped to the finest level actually present in
the file, so the default keeps all of them.
::

   #------------------- Filter Options -------------------------------------------------------
   filter_type = 1                            # DEF: 1; filter type as defined in PeleC. Only the box filter (1) is implemented
   base_fgr = 2                               # DEF: 2; filter-to-grid ratio on the base level, must be even
   same_fgr_all_levels = false                # DEF: false; keep the filter-to-grid ratio, rather than the absolute filter width, constant across levels

`base_fgr` is the filter-to-grid ratio on the base level and must be a positive
even number; the filter then spans ``base_fgr + 1`` cells with trapezoidal
weights (half weight on the two end cells), which sum to one. With
`same_fgr_all_levels` left at ``false`` the ratio is multiplied by the
refinement ratio on each finer level, so the filter keeps a constant *physical*
width over the whole AMR hierarchy. Setting it to ``true`` instead keeps the
ratio fixed, so the filter width shrinks with the cell size on finer levels.
::

   #------------------- Grid Options ---------------------------------------------------------
   max_grid_size = 32                         # DEF: grids of the input pltfile are used unchanged
   interp_type = 1                            # DEF: 1; interpolation used when FillPatching: 0 -> piecewise constant, 1 -> cell conservative linear

`max_grid_size` is optional here, unlike in :doc:`filterPlt`: if it is not set,
the ``BoxArray`` of the input plot file is used as is. Set it to re-box the data,
for instance to improve load balance. It must be at least as large as the filter
width on the finest level filtered, since each pass grows the boxes by half the
filter width in one direction. `interp_type` selects the interpolator used to
fill fine-level ghost cells from the underlying coarse level.

Notes
#############

Ghost cells at a non-periodic domain boundary are filled by first-order
extrapolation, so filtered values within half a filter width of such a boundary
are one-sided estimates. Periodic directions are detected from the plot file
geometry and wrap around as expected.

At coarse/fine interfaces each of the three passes fills fine-level ghost cells
by interpolating from the coarse level, which at that point already holds the
result of the preceding passes. Filtering and interpolation do not commute
exactly, so values within a filter width of a coarse/fine interface are close
to, but not bitwise identical to, what :doc:`filterPlt` produces with the same
settings. Away from level boundaries the two tools agree.

The coarse levels of the output are not averaged down from the finer levels;
each level is filtered from its own data, as in :doc:`filterPlt`.
