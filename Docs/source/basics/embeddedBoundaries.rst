.. highlight:: bash

.. _embedded_boundaries:

Embedded boundaries (EB)
************************

Plotfiles written by EB simulations (e.g. PeleLMeX or PeleC with
``USE_EB = TRUE``) contain cells that are fully or partially covered by the
geometry. Tools with EB support take the geometry into account so that
covered cells do not pollute their results. This page describes what is common
to all of them; the tool pages only list what is specific to each tool.

Building with EB support
########################

EB support is a compile-time option. Build the tool with ``USE_EB=TRUE``: ::

   make -j EBASE=grad DIM=3 USE_EB=TRUE

This compiles the AMReX EB sources, adds ``Src/EBUserDefined`` and
``Tools/EB`` to the include path and defines ``AMREX_USE_EB``. A tool built
without ``USE_EB=TRUE`` ignores the geometry completely and treats covered
cells as regular fluid cells.
``USE_EB`` can be combined with ``DIM=2``/``DIM=3`` and with the GPU back ends
(``USE_CUDA``, ``USE_HIP``, ``USE_SYCL``); see :ref:`eb_gpu_support` for what
actually runs on a GPU.

.. note::

   AMReX places EB and non-EB objects in the same build directory
   (``tmp_build_dir/o/3d.gnu.EXE``). Run ``make realclean`` when switching
   ``USE_EB`` on or off, otherwise stale objects from the other build are
   linked in.

Supported tools
###############

The tools use the geometry in one of two ways: they either rebuild it from the
input file, or they read the volume fraction ``volFrac`` that the simulation
wrote into the plotfile.

.. list-table::
   :header-rows: 1
   :widths: 22 48 15 15

   * - Tool
     - What changes with EB
     - Rebuilds geometry (``eb2.*``)
     - Needs ``volFrac`` in plotfile
   * - :doc:`../analysis/grad`
     - Gradients from an EB-aware operator (``MLEBABecLap``), evaluated at face
       centroids and averaged to cell centres.
     - yes
     - no
   * - :doc:`../analysis/curvature`
     - Progress variable is zeroed in covered cells; gradients, smoothing,
       mean/Gaussian curvature, strain rate and normal velocity use EB-aware
       operators; results are only computed where ``volFrac > 0``.
     - yes
     - no
   * - :doc:`../analysis/integral`
     - Every integrand is weighted by the cell volume fraction, so cut cells
       contribute their partial volume/area/length.
     - no
     - yes
   * - :doc:`../analysis/isosurface`
     - Fully covered cells (``volFrac == 0``) are masked out of the isosurface
       extraction.
     - no
     - yes

All other tools ignore EB, even when built with ``USE_EB=TRUE``.

Rebuilding the geometry
#######################

Tools that need the full EB geometry (cut-cell volume and area fractions, face
centroids) rebuild it at start-up, using the grids of the input plotfile. The
geometry must be the same one that produced the plotfile, so reuse the EB
inputs of the simulation.

::

   #------------------- EB geometry ---------------------------------------------------------
   eb2.geom_type = UserDefined            # Required: AMReX EB2 geometry type or "UserDefined"
   # eb2.max_level_generation = 2         # DEF: finestLevel; must equal finestLevel

`eb2.geom_type` selects the geometry and must always be given. Any geometry
type known to AMReX ``EB2::Build`` (for example ``box``, ``cylinder``,
``sphere``, ``plane`` or ``all_regular``) can be used, together with the
corresponding ``eb2.*`` parameters of that type (see the AMReX documentation).
An unknown type aborts the run. ``UserDefined`` calls the user-provided
geometry described below.

`eb2.max_level_generation` is the AMR level on which the geometry is generated;
the coarser levels are obtained by coarsening. It defaults to the finest level
that is processed (``finestLevel``) and must be equal to it; any other value
aborts the run, because a larger value has no geometry and a smaller one would
leave the finer levels without EB data.

For developers: a tool rebuilds the geometry with
``pele_analysis::buildEBFactories(geoms, grids, dmap, verbose)`` from
``Tools/EB/PeleAnalysis_EB.H``, which reads the inputs above and returns one
``EBFArrayBoxFactory`` per level.

User-defined geometry
=====================

For ``eb2.geom_type = UserDefined`` the tools call the function
``EBUserDefined(geom, required_coarsening_level, max_coarsening_level)`` from
``Src/EBUserDefined/PeleLMeX_EBUserDefined.H``. The folder follows the PeleLMeX
convention, so a PeleLMeX case can be reused directly:

#. Copy ``PeleLMeX_EBUserDefined.H`` and ``pelelmex_prob_parm.H`` from the
   PeleLMeX case into ``Src/EBUserDefined/``, replacing the shipped files.
#. Rebuild the tool with ``USE_EB=TRUE`` (``make realclean`` first).
#. Add the geometry parameters that ``EBUserDefined`` reads (typically
   ``prob.*``) to the tool input file, with the values used in the simulation.

The shipped files are an example (a prechamber, a connecting cylinder and a
main chamber) that reads ``prob.prechamber_side_x``,
``prob.prechamber_side_y``, ``prob.cylinder_radius`` and
``prob.cylinder_length``. They only serve as a template and must be replaced
for any other geometry.

.. _eb_gpu_support:

GPU support
###########

All EB-enabled tools compile with CUDA, HIP and SYCL, with and without EB (this
is checked in CI). Compiling is not the same as running correctly: the table
below is the result of a code review of how each tool accesses its data on a
GPU build, where AMReX keeps field data in device memory. The tools have not
been run on a GPU as part of this review.

.. list-table::
   :header-rows: 1
   :widths: 22 18 60

   * - Tool
     - Status on GPU
     - Notes
   * - ``grad``
     - Runs on GPU
     - All field operations are AMReX solvers, ``MultiFab`` operations or
       ``ParallelFor`` kernels.
   * - ``curvature``
     - Runs on GPU
     - As ``grad``. The file min/max reduction synchronises once per box, which
       is correct but slow.
   * - ``integral``
     - Wrong on multi-level data
     - The step that excludes coarse cells covered by a finer level is
       compiled out on GPU builds, so covered regions are counted twice when
       ``finestLevel > 0``. Single-level integrals are correct.
   * - ``isosurface``
     - Not supported
     - The data preparation runs on the GPU, but the surface extraction, the
       distance function and the output are host loops over device memory and
       would crash or read invalid data.
