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
(``USE_CUDA``, ``USE_HIP``, ``USE_SYCL``).

.. note::

   AMReX places EB and non-EB objects in the same build directory
   (``tmp_build_dir/o/3d.gnu.EXE``). Run ``make realclean`` when switching
   ``USE_EB`` on or off, otherwise stale objects from the other build are
   linked in.

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
