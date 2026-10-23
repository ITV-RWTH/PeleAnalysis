.. highlight:: bash

isosurface
**********

Covers the :doc:`isosurface </analysis/isosurface>` tool on periodic domains:
iso lines (2D) and iso surfaces (3D) that cross periodic boundaries, on uniform
and two-level grids.

**Location:** ``Tests/isosurface/``

Expected values
###############

All cases use the iso line :math:`c = 0.5` of
:math:`c = \cos(2\pi x)\cos(2\pi y)`. It has two congruent closed lines per unit
square. The reference length :math:`L_{ref}` is measured on the non-periodic
window :math:`[0.25, 0.75]^2`, which holds one of them alone. The periodic cases
shift the field by 10 cells, so their cell centres coincide with those of the
reference window, and one line reaches less than a cell across both periodic
boundaries. The measure must therefore be :math:`2 L_{ref}` in 2D and
:math:`0.5 L_{ref}` in 3D (the 2D lines extruded over :math:`z \in [0, 0.25]`),
to round-off. ``check_mef.py`` measures the MEF in the minimum image of the
periodic directions and also checks that every node lies inside the domain and
that no element edge is longer than 1.5 finest cell diagonals.

Test matrix
###########

.. list-table::
   :header-rows: 1
   :widths: 10 35 55

   * - #
     - Feature
     - Pass criterion
   * - Ref
     - Reference line, non-periodic
     - Within 0.2% of the exact length; nodes inside; no long edge
   * - per_2d
     - Iso lines across periodic boundaries (2D)
     - Measure :math:`2 L_{ref}` (MEF and ``computeArea``); nodes inside the
       domain; no long edge
   * - per_amr_2d
     - As per_2d on two levels
     - As per_2d
   * - per_3d
     - Iso surfaces across periodic boundaries (3D)
     - Measure :math:`0.5 L_{ref}` (MEF and ``computeArea``); nodes inside the
       domain; no long edge
   * - per_amr_3d
     - As per_3d on two levels
     - As per_3d
