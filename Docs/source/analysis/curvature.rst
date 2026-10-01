curvature
=========

Description
-----------

``curvature`` computes curvature-related quantities from AMReX plotfiles
based on a user-defined progress variable. The tool reads an input
plotfile, evaluates geometric and kinematic quantities, and writes a
new plotfile containing the computed fields.

Sign convention
^^^^^^^^^^^^^^^

The progress variable is :math:`c = (C - C_{min})/(C_{max} - C_{min})`, so
low :math:`c` is the unburnt side. The flame normal
:math:`\mathbf{n} = -\nabla c/|\nabla c|` points towards the unburnt gas, and
the mean curvature (:math:`\kappa = \tfrac{1}{2}\nabla\cdot\mathbf{n}` in 3D,
:math:`\kappa = \nabla\cdot\mathbf{n}` in 2D) is positive where the flame front
is convex towards the unburnt gas. The normal
velocity :math:`\mathbf{u}\cdot\mathbf{n}` is positive for flow towards the
unburnt gas. The convention is the same with and without embedded boundaries.

Usage
-----

.. code-block:: bash

   curvature curvature.inp


Input File
----------

Example input file:

.. code-block:: none

   infile = plt00000
   outfile = plt00000_curvature
   finestLevel = 0
   is_per = 1 1 0
   sym_dir = 0 0 0

   Aux_Variables = 2 3 4
   progressName = Y_H2
   progMin = 0.0
   progMax = 1.0
   useFileMinMax = false
   threshold_prog = false
   threshold_value = 0.01

   do_gaussCurv = true
   do_smooth = true
   smoothing_time = 1.0e-7

   do_strain = true
   getStrainTensor = true
   do_velnormal = true


Parameters
----------

``infile``
   Input AMReX plotfile.

``outfile``
   Output plotfile containing computed quantities.

``finestLevel``
   Finest AMR level to be processed.

``is_per``
   Periodicity in each direction (``1`` periodic, ``0`` not). Required; must
   match the simulation.

``sym_dir``
   Symmetry in each non-periodic direction (``1`` applies an odd-reflection
   boundary condition to the gradient operators, ``0`` a zero-gradient
   condition). Default: ``0`` in all directions.

``Aux_Variables``
   Names of variables copied unchanged from input to output plotfile.

``progressName``
   Name of the progress variable.

``progMin``, ``progMax``
   Limits applied to the progress variable.

``useFileMinMax``
   Use minimum and maximum values from the plotfile.

``threshold_prog``
   Enable clipping outside the flame front.

``threshold_value``
   Threshold used for clipping.

``do_gaussCurv``
   Compute Gaussian curvature.

``do_smooth``
   Apply smoothing to the progress variable.

``smoothing_time``
   Smoothing time scale.

``do_strain``
   Compute strain-related quantities.

``getStrainTensor``
   Output strain tensor components.

``do_velnormal``
   Compute normal velocity.

``n_files``
   Maximum number of binary files used to write the output plotfile data
   (AMReX ``VisMF::SetNOutFiles``). Lower this to reduce the number of files
   created for large parallel post-processing runs. AMReX clamps the value to
   the number of MPI ranks, so a serial run always writes a single data file.
   Default: the AMReX default.

Embedded boundaries
-------------------

When built with ``USE_EB=TRUE``, ``curvature`` rebuilds the EB geometry, sets
the progress variable to zero in covered cells, and computes all gradients,
the smoothing, the curvatures, the strain rate and the normal velocity with
EB-aware operators, only where the volume fraction is positive. This requires
the ``eb2.*`` geometry inputs of the simulation; see :ref:`embedded_boundaries`.
``curvature`` runs on GPU builds.
