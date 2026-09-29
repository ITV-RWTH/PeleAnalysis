.. highlight:: bash

plotDisplacementSpeed
*********************
Compute the flame displacement speed of a fuel species from an AMReX plot file
produced by a reacting flow solver such as PeleLMeX, split into its diffusive
and reactive contributions. The fuel mass fraction :math:`Y_F` is used as the
progress variable. Transport coefficients and the mean molecular weight are
evaluated with PelePhysics, the gradients and flux divergences are computed
with AMReX linear operators, and the result is written as a new AMReX plot
file.

The displacement speed is evaluated as

.. math::

   S_d = S_{d,Y} + S_{d,T} + S_{d,W} + S_{d,C},

with

.. math::

   S_{d,Y} = -\frac{\nabla\cdot\left(\mathcal{D}_Y \nabla Y_F\right)}{\rho\,|\nabla Y_F|}, \quad
   S_{d,T} = -\frac{\nabla\cdot\left(\mathcal{D}_T \nabla T\right)}{\rho\,|\nabla Y_F|}, \quad
   S_{d,W} = -\frac{\nabla\cdot\left(\mathcal{D}_W \nabla \overline{W}\right)}{\rho\,|\nabla Y_F|}, \quad
   S_{d,C} = -\frac{\dot{\omega}_F}{\rho\,|\nabla Y_F|},

where :math:`\rho` is the density, :math:`T` the temperature,
:math:`\overline{W}` the mean molecular weight and :math:`\dot{\omega}_F` the
fuel reaction rate taken from the plot file (``I_R(<fuel>)``). The coefficients
are built from the mixture-averaged diffusion coefficient :math:`D_F` of the
fuel and its molecular weight :math:`W_F`:

- Fickian term: :math:`\mathcal{D}_Y = D_F\,\overline{W}/W_F`
- Soret term: :math:`\mathcal{D}_T = 0.664\,\theta_F/T`, with the thermal
  diffusion coefficient :math:`\theta_F` (only with ``transport.use_soret = 1``,
  otherwise zero)
- Mean molecular weight term: :math:`\mathcal{D}_W = D_F\,Y_F/W_F` (only with
  ``use_wbar = 1``, otherwise zero)

Quantities are converted from the CGS units of PelePhysics to SI, so the plot
file is expected in SI units (as written by PeleLMeX). The displacement speed
is only evaluated where :math:`|\nabla Y_F| > 10^{-2}` (in 1/m); all other
cells are set to zero.

Usage: ::

   ./plotDisplacementSpeed3d.gnu.MPI.ex infile=<s> fuelName=<s> [options]

Example: ::

   ./plotDisplacementSpeed3d.gnu.MPI.ex ./InputsSamples/plotDisplacementSpeed.inp

.. note::

   The input plot file must contain ``Y(<species>)`` for all species of the
   compiled mechanism, ``temp``, ``density``, ``I_R(<fuel>)`` and ``rhoh``.
   The tool aborts if one of them is missing.

Tool Options
############
::

   #------------------- IO CONTROL -----------------------------------------------------------
   infile = plt00500                          # Input AMReX plot file
   outfile = plt00500_sd                      # DEF: <infile>_sd; Name of the output plot file
   n_files = 64                               # DEF: AMReX default; cap on the number of plotfile data files

`infile` is the AMReX plot file to read. `outfile` sets the name of the output
plot file and defaults to ``<infile>_sd``, where ``<infile>`` is the last path
component of the input. `n_files` caps the number of binary files used to write
the output plot file data (AMReX ``VisMF::SetNOutFiles``); AMReX clamps it to
the number of MPI ranks.
::

   #------------------- Species --------------------------------------------------------------
   fuelName = H2                              # Fuel species; must exist in the compiled mechanism

`fuelName` selects the species whose mass fraction is used as the progress
variable and whose displacement speed is computed. It must be given and must be
part of the compiled mechanism.
::

   #------------------- Transport contributions ----------------------------------------------
   transport.use_soret = 0                    # DEF: 0; include the Soret (thermal diffusion) term
   use_wbar = 1                               # DEF: 1; include the mean molecular weight term

`transport.use_soret` switches the Soret contribution :math:`S_{d,T}` on (``1``)
or off (``0``). `use_wbar` switches the mean molecular weight contribution
:math:`S_{d,W}` on or off. Disabled contributions are written as zero.
::

   #------------------- AMR Control ----------------------------------------------------------
   finestLevel = 2                            # DEF: finest level in file; Finest AMR level to process

`finestLevel` sets the finest AMR level to include. All levels from 0 up to and
including `finestLevel` are processed and written to the output.
::

   #------------------- Boundary conditions --------------------------------------------------
   is_per = 1 1 1                             # DEF: 1 1 1; periodicity in each direction
   sym_dir = 0 0 0                            # DEF: 0 0 0; symmetry in each non-periodic direction

`is_per` sets the periodicity per direction and should match the simulation.
In non-periodic directions, `sym_dir = 1` applies an odd-reflection (symmetry)
boundary condition to the gradient operators; otherwise a zero-gradient
(Neumann) condition is used.
::

   #------------------- Additional Flags -----------------------------------------------------
   verbose = 1                                # DEF: 0; print progress information

Embedded boundaries
###################
When built with ``USE_EB=TRUE``, the tool rebuilds the EB geometry and uses
EB-aware operators for all gradients and divergences; the displacement speed
is only written where the volume fraction is positive. This requires the
``eb2.*`` geometry inputs of the simulation; see :ref:`embedded_boundaries`.
The tool runs on GPU builds.

Output
######
The output plot file contains six variables:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Variable
     - Content
   * - ``Sd(<fuel>)_dY``
     - Fickian diffusion contribution :math:`S_{d,Y}`
   * - ``Sd(<fuel>)_dT``
     - Soret contribution :math:`S_{d,T}`
   * - ``Sd(<fuel>)_dW``
     - Mean molecular weight contribution :math:`S_{d,W}`
   * - ``Sd(<fuel>)_C``
     - Reaction contribution :math:`S_{d,C}`
   * - ``Sd(<fuel>)``
     - Total displacement speed :math:`S_d`
   * - ``Sd(<fuel>)rhoh``
     - :math:`S_d` multiplied by the plot file field ``rhoh``

Dependencies
############
This tool requires a PelePhysics-enabled build. The chemical mechanism, the
equation of state and the transport model are compiled in
(``Chemistry_Model``, ``Eos_Model``, ``Transport_Model`` in
``Src/ModelSpecificAnalysis/GNUmakefile``; the defaults are ``drm19``,
``Fuego`` and ``Simple``) and must match the simulation. SUNDIALS has to be
built once with ``make TPL`` in ``Src/ModelSpecificAnalysis``.
