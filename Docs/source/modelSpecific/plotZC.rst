plotZC
======

Description
-----------

``plotZC`` computes the Bilger mixture fraction ``Z`` together with two
progress-variable definitions ``CF`` and ``CP`` from an AMReX plotfile. The
mixture fraction is evaluated from the elemental (C, H, O) composition of the
mixture using the Bilger formulation and the PelePhysics equation of state,
so it is consistent with the compiled chemical mechanism.

For each cell the tool computes

.. math::

   Z = \frac{\beta - \beta_\mathrm{ox}}{\beta_\mathrm{fu} - \beta_\mathrm{ox}},
   \qquad
   \beta = \sum_n \beta_n Y_n

where :math:`\beta_n` are the species Bilger coupling factors built from the
elemental composition and atomic weights, and :math:`\beta_\mathrm{fu}`,
:math:`\beta_\mathrm{ox}` are the values in the pure fuel and oxidizer
streams. The fuel stream is taken as pure fuel, and the oxidizer stream as
air (:math:`Y_{O_2}=0.233`, :math:`Y_{N_2}=0.767`).

Two progress variables based on normalized mass fractions are also written:

.. math::

   C_F = 1 - \frac{Y_\mathrm{fuel}}{Y_{\mathrm{fuel},\max}},
   \qquad
   C_P = \frac{Y_\mathrm{prod}}{Y_{\mathrm{prod},\max}}

The normalization limits are taken by default from the min/max of the fuel and
product fields over all AMR levels, but can be overridden by the user.


Usage
-----

.. code-block:: bash

   plotZC infile=<s> [options]


Input File
----------

The program is controlled via command-line arguments or an input file with
the following structure:

.. code-block:: none

   infile      = plt000000
   fuelName    = H2
   productName = H2O
   finestLevel = 2


Parameters
~~~~~~~~~~

``infile``
   Input AMReX plotfile. It must contain the species mass fractions
   ``Y(<species>)`` for all species in the compiled mechanism.

``fuelName``
   Name of the fuel species used for the fuel-based progress variable and to
   define the pure-fuel stream. Default: ``H2``.

``productName``
   Name of the product species used for the product-based progress variable.
   Default: ``H2O``.

``finestLevel``
   Finest AMR level up to which the fields are computed.
   Default: plotfile finest level.

``YFminmax``
   Optional two-value override for the min/max of the fuel mass fraction used
   to normalize ``CF``. Default: computed from the data.

``YPminmax``
   Optional two-value override for the min/max of the product mass fraction
   used to normalize ``CP``. Default: computed from the data.

``is_per``
   Periodicity flags in each spatial direction (0: non-periodic, 1: periodic).
   Default: ``1 1 1``.


Output
------

A new AMReX plotfile named ``<infile>_ZC`` containing:

- ``Z``  — Bilger mixture fraction
- ``CF`` — fuel-based progress variable
- ``CP`` — product-based progress variable

The output inherits the domain geometry, coordinate system, and box
structure from the input plotfile.


Typical Applications
--------------------

- Mixture-fraction / progress-variable conditioning of reacting-flow data
- Flamelet and manifold post-processing
- Diagnostics of partially premixed and non-premixed combustion


Notes
-----

This tool requires a PelePhysics-enabled build. The chemical mechanism is
compiled in at build time and determines the species list and the elemental
composition used to build the Bilger factors; the tool cannot be used with a
generic AMReX build. The oxidizer stream is currently hard-coded as air
(:math:`Y_{O_2}=0.233`, :math:`Y_{N_2}=0.767`) and the output filename
``<infile>_ZC`` cannot be overridden.
