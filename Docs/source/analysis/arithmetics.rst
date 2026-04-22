arithmetics
===========

Description
-----------
``arithmetics`` applies a simple arithmetic operation (``+``, ``-``, ``*``, or ``/``)
to two field variables from an AMReX plotfile. The operator is selected at compile
time via a preprocessor define in ``arithmetics.cpp``. The tool reads a plotfile,
computes the result of ``A OP B``, and writes a new plotfile containing all original
variables plus the derived field.

Usage
-----
.. code-block:: bash

   arithmetics3d.gnu.MPI.OMP.ex infile=FILE inVarAName=NAME inVarBName=NAME outVarName=NAME [OPTIONS]

Parameters
~~~~~~~~~

``infile``
   Path to the AMReX plotfile.

``inVarAName``
   Name of operand A (must exist in the plotfile).

``inVarBName``
   Name of operand B (must exist in the plotfile).

``outVarName``
   Name of the resulting output variable.

``outfile``
   Output plotfile name.
   Default: ``<infile>_multiply``.

``finestLevel``
   Finest AMR level to process.
   Default: plotfile finest level.

``is_per``
   Periodicity flags in each spatial direction (0: non-periodic, 1: periodic).
   Default: ``1 1 1``.

``verbose``
   Enable verbose AMReX output if present.

Output
------
A new AMReX plotfile containing:

- all original input variables (copied unchanged)
- the derived field ``outVarName = A OP B``

Compile-time Operator
---------------------
The arithmetic operator is set before compilation in ``arithmetics.cpp``:

.. code-block:: cpp

   #define ARITHMETICS_OPERATION +   // or -, *, /

Recompile after changing the operator.

Typical Applications
--------------------
- Computing dimensionsless quantities
- Weighting a quantity by another quantity

Notes
-----
Both ``inVarAName`` and ``inVarBName`` must exist in the plotfile — the tool
will abort with an error message if either is not found. The tool is compatible
with MPI and OpenMP parallelism.
