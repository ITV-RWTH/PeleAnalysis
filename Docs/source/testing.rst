.. highlight:: bash

Testing
*******

Functional test suites live in the ``Tests/`` directory at the repository
root. Each subdirectory covers one tool and is self-contained: it provides
synthetic input generators, tool input files, an automated runner, and a
``TESTING.md`` with per-test expected outputs.

Running a suite
###############

::

   cd Tests/jpdf
   ./run_tests.sh               # compile executables then run all tests
   ./run_tests.sh --no-compile  # skip build, use existing executables in Src/

The runner prints a colour-coded PASS/FAIL line for each assertion and exits
with a non-zero status if any test fails. All artefacts land in
``Tests/<toolname>/testrun/`` and are excluded from version control.

Available suites
################

jpdf
~~~~

**Location:** ``Tests/jpdf/``

Covers the :doc:`analysis/jpdf` tool (Joint PDFs and 2D conditional means).

.. list-table::
   :header-rows: 1
   :widths: 5 30 65

   * - #
     - Feature
     - Pass criterion
   * - 1
     - All six output formats (plotfile, gnuplot, MATLAB, Tecplot, FAB, scatter)
     - All file types present; PDF sum = 1.0
   * - 2
     - 2D conditional mean (independent variable)
     - ``condMean_var_cond_on_*`` file created; all non-zero entries ≈ 0.5
   * - 3
     - Duplicate variable in ``condMean_vars``
     - Tool runs without error; deduplication confirmed in verbose output
   * - 4
     - ``useminmax`` range override + bin clamping
     - Header encodes overridden axis range; verbose reports ``v1g > 0``
   * - 5
     - ``do_conditioning=1`` (range filter)
     - PDF restricted to conditioned cells; sum = 1.0
   * - 6
     - ``do_conditioning=2`` (c(1−c) filter)
     - Tail bins of conditioning variable excluded; sum = 1.0
   * - 7
     - ``norm_cVal=1`` (normalised conditioning)
     - Only normalised-range cells contribute; sum = 1.0
   * - 8
     - Temporal averaging (``do_average=1``)
     - ``JPDFAverage*/`` created; averaged PDF matches per-file PDF exactly

Each test runs in both **2D** (tests 1, 2, 5, 8) and **3D** serial, and
tests 1, 5, 8 are additionally validated with **MPI** at 2 and 4 ranks.
MPI results are compared element-wise to serial (max diff < 10⁻¹⁰).

Known gaps
~~~~~~~~~~

- **AMR multi-level:** ``generateTestPlt`` produces single-level plotfiles;
  the ``finestLevel`` parameter is not exercised.
- **Slash-in-variable-name:** ``ProtectSlashes`` (``Y(OH)``-style names)
  cannot be tested with ``generateTestPlt`` due to ParmParse key
  constraints; requires a hand-crafted plotfile.

Adding a new suite
##################

1. Create ``Tests/<toolname>/`` with input files and a ``run_tests.sh``
   modelled on ``Tests/jpdf/run_tests.sh``.
2. Set ``SRC_DIR="$SCRIPT_DIR/../../Src"`` so the runner finds the built
   executables.
3. Write a ``TESTING.md`` documenting expected outputs for each test.
4. Add the tool to the table in ``Tests/README.rst`` and to this page.
