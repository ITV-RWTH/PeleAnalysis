.. highlight:: bash

stream
******

Covers the :doc:`partStream </analysis/partStream>` and
:doc:`streamBinTubeStats </analysis/streamBinTubeStats>` tools. Every case runs
the whole chain ``generateTestPlt`` → ``isosurface`` → ``partStream`` →
``streamBinTubeStats``.

**Location:** ``Tests/stream/``

Expected values
###############

The cases are built so that the answers are known exactly. With 32 cells on a
unit domain, ``hRK = 0.5`` and ``cSpace = 0``, each step advances by
``dx/2 = 1/64``, and ``Nsteps = 20`` gives 39 points, i.e. 38 intervals:

.. math::

   \ell = 38/64 = 0.59375

For the constant field ``q = 1`` the tube integral reported as ``q_volInt``
equals that stream length. The ``c = 0.5`` isosurface of a plane that also
wraps across the periodic boundary consists of two sheets of unit extent, so
the total area is 2 and the total stream-tube volume is
:math:`2\,\ell = 1.1875`. (The second sheet is placed at ``x = 0.5`` rather
than at the boundary by ``isosurface``; the expected values depend only on the
area.) Before streams
were allowed to cross periodic boundaries the same case gave 0.497 to 0.594 per
tube and 1.090 in total.

Test matrix
###########

.. list-table::
   :header-rows: 1
   :widths: 10 35 55

   * - #
     - Feature
     - Pass criterion
   * - T1
     - Streams crossing a periodic boundary (3D)
     - ``q_volInt = 0.59375`` and ``q_avg = 1`` in every element; total area 2;
       total volume 1.1875; no NaN or Inf; mapped ``q = 1`` along every stream;
       no particles lost
   * - T2
     - Vector field vanishing for ``x > 0.9``
     - Positions and statistics stay finite; no stream runs past ``x = 0.95``;
       the tubes that reach the vanishing field are shorter, none longer
   * - T4
     - Seeds between the last cell centre and a non-periodic wall
     - Mapped ``temp = 300`` and ``q = 1`` exactly, in serial and on 4 ranks
   * - T5
     - The same chain in 2D
     - ``q_volInt = 0.59375``; total length 2; total volume 1.1875; no NaN
   * - T3
     - Seeds outside the domain in a non-periodic direction
     - Exits non-zero, the message names the direction, no output written
   * - T6
     - ``writeStreamBin`` without a seed surface
     - Exits non-zero
   * - T7
     - Stream directory without a readable ``Header``
     - Exits non-zero, the message names the missing file
   * - T8
     - A surface node whose stream is missing
     - Exits non-zero, naming how many nodes and elements are affected
       (``break_stream_bin.py`` removes one record from a finished directory)
   * - T-MPI2, T-MPI4
     - ``partStream`` on 2 and 4 ranks
     - Reproduces the T1 numbers; the surface file matches the serial one to
       1e-12
   * - T-MPI4 abort
     - Seeds outside the domain, 4 ranks
     - All ranks exit non-zero within the timeout, i.e. no deadlock

Running
#######

::

   cd Tests/stream
   ./run_tests.sh               # compile + run
   ./run_tests.sh --no-compile  # reuse the executables in Src/

MPI tests are skipped when neither ``mpirun`` nor ``mpiexec`` is available, and
T8 is skipped without ``python3``. See ``Tests/stream/TESTING.md`` for the full
derivation of the expected values.
