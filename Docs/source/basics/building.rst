.. highlight:: bash

Building the tools
******************

Every tool is a separate executable, built with the GNU make system of AMReX.
The executable is written to the directory where ``make`` runs and is named
``<tool><DIM>d.<compiler>[.MPI][.OMP][.CUDA].ex``, e.g. ``integral3d.gnu.MPI.ex``.

AMReX and PelePhysics
#####################

The tools need AMReX and PelePhysics. A recursive clone provides both as
submodules, and the makefiles use them by default::

   git clone --recursive --shallow-submodules --single-branch git@github.com:ITV-RWTH/PeleAnalysis.git

To build against the checkouts of an existing Pele installation instead, set
``PELE_PHYSICS_HOME`` and ``AMREX_HOME`` to absolute paths, either in the
environment or on the ``make`` command line::

   export PELE_PHYSICS_HOME=/path/to/PelePhysics
   export AMREX_HOME=$PELE_PHYSICS_HOME/Submodules/amrex

Tools in ``Src``
################

Select the tool with ``EBASE`` and set the build options on the command line::

   cd Src
   make -j8 EBASE=integral DIM=3 USE_MPI=TRUE

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Option
     - Meaning
   * - ``EBASE``
     - Tool to build, the name of its source file without ``.cpp``
   * - ``DIM``
     - Spatial dimension of the plotfiles, ``2`` or ``3`` (default ``2``)
   * - ``USE_MPI``, ``USE_OMP``, ``USE_CUDA``
     - Parallel backends, ``TRUE`` or ``FALSE`` (default ``FALSE``). Which
       backends a tool really uses is listed on the :doc:`front page </index>`.
   * - ``COMP``
     - Compiler suite (default ``gnu``)
   * - ``DEBUG``
     - ``TRUE`` for a debug build (default ``FALSE``)

To build several tools in one go, pass them as targets to ``GNUmakefile.all``;
it knows the tools in its ``TARGETS`` list. It defaults to ``DIM = 3`` and to a
PelePhysics checkout next to PeleAnalysis; override these on the command line::

   make -j8 -f GNUmakefile.all integral isosurface progVar DIM=2

Some tools need additional libraries:

- ``decimateMEF``: build qslim first with ``make qslim EBASE=decimateMEF``.
- ``convert2hdf5``: HDF5, see :doc:`convert2hdf5 </analysis/convert2hdf5>`.
- ``optimalEstimatorTraining``, ``optimalEstimatorInfer``: libtorch, see
  :doc:`optimalEstimator </analysis/optimalEstimator>`.

Model-specific tools
####################

The tools in ``Src/ModelSpecificAnalysis`` evaluate thermochemistry and are
compiled with a chemical mechanism of PelePhysics. Build SUNDIALS once per
configuration with ``make TPL``, then the tool::

   cd Src/ModelSpecificAnalysis
   make TPL
   make -j8 EBASE=computeMixtureFraction Chemistry_Model=drm19

``Chemistry_Model`` is the name of a directory in ``PelePhysics/Mechanisms`` and
must be the mechanism of the simulation. Here ``DIM`` defaults to ``3``.

Cleaning up
###########

Object files are kept in ``tmp_build_dir/``. ``make cleanconfig`` removes the
objects and the executable of the current configuration (same ``EBASE``,
``DIM`` and options); ``make clean`` removes all executables and
``tmp_build_dir/``.
