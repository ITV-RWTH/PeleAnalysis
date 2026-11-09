
PeleAnalysis
============

This repository contains a collection of standalone routines for processing plotfiles created with the AMReX software framework for block-structured adaptive mesh refinement simulations.  Documentation is under development, but is
available at ::

   https://peleanalysis.readthedocs.io/en/latest/

AMReX is required for these tools, and is available https://github.com/AMReX-Codes/amrex

To clone this repository including the required submodules such as PelePhysics, AMReX, and qslim, use the ``--recursive`` option of the git clone command. The additional options ``--shallow-submodules`` and ``--single-branch`` help decreasing the size of your clone. ::

   git clone --recursive --shallow-submodules --single-branch git@github.com:ITV-RWTH/PeleAnalysis.git

How to build the tools, also against the PelePhysics and AMReX of an existing Pele installation, is described in the `building instructions <https://itv-rwth.github.io/PeleAnalysis/basics/building.html>`_.


Testing
^^^^^^^
Functional regression tests live in the ``Tests/`` directory.
Each subdirectory provides a self-contained ``run_tests.sh`` for one tool (not complete yet)::

   cd Tests/<toolname>
   ./run_tests.sh

See the `Testing documentation <https://itv-rwth.github.io/PeleAnalysis/testing.html>`_ for the full list of available suites and expected outputs.

Documentation
^^^^^^^^^^^^^
Documentation for the analysis routines exists in the Docs directory. To build the documentation::

    cd Docs
    make html

You can also view the compiled documentation on https://itv-rwth.github.io/PeleAnalysis/

Acknowledgment
^^^^^^^^^^^^^^
This research was supported by the Exascale Computing Project (ECP), Project
Number: 17-SC-20-SC, a collaborative effort of two DOE organizations -- the
Office of Science and the National Nuclear Security Administration --
responsible for the planning and preparation of a capable exascale ecosystem --
including software, applications, hardware, advanced system engineering, and
early testbed platforms -- to support the nation's exascale computing
imperative.
