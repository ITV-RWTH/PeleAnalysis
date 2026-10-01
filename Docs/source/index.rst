.. PeleAnalysis documentation master file

..
   NOTE FOR AI AGENTS / AUTOMATED TOOLS: The parallel-support tables
   on this page must be kept up to date. Whenever you add, remove, rename or
   modify a tool in Src/ or Src/ModelSpecificAnalysis/, or change how it uses
   OpenMP, MPI or GPU (AMReX ParallelFor, MFIter tiling, omp pragmas,
   DistributionMapping, reductions, AmrData::FillVar, host access to FAB
   data, ...), re-evaluate that tool against the definitions in the legend and
   update its row (and its Basis entry) in Docs/source/index.rst.

.. raw:: html

   <!--
   NOTE FOR AI AGENTS / AUTOMATED TOOLS: The parallel-support tables
   on this page must be kept up to date. Whenever you add, remove, rename or
   modify a tool in Src/ or Src/ModelSpecificAnalysis/, or change how it uses
   OpenMP, MPI or GPU (AMReX ParallelFor, MFIter tiling, omp pragmas,
   DistributionMapping, reductions, AmrData::FillVar, host access to FAB
   data, ...), re-evaluate that tool against the definitions in the legend and
   update its row (and its Basis entry) in Docs/source/index.rst.
   -->

PeleAnalysis documentation
==========================

PeleAnalysis is a collection of AMReX-based post-processing tools for Pele
plotfiles. The tables below summarize, for every tool in ``Src`` and
``Src/ModelSpecificAnalysis``, whether it can be run with OpenMP, MPI and
on GPUs. Click a tool name to open its documentation page.

Legend
------

A backend counts as supported only if the tool compiles with it, runs, gives
results identical to the serial build (up to floating-point round-off from
reordered reductions), **and** actually uses the backend for its main work.

* ✅ **Complete**: compiles, runs, gives identical results and the backend
  really parallelizes the main work.
* ➖ **Unused**: compiles, runs and gives correct results, but the backend does
  not parallelize the main work (e.g. serial loops, all work done on one MPI
  rank or repeated on every rank, host-only code in a GPU build).
* ❌ **No**: does not compile, crashes/aborts, or gives wrong or
  non-reproducible results with this backend.

The *Basis* column states how the entry was determined:

* **R** -- runtime-verified for OpenMP and MPI: serial, OpenMP (4 threads)
  and MPI (4 ranks) builds were run on synthetic 3D AMR plotfiles and the
  outputs compared against the serial result. The GPU column is always based
  on source-code review (no GPU hardware was available).
* **S** -- source-code review only (the tool needs external libraries such
  as HDF5 or libtorch, or input not producible in the test setup).

Tools without a link have no documentation page yet.

General tools (``Src``)
-----------------------

.. list-table:: Parallel support of the general tools
   :header-rows: 1
   :widths: 40 15 15 15 15

   * - Tool
     - OpenMP
     - MPI
     - GPU
     - Basis
   * - :doc:`amrToFE <analysis/amrToFE>`
     - ➖
     - ❌
     - ❌
     - R
   * - :doc:`arithmetics <analysis/arithmetics>`
     - ✅
     - ✅
     - ✅
     - R
   * - :doc:`avgPlotfiles <analysis/avgPlotfiles>`
     - ✅
     - ✅
     - ✅
     - R
   * - :doc:`avgToPlane <analysis/avgToPlane>`
     - ✅
     - ✅
     - ❌
     - R
   * - binMEF
     - ➖
     - ✅
     - ❌
     - R
   * - buildDistance
     - ❌
     - ❌
     - ❌
     - R
   * - :doc:`checkIso <analysis/checkIso>`
     - ➖
     - ➖
     - ❌
     - R
   * - :doc:`ciao2plt <analysis/ciao2plt>`
     - ➖
     - ✅
     - ❌
     - S
   * - combineMEF
     - ❌
     - ❌
     - ❌
     - R
   * - :doc:`combinePlts <analysis/combinePlts>`
     - ➖
     - ✅
     - ✅
     - R
   * - :doc:`conditionalMean <analysis/conditionalMean>`
     - ➖
     - ✅
     - ❌
     - R
   * - :doc:`convert2hdf5 <analysis/convert2hdf5>`
     - ➖
     - ✅
     - ✅
     - S
   * - :doc:`curvature <analysis/curvature>`
     - ➖
     - ✅
     - ❌
     - R
   * - decimateMEF
     - ➖
     - ➖
     - ➖
     - R
   * - :doc:`diffPlts <analysis/diffPlts>`
     - ✅
     - ✅
     - ✅
     - R
   * - :doc:`dumpFABslice <analysis/dumpFABslice>`
     - ➖
     - ➖
     - ❌
     - R
   * - :doc:`favreAvgPlotfiles <analysis/favreAvgPlotfiles>`
     - ✅
     - ✅
     - ✅
     - R
   * - :doc:`filterPlt <analysis/filterPlt>`
     - ✅
     - ✅
     - ✅
     - R
   * - :doc:`flattenAMRFile <analysis/flattenAMRFile>`
     - ✅
     - ✅
     - ✅
     - R
   * - :doc:`generateTestPlt <analysis/generateTestPlt>`
     - ➖
     - ✅
     - ❌
     - R
   * - :doc:`grad <analysis/grad>`
     - ✅
     - ✅
     - ❌
     - R
   * - :doc:`integral <analysis/integral>`
     - ➖
     - ✅
     - ❌
     - R
   * - isoMEF
     - ➖
     - ➖
     - ❌
     - R
   * - :doc:`isosurface <analysis/isosurface>`
     - ➖
     - ❌
     - ❌
     - R
   * - :doc:`jpdf <analysis/jpdf>`
     - ➖
     - ❌
     - ❌
     - R
   * - :doc:`makeImage <analysis/makeImage>`
     - ➖
     - ➖
     - ❌
     - R
   * - :doc:`makePlotfile <analysis/makePlotfile>`
     - ➖
     - ✅
     - ❌
     - S
   * - mergeMEF
     - ❌
     - ❌
     - ❌
     - R
   * - multMEF
     - ❌
     - ❌
     - ❌
     - R
   * - :doc:`optimalEstimatorInfer <analysis/optimalEstimator>`
     - ➖
     - ✅
     - ❌
     - S
   * - :doc:`optimalEstimatorTraining <analysis/optimalEstimator>`
     - ➖
     - ❌
     - ❌
     - S
   * - :doc:`partStream <analysis/partStream>`
     - ➖
     - ✅
     - ❌
     - R
   * - :doc:`progVar <analysis/progVar>`
     - ✅
     - ✅
     - ❌
     - R
   * - :doc:`qCriterion <analysis/qCriterion>`
     - ❌
     - ✅
     - ❌
     - R
   * - :doc:`regridPlt <analysis/regridPlt>`
     - ➖
     - ✅
     - ✅
     - R
   * - :doc:`rmsVel <analysis/rmsVel>`
     - ➖
     - ✅
     - ❌
     - R
   * - :doc:`scaleMEF <analysis/scaleMEF>`
     - ➖
     - ➖
     - ❌
     - R
   * - :doc:`sequentialFilterPlt <analysis/sequentialFilterPlt>`
     - ✅
     - ✅
     - ✅
     - R
   * - sliceMEF
     - ❌
     - ❌
     - ❌
     - R
   * - slicePlot
     - ➖
     - ❌
     - ❌
     - R
   * - smoothMEF
     - ❌
     - ❌
     - ❌
     - R
   * - :doc:`stream2plt <analysis/stream2plt>`
     - ➖
     - ❌
     - ❌
     - S
   * - :doc:`streamBinTubeStats <analysis/streamBinTubeStats>`
     - ❌
     - ❌
     - ➖
     - R
   * - :doc:`subPlt <analysis/subPlt>`
     - ➖
     - ✅
     - ✅
     - R
   * - :doc:`surfDATtoMEF <analysis/surfMEFtools>`
     - ➖
     - ➖
     - ❌
     - R
   * - :doc:`surfMEFtoDAT <analysis/surfMEFtools>`
     - ➖
     - ➖
     - ❌
     - R
   * - :doc:`template <basics/template>`
     - ✅
     - ✅
     - ❌
     - R
   * - :doc:`trimMEFgen <analysis/trimMEFgen>`
     - ➖
     - ➖
     - ❌
     - R
   * - :doc:`turbfile <analysis/turbfile>`
     - ➖
     - ❌
     - ❌
     - R

Model-specific tools (``Src/ModelSpecificAnalysis``)
----------------------------------------------------

.. list-table:: Parallel support of the model-specific tools
   :header-rows: 1
   :widths: 40 15 15 15 15

   * - Tool
     - OpenMP
     - MPI
     - GPU
     - Basis
   * - :doc:`computeMixtureFraction <modelSpecific/computeMixtureFraction>`
     - ✅
     - ✅
     - ❌
     - R
   * - :doc:`plotQPD <modelSpecific/plotQPD>`
     - ❌
     - ✅
     - ❌
     - R
   * - :doc:`plotTransportCoeff <modelSpecific/plotTransportCoeff>`
     - ✅
     - ✅
     - ❌
     - R
   * - :doc:`plotTYtoLe <modelSpecific/plotTYtoLe>`
     - ✅
     - ✅
     - ❌
     - R
   * - :doc:`plotXtoY <modelSpecific/plotXtoY>`
     - ➖
     - ✅
     - ❌
     - R
   * - :doc:`plotYtoX <modelSpecific/plotYtoX>`
     - ➖
     - ✅
     - ❌
     - R
   * - sCO2
     - ❌
     - ❌
     - ❌
     - R
   * - :doc:`testQPDtools <modelSpecific/testQPDtools>`
     - ➖
     - ➖
     - ➖
     - R
   * - :doc:`testTsolve <modelSpecific/testTsolve>`
     - ➖
     - ✅
     - ❌
     - R

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Basics:

   basics/data
   basics/template
   basics/analysis_util
   basics/embeddedBoundaries

.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Analysis:

   analysis/amrToFE
   analysis/arithmetics
   analysis/avgPlotfiles
   analysis/avgToPlane
   analysis/checkIso
   analysis/ciao2plt
   analysis/combinePlts
   analysis/conditionalMean
   analysis/convert2hdf5
   analysis/curvature
   analysis/diffPlts
   analysis/dumpFABslice
   analysis/favreAvgPlotfiles
   analysis/filterPlt
   analysis/flattenAMRFile
   analysis/generateTestPlt
   analysis/grad
   analysis/integral
   analysis/isosurface
   analysis/jpdf
   analysis/makeImage
   analysis/makePlotfile
   analysis/optimalEstimator
   analysis/partStream
   analysis/progVar
   analysis/qCriterion
   analysis/regridPlt
   analysis/rmsVel
   analysis/scaleMEF
   analysis/sequentialFilterPlt
   analysis/stream2plt
   analysis/streamBinTubeStats
   analysis/subPlt
   analysis/surfMEFtools
   analysis/trimMEFgen
   analysis/turbfile


.. toctree::
   :maxdepth: 2
   :hidden:
   :caption: Model Specific:

   modelSpecific/computeMixtureFraction
   modelSpecific/plotDisplacementSpeed
   modelSpecific/plotQPD
   modelSpecific/plotTransportCoeff
   modelSpecific/plotTYtoLe
   modelSpecific/plotXtoY
   modelSpecific/plotYtoX
   modelSpecific/testQPDtools
   modelSpecific/testTsolve


.. toctree::
   :maxdepth: 1
   :hidden:
   :caption: Testing:

   testing

Indices and tables
==================
* :ref:`genindex`
* :ref:`search`
