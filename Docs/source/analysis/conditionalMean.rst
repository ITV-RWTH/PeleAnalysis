.. highlight:: bash


conditionalMean
******************

The tool `conditionalMean` computes conditional averages of selected variables
from AMReX plotfiles based on a specified conditioning variable.
The tool is intended for post‑processing turbulent combustion datasets
and enables analysis such as conditional statistics with respect to
mixture fraction, temperature, or progress variables.

The executable reads AMReX plotfiles and bins data according to the
conditioning variable, computing mean values within each bin.

Usage: ::

   ./conditionalMean2d.gnu.MPI.ex infile=$(ls -d plt*) binComp=0 avgComps=1 2 [options]

Example: ::

  ./conditionalMean2d.gnu.MPI.ex ./InputSamples/conditionalMean.inp

Example Input File ``conditionalMean.inp``::

        #------------------- IO CONTROL -----------------------------------------------------------
        infile = plt00000 plt00001                 # Plot file list
        finestLevel = 0                            # DEF: finest level of plot file; Sets the 
                                           # finest level to read.
        doBin = true                               # [true, false], DEF: true; infile format
                                           # (true: ".plt", false: ".dat")
        outSuffix = "_conditionalMean"             # DEF: ""; Suffix to add to the pltfile name.
        aja = false                                # DEF: false; Put the header in a separate file                                            # for gnuplot/matlab.
        writeBinMinMax = false                     # DEF: False; Write min/max values for each bin.

        #------------------- Variables ------------------------------------------------------------
        binComp = 0                                # ID of variable to condition on.
        avgComps = 1 2 3                           # Variable IDs to average.
        binMin = 0.0                               # DEF: 0.0; Min value for bins.
        binMax = 1.0                               # DEF: 1.0; Max value for bins.
        bounds = -0.1 -0.1 -0.1 0.1 0.1 0.1        # DEF: all; Domain bounds
        nBins = 64                                 # DEF: 64. Number of bins in PDF.


Output
------

The output file ``CM_<binVar>.dat`` (written to the working directory)
contains per bin the bin center, the sum, the sum of squares, the mean and
the standard deviation of each requested variable, the sample weight ``N``
and the fraction ``p`` of all samples. Every cell is weighted by its volume
in units of a finest-level cell of the first plotfile, so ``N`` is the
binned volume in those units.

All plotfiles in ``infile`` are pooled into one set of statistics (e.g. a
time average), and ``bounds`` restricts the samples to a sub-domain of each
of them. Because sums and sums of squares are written, statistics of
different files or sub-domains can be added afterwards.

Typical Applications
--------------------

- Conditional temperature statistics
- Mixture fraction conditioned quantities
- Turbulence–chemistry interaction analysis
- Flame structure diagnostics

Notes
-----

The conditioning variable must exist in the plotfile.
Large plotfiles may require significant memory depending
on bin resolution.

