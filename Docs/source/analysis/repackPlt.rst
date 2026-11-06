.. highlight:: bash

repackPlt
*********

Rewrite an existing AMReX plot file with a user-chosen number of data files. PeleLMeX (and AMReX in general) writes plot file level data as ``Level_<l>/Cell_D_<nnnnn>`` files, by default one file per MPI rank per level (capped at 256). Large runs can therefore produce thousands of small files, which strains file-system quotas, archiving and transfers. ``repackPlt`` reads the level data of a plot file and writes it again with at most ``n_files`` data files per level.

Only the level data (``Level_*/<mf>_H`` and ``Level_*/<mf>_D_*``) is rewritten. Everything else (the top-level ``Header``, ``job_info``, particle directories, any other file or directory) is copied unchanged, so apart from the data-file layout the output is identical to the input. Cell values are bit-identical unless the precision is changed with ``float_type``. The VisMF header version and the precision of the input are kept.

Usage: ::

   mpirun -np <N> ./repackPlt3d.gnu.MPI.ex infile=<s> n_files=<n> [outfile=<s>] [options]

Example: ::

   mpirun -np 8 ./repackPlt3d.gnu.MPI.ex ./InputSamples/repackPlt.inp

.. note::

   AMReX caps the number of data files at the number of MPI ranks the tool runs on. A serial run always writes exactly one file per level. To *increase* the number of files (e.g. to split a single-file plot file into 8 files), run on at least ``n_files`` ranks. A warning is printed when ``n_files`` is larger than the rank count. Reducing the number of files works with any number of ranks.

Tool Options
#############

::

   #------------------- IO CONTROL -----------------------------------------------------------
   infile      = plt00500                     # Input AMReX plot file
   outfile     = plt00500_nf8                 # DEF: <infile>_nf<n_files>; Output plot file name

`infile` is the plot file directory to re-pack. `outfile` is the name of the new plot file; it defaults to ``<infile>_nf<n_files>`` (e.g. ``plt00500_nf8``) in the current directory. `outfile` must differ from `infile`; the input is never modified. If `outfile` already exists it is renamed to ``<outfile>.old.<id>`` (standard AMReX behaviour) before writing.

::

   #------------------- FILE LAYOUT ----------------------------------------------------------
   n_files     = 8                            # DEF: 1; Max number of data files per level (capped at MPI ranks)

`n_files` sets the maximum number of ``Cell_D`` files written per level (more precisely per MultiFab). With ``n_files`` smaller than the number of ranks, AMReX writes the data in sets of ranks that share a file. The default of 1 gives the smallest number of files.

::

   #------------------- PRECISION ------------------------------------------------------------
   float_type = 32                            # DEF: same as input; 32 or 64 bit output data

By default the output has the same precision as the input: plot files written in single precision (e.g. PeleLMeX with ``fab.format = NATIVE_32``) stay single precision. Setting `float_type` to ``32`` or ``64`` forces the precision of the written data. Converting 64-bit data to 32 bit loses precision.

::

   #------------------- Additional Flags -----------------------------------------------------
   verbose = 1                                # DEF: 0; Print progress information

With `verbose` set to 1, the tool reports each copied file and each re-packed MultiFab (number of boxes, components, precision and header version).

.. note::

   The tool works on one MultiFab at a time, so peak memory is about the size of the largest level divided over the ranks. Running with MPI parallelizes reading and writing; OpenMP and GPUs are not used.
