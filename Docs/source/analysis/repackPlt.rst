.. highlight:: bash

repackPlt
*********

Rewrite one or more existing AMReX plot files with a user-chosen number of data files. PeleLMeX (and AMReX in general) writes plot file level data as ``Level_<l>/Cell_D_<nnnnn>`` files, by default one file per MPI rank per level (capped at 256). Large runs can therefore produce thousands of small files, which strains file-system quotas, archiving and transfers. ``repackPlt`` reads the level data of each given plot file and writes it again with at most ``n_files`` data files per level. Each plot file is processed independently and written under its own name into an output directory, e.g. ``repack_nf8/plt00500``.

Only the level data (``Level_*/<mf>_H`` and ``Level_*/<mf>_D_*``) is rewritten. Everything else (the top-level ``Header``, ``job_info``, particle directories, any other file or directory) is copied unchanged, so apart from the data-file layout the output is identical to the input. Cell values are bit-identical unless the precision is changed with ``float_type``. The VisMF header version and the precision of the input are kept.

Usage: ::

   mpirun -np <N> ./repackPlt3d.gnu.MPI.ex infiles=<s> [<s> ...] n_files=<n> [outdir=<s>] [options]

Examples: ::

   mpirun -np 8 ./repackPlt3d.gnu.MPI.ex ./InputSamples/repackPlt.inp
   mpirun -np 8 ./repackPlt3d.gnu.MPI.ex infiles= $(ls -d plt0*) n_files=8

A shell glob is not expanded inside ``infiles=plt0*`` because the word starts with ``infiles=``. Either put a space after ``=`` and let the shell list the plot files (second example), or list them in an inputs file.

.. note::

   AMReX caps the number of data files at the number of MPI ranks the tool runs on. A serial run always writes exactly one file per level. To *increase* the number of files (e.g. to split a single-file plot file into 8 files), run on at least ``n_files`` ranks. A warning is printed when ``n_files`` is larger than the rank count. Reducing the number of files works with any number of ranks.

Tool Options
#############

::

   #------------------- IO CONTROL -----------------------------------------------------------
   infiles     = plt00500 plt01000 plt01500   # Input AMReX plot files (one or more)
   outdir      = repack_nf8                   # DEF: repack_nf<n_files>; Output directory

`infiles` lists one or more plot file directories to re-pack. Each one is written to ``<outdir>/<name>``, keeping its name. `outdir` defaults to ``repack_nf<n_files>`` (e.g. ``repack_nf8``) in the current directory and is created if needed. The inputs are never modified. If ``<outdir>/<name>`` already exists, it is renamed to ``<name>.old.<id>`` (standard AMReX behaviour) before writing.

All inputs are checked before anything is written. The tool aborts without writing any output if an input is not a plot file, if two inputs have the same name (e.g. ``run1/plt00500`` and ``run2/plt00500``, which would collide in `outdir`), or if `outdir` is the folder that holds the inputs.

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

The tool always prints one progress line per plot file. With `verbose` set to 1, it also reports each copied file and each re-packed MultiFab (number of boxes, components, precision and header version).

.. note::

   The tool works on one plot file and one MultiFab at a time, so peak memory is about the size of the largest level divided over the ranks. Running with MPI parallelizes reading and writing; OpenMP and GPUs are not used.
