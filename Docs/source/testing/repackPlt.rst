.. highlight:: bash

repackPlt
*********

Covers the :doc:`repackPlt </analysis/repackPlt>` tool (rewriting a plotfile
with a chosen number of data files).

**Location:** ``Tests/repackPlt/``

Data equality is checked box by box with the included ``compare_plt.py``,
which reads the VisMF headers and data files directly and does not depend on
any AMReX tool.

.. list-table::
   :header-rows: 1
   :widths: 5 35 60

   * - #
     - Feature
     - Pass criterion
   * - 1
     - Serial, default ``n_files=1`` and ``outdir``
     - ``repack_nf1/plt_repack`` has 1 ``Cell_D`` file per level; data
       identical; ``Header`` and extra files/directories copied byte for byte
   * - 2
     - VisMF header version 4 input
     - Output header version 4; data identical
   * - 3
     - 32-bit input
     - Output 32-bit; data identical
   * - 4
     - ``float_type=32`` on 64-bit input
     - Output 32-bit; values within 1e-5
   * - 5
     - Batch of 3 plotfiles into one ``outdir``
     - 3 plotfiles with unchanged names; 1 file per level; data identical
   * - 6--10
     - ``outdir`` = input folder, duplicate names, missing input in a batch,
       no ``infiles``, ``n_files=0``
     - Tool aborts before writing; inputs untouched; no output
   * - MPI
     - 4 ranks with ``n_files`` = 1, 2, 4; 2 ranks with ``n_files=16``;
       single-file input split into 4 files; batch of 3 with ``n_files=2``
     - Exactly the expected number of files per level (capped at the rank
       count, with a warning); data identical
