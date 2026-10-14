.. highlight:: bash

chemTableLookup
***************

Look up tabulated quantities of a chemistry table in OpenFOAM format, as written by the
``chemtable_generator`` (``tableProperties`` plus one ``<name>_prem_table`` file per quantity), at
the state given by plotfile variables. A typical use is the a priori evaluation of a tabulated LES
model: the Favre-filtered DNS state (e.g. progress variable, its sub-filter variance and the
mixture fraction) is passed in ``coordVars`` and the table quantities are returned at every cell.

* The table parameters and their axes are read from the ``params`` list of ``tableProperties``;
  ``coordVars`` names one plotfile variable per parameter, in that order. Up to four parameters
  are supported, and the axes need not be uniform.
* The values are interpolated multilinearly. States outside the table range are clamped to its
  boundary, and the fraction of cells below and above the range is reported per parameter and
  level.
* The table is read on the I/O rank and broadcast; on GPUs it is copied to the device.

The output holds the ``keepVars`` followed by the tables, named ``<outPrefix><table>``.

Usage: ::

  ./chemTableLookup3d.gnu.MPI.ex infile=<s> tableDir=<s> tables=<s> ... coordVars=<s> ... [options]

Example: ::

  ./chemTableLookup3d.gnu.MPI.ex ./InputSamples/chemTableLookup.inp

Example Input File ``chemTableLookup.inp``::

        # Example input file for chemTableLookup
        
        #------------------- IO CONTROL -----------------------------------------------------------
        infile = plt00000_favre                   # input pltfile with the lookup state
        outfile = plt00000_table                  # DEF: <infile>_table, output pltfile
        n_files = 64                              # DEF: AMReX default; max number of plotfile data
                                                  # files written (VisMF), clamped to the MPI rank count
        
        #------------------- Table ----------------------------------------------------------------
        tableDir = ./chemTables/table_H2O_filt    # directory with tableProperties and the tables
        tableSuffix = _prem_table                 # DEF: _prem_table; table file = <name><tableSuffix>
        tables = DC Dz DzSRC T                    # tables to look up
        
        #------------------- Operation control ----------------------------------------------------
        coordVars = Ct varC Zt                    # one plotfile variable per table parameter, in the
                                                  # order of 'params' in tableProperties
        outPrefix = tab_                          # DEF: none; prefix of the output variable names
        keepVars = Ct varC Zt                     # DEF: none; variables copied to the output
        finestLevel = 1000                        # DEF: finest level; max level to process
