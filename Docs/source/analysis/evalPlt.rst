.. highlight:: bash

evalPlt
*******

Evaluate user-defined expressions of plotfile variables and write the results to a new plotfile.
The expressions are parsed at run time by the AMReX parser (``amrex::Parser``), so arbitrary
combinations of variables, named constants and the usual functions (``sqrt``, ``exp``, ``log``,
``pow``, ``min``, ``max``, ``abs``, ``sin``, ``tanh``, ...) can be formed in one call, replacing chains
of :doc:`arithmetics` calls.

* Several input plotfiles on the same domain may be given; a variable is read from the first file
  that contains it. Files with different grids at a level are interpolated onto the grids of the
  first file.
* The expressions in ``outVars`` are evaluated in the given order, and each may use the results of
  the earlier ones.
* Plotfile variables whose names are identifiers (``density``, ``temp``, ``Ct_gx``) are used by
  name. Names such as ``Y(H2)`` are given a symbol with ``alias = "YH2:Y(H2)"``.
* ``constants`` defines named constants; ``pi`` is predefined. The cell-center coordinates are
  available as ``x``, ``y``, ``z``, unless a variable of that name exists.
* An expression may reference at most 16 distinct variables.

The output holds the ``keepVars`` followed by the computed variables (or the subset given in
``writeVars``). Unquoted expressions with blanks on the command line are joined back together.

Usage: ::

  ./evalPlt3d.gnu.MPI.ex infiles=<s> [<s> ...] outVars=<s> [<s> ...] expr.<name>=<s> [options]

Example: ::

  ./evalPlt3d.gnu.MPI.ex ./InputSamples/evalPlt.inp
  ./evalPlt3d.gnu.MPI.ex infiles=plt00000 outVars=rhoT expr.rhoT="density*temp"

Example Input File ``evalPlt.inp``::

        # Example input file for evalPlt
        
        #------------------- IO CONTROL -----------------------------------------------------------
        infiles = plt00000 plt00000_D             # plotfile(s) on the same domain; a variable is read
                                                  # from the first file that contains it
        outfile = plt00000_eval                   # DEF: <infiles[0]>_eval, output pltfile
        n_files = 64                              # DEF: AMReX default; max number of plotfile data
                                                  # files written (VisMF), clamped to the MPI rank count
        
        #------------------- Operation control ----------------------------------------------------
        alias = "YH2:Y(H2)" "YO2:Y(O2)" "rhoDH2:rhoD(H2)"  # DEF: none; symbol:variable, for variable
                                                  # names that are not identifiers. Names that are
                                                  # identifiers (density, temp, ...) need no alias
        constants = "nu:7.936" "YO2air:0.2329"    # DEF: none (pi is predefined); named constants
        keepVars = density                        # DEF: none; input variables copied to the output
        outVars = Z rhoZ DZ                       # variables to compute, evaluated in this order;
                                                  # each may use the earlier ones
        expr.Z = "(nu*YH2 + YO2air - YO2)/(nu + YO2air)"
        expr.rhoZ = "density*Z"
        expr.DZ = "max(rhoDH2, 0.0)"
        #writeVars = Z rhoZ                       # DEF: all outVars; subset of outVars to write
        finestLevel = 1000                        # DEF: finest level; max level to process
        # The cell-center coordinates are available as x, y, z (unless a variable has that name)
