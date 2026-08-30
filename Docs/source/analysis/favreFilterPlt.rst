.. highlight:: bash

Apply Favre-filtering to a Plot file
************************************

The tools are executed in the following order:

1. FavrePreMultiplyPlt - multiple all the existing fields with density
2. filterPlt - apply a filter to a Plot file
3. FavrePostDividePlt - divide all the field by filtered density

The present set of tools utilizes the PelePhysics PltFileManager with the preprocessor, FavrePreMultiplyPlt
and postprocessor, FavrePostDividePlt to produce favre filtered fields in plot files
To compile, it is necessary to define the AMREX_HOME and PELE_PHYSICS_HOME
variables in the GNUmakefile. For multilevel plot files, filtering can either
be done with a constant absolute filter width, or a constant filter to grid
ratio across levels. The filter to grid ratio on the base level must always be
even. Consult `PelePhysics <https://amrex-combustion.github.io/PelePhysics/Utility.html#filter>`_
to see different filter types offered.

Usage: ::

   ./FavrePreMultiplyPlt.gnu.MPI.ex infile=plt00000 outfile=plt00000_DensityWt [options]
   ./filterPlt3d.gnu.MPI.ex infile=plt00000_DensityWt outfile=plt00000_DensityWt_Filtered [options]
   ./FavrePostDividePlt.gnu.MPI.ex infile=plt00000_Filtered outfile=plt00000_Favre_Filtered [options]

Help: ::

   ./FavrePreMultiplyPlt.gnu.MPI.ex help=true
   ./filterPlt3d.gnu.MPI.ex help=true
   ./FavrePostDividePlt.gnu.MPI.ex help=true

Example: ::
   
   ./FavrePreMultiplyPlt.gnu.MPI.ex ./InputSamples/FavrePreMultiply.inp
   ./filterPlt3d.gnu.MPI.ex ./InputSamples/filterPlt.inp
   ./FavrePostDividePlt.gnu.MPI.ex ./InputSamples/FavrePostDivide.inp

Example Input File ``FavrePreMultiply.inp``::

        #------------------- IO CONTROL -----------------------------------------------------------
        infiles = plt00000                        # Input pltfile
        outfile = plt00000_DensityWt              # Name of output file
    
        #------------------- Operation control ----------------------------------------------------
        variables = density temp HeatRelease      # DEF: all possible, list of variable names (density is required)

Example Input File ``filterPlt.inp``::

        #------------------- IO CONTROL -----------------------------------------------------------
        infiles = plt00000_DensityWt              # pltfiles to average
        outfile = plt00000_Filtered		  # DEF: plt_averaged, Name of output file
        
        #------------------- Operation control ----------------------------------------------------
        variables = temp HeatRelease              # DEF: all possible, list of variable names to average
        max_filter_level = 4			  # DEF: 1000, max level to consider for filtering
        filter_type = 1				  # DEF: 1, filter type as defined in PeleC (1->box, 2->Gaussian, etc)
        base_fgr = 2				  # DEF: 2,is the desired filter to grid ratio on the base level, must be even
        same_fgr_all_levels = false		  # DEF: false, if true the same filter to grid ratio is kept on all levels (rather than absolute filter width)

        max_grid_size = 32			  # DEF: 32, output max_grid_size. 
        interp_type = 1				  # DEF: 1, determines the type of interpolation when FillPatching: 0 -> piecewise constant, 1 -> cell cons linear

Example Input File ``FavrePostDivide.inp``::

        #------------------- IO CONTROL -----------------------------------------------------------
        infiles = plt00000_Filtered               # Input pltfile
        outfile = plt00000_Favre_Filtered         # Name of output file

        #------------------- Operation control ----------------------------------------------------
        variables = density temp HeatRelease      # DEF: all possible, list of variable names (density is required)
