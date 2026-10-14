.. highlight:: bash

divergence
**********

Compute the divergence of one or more vector fields of a plotfile. The derivatives are evaluated
as in :doc:`grad`: the face-centered gradient of every component is obtained from the MLMG Poisson
operator (which handles the coarse-fine interfaces of AMR data and the domain boundary conditions)
and averaged to the cell centers, i.e. second-order central differences in the interior. At
non-periodic boundaries a Neumann condition is applied, or ``reflect_odd`` where ``sym_dir`` is set.

The components of every vector are listed in ``vectorVars``, ``AMREX_SPACEDIM`` per vector, and
``outVars`` names the divergence of each vector. The output holds the ``keepVars`` followed by the
divergences.

Usage: ::

  ./divergence3d.gnu.MPI.ex infile=<s> vectorVars=<s> ... outVars=<s> ... [options]

Example: ::

  ./divergence3d.gnu.MPI.ex ./InputSamples/divergence.inp

Example Input File ``divergence.inp``::

        # Example input file for divergence
        
        #------------------- IO CONTROL -----------------------------------------------------------
        infile = plt00000                         # input pltfile
        outfile = plt00000_div                    # DEF: <infile>_div, output pltfile
        n_files = 64                              # DEF: AMReX default; max number of plotfile data
                                                  # files written (VisMF), clamped to the MPI rank count
        
        #------------------- Operation control ----------------------------------------------------
        vectorVars = x_velocity y_velocity z_velocity  # components of the vector(s), AMREX_SPACEDIM
                                                  # per vector; several vectors may be listed
        outVars = divU                            # one output name per vector
        keepVars = density                        # DEF: none; variables copied to the output
        finestLevel = 1000                        # DEF: finest level; max level to process
        sym_dir = 0 0 0                           # DEF: 0 0 0, sets case symmetry
        is_per = 1 1 0                            # DEF: 0 0 0, sets case periodicity
