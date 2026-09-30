#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil_C.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>
#include <AMReX_VisMF.H>

#include <AMReX_MLMG.H>

// For min/max
#include <AMReX_Reduce.H>

#ifdef AMREX_USE_EB
#include <AMReX_MLEBABecLap.H>
#include <AMReX_EBMultiFabUtil.H>
#include <PeleAnalysis_EB.H>
#else
#include <AMReX_MLPoisson.H>
#include <AMReX_MLABecLaplacian.H>
#endif

using namespace amrex;

static void
print_usage(int, char* argv[])
{
  std::cerr << "Usage:\n"
            << "  " << argv[0]
            << " infile=FILE progressName=NAME is_per=\"X Y [Z]\" [OPTIONS]\n\n"

            << "Required arguments:\n"
            << "  infile=FILE        AMReX plotfile\n"
            << "  progressName=NAME  Variable used to build the progress "
               "variable\n"
            << "  is_per=X Y [Z]     Periodicity in each direction (0 or 1)\n\n"

            << "Options:\n"
            << "  -h, --help         Show this help message\n\n"
            << "Visit PeleAnalysis/Src/InputSamples for examples or refer to "
            << "the documentation.\n";

  std::exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile, std::string("/"));
  return tokens[tokens.size() - 1];
}

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  {
    if (argc < 2) {
      print_usage(argc, argv);
    } else if (
      (std::strcmp(argv[1], "-h") == 0) ||
      (std::strcmp(argv[1], "--help") == 0)) {
      print_usage(argc, argv);
    }

    // ---------------------------------------------------------------------
    // Set defaults input values
    // ---------------------------------------------------------------------
    std::string progressName = "temp";
    Real progMin = 1.0e20;
    Real progMax = -1.0e20;
    int finestLevel = 1000;
    int verbose = 0;
    int floorIt = 0;
    int useFileMinMax = 1;
    bool do_strain = false;
    bool do_gaussCurv = false;
    bool getStrainTensor = false;
    bool do_velnormal = false;
    bool do_threshold = false;
    Real threshold = 0.0001;
    bool do_smooth = false;
    Real smooth_time = 1.0e-7;
    int nAuxVar = 0;
    int n_files = amrex::VisMF::GetNOutFiles();

    // ---------------------------------------------------------------------
    // ParmParse
    // ---------------------------------------------------------------------
    ParmParse pp;

    // IO
    pp.query("verbose", verbose);
    pp.query(
      "n_files", n_files); // Cap on the number of plotfile data files (VisMF)
    std::string plotFileName;
    pp.get("infile", plotFileName);
    std::string outfile(getFileRoot(plotFileName) + "_K");
    pp.query("outfile", outfile);
    pp.query("finestLevel", finestLevel);
    pp.query("do_gaussCurv", do_gaussCurv);

    // Progress variable
    pp.query("progressName", progressName);
    pp.query("progMin", progMin);
    pp.query("progMax", progMax);
    pp.query("floorIt", floorIt);
    pp.query("useFileMinMax", useFileMinMax);

    // Clip results outside the flame front ? ( C ~ 0 or C ~ 1)
    pp.query("threshold_prog", do_threshold);
    pp.query("threshold_value", threshold);

    // Progress variable smoothing
    pp.query("do_smooth", do_smooth);
    pp.query("smoothing_time", smooth_time);

    // Pertaining to the velocity computation
    pp.query("do_strain", do_strain);
    if (do_strain) {
      pp.query("getStrainTensor", getStrainTensor);
    }
    pp.query("do_velnormal", do_velnormal);

    // Auxiliary variables
    nAuxVar = pp.countval("Aux_Variables");
    Vector<std::string> AuxVar(nAuxVar);
    for (int ivar = 0; ivar < nAuxVar; ++ivar) {
      pp.get("Aux_Variables", AuxVar[ivar], ivar);
    }

    Print() << "infile = " << plotFileName << "\n";
    Print() << "reading plt file = " << plotFileName << "\n";

    PlotFileData pf(plotFileName);

    // Plotfile global infos
    finestLevel = std::min(finestLevel, pf.finestLevel());
    int Nlev = finestLevel + 1;
    const Vector<std::string>& plotVarNames = pf.varNames();
    RealBox rb(&(pf.probLo()[0]), &(pf.probHi()[0]));

    auto it = std::find(plotVarNames.begin(), plotVarNames.end(), progressName);
    if (it == plotVarNames.end()) {
      amrex::Abort("Wrong progress variable name: " + progressName);
    }
    int idC = std::distance(plotVarNames.begin(), it);

    // ---------------------------------------------------------------------
    // Variables index management
    // ---------------------------------------------------------------------
    const int idCst = 0;
    const int idVst = idCst + 1;
    int nCompIn = idVst;

    Vector<std::string> inVarNames(nCompIn);
    inVarNames[idCst] = plotVarNames[idC];

    if (do_strain) {
      nCompIn += AMREX_SPACEDIM;
      inVarNames.resize(nCompIn);
      inVarNames[idVst + 0] = "x_velocity";
      inVarNames[idVst + 1] = "y_velocity";
#if AMREX_SPACEDIM == 3
      inVarNames[idVst + 2] = "z_velocity";
#endif
    }

    if (nAuxVar > 0) {
      inVarNames.resize(nCompIn + nAuxVar);
      for (int ivar = 0; ivar < nAuxVar; ++ivar) {
        auto it =
          std::find(plotVarNames.begin(), plotVarNames.end(), AuxVar[ivar]);
        if (it == plotVarNames.end()) {
          amrex::Abort("Unknown auxiliary variable name: " + AuxVar[ivar]);
        }
        inVarNames[nCompIn] = AuxVar[ivar];
        nCompIn++;
      }
    }

    Vector<int> destFillComps(nCompIn);
    for (int i = 0; i < nCompIn; ++i) {
      destFillComps[i] = i;
    }

    // Start appending variables that we will compute
    const int idProg = nCompIn;      // progress variable
    const int idSmProg = idProg + 1; // smoothed progress variable
    const int idKm = idSmProg + 1;   // Mean curvature
    const int idN = idKm + 1;        // Flame normal components
    int idSR = -1;                   // Strain rate
    int nCompOut = 0;                // Total number of output variables

#if AMREX_SPACEDIM > 2
    const int idKg = idN + AMREX_SPACEDIM; // Gaussian curvature
    if (do_strain) {
      idSR = idKg + 1;
      nCompOut = idSR + 1;
    } else {
      nCompOut = idKg + 1;
    }
#else
    if (do_strain) {
      idSR = idN + AMREX_SPACEDIM;
      nCompOut = idSR + 1;
    } else {
      nCompOut = idN + AMREX_SPACEDIM;
    }
#endif

    int idROST = -1;
    if (getStrainTensor) {
      idROST = nCompOut;
      nCompOut = idROST + AMREX_SPACEDIM * AMREX_SPACEDIM; // Rate-of-strain
    }

    int idVelNormal = -1;
    if (do_velnormal) {
      idVelNormal = nCompOut;
      nCompOut += 1;
    }

    if (verbose) {
      Print() << "Will read the following variables: ";
      for (int i = 0; i < nCompIn; ++i) {
        Print() << " " << inVarNames[i];
      }
      Print() << '\n';
      Print() << "States out will be those plus: " << '\n';
      Print() << "   idProg:   " << idProg << '\n';
      Print() << "   idSmProg: " << idSmProg << '\n';
      Print() << "   idKm:     " << idKm << '\n';
#if AMREX_SPACEDIM > 2
      Print() << "   idKg:     " << idKg << '\n';
#endif
      if (do_strain) {
        Print() << "   idSR:     " << idSR << '\n';
      }
    }

    // Check symmetry/periodicity in given coordinate direction
    Vector<int> sym_dir(AMREX_SPACEDIM, 0);
    pp.queryarr("sym_dir", sym_dir, 0, AMREX_SPACEDIM);

    // Periodicity has no safe default: it must match the simulation
    Vector<int> is_per(AMREX_SPACEDIM, 0);
    pp.getarr("is_per", is_per, 0, AMREX_SPACEDIM);

    Print() << "Periodicity assumed for this case: ";
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      Print() << is_per[i] << " ";
    }
    Print() << "\n";

    int coord = 0;

    // ---------------------------------------------------------------------
    // Let's start the real work
    // ---------------------------------------------------------------------
    Vector<MultiFab*> state(Nlev);
    Vector<MultiFab*> flame_normal(Nlev);
    Vector<MultiFab*> cell_normal(Nlev);
    Vector<Geometry*> geoms(Nlev);
    Vector<Geometry> geomsOP(Nlev);
    Vector<BoxArray> grids(Nlev);
    Vector<DistributionMapping> dmap(Nlev);
    const int nGrow = 2;

    // Read the required data from pltfile
    FArrayBox tmp;
    for (int lev = 0; lev < Nlev; ++lev) {
      const BoxArray ba = pf.boxArray(lev);
      grids[lev] = ba;
      dmap[lev] = pf.DistributionMap(lev);
      geoms[lev] = new Geometry(pf.probDomain(lev), &rb, coord, &(is_per[0]));
      geomsOP[lev] = *geoms[lev];

      state[lev] = new MultiFab(ba, dmap[lev], nCompOut, nGrow);
      state[lev]->setVal(0.0);
      flame_normal[lev] = new MultiFab(ba, dmap[lev], AMREX_SPACEDIM, nGrow);
      flame_normal[lev]->setVal(0.0);
      cell_normal[lev] = new MultiFab(ba, dmap[lev], AMREX_SPACEDIM, nGrow);
      cell_normal[lev]->setVal(0.0);

      // Get input state data
      if (verbose)
        Print() << "Reading data for level " << lev << "\n";
      for (int n = 0; n < inVarNames.size(); ++n) {
        const MultiFab& src = pf.get(lev, inVarNames[n]);
        MultiFab::Copy(*state[lev], src, 0, n, 1, 0);
      }
    }

    //---------------------------------------------------------------------------
    // Building EB
    //---------------------------------------------------------------------------

#ifdef AMREX_USE_EB

    Vector<std::unique_ptr<EBFArrayBoxFactory>> eb_factory =
      pele_analysis::buildEBFactories(geomsOP, grids, dmap, verbose);

#endif // AMREX_USE_EB

    //---------------------------------------------------------------------------
    // Finding min and max values in the domain
    //---------------------------------------------------------------------------
    if (useFileMinMax || floorIt) {
      if (useFileMinMax) {

        if (verbose)
          Print() << "Getting Minimum and Maximum values..." << std::endl;

        Real global_min = 1.0e20;
        Real global_max = -1.0e20;

        for (int lev = 0; lev < Nlev; ++lev) {
          if (verbose)
            Print() << "Starting MFIter to find min/max for Level " << lev
                    << std::endl;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          {
            Real thread_min = 1.0e20;
            Real thread_max = -1.0e20;

            for (MFIter mfi(*state[lev], TilingIfNotGPU()); mfi.isValid();
                 ++mfi) {
              const Box& bx = mfi.tilebox();
              auto const& varBox = state[lev]->array(mfi, idCst);
#ifdef AMREX_USE_EB
              auto const& volFracBox =
                eb_factory[lev]->getVolFrac()[mfi].array();
#endif

              // Min/Max inside of tile (should be GPU safe)
              ReduceOps<ReduceOpMin, ReduceOpMax> reduce_op;
              ReduceData<Real, Real> reduce_data(reduce_op);
              using ReduceTuple = amrex::GpuTuple<Real, Real>;

              reduce_op.eval(
                bx, reduce_data,
                [=] AMREX_GPU_DEVICE(
                  int i, int j, int k) noexcept -> ReduceTuple {
#ifdef AMREX_USE_EB
                  if (volFracBox(i, j, k) <= 0.0)
                    return {1.0e20, -1.0e20};
#endif
                  Real val = varBox(i, j, k, 0);
                  return {val, val};
                });

              auto [tile_min, tile_max] =
                reduce_data.value(); // Result of tile reduction
              thread_min = std::min(thread_min, tile_min);
              thread_max = std::max(thread_max, tile_max);
            } // MFIter

#ifdef AMREX_USE_OMP
#pragma omp critical
#endif
            {
              global_min = std::min(global_min, thread_min);
              global_max = std::max(global_max, thread_max);
            }
          } // End OMP parallel region
        } // End Level loop

        // Global reduction across all ranks
        amrex::ParallelDescriptor::ReduceRealMin(global_min);
        amrex::ParallelDescriptor::ReduceRealMax(global_max);

        progMin = global_min;
        progMax = global_max;
      }

      Print() << "progressName = " << progressName << " at index: " << idC
              << "\n";
      Print() << "useFileMinMax = " << useFileMinMax << "\n";
      Print() << "Min/Max = " << progMin << " / " << progMax << "\n";

      ParallelDescriptor::Barrier();
    }

    if (progMin >= progMax) {
      amrex::Abort("progMin must be less than progMax");
    }

    //---------------------------------------------------------------------------
    // Building progress Variable
    //---------------------------------------------------------------------------

    // Build progress variable from state at idCst, put into idProg
    for (int lev = 0; lev < Nlev; ++lev) {
      MultiFab StateVar(*state[lev], amrex::make_alias, idCst, 1);
      MultiFab ProgressVar(*state[lev], amrex::make_alias, idProg, 1);
      for (MFIter mfi(*state[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const auto& StateVarFab = StateVar.array(mfi);
        const auto& ProgVarFab = ProgressVar.array(mfi);

// Get volume fraction if an EB is used
#ifdef AMREX_USE_EB
        auto const& volFracBox = eb_factory[lev]->getVolFrac()[mfi].array();
#endif

        Real invdenom = 1.0 / (progMax - progMin);

        amrex::ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

// Making sure that progress variable outside of domain is set to 0
#ifdef AMREX_USE_EB
            if (volFracBox(i, j, k) > 0) {
              ProgVarFab(i, j, k) = (StateVarFab(i, j, k) - progMin) * invdenom;
            } else {
              ProgVarFab(i, j, k) = 0;
            }
#else
                ProgVarFab(i,j,k) = ( StateVarFab(i,j,k) - progMin ) * invdenom;
#endif
          });
      }
      ProgressVar.FillBoundary(geoms[lev]->periodicity());

      if (verbose)
        Print() << "Progress variable computed for level " << lev << "\n";

    } // End lev loop

    //---------------------------------------------------------------------------
    // Optional smoothing
    //---------------------------------------------------------------------------

    if (do_smooth) {
      if (verbose)
        Print() << "Start smoothing of ProgVar!" << std::endl;
      // Set a composite ABec solve to smooth the progress variable
      // Try solving c^{n+1} - ∆t \nabla \cdot b_i \nabla c^{n+1} = c^{n)
      // In ABec terms  (\alpha A - \beta \nabla \cdot B_i \nabla) \phi = rhs
      // \alpha = 1, A = I
      // \beta = ∆t, b_i = ?? let's start with 1, switch to a D_c later if need
      // be. Adapt ∆t accordingly ... rhs = c^{n}

      LPInfo info;
      info.setAgglomeration(1);
      info.setConsolidation(1);
      info.setMetricTerm(false);

#ifdef AMREX_USE_EB
      MLEBABecLap mlabec(
        geomsOP, grids, dmap, info, amrex::GetVecOfConstPtrs(eb_factory));
#else
      MLABecLaplacian mlabec(geomsOP, grids, dmap, info);
#endif

      mlabec.setMaxOrder(4);

      const Real tol_rel = 1.e-12; // 1.e-12
      const Real tol_abs = 1.e-12;

      // Problem with Periodic or Neumann BC
      std::array<LinOpBCType, AMREX_SPACEDIM> lo_bc;
      std::array<LinOpBCType, AMREX_SPACEDIM> hi_bc;
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
        if (is_per[idim] == 1) {
          lo_bc[idim] = hi_bc[idim] = LinOpBCType::Periodic;
        } else {
          lo_bc[idim] = hi_bc[idim] = LinOpBCType::Neumann;
        }
      }
      mlabec.setDomainBC(lo_bc, hi_bc);

      for (int lev = 0; lev < Nlev; ++lev) {
        // for problem with homogeneous Neumann BC, we need to pass nullptr
        mlabec.setLevelBC(lev, nullptr);
      }
      Real alpha = 1.0;
      Real beta = smooth_time;
      mlabec.setScalars(alpha, beta);

      for (int lev = 0; lev < Nlev; ++lev) {
        // Set A = I  at each level
        MultiFab acoef(
          state[lev]->boxArray(), state[lev]->DistributionMap(), 1, 0);
        acoef.setVal(1.0);
        mlabec.setACoeffs(lev, acoef);

        // Set b_i = 1.0  at each level
        Array<MultiFab, AMREX_SPACEDIM> face_bcoef;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
          const BoxArray& ba = amrex::convert(
            state[lev]->boxArray(), IntVect::TheDimensionVector(idim));

          face_bcoef[idim].define(ba, state[lev]->DistributionMap(), 1, 0);
          face_bcoef[idim].setVal(1); // was 1.0 1e-30
        }
        mlabec.setBCoeffs(lev, amrex::GetArrOfConstPtrs(face_bcoef));
      }

      if (verbose)
        Print() << "Building Solution and RHS MultiFabs!" << std::endl;
      Vector<MultiFab> solution(Nlev);
      Vector<MultiFab> rhs(Nlev);
      for (int lev = 0; lev < Nlev; ++lev) {
        rhs[lev].define(grids[lev], dmap[lev], 1, 0);
        solution[lev].define(grids[lev], dmap[lev], 1, nGrow);
        MultiFab::Copy(rhs[lev], *state[lev], idProg, 0, 1, 0);
        solution[lev].setVal(0.0);
      }

      if (verbose)
        Print() << "Building MLMG for Smoothing!" << std::endl;
      MLMG mlmg(mlabec);
      mlmg.setMaxIter(100);
      mlmg.setVerbose(1);
      mlmg.solve(
        GetVecOfPtrs(solution), GetVecOfConstPtrs(rhs), tol_rel, tol_abs);

      for (int lev = 0; lev < Nlev; ++lev) {
        MultiFab::Copy(*state[lev], solution[lev], 0, idSmProg, 1, nGrow);
        state[lev]->FillBoundary(idSmProg, 1, geoms[lev]->periodicity());
      }
      if (verbose)
        Print() << "Progress variable smoothed successfully \n";
    } // do_smooth

    int idprogvar = do_smooth ? idSmProg : idProg;

    //---------------------------------------------------------------------------
    // Compute curvature
    //---------------------------------------------------------------------------

#ifdef AMREX_USE_EB
    LPInfo info_apply;
    info_apply.setMaxCoarseningLevel(0);
    // info_apply.setAgglomeration(1);
    // info_apply.setConsolidation(1);
    info_apply.setMetricTerm(false);
#else
    LPInfo info;
    info.setAgglomeration(1);
    info.setConsolidation(1);
    info.setMetricTerm(false);
    info.setMaxCoarseningLevel(0);
#endif

    for (int lev = 0; lev < Nlev; ++lev) {
      const BoxArray ba = pf.boxArray(lev);

      if (verbose)
        Print() << "Starting computation of mean curvature on level " << lev
                << "\n";

// Get face gradients of progress variable
#ifdef AMREX_USE_EB
      MLEBABecLap poisson(
        {*geoms[lev]}, {ba}, {dmap[lev]}, info_apply, {eb_factory[lev].get()});
      poisson.setVerbose(4);
      poisson.setMaxOrder(4);

      // Poisson like solver for EB
      poisson.setScalars(0.0, 1.0);
      poisson.setBCoeffs(0, -1.0); // Important to set only for lev 0 as solver
                                   // is initialised with only 1 level
#else
      //  Compute curvature using LinearOperators
      // Set-up Poisson Linear Solver
      MLPoisson poisson({*geoms[lev]}, {ba}, {dmap[lev]}, info);
      poisson.setMaxOrder(4);
#endif

      std::array<LinOpBCType, AMREX_SPACEDIM> lo_bc;
      std::array<LinOpBCType, AMREX_SPACEDIM> hi_bc;
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
        if (is_per[idim] == 1) {
          lo_bc[idim] = hi_bc[idim] = LinOpBCType::Periodic;
        } else {
          if (sym_dir[idim] == 1) {
            lo_bc[idim] = hi_bc[idim] = LinOpBCType::reflect_odd;
          } else {
            lo_bc[idim] = hi_bc[idim] = LinOpBCType::Neumann;
          }
        }
      }
      poisson.setDomainBC(lo_bc, hi_bc);

      // Coarse BC data must outlive setLevelBC, which reads it
      MultiFab ProgVarCoarse;
      if (lev > 0) {
        ProgVarCoarse.define(
          state[lev - 1]->boxArray(), state[lev - 1]->DistributionMap(), 1,
          state[lev - 1]->nGrow());
        MultiFab::Copy(ProgVarCoarse, *state[lev - 1], idprogvar, 0, 1, nGrow);
        poisson.setCoarseFineBC(&ProgVarCoarse, pf.refRatio(lev - 1));
      }

      MultiFab ProgVar(ba, dmap[lev], 1, nGrow);
      MultiFab::Copy(ProgVar, *state[lev], idprogvar, 0, 1, nGrow);

      poisson.setLevelBC(0, &ProgVar);

      MLMG mlmg(poisson);

      std::array<MultiFab, AMREX_SPACEDIM> face_gradient;
#ifdef AMREX_USE_EB
      AMREX_D_TERM(face_gradient[0].define(
        convert(ba, IntVect::TheDimensionVector(0)), dmap[lev], 1, 0, MFInfo(),
        *eb_factory[lev]);
                   , face_gradient[1].define(
                       convert(ba, IntVect::TheDimensionVector(1)), dmap[lev],
                       1, 0, MFInfo(), *eb_factory[lev]);
                   , face_gradient[2].define(
                       convert(ba, IntVect::TheDimensionVector(2)), dmap[lev],
                       1, 0, MFInfo(), *eb_factory[lev]););

      mlmg.getFluxes(
        {amrex::GetArrOfPtrs(face_gradient)}, {&ProgVar},
        MLMG::Location::FaceCentroid);
#else
      AMREX_D_TERM(
        face_gradient[0].define(
          convert(ba, IntVect::TheDimensionVector(0)), dmap[lev], 1, 0);
        , face_gradient[1].define(
            convert(ba, IntVect::TheDimensionVector(1)), dmap[lev], 1, 0);
        , face_gradient[2].define(
            convert(ba, IntVect::TheDimensionVector(2)), dmap[lev], 1, 0););

      mlmg.getFluxes(
        {amrex::GetArrOfPtrs(face_gradient)}, {&ProgVar},
        MLMG::Location::FaceCenter);

#endif

      // New Multifab cell avg gradient
      MultiFab cellavg_gradient(ba, dmap[lev], AMREX_SPACEDIM, 0);
      cellavg_gradient.setVal(0.0);

#ifdef AMREX_USE_EB
      // MLEBABecLap with beta = 1, B = -1 returns +grad(C) as flux, whereas
      // MLPoisson returns -grad(C); only the latter needs the sign flip.
      EB_average_face_to_cellcenter(
        cellavg_gradient, 0, amrex::GetArrOfConstPtrs(face_gradient));
#else
      average_face_to_cellcenter(
        cellavg_gradient, 0, amrex::GetArrOfConstPtrs(face_gradient));
      cellavg_gradient.mult(-1.0, 0, AMREX_SPACEDIM);
#endif

      // Compute ||\nabla c||

      MultiFab cellnorm_gradient(ba, dmap[lev], 1, 1);
      cellnorm_gradient.setVal(0.0);

      for (MFIter mfi(cellnorm_gradient, TilingIfNotGPU()); mfi.isValid();
           ++mfi) {
        const Box& bx = mfi.validbox();
        AMREX_D_TERM(auto const& Cx = cellavg_gradient.array(mfi, 0);
                     , auto const& Cy = cellavg_gradient.array(mfi, 1);
                     , auto const& Cz = cellavg_gradient.array(mfi, 2););
        const auto& normgrad = cellnorm_gradient.array(mfi);
        const auto& progvar = ProgVar.array(mfi);

#ifdef AMREX_USE_EB
        auto const& volFracBox = eb_factory[lev]->getVolFrac()[mfi].array();
#endif

        amrex::ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

#ifdef AMREX_USE_EB
            if (volFracBox(i, j, k) > 0) {
              normgrad(i, j, k) = std::max(
                1e-14,
                std::sqrt(AMREX_D_TERM(
                  std::pow(Cx(i, j, k), 2.0), +std::pow(Cy(i, j, k), 2.0),
                  +std::pow(Cz(i, j, k), 2.0))));
              normgrad(i, j, k) = -normgrad(i, j, k);
            }
#else
                normgrad(i,j,k) = std::max(1e-14, std::sqrt( AMREX_D_TERM (   std::pow(Cx(i,j,k),2.0),
                                                                            + std::pow(Cy(i,j,k),2.0),
                                                                            + std::pow(Cz(i,j,k),2.0)) ) ) ;

                normgrad(i,j,k) = -normgrad(i,j,k);
#endif
          });
      }
      cellnorm_gradient.FillBoundary(geoms[lev]->periodicity());

      // Get the cell centered flame normal N_i (pointing toward fresh gases)
      MultiFab::Copy(
        *cell_normal[lev], cellavg_gradient, 0, 0, AMREX_SPACEDIM, 0);
      cell_normal[lev]->FillBoundary(
        0, AMREX_SPACEDIM, geoms[lev]->periodicity());

      if (verbose)
        Print() << "Done with flame normal on level " << lev << "\n";

      //      At this point, I got the flame normal from clean gradients
      //      provided by the MLMG Now try to compute the divergence of the
      //      flame normal using another linear solver: 1 -> get the clean face
      //      gradient of flame normal components from MLMG: d N_i/d x_i 2 ->
      //      manually get the divergence from the gradients: ∑_i d N_i/d x_i

      //      Copy into level aware flame_normal MF and fill same level ghost
      //      cells on normal
      MultiFab::Copy(
        *flame_normal[lev], *cell_normal[lev], 0, 0, AMREX_SPACEDIM, 0);

      // Flame normal n = - grad(C)/norm(grad(C). norm(grad(C)) is only zero in
      // EB-covered cells (it is floored at 1e-14 elsewhere); the normal is set
      // to zero there instead of dividing 0/0.
      for (MFIter mfi(*flame_normal[lev], TilingIfNotGPU()); mfi.isValid();
           ++mfi) {
        const Box& bx = mfi.tilebox();
        const auto& normal = flame_normal[lev]->array(mfi);
        const auto& normgrad = cellnorm_gradient.const_array(mfi);
        amrex::ParallelFor(
          bx, AMREX_SPACEDIM,
          [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            normal(i, j, k, n) = (normgrad(i, j, k) != 0.0)
                                   ? normal(i, j, k, n) / normgrad(i, j, k)
                                   : 0.0;
          });
      }
      flame_normal[lev]->FillBoundary(
        0, AMREX_SPACEDIM, geoms[lev]->periodicity());

      //      Define curvature
      MultiFab Curv(ba, dmap[lev], 1, 0);
      Curv.setVal(0.0);

      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {

#ifdef AMREX_USE_EB
        MLEBABecLap poisson2(
          {*geoms[lev]}, {ba}, {dmap[lev]}, info_apply,
          {eb_factory[lev].get()});
        poisson2.setVerbose(4);
        poisson2.setMaxOrder(4);

        // Poisson like solver for EB
        poisson2.setScalars(0.0, 1.0);
        poisson2.setBCoeffs(0, -1.0);

#else
        MLPoisson poisson2({*geoms[lev]}, {ba}, {dmap[lev]}, info);
        poisson2.setMaxOrder(4);
#endif

        poisson2.setDomainBC(lo_bc, hi_bc);
        // Coarse BC data must outlive setLevelBC, which reads it
        MultiFab FlameNormalIdimCoarse;
        if (lev > 0) {
          FlameNormalIdimCoarse.define(
            flame_normal[lev - 1]->boxArray(),
            flame_normal[lev - 1]->DistributionMap(), 1, 0);
          MultiFab::Copy(
            FlameNormalIdimCoarse, *flame_normal[lev - 1], idim, 0, 1, 0);
          poisson2.setCoarseFineBC(
            &FlameNormalIdimCoarse, pf.refRatio(lev - 1));
        }

        MultiFab FlameNormalIdim(ba, dmap[lev], 1, 1);
        MultiFab::Copy(FlameNormalIdim, *flame_normal[lev], idim, 0, 1, 1);
        poisson2.setLevelBC(0, &FlameNormalIdim);

        MLMG mlmg2(poisson2);

        // Get the fluxes : d N_i / d x_j   , i = idim, j = 0, 1 (,2)
        std::array<MultiFab, AMREX_SPACEDIM> faceg;

#ifdef AMREX_USE_EB
        AMREX_D_TERM(faceg[0].define(
          convert(ba, IntVect::TheDimensionVector(0)), dmap[lev], 1, 0,
          MFInfo(), *eb_factory[lev]);
                     , faceg[1].define(
                         convert(ba, IntVect::TheDimensionVector(1)), dmap[lev],
                         1, 0, MFInfo(), *eb_factory[lev]);
                     , faceg[2].define(
                         convert(ba, IntVect::TheDimensionVector(2)), dmap[lev],
                         1, 0, MFInfo(), *eb_factory[lev]););

        mlmg2.getFluxes(
          {amrex::GetArrOfPtrs(faceg)}, {&FlameNormalIdim},
          MLMG::Location::FaceCentroid);
#else
        AMREX_D_TERM(
          faceg[0].define(
            convert(ba, IntVect::TheDimensionVector(0)), dmap[lev], 1, 0);
          , faceg[1].define(
              convert(ba, IntVect::TheDimensionVector(1)), dmap[lev], 1, 0);
          , faceg[2].define(
              convert(ba, IntVect::TheDimensionVector(2)), dmap[lev], 1, 0););

        mlmg2.getFluxes(
          {amrex::GetArrOfPtrs(faceg)}, {&FlameNormalIdim},
          MLMG::Location::FaceCenter);
#endif

        // Get cell centered d N_i / d x_y
        MultiFab cell_avgg(ba, dmap[lev], AMREX_SPACEDIM, 0);
        cell_avgg.setVal(0.0);

#ifdef AMREX_USE_EB
        EB_average_face_to_cellcenter(
          cell_avgg, 0, amrex::GetArrOfConstPtrs(faceg));
#else
        average_face_to_cellcenter(
          cell_avgg, 0, amrex::GetArrOfConstPtrs(faceg));
        cell_avgg.mult(-1.0, 0, AMREX_SPACEDIM);
#endif

        // Add d N_i / d x_i to curvature
        MultiFab::Add(
          Curv, cell_avgg, idim, 0, 1,
          0); // Adds cell_avgg to Curv starting from location idim to 0 in Curv
              // and only adds 1 comp + 0 ghost cells
      }

#if AMREX_SPACEDIM == 3
      // Mean curvature : 0.5 * \div \cdot n
      // TODO: I only need to do that in 3D ... right ?
      Curv.mult(0.5, 0, 1);
#endif

      // Clip curvature & flame normal for c < threshold or c > 1.0-threshold
      for (MFIter mfi(Curv, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        const auto& CurvFab = Curv.array(mfi);
        AMREX_D_TERM(const auto& FnormXFab = flame_normal[lev]->array(mfi, 0);
                     , const auto& FnormYFab = flame_normal[lev]->array(mfi, 1);
                     ,
                     const auto& FnormZFab = flame_normal[lev]->array(mfi, 2););
        const auto& progvar = ProgVar.array(mfi);
        amrex::ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (
              do_threshold && (progvar(i, j, k) < threshold ||
                               progvar(i, j, k) > 1.0 - threshold)) {
              CurvFab(i, j, k) = 0.0;
              AMREX_D_TERM(FnormXFab(i, j, k) = 0.0;, FnormYFab(i, j, k) = 0.0;
                           , FnormZFab(i, j, k) = 0.0;);
            }
          });
      }

      MultiFab::Copy(*state[lev], Curv, 0, idKm, 1, 0);
      MultiFab::Copy(
        *state[lev], *flame_normal[lev], 0, idN, AMREX_SPACEDIM, 0);

      if (verbose)
        Print() << "Mean curvature has been computed on level " << lev << "\n";

        //---------------------------------------------------------------------------
        // Compute gaussian curvature
        //---------------------------------------------------------------------------

        // Now work on the gaussian curvature: only if 3D and required
#if AMREX_SPACEDIM == 3
      if (do_gaussCurv) {

        // Start by getting the Hessian of the progress variable
        MultiFab Hessian(ba, dmap[lev], 9, 0);

        // Compute grad of grad in each dim and store in Hessian
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {

#ifdef AMREX_USE_EB
          MLEBABecLap poisson2(
            {*geoms[lev]}, {ba}, {dmap[lev]}, info_apply,
            {eb_factory[lev].get()});
          poisson2.setVerbose(4);
          poisson2.setMaxOrder(4);

          // Poisson like solver for EB
          poisson2.setScalars(0.0, 1.0);
          poisson2.setBCoeffs(0, -1.0);

#else
          MLPoisson poisson2({*geoms[lev]}, {ba}, {dmap[lev]}, info);
          poisson2.setMaxOrder(4);
#endif

          poisson2.setDomainBC(lo_bc, hi_bc);
          // Coarse BC data must outlive setLevelBC, which reads it
          MultiFab gradIdimCoarse;
          if (lev > 0) {
            gradIdimCoarse.define(
              cell_normal[lev - 1]->boxArray(),
              cell_normal[lev - 1]->DistributionMap(), 1, 0);
            MultiFab::Copy(
              gradIdimCoarse, *cell_normal[lev - 1], idim, 0, 1, 0);
            poisson2.setCoarseFineBC(&gradIdimCoarse, pf.refRatio(lev - 1));
          }
          MultiFab gradIdim(ba, dmap[lev], 1, 1);
          MultiFab::Copy(gradIdim, *cell_normal[lev], idim, 0, 1, 1);
          poisson2.setLevelBC(0, &gradIdim);

          MLMG mlmg2(poisson2);

          std::array<MultiFab, AMREX_SPACEDIM> faceg;

#ifdef AMREX_USE_EB
          AMREX_D_TERM(faceg[0].define(
            convert(ba, IntVect::TheDimensionVector(0)), dmap[lev], 1, 0,
            MFInfo(), *eb_factory[lev]);
                       , faceg[1].define(
                           convert(ba, IntVect::TheDimensionVector(1)),
                           dmap[lev], 1, 0, MFInfo(), *eb_factory[lev]);
                       , faceg[2].define(
                           convert(ba, IntVect::TheDimensionVector(2)),
                           dmap[lev], 1, 0, MFInfo(), *eb_factory[lev]););

          mlmg2.getFluxes(
            {amrex::GetArrOfPtrs(faceg)}, {&gradIdim},
            MLMG::Location::FaceCentroid);

#else
          AMREX_D_TERM(
            faceg[0].define(
              convert(ba, IntVect::TheDimensionVector(0)), dmap[lev], 1, 0);
            , faceg[1].define(
                convert(ba, IntVect::TheDimensionVector(1)), dmap[lev], 1, 0);
            , faceg[2].define(
                convert(ba, IntVect::TheDimensionVector(2)), dmap[lev], 1, 0););
          mlmg2.getFluxes({amrex::GetArrOfPtrs(faceg)}, {&gradIdim});
#endif

          // Get cell centered d C / d idim _x_y_z
          MultiFab cell_avgg(ba, dmap[lev], AMREX_SPACEDIM, 0);

#ifdef AMREX_USE_EB
          EB_average_face_to_cellcenter(
            cell_avgg, 0, amrex::GetArrOfConstPtrs(faceg));
#else
          average_face_to_cellcenter(
            cell_avgg, 0, amrex::GetArrOfConstPtrs(faceg));
          cell_avgg.mult(
            -1.0, 0, AMREX_SPACEDIM); // Only necessary for non EB cases
#endif

          // Copy in Hessian
          MultiFab::Copy(Hessian, cell_avgg, 0, (idim) * 3, AMREX_SPACEDIM, 0);
        }

        // Get the adjugate of the Hessian
        MultiFab AdjHessian(ba, dmap[lev], 9, 0);

        for (MFIter mfi(Hessian, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.validbox();
          const auto& HxiFab = Hessian.array(mfi, 0); // Cxx, Cxy, Cxz
          const auto& HyiFab = Hessian.array(mfi, 3); // Cyx, Cyy, Cyz
          const auto& HziFab = Hessian.array(mfi, 6); // Czx, Czy, Czz
          const auto& AdjHxiFab = AdjHessian.array(mfi, 0);
          const auto& AdjHyiFab = AdjHessian.array(mfi, 3);
          const auto& AdjHziFab = AdjHessian.array(mfi, 6);
          amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              AdjHxiFab(i, j, k, 0) = HyiFab(i, j, k, 1) * HziFab(i, j, k, 2) -
                                      HziFab(i, j, k, 1) * HyiFab(i, j, k, 2);
              AdjHyiFab(i, j, k, 0) = HyiFab(i, j, k, 2) * HziFab(i, j, k, 0) -
                                      HziFab(i, j, k, 2) * HyiFab(i, j, k, 0);
              AdjHziFab(i, j, k, 0) = HyiFab(i, j, k, 0) * HziFab(i, j, k, 1) -
                                      HziFab(i, j, k, 0) * HyiFab(i, j, k, 1);
              AdjHxiFab(i, j, k, 1) = HxiFab(i, j, k, 2) * HziFab(i, j, k, 1) -
                                      HziFab(i, j, k, 2) * HxiFab(i, j, k, 1);
              AdjHyiFab(i, j, k, 1) = HxiFab(i, j, k, 0) * HziFab(i, j, k, 2) -
                                      HziFab(i, j, k, 0) * HxiFab(i, j, k, 2);
              AdjHziFab(i, j, k, 1) = HxiFab(i, j, k, 1) * HziFab(i, j, k, 0) -
                                      HziFab(i, j, k, 1) * HxiFab(i, j, k, 0);
              AdjHxiFab(i, j, k, 2) = HxiFab(i, j, k, 1) * HyiFab(i, j, k, 2) -
                                      HyiFab(i, j, k, 1) * HxiFab(i, j, k, 2);
              AdjHyiFab(i, j, k, 2) = HxiFab(i, j, k, 2) * HyiFab(i, j, k, 0) -
                                      HyiFab(i, j, k, 2) * HxiFab(i, j, k, 0);
              AdjHziFab(i, j, k, 2) = HxiFab(i, j, k, 0) * HyiFab(i, j, k, 1) -
                                      HyiFab(i, j, k, 0) * HxiFab(i, j, k, 1);
            });
        }

        // Now get the gausian curvature
        MultiFab gCurv(ba, dmap[lev], 1, 0);
        for (MFIter mfi(gCurv, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.validbox();
          const auto& progvar = ProgVar.array(mfi);
          const auto& gCurvFab = gCurv.array(mfi);
          const auto& CgradNorm = cellnorm_gradient.array(mfi);
          const auto& CxFab = cellavg_gradient.array(mfi, 0);
          const auto& CyFab = cellavg_gradient.array(mfi, 1);
          const auto& CzFab = cellavg_gradient.array(mfi, 2);
          const auto& AdjHxiFab = AdjHessian.array(mfi, 0);
          const auto& AdjHyiFab = AdjHessian.array(mfi, 3);
          const auto& AdjHziFab = AdjHessian.array(mfi, 6);
          amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              gCurvFab(i, j, k) =
                (CxFab(i, j, k) * (AdjHxiFab(i, j, k, 0) * CxFab(i, j, k) +
                                   AdjHxiFab(i, j, k, 1) * CyFab(i, j, k) +
                                   AdjHxiFab(i, j, k, 2) * CzFab(i, j, k)) +
                 CyFab(i, j, k) * (AdjHyiFab(i, j, k, 0) * CxFab(i, j, k) +
                                   AdjHyiFab(i, j, k, 1) * CyFab(i, j, k) +
                                   AdjHyiFab(i, j, k, 2) * CzFab(i, j, k)) +
                 CzFab(i, j, k) * (AdjHziFab(i, j, k, 0) * CxFab(i, j, k) +
                                   AdjHziFab(i, j, k, 1) * CyFab(i, j, k) +
                                   AdjHziFab(i, j, k, 2) * CzFab(i, j, k)));
              // norm(grad(C)) is zero only in EB-covered cells
              gCurvFab(i, j, k) =
                (CgradNorm(i, j, k) != 0.0)
                  ? gCurvFab(i, j, k) / std::pow(CgradNorm(i, j, k), 4.0)
                  : 0.0;
              if (
                do_threshold && (progvar(i, j, k) < threshold ||
                                 progvar(i, j, k) > 1.0 - threshold)) {
                gCurvFab(i, j, k) = 0.0;
              }
            });
        }
        MultiFab::Copy(*state[lev], gCurv, 0, idKg, 1, 0);
        if (verbose)
          Print() << "Gaussian curvature has been computed on level " << lev
                  << "\n";
      }
#endif

      //---------------------------------------------------------------------------
      // Compute strain
      //---------------------------------------------------------------------------

      if (do_strain) {
        // Strain rate -nn:\nabla u + \nabla \cdot u
        if (verbose)
          Print() << "Starting computation of tangential strain rate for "
                  << lev << "\n";

        // Start by building the strain tensor
        MultiFab StrainT(ba, dmap[lev], AMREX_SPACEDIM * AMREX_SPACEDIM, 0);
        StrainT.setVal(0.0);

        // Compute strain tensor with a MLMG in each direction
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {

#ifdef AMREX_USE_EB
          MLEBABecLap poisson2(
            {*geoms[lev]}, {ba}, {dmap[lev]}, info_apply,
            {eb_factory[lev].get()});
          poisson2.setVerbose(4);
          poisson2.setMaxOrder(4);

          // Poisson like solver for EB
          poisson2.setScalars(0.0, 1.0);
          poisson2.setBCoeffs(0, -1.0);

#else
          MLPoisson poisson2({*geoms[lev]}, {ba}, {dmap[lev]}, info);
          poisson2.setMaxOrder(4);
#endif

          poisson2.setDomainBC(lo_bc, hi_bc);

          // Coarse BC data must outlive setLevelBC, which reads it
          MultiFab velIdimCoarse;
          if (lev > 0) {
            velIdimCoarse.define(
              state[lev - 1]->boxArray(), state[lev - 1]->DistributionMap(), 1,
              0);
            MultiFab::Copy(
              velIdimCoarse, *state[lev - 1], idVst + idim, 0, 1, 0);
            poisson2.setCoarseFineBC(&velIdimCoarse, pf.refRatio(lev - 1));
          }

          MultiFab velIdim(ba, dmap[lev], 1, 1);
          MultiFab::Copy(velIdim, *state[lev], idVst + idim, 0, 1, 1);
          poisson2.setLevelBC(0, &velIdim);

          MLMG mlmg2(poisson2);

          std::array<MultiFab, AMREX_SPACEDIM> faceg;

#ifdef AMREX_USE_EB
          AMREX_D_TERM(faceg[0].define(
            convert(ba, IntVect::TheDimensionVector(0)), dmap[lev], 1, 0,
            MFInfo(), *eb_factory[lev]);
                       , faceg[1].define(
                           convert(ba, IntVect::TheDimensionVector(1)),
                           dmap[lev], 1, 0, MFInfo(), *eb_factory[lev]);
                       , faceg[2].define(
                           convert(ba, IntVect::TheDimensionVector(2)),
                           dmap[lev], 1, 0, MFInfo(), *eb_factory[lev]););

          mlmg2.getFluxes(
            {amrex::GetArrOfPtrs(faceg)}, {&velIdim},
            MLMG::Location::FaceCentroid);

          // Get cell centered d u_idim / d _x_y(_z)
          MultiFab cell_avgg(
            ba, dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), *eb_factory[lev]);
          cell_avgg.setVal(0.0);
          EB_average_face_to_cellcenter(
            cell_avgg, 0, amrex::GetArrOfConstPtrs(faceg));
#else
          AMREX_D_TERM(
            faceg[0].define(
              convert(ba, IntVect::TheDimensionVector(0)), dmap[lev], 1, 0);
            , faceg[1].define(
                convert(ba, IntVect::TheDimensionVector(1)), dmap[lev], 1, 0);
            , faceg[2].define(
                convert(ba, IntVect::TheDimensionVector(2)), dmap[lev], 1, 0););

          mlmg2.getFluxes({amrex::GetArrOfPtrs(faceg)}, {&velIdim});

          // Get cell centered d u_idim / d _x_y(_z)
          MultiFab cell_avgg(ba, dmap[lev], AMREX_SPACEDIM, 0);
          cell_avgg.setVal(0.0);
          average_face_to_cellcenter(
            cell_avgg, 0, amrex::GetArrOfConstPtrs(faceg));
          cell_avgg.mult(-1.0, 0, AMREX_SPACEDIM);

#endif

          // Copy in strain tensor
          MultiFab::Copy(
            StrainT, cell_avgg, 0, (idim)*AMREX_SPACEDIM, AMREX_SPACEDIM, 0);
        }

        // Gather the components of strain rate
        MultiFab strainrate(ba, dmap[lev], 1, 0);
        strainrate.setVal(0.0);

        for (MFIter mfi(strainrate, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.validbox();
          const auto& progvar = ProgVar.array(mfi);
          const auto& srFab = strainrate.array(mfi);
          AMREX_D_TERM(const auto& NxFab = flame_normal[lev]->array(mfi, 0);
                       , const auto& NyFab = flame_normal[lev]->array(mfi, 1);
                       , const auto& NzFab = flame_normal[lev]->array(mfi, 2););
          AMREX_D_TERM(
            const auto& gradUx = StrainT.array(mfi, 0);
            , const auto& gradUy = StrainT.array(mfi, AMREX_SPACEDIM);
            , const auto& gradUz = StrainT.array(mfi, AMREX_SPACEDIM * 2););
          amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              srFab(i, j, k) = AMREX_D_TERM(
                AMREX_D_TERM(
                  -gradUx(i, j, k, 0) * NxFab(i, j, k) * NxFab(i, j, k),
                  -gradUx(i, j, k, 1) * NxFab(i, j, k) * NyFab(i, j, k),
                  -gradUx(i, j, k, 2) * NxFab(i, j, k) * NzFab(i, j, k)),
                AMREX_D_TERM(
                  -gradUy(i, j, k, 0) * NyFab(i, j, k) * NxFab(i, j, k),
                  -gradUy(i, j, k, 1) * NyFab(i, j, k) * NyFab(i, j, k),
                  -gradUy(i, j, k, 2) * NyFab(i, j, k) * NzFab(i, j, k)),
                AMREX_D_TERM(
                  -gradUz(i, j, k, 0) * NzFab(i, j, k) * NxFab(i, j, k),
                  -gradUz(i, j, k, 1) * NzFab(i, j, k) * NyFab(i, j, k),
                  -gradUz(i, j, k, 2) * NzFab(i, j, k) *
                    NzFab(i, j, k))); //-nn:\nabla u
              srFab(i, j, k) += AMREX_D_TERM(
                +gradUx(i, j, k, 0), +gradUy(i, j, k, 1),
                +gradUz(i, j, k, 2)); // + \nabla \cdot u
            });
        }

        MultiFab::Copy(*state[lev], strainrate, 0, idSR, 1, 0);
        if (verbose)
          Print() << "Tangential strain rate has been computed on level " << lev
                  << "\n";

        // Store the strain rate tensor if required
        if (getStrainTensor) {
          MultiFab::Copy(
            *state[lev], StrainT, 0, idROST, AMREX_SPACEDIM * AMREX_SPACEDIM,
            0);
        }
      }

      if (do_velnormal) {
        // Velocity normal to the flame front
        MultiFab velNormal(ba, dmap[lev], 1, 0);

        for (MFIter mfi(velNormal, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.validbox();
          const auto& progvar = ProgVar.array(mfi);
          const auto& velNormFab = velNormal.array(mfi);
          AMREX_D_TERM(const auto& NxFab = flame_normal[lev]->array(mfi, 0);
                       , const auto& NyFab = flame_normal[lev]->array(mfi, 1);
                       , const auto& NzFab = flame_normal[lev]->array(mfi, 2););
          AMREX_D_TERM(const auto& Ux = state[lev]->array(mfi, idVst + 0);
                       , const auto& Uy = state[lev]->array(mfi, idVst + 1);
                       , const auto& Uz = state[lev]->array(mfi, idVst + 2););
          amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              velNormFab(i, j, k) = AMREX_D_TERM(
                +Ux(i, j, k) * NxFab(i, j, k), +Uy(i, j, k) * NyFab(i, j, k),
                +Uz(i, j, k) * NzFab(i, j, k));
              if (
                do_threshold && (progvar(i, j, k) < threshold ||
                                 progvar(i, j, k) > 1.0 - threshold)) {
                velNormFab(i, j, k) = 0.0;
              }
            });
        }
        MultiFab::Copy(*state[lev], velNormal, 0, idVelNormal, 1, 0);
        if (verbose)
          Print()
            << "Flow velocity normal to the flame has been computed on level "
            << lev << "\n";
      }

    } // level-loop

    // ---------------------------------------------------------------------
    // Set-up the output
    // ---------------------------------------------------------------------
    Vector<std::string> nnames(nCompOut);

    // Plot file variables
    for (int i = 0; i < nCompIn; ++i) {
      nnames[i] = inVarNames[i];
    }

    // Computed variables
    nnames[idProg] = "Progress";
    nnames[idSmProg] = "SmoothedProgress";
    nnames[idKm] = "MeanCurvature_" + progressName;
    nnames[idN] = "FlameNormalX_" + progressName;
    nnames[idN + 1] = "FlameNormalY_" + progressName;
#if AMREX_SPACEDIM == 3
    nnames[idN + 2] = "FlameNormalZ_" + progressName;
    nnames[idKg] = "GaussianCurvature_" + progressName;
#endif
    if (do_strain)
      nnames[idSR] = "StrainRate_" + progressName;

    if (getStrainTensor) {
      std::string dirChar[3] = {"x", "y", "z"};
      for (int i = 0; i < AMREX_SPACEDIM * AMREX_SPACEDIM; ++i) {
        int n = i % AMREX_SPACEDIM;
        int m = i / AMREX_SPACEDIM;
        nnames[idROST + i] = "ROST_dU" + dirChar[m] + "d" + dirChar[n];
      }
    }

    if (do_velnormal) {
      nnames[idVelNormal] = "VelFlameNormal";
    }

    // Write to plotfile
    // Remove GC from outstate
    Vector<MultiFab*> ostate(Nlev);
    for (int lev = 0; lev < Nlev; ++lev) {
      const BoxArray ba = state[lev]->boxArray();
      ostate[lev] = new MultiFab(ba, dmap[lev], nCompOut, 0);
      MultiFab::Copy(*ostate[lev], *state[lev], 0, 0, nCompOut, 0);
    }
    Real time = pf.time();
    Print() << "Writing new data to " << outfile << " , for TS: " << time
            << "\n";
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev - 1, IntVect(2));
    for (int lev = 0; lev < Nlev - 1; ++lev) {
      refRatios[lev] = IntVect(pf.refRatio(lev));
    }
    VisMF::SetNOutFiles(n_files);
    amrex::WriteMultiLevelPlotfile(
      outfile, Nlev, GetVecOfConstPtrs(ostate), nnames, geomsOP, time, isteps,
      refRatios);
  }
  ParallelDescriptor::Barrier();
  amrex::Finalize();
  return 0;
}
