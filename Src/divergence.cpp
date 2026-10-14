#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>

using namespace amrex;

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "Computes the divergence of one or more vector fields of a plotfile.\n"
    << "Derivatives are evaluated as in the grad tool (face-centered "
       "gradients\n"
    << "from the MLMG Poisson operator, averaged to cell centers).\n\n"
    << "Usage:\n"
    << "  " << argv[0]
    << " infile=FILE vectorVars=NAMES outVars=NAMES [OPTIONS]\n\n"

    << "Required arguments:\n"
    << "  infile=FILE              AMReX plotfile\n"
    << "  vectorVars=NAMES         components of the vectors, "
    << AMREX_SPACEDIM << " per vector\n"
    << "                           (e.g. Fx Fy Fz Gx Gy Gz)\n"
    << "  outVars=NAMES            one output name per vector\n\n"

    << "Options:\n"
    << "  outfile=FILE             output plotfile [DEF: <infile>_div]\n"
    << "  keepVars=NAMES           variables copied unchanged to the output\n"
    << "  finestLevel=N            finest AMR level to process [DEF: "
       "finest]\n"
    << "  is_per=I I I             periodicity per direction [DEF: 0 0 0]\n"
    << "  sym_dir=I I I            reflect_odd boundaries per direction "
       "[DEF: 0 0 0]\n"
    << "  n_files=N                max number of plotfile data files\n"
    << "  -h, --help               show this help message\n\n"

    << "Visit PeleAnalysis/Src/InputsSamples for examples or refer to "
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

    ParmParse pp;

    std::string infile;
    pp.get("infile", infile);
    std::string outfile(getFileRoot(infile) + "_div");
    pp.query("outfile", outfile);

    // Initialize DataService
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if (!dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = 1000;
    pp.query("finestLevel", finestLevel);
    finestLevel = std::min(finestLevel, amrData.FinestLevel());
    const int Nlev = finestLevel + 1;
    RealBox rb(&(amrData.ProbLo()[0]), &(amrData.ProbHi()[0]));

    // ---------------------------------------------------------------------
    // Variables index management
    // ---------------------------------------------------------------------
    const int nVecComp = pp.countval("vectorVars");
    if (nVecComp < AMREX_SPACEDIM || nVecComp % AMREX_SPACEDIM != 0) {
      Abort(
        "vectorVars must hold " + std::to_string(AMREX_SPACEDIM) +
        " components per vector");
    }
    const int nVec = nVecComp / AMREX_SPACEDIM;
    Vector<std::string> vecVars(nVecComp);
    pp.getarr("vectorVars", vecVars, 0, nVecComp);

    if (pp.countval("outVars") != nVec) {
      Abort("outVars must hold one name per vector in vectorVars");
    }
    Vector<std::string> outVars(nVec);
    pp.getarr("outVars", outVars, 0, nVec);

    const int nKeep = pp.countval("keepVars");
    Vector<std::string> keepVars(nKeep);
    if (nKeep > 0) {
      pp.getarr("keepVars", keepVars, 0, nKeep);
    }

    // State: vector components first, then the carried variables
    Vector<std::string> inVarNames = vecVars;
    for (const auto& var : keepVars) {
      inVarNames.push_back(var);
    }
    for (const auto& var : inVarNames) {
      if (amrData.StateNumber(var) < 0) {
        Abort("Unknown variable name: " + var);
      }
    }
    const int nCompIn = inVarNames.size();
    Vector<int> destFillComps(nCompIn);
    for (int i = 0; i < nCompIn; ++i) {
      destFillComps[i] = i;
    }

    // Output: carried variables, then one divergence per vector
    const int idDiv = nKeep;
    const int nCompOut = nKeep + nVec;

    // Check symmetry/periodicity in given coordinate direction
    Vector<int> sym_dir(AMREX_SPACEDIM, 0);
    pp.queryarr("sym_dir", sym_dir, 0, AMREX_SPACEDIM);

    Vector<int> is_per(AMREX_SPACEDIM, 0);
    pp.queryarr("is_per", is_per, 0, AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      Print() << is_per[idim] << " ";
    }
    Print() << "\n";

    int coord = 0;

    // ---------------------------------------------------------------------
    // Let's start the real work
    // ---------------------------------------------------------------------
    Vector<MultiFab> state(Nlev);
    Vector<MultiFab> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);
    Vector<BoxArray> grids(Nlev);
    Vector<DistributionMapping> dmap(Nlev);
    const int nGrow = 1;

    // Read data on all the levels
    for (int lev = 0; lev < Nlev; ++lev) {
      const BoxArray ba = amrData.boxArray(lev);
      grids[lev] = ba;
      dmap[lev] = DistributionMapping(ba);
      geoms[lev] =
        Geometry(amrData.ProbDomain()[lev], &rb, coord, &(is_per[0]));
      state[lev].define(grids[lev], dmap[lev], nCompIn, nGrow);
      outdata[lev].define(grids[lev], dmap[lev], nCompOut, 0);

      Print() << "Reading data for level: " << lev << std::endl;
      amrData.FillVar(state[lev], lev, inVarNames, destFillComps);
      state[lev].FillBoundary(0, nVecComp, geoms[lev].periodicity());

      for (int n = 0; n < nKeep; ++n) {
        MultiFab::Copy(outdata[lev], state[lev], nVecComp + n, n, 1, 0);
      }
      outdata[lev].setVal(0.0, idDiv, nVec);
    }

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

    // One component at a time: d(F_d)/dx_d is the d-th face gradient of F_d,
    // averaged to the cell center, and is added to the divergence
    for (int v = 0; v < nVec; ++v) {
      Print() << "Computing " << outVars[v] << std::endl;
      for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const int comp = v * AMREX_SPACEDIM + d;

        // Get face-centered gradients from MLMG
        LPInfo info;
        info.setAgglomeration(1);
        info.setConsolidation(1);
        info.setMetricTerm(false);
        info.setMaxCoarseningLevel(0);
        MLPoisson poisson({geoms}, {grids}, {dmap}, info);
        poisson.setMaxOrder(4);
        poisson.setDomainBC(lo_bc, hi_bc);

        // Need to apply the operator to ensure CF consistency with composite
        // solve
        Vector<Array<MultiFab, AMREX_SPACEDIM>> grad(Nlev);
        Vector<std::unique_ptr<MultiFab>> phi;
        Vector<MultiFab> laps;
        for (int lev = 0; lev < Nlev; ++lev) {
          for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
            grad[lev][idim].define(
              amrex::convert(grids[lev], IntVect::TheDimensionVector(idim)),
              dmap[lev], 1, 0);
          }
          phi.push_back(
            std::make_unique<MultiFab>(state[lev], amrex::make_alias, comp, 1));
          poisson.setLevelBC(lev, phi[lev].get());
          laps.emplace_back(grids[lev], dmap[lev], 1, 1);
        }

        MLMG mlmg(poisson);
        mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
        mlmg.getFluxes(
          GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi),
          MLMG::Location::FaceCenter);

        // The MLMG fluxes are -grad(phi)
        const Dim3 shift = IntVect::TheDimensionVector(d).dim3();
        for (int lev = 0; lev < Nlev; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
          for (MFIter mfi(outdata[lev], TilingIfNotGPU()); mfi.isValid();
               ++mfi) {
            const Box& bx = mfi.tilebox();
            auto const& flux = grad[lev][d].const_array(mfi);
            auto const& div = outdata[lev].array(mfi, idDiv + v);
            amrex::ParallelFor(
              bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                div(i, j, k) -=
                  0.5 *
                  (flux(i, j, k) + flux(i + shift.x, j + shift.y, k + shift.z));
              });
          }
        }
      }
    }

    // ---------------------------------------------------------------------
    // Write the results
    // ---------------------------------------------------------------------
    Vector<std::string> nnames = keepVars;
    for (const auto& var : outVars) {
      nnames.push_back(var);
    }

    // Cap the number of plotfile data files via the n_files option (AMReX)
    int n_files = amrex::VisMF::GetNOutFiles();
    pp.query("n_files", n_files);
    amrex::VisMF::SetNOutFiles(n_files);

    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev - 1);
    for (int lev = 0; lev < Nlev - 1; ++lev) {
      int r = amrData.RefRatio()[lev];
      refRatios[lev] = IntVect(AMREX_D_DECL(r, r, r));
    }
    amrex::WriteMultiLevelPlotfile(
      outfile, Nlev, GetVecOfConstPtrs(outdata), nnames, geoms, amrData.Time(),
      isteps, refRatios);
  }
  amrex::Finalize();
  return 0;
}
