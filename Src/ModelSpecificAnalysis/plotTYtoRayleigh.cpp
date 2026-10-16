#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_Array.H>

#include <PelePhysics.H>
#include "RayleighCrossSections.H"

using namespace amrex;

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "Usage:\n"
    << "  " << argv[0] << " infile=FILE lambda=val [OPTIONS]\n\n"

    << "Computes the molar-weighted mixture Rayleigh cross-section\n"
    << "  sigma_mix = sum_k X_k sigma_k(nu)\n"
    << "from the Y(*) fields of a plotfile.\n\n"

    << "Required arguments:\n"
    << "  infile=FILE        AMReX plotfile to read\n"
    << "  lambda=val      Laser wavelength [nm] \n\n"

    << "Options:\n"
    << "  outfile=FILE       Output plotfile (default: <infile>_Rayleigh)\n"
    << "  Aux_Variables=\"...\" Variables copied unchanged to the output\n"
    << "  n_files=N          Max number of plotfile data files\n"
    << "  finestLevel=N      Cap AMR levels (default: all levels in file)\n"
    << "  is_per=\"1 1 1\"     Periodicity flags per dimension (default: all "
       "1)\n"
    << "  verbose            Enable verbose plotfile I/O\n"
    << "  -h, --help         Show this help message\n\n"

    << "Visit PeleAnalysis/Src/InputSamples for examples or refer to "
    << "the documentation.\n";

  std::exit(1);
}

static std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile, std::string("/"));
  return tokens[tokens.size() - 1];
}

int
main(int argc, char* argv[])
{
  Initialize(argc, argv);
  {
    if (
      argc < 2 || std::strcmp(argv[1], "-h") == 0 ||
      std::strcmp(argv[1], "--help") == 0) {
      print_usage(argc, argv);
    }
    ParmParse pp;

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    std::string infile;
    pp.get("infile", infile);

    // Laser wavelength [nm]
    Real lambda_nm = 0.0;
    pp.get("lambda", lambda_nm);
    if (lambda_nm <= 0.0) {
      amrex::Abort("plotYTtoRayleigh: lambda must be positive [nm]");
    }

    Vector<std::string> auxVars;
    const int nAux = pp.countval("Aux_Variables");
    if (nAux > 0) {
      pp.getarr("Aux_Variables", auxVars, 0, nAux);
    }

    Vector<int> is_per(AMREX_SPACEDIM, 1);
    pp.queryarr("is_per", is_per, 0, AMREX_SPACEDIM);

    // -------------------------------------------------------------------------
    // Open plotfile
    // -------------------------------------------------------------------------
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if (!dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel", finestLevel);
    finestLevel = std::min(finestLevel, amrData.FinestLevel());
    const int Nlev = finestLevel + 1;

    // -------------------------------------------------------------------------
    // Per-species cross-sections at this wavelength (host, once)
    // -------------------------------------------------------------------------
    Vector<std::string> spec_names;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
      spec_names);

    GpuArray<Real, NUM_SPECIES> sigma;
    rayleigh::cross_sections(lambda_nm, sigma.data());
    Print() << "Rayleigh cross-sections [m^2] at lambda = " << lambda_nm << " nm):\n";
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Print() << "  " << spec_names[n] << ": " << sigma[n]
              << (sigma[n] == 0.0 ? "  (no data)" : "") << "\n";
    }

    // -------------------------------------------------------------------------
    // Input layout: Y(0..NUM_SPECIES-1), temp, aux
    // -------------------------------------------------------------------------
    constexpr int idYin = 0;
    constexpr int idTin = NUM_SPECIES;
    constexpr int idAuxIn = NUM_SPECIES + 1;
    const int nCompIn = NUM_SPECIES + 1 + nAux;

    Vector<std::string> inNames;
    for (int n = 0; n < NUM_SPECIES; ++n) {
      inNames.push_back("Y(" + spec_names[n] + ")");
    }
    inNames.push_back("temp");
    for (const auto& v : auxVars) {
      inNames.push_back(v);
    }
    Vector<int> destFillComps(nCompIn);
    for (int n = 0; n < nCompIn; ++n) {
      if (amrData.StateNumber(inNames[n]) < 0) {
        amrex::Abort("plotTYtoRayleigh: cannot find " + inNames[n] + " in " + infile);
      }
      destFillComps[n] = n;
    }

    // -------------------------------------------------------------------------
    // Output layout: RayleighSigma, Temp, aux
    // -------------------------------------------------------------------------
    constexpr int idSigOut = 0;
    constexpr int idIOut = 1;
    constexpr int idLogSigOut = 2;
    constexpr int idLogIOut = 3;
    constexpr int idTOut = 4;
    constexpr int idLogTOut = 5;
    constexpr int idAuxOut = 6;
    const int nCompOut = idAuxOut + nAux;

    Vector<std::string> outNames = {"RayleighSigma", "RayleighIntensity","log(RayleighSigma)", "log(RayleighIntensity)","temp","log(temp)"};
    for (const auto& v : auxVars) {
      outNames.push_back(v);
    }

    // -------------------------------------------------------------------------
    // Geometry
    // -------------------------------------------------------------------------
    Vector<Geometry> geoms(Nlev);
    RealBox real_box(
      {AMREX_D_DECL(
        amrData.ProbLo()[0], amrData.ProbLo()[1], amrData.ProbLo()[2])},
      {AMREX_D_DECL(
        amrData.ProbHi()[0], amrData.ProbHi()[1], amrData.ProbHi()[2])});
    Array<int, AMREX_SPACEDIM> periodic{
      AMREX_D_DECL(is_per[0], is_per[1], is_per[2])};
    geoms[0] = Geometry(
      amrData.ProbDomain()[0], real_box, amrData.CoordSys(), periodic);
    Vector<IntVect> refRatios(Nlev - 1);
    for (int lev = 1; lev < Nlev; ++lev) {
      refRatios[lev - 1] = IntVect(amrData.RefRatio()[lev - 1]);
      geoms[lev] = amrex::refine(geoms[lev - 1], refRatios[lev - 1]);
    }

    // -------------------------------------------------------------------------
    // Process
    // -------------------------------------------------------------------------
    auto eos = pele::physics::PhysicsType::eos();
    Vector<MultiFab> outmf(Nlev);

    for (int lev = 0; lev < Nlev; ++lev) {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      MultiFab indata(ba, dm, nCompIn, 0);
      outmf[lev].define(ba, dm, nCompOut, 0);

      Print() << "Reading data for level " << lev << "\n";
      amrData.FillVar(indata, lev, inNames, destFillComps);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(indata, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto const& in = indata.const_array(mfi);
        auto const& out = outmf[lev].array(mfi);
        amrex::ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Real Yl[NUM_SPECIES];
            Real Xl[NUM_SPECIES];
            for (int n = 0; n < NUM_SPECIES; ++n) {
              Yl[n] = in(i, j, k, idYin + n);
            }
            eos.Y2X(Yl, Xl);
            Real sig_mix = 0.0;
            for (int n = 0; n < NUM_SPECIES; ++n) {
              sig_mix += Xl[n] * sigma[n];
            }
	    const Real T = in(i,j,k,idTin);
            out(i, j, k, idSigOut) = sig_mix;
	    out(i,j,k,idLogSigOut) = std::log(sig_mix);
	    out(i,j,k,idTOut) = T;
	    out(i,j,k,idLogTOut) = std::log(T);
	    const Real I = sig_mix/T;
            out(i, j, k, idIOut) = I;
	    out(i,j,k,idLogIOut) = std::log(I);
          });
      }

      // Copy auxiliary variables unchanged
      if (nAux > 0) {
        MultiFab::Copy(outmf[lev], indata, idAuxIn, idAuxOut, nAux, 0);
      }
      Print() << "Derive finished for level " << lev << "\n";
    }

    // -------------------------------------------------------------------------
    // Write output
    // -------------------------------------------------------------------------
    std::string outfile = getFileRoot(infile) + "_Rayleigh";
    pp.query("outfile", outfile);

    // Cap the number of plotfile data files via the n_files option (AMReX)
    int n_files = VisMF::GetNOutFiles();
    pp.query("n_files", n_files);
    VisMF::SetNOutFiles(n_files);

    Print() << "Writing new data to " << outfile << "\n";
    Vector<int> isteps(Nlev, 0);
    WriteMultiLevelPlotfile(
      outfile, Nlev, GetVecOfConstPtrs(outmf), outNames, geoms, amrData.Time(),
      isteps, refRatios);
  }
  Finalize();
  return 0;
}
