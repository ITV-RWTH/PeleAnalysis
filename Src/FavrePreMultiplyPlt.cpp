#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

static void
print_usage(int, char* argv[])
{
  std::cerr << "usage:\n"
            << "  " << argv[0]
            << " infile=IN_PLTFILE outfile=OUT_PLTFILE [OPTIONS]\n\n"

            << "Required arguments:\n"
            << "  infile=IN_PLTFILE         AMReX plotfile\n"
            << "  outfile=OUT_PLTFILE       Output plotfile with fields "
               "multiplied by density\n\n"

            << "Options:\n"
            << "  -h, --help         Show this help message\n\n"

            << "Visit PeleAnalysis/Src/InputSamples for examples or refer to "
            << "the documentation.\n";

  exit(1);
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
  Initialize(argc, argv);
  {
    if (argc < 2)
      print_usage(argc, argv);

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc, argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    // get infile and outfile name
    std::string plotFileName;
    pp.get("infile", plotFileName);
    std::string outFileName;
    pp.get("outfile", outFileName);

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    // setting up reading for pltfile
    DataServices dataServices(plotFileName, fileType);

    // check if plotfile is good
    if (!dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }

    AmrData& amrData = dataServices.AmrDataRef();

    // getting the finest level (or whatever user sets)
    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel", finestLevel);
    int Nlev = finestLevel + 1;

    // setting up periodicity
    Vector<int> is_per(AMREX_SPACEDIM, 1);
    pp.queryarr("is_per", is_per, 0, AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      Print() << is_per[idim] << " ";
    }

    Print() << "\n";

    RealBox rb(&(amrData.ProbLo()[0]), &(amrData.ProbHi()[0]));

    // setting up pltfile boxes
    Vector<MultiFab> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);

    int nGrow = 0;

    Vector<std::string> inNames;

    int nvar = pp.countval("variables");
    if (nvar == 0) {
      Abort("You must specify variables= including density");
    }

    pp.getarr("variables", inNames);

    // Validate variables exist
    Vector<std::string> plotVars = amrData.PlotVarNames();
    for (const auto& v : inNames) {
      if (std::find(plotVars.begin(), plotVars.end(), v) == plotVars.end()) {
        Abort("Variable '" + v + "' not found in plotfile");
      }
    }

    Vector<std::string> outNames = inNames;
    int nComp = inNames.size();

    Vector<int> destFillComps(nComp);
    for (int i = 0; i < nComp; ++i) {
      destFillComps[i] = i;
    }

    int ID_rho = -1;

    for (int i = 0; i < nComp; ++i) {
      destFillComps[i] = i;
      if (inNames[i] == "density")
        ID_rho = i;
    }

    // Check if density is written in the infile
    if (ID_rho < 0) {
      Abort("density not found in plot file!");
    }

    for (int lev = 0; lev < Nlev; ++lev) {

      BoxArray ba = amrData.boxArray(lev);
      ba.maxSize(32); // or 64 depending on memory
      const DistributionMapping dm(ba);

      outdata[lev] = MultiFab(ba, dm, nComp, nGrow);
      MultiFab indata(ba, dm, nComp, nGrow);
      // outdata[lev].setVal(0.0);

      int coord = 0;
      geoms[lev] =
        Geometry(amrData.ProbDomain()[lev], &rb, coord, &(is_per[0]));

      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata, lev, inNames, destFillComps);
      Print() << "Data has been read for level " << lev << std::endl;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(indata, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto const& out_a = outdata[lev].array(mfi);
        auto const& in_a = indata.array(mfi);
        amrex::ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            // Do something funky
            for (int n = 0; n < nComp; n++) {
              if (n == ID_rho) {
                out_a(i, j, k, n) = in_a(i, j, k, n);
              } else {
                out_a(i, j, k, n) = in_a(i, j, k, ID_rho) * in_a(i, j, k, n);
              }
            }
          });
      }
      amrex::Gpu::synchronize();
    }

    // write pltfile
    Print() << "Writing new data to " << outFileName << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev - 1, {AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(
      outFileName, Nlev, GetVecOfConstPtrs(outdata), outNames, geoms, 0.0,
      isteps, refRatios);
  }
  Finalize();
  return 0;
}
