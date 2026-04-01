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
  std::cerr << "Usage:\n"
            << "  " << argv[0] << " infile=FILE [OPTIONS]\n\n"

            << "Required arguments:\n"
            << "  infile=FILE        AMReX plotfile\n\n"
            << "  outfile=FILE       AMReX plotfile [DEF=\%infile_multiply]\n\n"
	    << "  inVarAName=NAME\n\n"
	    << "  inVarBName=NAME\n\n"
	    << "  outVarName=NAME\n\n"

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
  Initialize(argc, argv);
  {
    if (argc < 2) {
      print_usage(argc, argv);
    } else if (
      (std::strcmp(argv[1], "-h") == 0) ||
      (std::strcmp(argv[1], "--help") == 0)) {
      print_usage(argc, argv);
    }
    ParmParse pp;

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    std::string infileName;
    pp.get("infile", infileName);
    std::string inVarAName;
    std::string outfileName(getFileRoot(infileName) + "_multiply");
    pp.query("outfile", outfileName);
    pp.get("inVarAName", inVarAName);
    std::string inVarBName;
    pp.get("inVarBName", inVarBName);
    std::string outVarName;
    pp.get("outVarName", outVarName);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(infileName, fileType);
    if (!dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel", finestLevel);
    int Nlev = finestLevel + 1;

    Vector<int> is_per(AMREX_SPACEDIM, 1);
    pp.queryarr("is_per", is_per, 0, AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      Print() << is_per[idim] << " ";
    }

    RealBox rb(&(amrData.ProbLo()[0]), &(amrData.ProbHi()[0]));

    Vector<MultiFab> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);

    int nGrow = 0;
    int nCompIn = amrData.NComp();
    int nCompOut = nCompIn + 1;
    int outVar_id = nCompOut - 1;

    Vector<std::string> inNames = amrData.PlotVarNames();
    auto id = std::find(inNames.begin(), inNames.begin(), inVarAName);
    int inVarA_id; 
    if (id != inNames.end()) {
	inVarA_id = std::distance(inNames.begin(), id);
    } else {
	Abort("Variable " + inVarAName + " not found in file " + infileName + "!");
    } 
    id = std::find(inNames.begin(), inNames.begin(), inVarBName);
    int inVarB_id;
    if (id != inNames.end()) {
	inVarB_id = std::distance(inNames.begin(), id);
    } else {
	Abort("Variable " + inVarBName + " not found in file " + infileName + "!");
    } 
    Vector<std::string> outNames = amrData.PlotVarNames();
    outNames.push_back(outVarName);

    Vector<int> destFillComps(nCompOut);
    for (int i = 0; i < nCompOut; ++i)
      destFillComps[i] = i;

    for (int lev = 0; lev < Nlev; ++lev) {

      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);

      outdata[lev] = MultiFab(ba, dm, nCompOut, nGrow);
      MultiFab indata(ba, dm, nCompIn, nGrow);

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
        const Box& bx = mfi.tilebox();
        auto const& out_a = outdata[lev].array(mfi);
        auto const& in_a = indata.array(mfi);
	amrex::ParallelFor(
          bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            out_a(i, j, k, outVar_id) = in_a(i, j, k, inVarA_id) * in_a(i, j, k, inVarB_id);
          });
      }
    }

    Print() << "Writing new data to " << outfileName << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev - 1, {AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(
      outfileName, Nlev, GetVecOfConstPtrs(outdata), outNames, geoms, 0.0, isteps,
      refRatios);
  }
  Finalize();
  return 0;
}
