#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_FillPatchUtil.H>
#include <Filter.H>
#include <PltFileManager.H>
#include <PltFileManagerBCFill.H>

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "Usage:\n"
    << "  " << argv[0] << " infile=FILE [OPTIONS]\n\n"

    << "Required arguments:\n"
    << "  infile=FILE        AMReX plotfile\n\n"

    << "Options:\n"
    << "  variables=LIST     Space-separated names of the variables to filter\n"
    << "                     (default: every variable in the file)\n"
    << "  -h, --help         Show this help message\n\n"
    << "Visit PeleAnalysis/Src/InputsSamples for examples or refer to "
    << "the documentation.\n";

  std::exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = amrex::Tokenize(infile, std::string("/"));
  return tokens[tokens.size() - 1];
}

int
main(int argc, char** argv)
{

  amrex::Initialize(argc, argv);

  if (argc < 2) {
    print_usage(argc, argv);
  } else if (
    (std::strcmp(argv[1], "-h") == 0) ||
    (std::strcmp(argv[1], "--help") == 0)) {
    print_usage(argc, argv);
  }

  amrex::ParmParse pp;

  std::string infile = "";
  int finestLevel = 1000;
  amrex::Vector<Filter> les_filter;
  int les_filter_type = 1;
  int les_filter_fgr = 2;
  bool same_fgr_all_levels = false;
  int interp_type = 1;

  pp.get("infile", infile);
  pp.query("max_filter_level", finestLevel);
  pp.query("filter_type", les_filter_type);
  pp.query("base_fgr", les_filter_fgr); // filter to grid ratio
  pp.query("same_fgr_all_levels", same_fgr_all_levels);
  int max_grid_size = 32;
  pp.query("max_grid_size", max_grid_size);
  pp.query("interp_type", interp_type);

  // use PltFileManager to load data
  amrex::Vector<pele::physics::pltfilemanager::PltFileManager*> plt_file_data(
    1);
  plt_file_data[0] = new pele::physics::pltfilemanager::PltFileManager(infile);
  // Plotfile global infos
  int Nlev = std::min(finestLevel + 1, plt_file_data[0]->getNlev());

  // Variable names
  int ncomp_filter;
  amrex::Vector<std::string> variableNames;
  amrex::Vector<int> var_idxs;
  int nvar = pp.countval("variables");
  bool all_vars = nvar > 0 ? false : true;
  if (!all_vars) {
    ncomp_filter = nvar;
    pp.getarr("variables", variableNames);
    const amrex::Vector<std::string>& plotVarNames =
      plt_file_data[0]->getVariableList();
    for (int var = 0; var < nvar; ++var) {
      int pvar;
      for (pvar = 0; pvar < plotVarNames.size(); ++pvar) {
        if (variableNames[var] == plotVarNames[pvar]) {
          var_idxs.push_back(pvar);
          break;
        }
      }
      if (pvar == plotVarNames.size()) {
        amrex::Abort("Variable '" + variableNames[var] + "' not found in file");
      }
    }
  } else {
    variableNames = plt_file_data[0]->getVariableList();
    ncomp_filter = variableNames.size();
  }

  amrex::Vector<amrex::MultiFab> indata(Nlev);
  amrex::Vector<amrex::MultiFab> outdata(Nlev);
  amrex::Vector<amrex::Geometry> level_geometries;
  amrex::Vector<amrex::IntVect> ref_ratio;

  // Load the data from the Pltfile
  amrex::Print() << "Reading data..." << std::endl;
  int les_filter_fgr_lev = les_filter_fgr;
  for (int lev = 0; lev < Nlev; ++lev) {
    amrex::Print() << "on level " << lev << std::endl;

    // Initialize filter stuff
    if (lev > 0) {
      ref_ratio.emplace_back(plt_file_data[0]->getRefRatio(lev - 1));
    }
    if ((!same_fgr_all_levels) && (lev > 0)) {
      les_filter_fgr_lev *= plt_file_data[0]->getRefRatio(lev - 1);
    }

    les_filter.push_back(Filter(les_filter_type, les_filter_fgr_lev));
    int nGrowF = les_filter[lev].get_filter_ngrow();
    int nGrow = 0;

    amrex::BoxArray ba(plt_file_data[0]->getGrid(lev));
    ba.maxSize(max_grid_size);
    const amrex::DistributionMapping dm = amrex::DistributionMapping(ba);
    level_geometries.push_back(plt_file_data[0]->getGeom(lev));
    indata[lev].define(ba, dm, ncomp_filter, nGrowF);
    outdata[lev].define(ba, dm, ncomp_filter, nGrow);

    if (all_vars) {
      plt_file_data[0]->fillPatchFromPlt(
        lev, level_geometries[lev], 0, 0, ncomp_filter, indata[lev],
        interp_type);
    } else {
      for (int var = 0; var < nvar; ++var) {
        plt_file_data[0]->fillPatchFromPlt(
          lev, level_geometries[lev], var_idxs[var], var, 1, indata[lev],
          interp_type);
      }
    }
  }
  amrex::Print() << "Done!" << std::endl;

  // fillPatchFromPlt doesn't fill ghost cells, so fill those now
  // Note: domain boundary cells will be FOExtraped because we don't know
  // anything better to do
  amrex::Print() << "FillPatching data..." << std::endl;
  amrex::Vector<amrex::BCRec> dummyBCRec(ncomp_filter);
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    if (level_geometries[0].isPeriodic(idim)) {
      for (int n = 0; n < ncomp_filter; n++) {
        dummyBCRec[n].setLo(idim, amrex::BCType::int_dir);
        dummyBCRec[n].setHi(idim, amrex::BCType::int_dir);
      }
    } else {
      for (int n = 0; n < ncomp_filter; n++) {
        dummyBCRec[n].setLo(idim, amrex::BCType::foextrap);
        dummyBCRec[n].setHi(idim, amrex::BCType::foextrap);
      }
    }
  }
  amrex::InterpBase* mapper = &amrex::mf_pc_interp;
  if (interp_type == 1) {
    mapper = &amrex::mf_cell_cons_interp;
  }
  amrex::Vector<amrex::Vector<amrex::MultiFab*>> smf(Nlev);
  amrex::Vector<amrex::Vector<amrex::Real>> stime(Nlev, {0.0});
  amrex::Vector<amrex::PhysBCFunct<
    amrex::GpuBndryFuncFab<pele::physics::pltfilemanager::FillExtDirDummy>>>
    bndry_func;
  for (int lev = 0; lev < Nlev; ++lev) {
    smf[lev] = {&indata[lev]};
    bndry_func.emplace_back(
      level_geometries[lev], dummyBCRec,
      pele::physics::pltfilemanager::FillExtDirDummy{});
  }
  for (int lev = 0; lev < Nlev; ++lev) {
    amrex::Print() << "on level " << lev << std::endl;
    // Unlike FillPatchTwoLevels, this also uses levels below lev-1 where needed
    amrex::FillPatchNLevels(
      indata[lev], lev, indata[lev].nGrowVect(), 0.0, smf, stime, 0, 0,
      ncomp_filter, level_geometries, bndry_func, 0, ref_ratio, mapper,
      dummyBCRec, 0);
  }
  amrex::Print() << "Done!" << std::endl;

  amrex::Print() << "Filtering data..." << std::endl;
  for (int lev = 0; lev < Nlev; ++lev) {
    amrex::Print() << "on level " << lev << std::endl;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(indata[lev], amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

      amrex::FArrayBox& fab_in = (indata[lev])[mfi];
      amrex::FArrayBox& fab_out = (outdata[lev])[mfi];
      const amrex::Box& box = mfi.tilebox();

      les_filter[lev].apply_filter(box, fab_in, fab_out);
    }
  }
  amrex::Print() << "Done!" << std::endl;

  amrex::Print() << "Saving filtered data..." << std::endl;
  std::string outfile(getFileRoot(infile) + "_filtered");

  // Cap the number of plotfile data files via the n_files option (AMReX)
  int n_files = amrex::VisMF::GetNOutFiles();
  pp.query("n_files", n_files);
  amrex::VisMF::SetNOutFiles(n_files);

  amrex::WriteMultiLevelPlotfile(
    outfile, Nlev, amrex::GetVecOfConstPtrs(outdata), variableNames,
    level_geometries, plt_file_data[0]->getTime(), amrex::Vector<int>(Nlev, 0),
    ref_ratio);
  amrex::Print() << "Done!" << std::endl;

  amrex::Finalize();
  return 0;
}
