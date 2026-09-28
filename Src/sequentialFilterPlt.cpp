#include <algorithm>
#include <cstring>
#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_FillPatchUtil.H>
#include <PltFileManager.H>
#include <PltFileManagerBCFill.H>

namespace {

using BndryFunc = amrex::PhysBCFunct<
  amrex::GpuBndryFuncFab<pele::physics::pltfilemanager::FillExtDirDummy>>;

void
print_usage(int, char* argv[])
{
  std::cerr
    << "Applies a box filter to a plotfile as three successive 1D passes (x,\n"
    << "then y, then z), so that only a one-directional grow region is ever\n"
    << "allocated. Use this instead of filterPlt when the 3D grow region of a\n"
    << "single-pass filter does not fit in memory. Only the box filter is\n"
    << "implemented; must be compiled with DIM = 3.\n\n"

    << "Usage:\n"
    << "  " << argv[0] << " infile=FILE [OPTIONS]\n\n"

    << "Required arguments:\n"
    << "  infile=FILE              AMReX plotfile\n\n"

    << "Options:\n"
    << "  outfile=FILE             Output plotfile (default: <infile>_filtered)\n"
    << "  variables=LIST           Space-separated names of the variables to\n"
    << "                           filter (default: every variable in the file)\n"
    << "  base_fgr=N               Even filter-to-grid ratio on level 0 (default: 2)\n"
    << "  same_fgr_all_levels=BOOL Keep the ratio instead of the physical width\n"
    << "                           on finer levels (default: false)\n"
    << "  max_filter_level=N       Finest level to filter (default: all)\n"
    << "  max_grid_size=N          Re-box the input grids (default: unchanged)\n"
    << "  interp_type=N            0: piecewise constant, 1: conservative linear\n"
    << "                           (default: 1)\n"
    << "  n_files=N                Max number of plotfile data files\n"
    << "  -h, --help               Show this help message\n\n"

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

amrex::IntVect
ghostsAlong(int dir, int ngrow)
{
  amrex::IntVect ng(0);
  if (dir < AMREX_SPACEDIM) {
    ng[dir] = ngrow;
  }
  return ng;
}

void
fillGhostCells(
  amrex::Vector<amrex::MultiFab>& mf,
  const amrex::Vector<amrex::Geometry>& geom,
  const amrex::Vector<amrex::IntVect>& ref_ratio,
  const amrex::Vector<amrex::BCRec>& bcr,
  amrex::InterpBase* mapper)
{
  const int nlev = static_cast<int>(mf.size());
  amrex::Vector<amrex::Vector<amrex::MultiFab*>> smf(nlev);
  amrex::Vector<amrex::Vector<amrex::Real>> stime(nlev, {0.0});
  amrex::Vector<BndryFunc> bndry_func;
  for (int lev = 0; lev < nlev; ++lev) {
    smf[lev] = {&mf[lev]};
    bndry_func.emplace_back(
      geom[lev], bcr, pele::physics::pltfilemanager::FillExtDirDummy{});
  }
  // Unlike FillPatchTwoLevels, this also uses levels below lev-1 where needed
  for (int lev = 0; lev < nlev; ++lev) {
    amrex::FillPatchNLevels(
      mf[lev], lev, mf[lev].nGrowVect(), 0.0, smf, stime, 0, 0,
      mf[lev].nComp(), geom, bndry_func, 0, ref_ratio, mapper, bcr, 0);
  }
}

void
filterAlong(
  int dir,
  const amrex::Vector<amrex::MultiFab>& src,
  amrex::Vector<amrex::MultiFab>& dst,
  const amrex::Vector<amrex::Vector<amrex::Real>>& weights)
{
  for (int lev = 0; lev < static_cast<int>(src.size()); ++lev) {
    const int nweights = static_cast<int>(weights[lev].size());
    const int ng = nweights / 2;
    amrex::Gpu::DeviceVector<amrex::Real> d_weights(nweights);
    amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, weights[lev].begin(), weights[lev].end(),
      d_weights.begin());
    const amrex::Real* w = d_weights.data();
    const amrex::IntVect e = amrex::IntVect::TheDimensionVector(dir);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(dst[lev], amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& bx = mfi.tilebox();
      amrex::Array4<const amrex::Real> in = src[lev].const_array(mfi);
      amrex::Array4<amrex::Real> out = dst[lev].array(mfi);
      amrex::ParallelFor(
        bx, dst[lev].nComp(),
        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          amrex::Real sum = 0.0;
          for (int m = -ng; m <= ng; ++m) {
            sum += w[ng + m] * in(i + m * e[0], j + m * e[1], k + m * e[2], n);
          }
          out(i, j, k, n) = sum;
        });
    }
    // d_weights is freed at the end of this iteration; the kernels must finish.
    amrex::Gpu::streamSynchronize();
  }
}

} // namespace

int
main(int argc, char** argv)
{
  amrex::Initialize(argc, argv);
  static_assert(AMREX_SPACEDIM == 3, "This tool must be compiled with DIM = 3");

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
  int les_filter_type = 1;
  int les_filter_fgr = 2;
  bool same_fgr_all_levels = false;
  int interp_type = 1;
  int max_grid_size = -1;
  pp.get("infile", infile);
  pp.query("max_filter_level", finestLevel);
  pp.query("filter_type", les_filter_type);
  // Fail loudly rather than silently box-filtering for other PeleC types
  if (les_filter_type != 1) {
    amrex::Abort("sequentialFilterPlt only implements the box filter "
                 "(filter_type = 1); use filterPlt for the other types");
  }
  pp.query("base_fgr", les_filter_fgr); // filter to grid ratio
  AMREX_ALWAYS_ASSERT(les_filter_fgr > 0 && les_filter_fgr % 2 == 0);
  pp.query("same_fgr_all_levels", same_fgr_all_levels);
  pp.query("interp_type", interp_type);
  pp.query("max_grid_size", max_grid_size);

  std::string outfile(getFileRoot(infile) + "_filtered");
  pp.query("outfile", outfile);

  // Cap the number of plotfile data files via the n_files option (AMReX)
  int n_files = amrex::VisMF::GetNOutFiles();
  pp.query("n_files", n_files);
  amrex::VisMF::SetNOutFiles(n_files);

  // use PltFileManager to load data
  auto plt_file_data =
    std::make_unique<pele::physics::pltfilemanager::PltFileManager>(infile);
  // Plotfile global infos
  int Nlev = std::min(finestLevel + 1, plt_file_data->getNlev());

  // Variable names
  amrex::Vector<std::string> variableNames;
  amrex::Vector<int> var_idxs;
  const bool all_vars = pp.countval("variables") == 0;
  if (!all_vars) {
    pp.getarr("variables", variableNames);
    const amrex::Vector<std::string>& plotVarNames =
      plt_file_data->getVariableList();
    for (const auto& name : variableNames) {
      auto it = std::find(plotVarNames.begin(), plotVarNames.end(), name);
      if (it == plotVarNames.end()) {
        amrex::Abort("Variable '" + name + "' not found in file");
      }
      var_idxs.push_back(static_cast<int>(it - plotVarNames.begin()));
    }
  } else {
    variableNames = plt_file_data->getVariableList();
  }
  const int ncomp_filter = static_cast<int>(variableNames.size());

  amrex::Vector<amrex::MultiFab> src(Nlev);
  amrex::Vector<amrex::MultiFab> dst(Nlev);
  amrex::Vector<amrex::Geometry> level_geometries;
  amrex::Vector<amrex::IntVect> ref_ratio;
  amrex::Vector<int> nGrow(Nlev);
  amrex::Vector<amrex::Vector<amrex::Real>> weights(Nlev);

  // Load the data from the Pltfile
  amrex::Print() << "Reading data..." << std::endl;
  int les_filter_fgr_lev = les_filter_fgr;
  for (int lev = 0; lev < Nlev; ++lev) {
    amrex::Print() << "on level " << lev << std::endl;

    if (lev > 0) {
      ref_ratio.emplace_back(plt_file_data->getRefRatio(lev - 1));
      if (!same_fgr_all_levels) {
        les_filter_fgr_lev *= plt_file_data->getRefRatio(lev - 1);
      }
    }
    // Trapezoidal box filter: fgr+1 points, half weight on the two end points
    nGrow[lev] = les_filter_fgr_lev / 2;
    weights[lev].assign(les_filter_fgr_lev + 1, 1.0 / les_filter_fgr_lev);
    weights[lev].front() = weights[lev].back() = 0.5 / les_filter_fgr_lev;

    amrex::BoxArray ba(plt_file_data->getGrid(lev));
    if (max_grid_size > 0) {
      ba.maxSize(max_grid_size);
    }
    const amrex::DistributionMapping dm(ba);

    amrex::Print() << "Number of boxes on level " << lev << " " << ba.size()
                   << std::endl;
    level_geometries.push_back(plt_file_data->getGeom(lev));
    src[lev].define(ba, dm, ncomp_filter, ghostsAlong(0, nGrow[lev]));
    dst[lev].define(ba, dm, ncomp_filter, ghostsAlong(1, nGrow[lev]));

    if (all_vars) {
      plt_file_data->fillPatchFromPlt(
        lev, level_geometries[lev], 0, 0, ncomp_filter, src[lev], interp_type);
    } else {
      for (int var = 0; var < ncomp_filter; ++var) {
        plt_file_data->fillPatchFromPlt(
          lev, level_geometries[lev], var_idxs[var], var, 1, src[lev],
          interp_type);
      }
    }
  }
  amrex::Print() << "Done!" << std::endl;

  // Non-periodic domain boundaries are FOExtraped for lack of anything better
  amrex::Vector<amrex::BCRec> dummyBCRec(ncomp_filter);
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    const int bctype = level_geometries[0].isPeriodic(idim)
                         ? amrex::BCType::int_dir
                         : amrex::BCType::foextrap;
    for (int n = 0; n < ncomp_filter; n++) {
      dummyBCRec[n].setLo(idim, bctype);
      dummyBCRec[n].setHi(idim, bctype);
    }
  }
  amrex::InterpBase* mapper = &amrex::mf_pc_interp;
  if (interp_type == 1) {
    mapper = &amrex::mf_cell_cons_interp;
  }

  // src is grown along dir, dst along dir+1; after the swap dst is regrown
  const char* dir_names = "xyz";
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    amrex::Print() << "Filtering in " << dir_names[dir] << "..." << std::endl;
    fillGhostCells(src, level_geometries, ref_ratio, dummyBCRec, mapper);
    filterAlong(dir, src, dst, weights);
    std::swap(src, dst);
    if (dir + 1 < AMREX_SPACEDIM) {
      for (int lev = 0; lev < Nlev; ++lev) {
        const amrex::BoxArray ba = dst[lev].boxArray();
        const amrex::DistributionMapping dm = dst[lev].DistributionMap();
        dst[lev].define(
          ba, dm, ncomp_filter, ghostsAlong(dir + 2, nGrow[lev]));
      }
    }
  }
  amrex::Print() << "Done!" << std::endl;

  amrex::Print() << "Saving filtered data..." << std::endl;
  amrex::WriteMultiLevelPlotfile(
    outfile, Nlev, amrex::GetVecOfConstPtrs(src), variableNames,
    level_geometries, plt_file_data->getTime(), amrex::Vector<int>(Nlev, 0),
    ref_ratio);
  amrex::Print() << "Done!" << std::endl;

  amrex::Finalize();
  return 0;
}
