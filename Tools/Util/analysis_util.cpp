#include <analysis_util.H>

namespace analysis_util {

int nlev = -1;
amrex::Vector<amrex::BoxArray> grids;
amrex::Vector<amrex::DistributionMapping> dmap;
amrex::Vector<int> ref_ratio;
amrex::Vector<amrex::Geometry> geoms;
amrex::IntVect periodicity;
bool initialized = false;

void
init(
  const amrex::Vector<amrex::MultiFab>& a_mf,
  const amrex::Vector<amrex::Geometry>& a_geoms)
{
  // initialise
  if (nlev == -1) {
    nlev = a_mf.size();
  }
  if (grids.empty()) {
    grids.reserve(nlev);
    for (int lev = 0; lev < nlev; ++lev) {
      grids.emplace_back(a_mf[lev].boxArray());
    }
  }
  if (dmap.empty()) {
    dmap.reserve(nlev);
    for (int lev = 0; lev < nlev; ++lev) {
      dmap.emplace_back(a_mf[lev].DistributionMap());
    }
  }
  if (ref_ratio.empty()) {
    ref_ratio.assign(nlev, 2);
  }
  if (geoms.empty()) {
    geoms = a_geoms;
  }
  periodicity = a_geoms[0].periodicity().intVect();

  initialized = true;
}

amrex::Vector<std::unique_ptr<amrex::iMultiFab>>
get_covered_mf()
{
  AMREX_ALWAYS_ASSERT(initialized);

  amrex::Vector<std::unique_ptr<amrex::iMultiFab>> mask_mf(nlev);
  for (int lev = 0; lev < nlev; ++lev) {
    mask_mf[lev] =
      std::make_unique<amrex::iMultiFab>(grids[lev], dmap[lev], 1, 0);
    mask_mf[lev]->setVal(1);
  }
  const int finest_level = nlev - 1;
  for (int lev = 0; lev < finest_level; ++lev) {
    amrex::BoxArray baf = grids[lev + 1];
    baf.coarsen(ref_ratio[lev]);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {
      std::vector<std::pair<int, amrex::Box>> isects;
      for (amrex::MFIter mfi(*mask_mf[lev], amrex::TilingIfNotGPU());
           mfi.isValid(); ++mfi) {
        auto const& mask = mask_mf[lev]->array(mfi);
        baf.intersections(grids[lev][mfi.index()], isects);
        for (const auto& is : isects) {
          amrex::ParallelFor(
            is.second, [mask] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              mask(i, j, k) = 0;
            });
        }
      }
    }
  }
  return mask_mf;
}

} // namespace analysis_util
