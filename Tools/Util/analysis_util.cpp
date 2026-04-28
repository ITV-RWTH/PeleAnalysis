#include <analysis_util.H>

namespace analysis_util {

amrex::Vector<std::unique_ptr<amrex::iMultiFab>>
get_covered_mf(
  const amrex::Vector<amrex::MultiFab>& mf,
  const amrex::Vector<int>&             ref_ratios)
{
  const int nlev = static_cast<int>(mf.size());
  amrex::Vector<std::unique_ptr<amrex::iMultiFab>> mask_mf(nlev);
  for (int lev = 0; lev < nlev; ++lev) {
    mask_mf[lev] = std::make_unique<amrex::iMultiFab>(
      mf[lev].boxArray(), mf[lev].DistributionMap(), 1, 0);
    mask_mf[lev]->setVal(1);
  }
  const int finest_level = nlev - 1;
  for (int lev = 0; lev < finest_level; ++lev) {
    amrex::BoxArray baf = mf[lev + 1].boxArray();
    baf.coarsen(ref_ratios[lev]);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    {
      std::vector<std::pair<int, amrex::Box>> isects;
      for (amrex::MFIter mfi(*mask_mf[lev], amrex::TilingIfNotGPU());
           mfi.isValid(); ++mfi) {
        auto const& mask = mask_mf[lev]->array(mfi);
        baf.intersections(mf[lev].boxArray()[mfi.index()], isects);
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
