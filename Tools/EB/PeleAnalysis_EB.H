#ifndef PELEANALYSIS_EB_H
#define PELEANALYSIS_EB_H

#ifdef AMREX_USE_EB

#include <memory>
#include <string>

#include <AMReX_EB2.H>
#include <AMReX_EBFabFactory.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <pelelmex_prob_parm.H>
#include <PeleLMeX_EBUserDefined.H>

namespace pele_analysis {

// Build the EB geometry from the eb2.* inputs and return one EB factory per
// level. geoms, grids and dmap hold levels 0..finest level processed.
//
// eb2.geom_type selects the geometry: "UserDefined" calls EBUserDefined() from
// Src/EBUserDefined, anything else is passed to EB2::Build (an unknown type
// aborts there). The geometry is generated on eb2.max_level_generation and
// coarsened to the coarser levels, so it must be the finest level processed:
// a larger value has no geometry, a smaller one leaves the finer levels
// without EB data.
inline amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>>
buildEBFactories(
  const amrex::Vector<amrex::Geometry>& geoms,
  const amrex::Vector<amrex::BoxArray>& grids,
  const amrex::Vector<amrex::DistributionMapping>& dmap,
  int verbose = 0)
{
  BL_PROFILE("PeleAnalysis::buildEBGeometry()");
  if (verbose > 0) {
    amrex::Print() << "Start building EB!" << std::endl;
  }

  const int nlev = static_cast<int>(geoms.size());
  const int finestLevel = nlev - 1;
  const int max_coarsening_level = 100;
  const int req_coarsening_level = finestLevel;

  amrex::ParmParse ppeb2("eb2");
  std::string geom_type;
  ppeb2.get("geom_type", geom_type);

  int max_lvl_eb = finestLevel;
  ppeb2.query("max_level_generation", max_lvl_eb);
  if (max_lvl_eb != finestLevel) {
    amrex::Abort(
      "eb2.max_level_generation = " + std::to_string(max_lvl_eb) +
      " must equal the finest level processed (" + std::to_string(finestLevel) +
      ")");
  }

  if (geom_type == "UserDefined") {
    EBUserDefined(
      geoms[max_lvl_eb], req_coarsening_level, max_coarsening_level);
  } else {
    amrex::EB2::Build(
      geoms[max_lvl_eb], req_coarsening_level, max_coarsening_level);
  }

  if (verbose > 0) {
    amrex::Print() << "Setting up EBFactory!" << std::endl;
  }
  amrex::Vector<std::unique_ptr<amrex::EBFArrayBoxFactory>> eb_factory(nlev);
  const amrex::EB2::IndexSpace& eb_is = amrex::EB2::IndexSpace::top();
  for (int lev = 0; lev < nlev; ++lev) {
    const amrex::EB2::Level& eb_level = eb_is.getLevel(geoms[lev]);
    eb_factory[lev] = std::make_unique<amrex::EBFArrayBoxFactory>(
      eb_level, geoms[lev], grids[lev], dmap[lev], amrex::Vector<int>{2, 2, 2},
      amrex::EBSupport::full);
  }
  return eb_factory;
}

} // namespace pele_analysis

#endif // AMREX_USE_EB
#endif // PELEANALYSIS_EB_H
