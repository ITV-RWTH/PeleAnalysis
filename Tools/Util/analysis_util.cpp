#include <analysis_util.H>

namespace analysis_util {

  int nlev;
  amrex::Vector<amrex::BoxArray> grids;
  amrex::Vector<amrex::DistributionMapping> dmap;
  amrex::Vector<int> ref_ratio;
  amrex::Vector<amrex::Geometry> geoms;
  amrex::IntVect periodicity;
  bool initialized = false;
  
  void init(const amrex::Vector<amrex::MultiFab>& a_mf, const amrex::Vector<amrex::Geometry>& a_geoms) {
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
      ref_ratio.assign(nlev,2);
    }
    if (geoms.empty()) {
      geoms = a_geoms;
    }
    periodicity = a_geoms[0].periodicity().intVect();    
    
    initialized = true;
  }

}
