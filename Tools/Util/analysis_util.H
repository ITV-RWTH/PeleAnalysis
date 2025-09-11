#ifndef analysis_util_H
#define analysis_util_H

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Vector.H>

namespace analysis_util {

  // --- Global variables ---
  extern int nlev;
  extern amrex::Vector<amrex::BoxArray> grids;
  extern amrex::Vector<amrex::DistributionMapping> dmap;
  extern amrex::Vector<int> ref_ratio;
  extern amrex::Vector<amrex::Geometry> geoms;
  extern amrex::IntVect periodicity;
  extern bool initialized;
  
  // --- Functions ---
  // Fill namespace variables
  void init(const amrex::Vector<amrex::MultiFab>& a_mf, const amrex::Vector<amrex::Geometry>& a_geoms);

  // --- Integration tools ---
  // Return a vector of iMultiFabs representing covered regions
  [[nodiscard]]
  amrex::Vector<std::unique_ptr<amrex::iMultiFab>> get_covered_mf();
  
  // Integrate MultiFabs across specified axes
  [[nodiscard]]
  std::unique_ptr<amrex::Gpu::ManagedVector<amrex::Real>> integrate(const amrex::Vector<amrex::MultiFab>& a_mf,int scomp,int ncomp,const amrex::Vector<int>& axes_to_integrate);

  // --- Gradient tools ---

  void gradient(const amrex::Vector<amrex::MultiFab> a_mf, const int scomp, const int ncomp, amrex::Vector<amrex::MultiFab> grad_mf);

  // --- I/O tools ---
  
} // namespace analysis_util

#endif // analysis_util_H
