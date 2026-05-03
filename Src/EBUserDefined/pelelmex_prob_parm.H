#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>

#ifdef USE_HIT
#include <AMReX_GpuMemory.H>
#endif

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean   = 101325.0_rt;
#if AMREX_USE_EB
  amrex::Real prechamber_side_x;
  amrex::Real prechamber_side_y;
  #if USE_DI
    // Placeholder
  #else
    amrex::Real kernel_diameter = 1e-3;
  #endif
#else
  amrex::Real jet_diameter  = 0.016;
  amrex::Real flame_height  = 0.1;
  amrex::Real wall_temp = 298.0;
#endif

#if USE_HIT
int input_resolution = 0;
  amrex::Real urms0 = 1.0;
  amrex::Real uin_norm = 1.0;

  amrex::Real Linput = 0.0;

  amrex::Real* d_uinput = nullptr;
  amrex::Real* d_vinput = nullptr;
  amrex::Real* d_winput = nullptr;
  amrex::Real* d_xarray = nullptr;
  amrex::Real* d_xdiff = nullptr;
#endif

};
#endif


// Struct "ProbParm" with the user-defined variables for the usage in e.g. pelelmex_prob.H -> Allows the change of these run-time parameter without needing to compile again? 