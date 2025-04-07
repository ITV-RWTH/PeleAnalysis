#ifndef _EBUSERDEFINED_H
#define _EBUSERDEFINED_H

using namespace amrex;

#ifdef AMREX_USE_EB
#include <AMReX_ParmParse.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
void
EBUserDefined(
  const Geometry& geom,
  const int required_coarsening_level,
  const int max_coarsening_level)
{
  
  // ParmParse your geometry parameters
  // --> fixed geometry of the challenge problem
  ParmParse pp("prob");

  //RealArray prechamberHi;
  //GpuArray<Real,AMREX_SPACEDIM> prechamberHi, prechamberLo;
  //Vector<Real> prechamberHi(AMREX_SPACEDIM);
  //Vector<Real> prechamberLo(AMREX_SPACEDIM);
  //pp.getarr("prechamber_hi",prechamberHi);
  //pp.getarr("prechamber_lo",prechamberLo);
  
  Real prechamber_side_x, prechamber_side_y,cylinderRadius, cylinderLength;
  pp.get("cylinder_radius",cylinderRadius);
  pp.get("cylinder_length",cylinderLength);
  pp.get("prechamber_side_x",prechamber_side_x);
  pp.get("prechamber_side_y",prechamber_side_y);
  RealArray prechamberLo,prechamberHi,cylinderCentre,chamberLo,chamberHi;
  AMREX_D_TERM(prechamberLo[0] = geom.ProbLo(0);,prechamberLo[1] = -prechamber_side_y/2.0;,prechamberLo[2] = -prechamber_side_y/2.0;); // Macro defines variables based on dimension  -> Raik question: Third dimension of prechamber?
  Print() << prechamberLo << std::endl; // Just debugging or function?
  AMREX_D_TERM(prechamberHi[0] = geom.ProbLo(0)+prechamber_side_x;,prechamberHi[1] = prechamber_side_y/2.0;,prechamberHi[2] = prechamber_side_y/2.0;);
  Print() << prechamberHi << std::endl;
  AMREX_D_TERM(cylinderCentre[0] = prechamberHi[0]+cylinderLength/2.0;,cylinderCentre[1] = 0.0;,cylinderCentre[2] = 0.0;);
  Print() << cylinderCentre << std::endl;
  AMREX_D_TERM(chamberLo[0] = geom.ProbLo(0) + prechamber_side_x + cylinderLength;,chamberLo[1] = geom.ProbLo(1);,chamberLo[2] = geom.ProbLo(2););
  AMREX_D_TERM(chamberHi[0] = geom.ProbHi(0);,chamberHi[1] = geom.ProbHi(1);,chamberHi[2] = geom.ProbHi(2););

  //RealArray endEB(AMREX_D_DECL(chamberstart,geom.ProbHi()[1],geom.ProbHi()[2]));
  //EB2::BoxIF walls(geom.ProbLo(),endEB,false);

  // Initialisation of the different shapes see https://amrex-codes.github.io/amrex/doxygen/index.html
  EB2::BoxIF prechamber(prechamberLo, prechamberHi,false);
  EB2::CylinderIF cylinder(cylinderRadius,cylinderLength,0,cylinderCentre,false);
  EB2::BoxIF chamber(chamberLo,chamberHi,false);

  auto domaininv = EB2::makeUnion(prechamber,cylinder,chamber); // Combines the shapes to one domain (automatic detects typ of variable)
  auto domain = EB2::makeComplement(domaininv); // Complement of an object. E.g. a sphere with fluid on outside becomes a sphere with fluid inside.
  
  // Build your geometry shop using EB2::makeShop
  auto gshop = EB2::makeShop(domain); // "Shop item" we can reuse? 

  // Build geom using EB2::Build
  EB2::Build(gshop, geom, required_coarsening_level, max_coarsening_level);
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setEBState(
  const amrex::Real* /*xEBface[AMREX_SPACEDIM]*/,
  amrex::Real* /*s_ext[NVAR]*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& /*prob_parm*/)
{
  // Provide the state on the EB centroid provided
  // by xEBface
  // Not the entire has to be specified. For instance, if using isothermal
  // boundary, only s_ext[TEMP] will be used.
  // Aborting here (message will not show on GPUs)
  Abort("Default setEBState called. Copy and edit this file for your needs");
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setEBType(
  const amrex::Real* /*xEBface[AMREX_SPACEDIM]*/,
  amrex::Real& /*EBflagType*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& /*prob_parm*/)
{
  // Provide a factor applied on the EB face diffusion
  // coefficient: 0.0 will turn off the diffusion (homogeneous Neumann BC)
  // while 1.0 will let it unchanged. Intermediate values can be specified
  // for partially Dirichlet BCs
  Abort("Default setEBType called. Copy and edit this file for your needs");
}
#endif
#endif
