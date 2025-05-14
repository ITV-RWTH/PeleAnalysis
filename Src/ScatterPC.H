#ifndef SCATTER_PARTICLE_CONTAINER_H_
#define SCATTER_PARTICLE_CONTAINER_H_

#include <AMReX_Particles.H>

typedef amrex::Array<amrex::Real,AMREX_SPACEDIM> dim3;
typedef amrex::Array<int,AMREX_SPACEDIM> int3;

struct RealData
{
  enum {
    AMREX_D_DECL(xloc = 0,
                 yloc,
                 zloc),
  };
};

class ScatterParticleContainer
  : public amrex::ParticleContainer<AMREX_SPACEDIM, 0, 0, 0>
{
public:

  using MyParIter = amrex::ParIter<AMREX_SPACEDIM, 0, 0, 0>;

  ScatterParticleContainer(const int                                        a_nPts,
			   const int                                       a_nPartsPerLoc,
                           const amrex::Vector<amrex::Geometry>            & a_geoms,
                           const amrex::Vector<amrex::DistributionMapping> & a_dmaps,
                           const amrex::Vector<amrex::BoxArray>            & a_bas,
                           const amrex::Vector<int>                        & a_rrs,
			   const int                                       a_nVar,
			   const amrex::IntVect                            a_is_per,
			   const amrex::Vector<std::string>                & a_outVarNames,
			   const int                                       a_plane);

  void InitParticles(const amrex::Vector<dim3>& a_locs);

  void SetParticleLocation(int a_streamPos);

  void ComputeNextLocation(int                                    a_fromLoc,
                           const amrex::Real                            a_stepSize);

  void InterpDataAtLocation(int                                    a_fromLoc,
			    const amrex::Vector<amrex::MultiFab> & a_vectorField,
			    const amrex::Vector<int>               a_order);
  
  
  void WriteScatterAsBinary(const std::string& outfile,
			   amrex::Vector<int>& faceData);
  
  void WriteScatterAsTecplot(const std::string& outfile);
  
protected:
  int Nlev;
  int nPts;
  int nPartsPerLoc;
  int sizeOfRealData;
  int pcomp;
  int fcomp;
  int plane;
  amrex::IntVect is_per;
  amrex::Vector<std::string> outVarNames;
  amrex::Vector<dim3> directions;
};

#endif
