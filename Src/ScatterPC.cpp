#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Utility.H>
#include <AMReX_MultiFab.H>

#include "ScatterPC.H"

using namespace amrex;

ScatterParticleContainer::
ScatterParticleContainer(const int                                 a_nPts,
			 const int                        a_nPartsPerLoc, 			 
			 const Vector<Geometry>            & a_geoms,
			 const Vector<DistributionMapping> & a_dmaps,
			 const Vector<BoxArray>            & a_bas,
			 const Vector<int>                 & a_rrs,
			 const int                         a_nVar,
			 const IntVect                     a_is_per,
			 const Vector<std::string>         & a_outVarNames,
			 const int                         a_plane)
: ParticleContainer<AMREX_SPACEDIM, 0, 0, 0> (a_geoms, a_dmaps, a_bas, a_rrs)
{
  Nlev = a_geoms.size();
  nPts = a_nPts;
  nPartsPerLoc = a_nPartsPerLoc;
  fcomp = a_nVar;
  pcomp = fcomp + 1 + AMREX_SPACEDIM;
  sizeOfRealData = pcomp * nPts;
  is_per = a_is_per;
  outVarNames = a_outVarNames;
  ResizeRuntimeRealComp(sizeOfRealData,true);
  directions.resize(nPartsPerLoc);
  plane = a_plane;
}


static void fibonacci_lattice(Vector<dim3>& points, const int nPts) {
  constexpr Real pi = 3.14159265358979323846;
  constexpr Real phi = pi * (3.0-std::sqrt(5.0));
  if (nPts > 1) {
    for (int i = 0; i < nPts; ++i) {
      Real y = 1 - (i / double(nPts - 1)) * 2;
      Real radius = sqrt(1 - y * y);
      
      Real theta = phi * i;
      
      // Convert spherical to Cartesian coordinates
      Real x = radius * std::cos(theta);
      Real z = radius * std::sin(theta);
      points[i] = {AMREX_D_DECL(x,y,z)};
    }
  } else {
    points[0] = {AMREX_D_DECL(1,0,0)};    
  }
}

static void planar_lattice(Vector<dim3>& points, const int nPts ,const int plane) {
  constexpr Real twopi = 2*3.14159265358979323846;
  const Real phi = twopi / double(nPts);
  for (int i = 0; i < nPts; ++i) {
    points[i][plane] = 0;
    points[i][plane+1 % AMREX_SPACEDIM] = std::sin(i*phi);
#if AMREX_SPACEDIM == 3
    points[i][plane+2 % AMREX_SPACEDIM] = std::cos(i*phi);
#endif
    //Print() << points[i] << std::endl;
  }
}


void
ScatterParticleContainer::
InitParticles (const Vector<dim3>& locs)
{
  BL_PROFILE("ScatterParticleContainer::InitParticles");

  const int nProc = ParallelDescriptor::NProcs();
  const int myProc = ParallelDescriptor::MyProc();
  const int nLocs = locs.size();
  auto& particle_tile = DefineAndReturnParticleTile(0, myProc, 0);
  //lets just initialise on a random processor and redistribute afterwards 
  //trying to find the right processor here doesn't work well with particles near boundaries

  if (plane < 0) {
    fibonacci_lattice(directions, nPartsPerLoc);	
  } else {
    planar_lattice(directions,nPartsPerLoc,plane);
  }
  for (int n = 0; n < nLocs; ++n) {
    if (n%nProc == myProc) {
    //we want multiple particles starting in the same place
      for (int k = 0; k < nPartsPerLoc; ++k) {
	const int p_id = nPartsPerLoc*n+k+1; 
	ParticleType p;
	p.id() = p_id;
	p.cpu() = myProc;
	//Initialise at same position	
	AMREX_D_EXPR(p.pos(0) = locs[n][0],p.pos(1) = locs[n][1],p.pos(2) = locs[n][2]);	
	//each particle gets given a vector to go in
	AMREX_D_EXPR(p.rdata(0) = directions[k][0],p.rdata(1) = directions[k][1],p.rdata(2) = directions[k][2]);
	//Add each particle to the tile and allocate real space in order
	particle_tile.push_back(p);
	for (int i = 0; i<NumRuntimeRealComps(); ++i) {
	  particle_tile.push_back_real(i, i < AMREX_SPACEDIM ? locs[n][i] : 0.0);
	}
	//std::cout << "NReal " << p.NReal << std::endl;
      }
    }
  }
  Redistribute();
}


void
ScatterParticleContainer::
SetParticleLocation(const int a_Loc)
{
  BL_PROFILE("ScatterParticleContainer::SetParticleLocation");
  //offset to find coordinates
  int offset = pcomp * a_Loc;// + AMREX_SPACEDIM;
  dim3 newpos;
  //std::cout << "Setting particle location at position " << a_Loc << std::endl;
  for (int lev = 0; lev < Nlev; ++lev)
  {
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();

      for (size_t pindex=0; pindex<aos.size(); ++pindex)
      {
        ParticleType& p = aos[pindex];
        if (p.id() > 0)
        {
          newpos = {AMREX_D_DECL(soa.GetRealData(offset + RealData::xloc)[pindex],
                                 soa.GetRealData(offset + RealData::yloc)[pindex],
                                 soa.GetRealData(offset + RealData::zloc)[pindex])};
	  //Print() << "New position of particle " << pindex << " :" << newpos << std::endl;
	  AMREX_D_EXPR(p.pos(0) = newpos[0], p.pos(1) = newpos[1], p.pos(2) = newpos[2]);
	  //Print() << "Set position of particle " << pindex << " :" << p.pos(0) << " " << p.pos(1) << " " << p.pos(2) << std::endl;
        }
      }
    }
  }
  Redistribute();
}


void
ScatterParticleContainer::
ComputeNextLocation(const int                      a_fromLoc,
                    const Real                     a_stepSize)
{
  BL_PROFILE("ScatterParticleContainer::ComputeNextLocation");
  //Print() << "about to set location" << std::endl;
  SetParticleLocation(a_fromLoc);
  const Real dxFine = Geom(Nlev - 1).CellSize()[0];
  const Real dr = dxFine*a_stepSize;
  const int new_loc_id = a_fromLoc + 1;
  const int offset = pcomp * a_fromLoc;// + AMREX_SPACEDIM;
  const int offset_new = pcomp * new_loc_id;
  for (int lev = 0; lev < Nlev; ++lev)
  {
    //const auto& geom = Geom(lev);
    //const auto& dx = geom.CellSize();
    //const auto& plo = geom.ProbLo(); 
    //const auto& phi = geom.ProbHi(); 

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();
      
      for (size_t pindex=0; pindex<aos.size(); ++pindex) {
        
	ParticleType& p = aos[pindex];
	//Print() << p.pos(0) << " " << p.pos(1) << " " << p.pos(2) << std::endl;
	dim3 x = {AMREX_D_DECL(soa.GetRealData(offset + RealData::xloc)[pindex], soa.GetRealData(offset + RealData::yloc)[pindex], soa.GetRealData(offset + RealData::zloc)[pindex])};
	dim3 dir = {AMREX_D_DECL(p.rdata(0), p.rdata(1), p.rdata(2))};
	//Print() << "aLoc " << a_fromLoc << std::endl;
	//Print() << "pindex " << pindex << std::endl;
	//Print() << "current loc: " << x << std::endl;
	//Print() << "going to move in dir " << dir << " with size " << dr << std::endl;
	if (p.id() > 0){
	  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
	    x[idim] += dir[idim]*dr;
	  }
        }
	//Print() << "next loc: " <<  x << std::endl;
	
	//std::cout << "next location " << x[0] << " " << x[1] << " " << x[2] << std::endl;
	AMREX_D_EXPR(soa.GetRealData(offset_new + RealData::xloc)[pindex] = x[0],
                     soa.GetRealData(offset_new + RealData::yloc)[pindex] = x[1],
                     soa.GetRealData(offset_new + RealData::zloc)[pindex] = x[2]);
	//also put the distance we are
	soa.GetRealData(offset_new + AMREX_SPACEDIM)[pindex] = (a_fromLoc+1)*dr;
      }
    }
  }
}


static bool ntrpv(const dim3 x,const FArrayBox& gfab,
                  const Real* dx,const Real* plo,const Real* phi,Vector<Real>& u, int nComp)
{
  int3 b;
  dim3 n;

  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    b[d] = (int)std::floor( (x[d] - plo[d]) / dx[d] - 0.5 );
    n[d] = ( x[d] - ( (b[d] + 0.5 ) * dx[d] + plo[d] ) )/dx[d];
    n[d] = std::max(0., std::min(1.,n[d]));
  }

  const auto& gbx = gfab.box();
  const auto& glo = gbx.smallEnd();
  const auto& ghi = gbx.bigEnd();
  for (int d=0; d<AMREX_SPACEDIM; ++d) {
    if (b[d] < glo[d] ||  b[d] > ghi[d]-1) {
      std::cout << "dir: " << d << std::endl;
      std::cout << "d,b,glo,ghi: " << d << " " << b[d] << " " << glo[d] << " " << ghi[d] << std::endl;
      std::cout << "x,plo,phi " << x[d] << " " << plo[d] << " " << phi[d] << std::endl;
      std::cout << "boxlo,boxhi " << plo[d]+glo[d]*dx[d] << " " << plo[d]+(ghi[d]+1)*dx[d] << std::endl;
      return false;
    }
  }
  const auto& g = gfab.array();
  for (int i=0; i<nComp; ++i) {
#if AMREX_SPACEDIM == 2
    u[i] =
      +   n[0]   *   n[1]   * g(b[0]+1,b[1]+1,0,i)
      +   n[0]   * (1-n[1]) * g(b[0]+1,b[1]  ,0,i)
      + (1-n[0]) *   n[1]   * g(b[0]  ,b[1]+1,0,i)
      + (1-n[0]) * (1-n[1]) * g(b[0]  ,b[1]  ,0,i);
#else
    u[i] =
      +    n[0]   *    n[1]  *    n[2]  * g(b[0]+1,b[1]+1,b[2]+1,i)
      +    n[0]   * (1-n[1]) *    n[2]  * g(b[0]+1,b[1]  ,b[2]+1,i)
      +    n[0]   *    n[1]  * (1-n[2]) * g(b[0]+1,b[1]+1,b[2]  ,i)
      +    n[0]   * (1-n[1]) * (1-n[2]) * g(b[0]+1,b[1]  ,b[2]  ,i)
      +  (1-n[0]) *    n[1]  *    n[2]  * g(b[0]  ,b[1]+1,b[2]+1,i)
      +  (1-n[0]) * (1-n[1]) *    n[2]  * g(b[0]  ,b[1]  ,b[2]+1,i)
      +  (1-n[0]) *    n[1]  * (1-n[2]) * g(b[0]  ,b[1]+1,b[2]  ,i)
      +  (1-n[0]) * (1-n[1]) * (1-n[2]) * g(b[0]  ,b[1]  ,b[2]  ,i);
#endif
  }
  return true;
}

void
ScatterParticleContainer::
InterpDataAtLocation(int                      a_fromLoc,
		     const Vector<MultiFab> & a_data,
		     const Vector<int>        a_order)
{
  BL_PROFILE("ScatterParticleContainer::InterpDataAtLocation");

  SetParticleLocation(a_fromLoc);

  int offset = pcomp * a_fromLoc + 1 + AMREX_SPACEDIM; // components on particle

  for (int lev = 0; lev < Nlev; ++lev)
  {
    const auto& geom = Geom(lev);
    const auto& dx = geom.CellSize();
    const auto& plo = geom.ProbLo();
    const auto& phi = geom.ProbHi();

    dim3 Lx;
    for (int iComp=0; iComp<AMREX_SPACEDIM; iComp++)
      Lx[iComp] = phi[iComp]-plo[iComp];

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();
      const FArrayBox& data = a_data[lev][pti];

      for (size_t pindex=0; pindex<aos.size(); ++pindex)
      {
        ParticleType& p = aos[pindex];
	if (p.id()>0) {
	  // where's the particle?
	  dim3 x = {AMREX_D_DECL(p.pos(0), p.pos(1), p.pos(2))};
	  Vector<Real> ntrpvOut(fcomp); // components in infile

	  // interpolate all data to particle location
	  ntrpv(x,data,dx,plo,phi,ntrpvOut,fcomp); // components in infile

	  // copy the interpolated data to the particle
	  // first DIM components are particle location
	  // next DIM components are stream location w/o adjusting for periodicity
	  // then we have the interpolated data we want

	  for (int iComp=0; iComp<fcomp; ++iComp) {
	    const int idxOnPart = offset + iComp;	    
	    if (a_fromLoc == 0) {
	      soa.GetRealData(idxOnPart)[pindex] = ntrpvOut[iComp];
	    } else {
	      soa.GetRealData(idxOnPart)[pindex] = std::pow(ntrpvOut[iComp]-soa.GetRealData(iComp + AMREX_SPACEDIM + 1)[pindex],a_order[iComp]);
	      //if we've got to the end overwrite the start point
	      if (a_fromLoc == nPts-1) {
		soa.GetRealData(1+AMREX_SPACEDIM+iComp)[pindex] = 0.0;
	      }
	    }
	    //std::cout << "pidx: " << pindex << ",pt: " << a_fromLoc << ",comp: " << iComp << ",val: " << soa.GetRealData(idxOnPart)[pindex] << std::endl;
	    
	    //subtract centre data
	  }
	  // let's figure out the location on the stream w/o adjusting for periodicity
	  /*if (a_fromLoc==0) { // nothing to do on the surface; just take a copy of the location
	    for (int d=0; d<AMREX_SPACEDIM; ++d) {
	      int idx = offset + d;
	      soa.GetRealData(idx+AMREX_SPACEDIM)[pindex] = soa.GetRealData(idx)[pindex];
	    }
	    } */
	  /*
	  if (a_fromLoc != 0) { // set new location adjusting for periodicity
	    // calculate the change in position delta
	    // add to the old position
	    // adjust delta if it's affected by periodicity (i.e. delta too big)
	    for (int d=0; d<AMREX_SPACEDIM; ++d) {
	      int  idx    = offset + d;
	      int  idxOld = idx - pcomp;
	      Real xnew   = soa.GetRealData(idx   )[pindex];
	      Real xold   = soa.GetRealData(idxOld)[pindex];
	      Real delta  = xnew-xold;
	      if (fabs(delta)>dx[d]) { // has been adjusted for periodicity
		//printf("%i %i %e %e %e",pindex,d,xnew,xold,delta);
		if (delta<0.) delta+=Lx[d];
		else          delta-=Lx[d];
		//printf(" --> %e\n",delta);
	      }
	      // store periodicity-adjusted copy of location at idx+SPACEDIM
	      Real sold   = soa.GetRealData(idxOld)[pindex];
	      Real snew   = sold+delta;
	      soa.GetRealData(idx)[pindex] = snew;
	    }
	  }
	  */
	}
      }
    }
  }
}


void
ScatterParticleContainer::
WriteScatterAsTecplot(const std::string& outFile)
{
  // Set location to first point on stream to guarantee partner line is local
  SetParticleLocation(0);
  
  // Create a folder and have each processor write their own data, one file per streamline
  auto myProc = ParallelDescriptor::MyProc();

  if (!amrex::UtilCreateDirectory(outFile, 0755))
    amrex::CreateDirectoryFailed(outFile);
  ParallelDescriptor::Barrier();

  bool will_write = false;
  for (int lev = 0; lev < Nlev && !will_write; ++lev)
  {
    for (MyParIter pti(*this, lev); pti.isValid() && !will_write; ++pti)
    {
      auto& aos = pti.GetArrayOfStructs();

      for (size_t pindex=0; pindex<aos.size() && !will_write; ++pindex)
      {
        ParticleType& p = aos[pindex];
        will_write |= (p.id() > 0);
      }
    }
  }

  if (will_write)
  {
    std::string fileName = outFile + "/str_";
    fileName = Concatenate(fileName,myProc) + ".dat";
    std::ofstream ofs(fileName.c_str());
    ofs << "VARIABLES = \"";
    for (int iComp=0; iComp<pcomp-1; ++iComp)
      ofs << outVarNames[iComp] << "\" \"";
    ofs << outVarNames[pcomp - 1] << "\"\n";
    
    for (int lev = 0; lev < Nlev; ++lev)
    {
      for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
      {
        auto& aos = pti.GetArrayOfStructs();
        auto& soa = pti.GetStructOfArrays();

        for (size_t pindex=0; pindex<aos.size(); ++pindex)
        {
          ofs << "ZONE I=1 J=" << nPts << " K=1 FORMAT=POINT\n";
	  for (int j=0; j<nPts; ++j)
	    {
	      // by including the spacedim offset, we use locations w/o periodicity adjustments
	      int offset = j*pcomp;
	      //std::cout << offset << std::endl;
	      for (int iComp=0; iComp<pcomp; ++iComp) {
		ofs << soa.GetRealData(offset + iComp)[pindex] << " ";
	      }
	      ofs << '\n';
	      
	    }
        }
      }
    }
    ofs.close();
  }
}


void
ScatterParticleContainer::
WriteScatterAsBinary(const std::string& outFile,
		    Vector<int>& faceData)
{
  // Set location to first point on stream to guarantee partner line is local
  //SetParticleLocation(0);

  // Create a folder and have each processor write their own data, one file per streamline
  auto myProc = ParallelDescriptor::MyProc();
  auto nProcs = ParallelDescriptor::NProcs();

  if (!amrex::UtilCreateDirectory(outFile, 0755))
    amrex::CreateDirectoryFailed(outFile);
  ParallelDescriptor::Barrier();

  bool will_write = false;
  for (int lev = 0; lev < Nlev && !will_write; ++lev)  {
    //Print() << will_write ? 1 : 0 << std::endl;
    for (MyParIter pti(*this, lev); pti.isValid() && !will_write; ++pti) {
      auto& aos = pti.GetArrayOfStructs();	
      for (size_t pindex=0; pindex<aos.size() && !will_write; ++pindex) {
	ParticleType& p = aos[pindex];
	will_write |= (p.id() > 0);
      }	
    }
  }
  
  // Need to count the total number of streams to be written
  // by all ptiters on all levels on this processor
  
  int nLines = 0;
    
  for (int lev = 0; lev < Nlev; ++lev) {
    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
      auto& aos = pti.GetArrayOfStructs();
      for (size_t pindex=0; pindex<aos.size(); ++pindex) {
	ParticleType& p = aos[pindex];
	if ( p.id()>0 ) {
	  nLines++;
	}
      }
    }
  }

  // write to a binary file
  std::string rootName = outFile + "/str_";
  std::string fileName = Concatenate(rootName,myProc) + ".bin";
  //std::string headName = Concatenate(rootName,myProc) + ".head";
  FILE *file=fopen(fileName.c_str(),"w");
  //FILE *head=fopen(headName.c_str(),"w");
  // total number of streams in file
  fwrite(&(nLines),sizeof(int),1,file);
  
  //int minId=100000000;
  //int maxId=-minId;
  for (int lev = 0; lev < Nlev; ++lev) {

    for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
      auto& aos = pti.GetArrayOfStructs();
      auto& soa = pti.GetStructOfArrays();
      
      size_t aosSize = aos.size();
      // construct the pair mapping for this pti
      /*std::map<int,int> pid_to_pindex;
      for (size_t pindex=0; pindex<aosSize; ++pindex) {
      ParticleType& p = aos[pindex];
	int pId = p.id();
	if (pId>0) {
	  pid_to_pindex[pId]=pindex;
	}
	}
      */
      
      for (size_t pindex=0; pindex<aosSize; ++pindex) {
	ParticleType& p = aos[pindex];	  
	int pId         = p.id();
	
	if (pId>0) {
	  // write info about this stream and its pair
	  
	  //minId=min(minId,pId);
	  //maxId=max(maxId,pId);
	  pId -= 1;
	  //Print() << pId << std::endl;
	  fwrite(&(pId),sizeof(int),1,file); // pair id

	  // if we write in the order paricle->component->position on surface,
	  // then end up with a single-component stream together in memory
	  // by including spacedim in offset, we use locations w/o periodicity adjustments
	  
	  for (int iComp=0; iComp<pcomp; ++iComp) {
	    for (int j=0; j<nPts; ++j) {
	      int offset = j*pcomp + iComp;

	      //std::cout << "pidx: " << pindex << ",pt: " << j << ", comp: " << iComp << ",val: " << soa.GetRealData(offset)[pindex] << std::endl; 
	      //Print() << "iLine " << pId << std::endl;
	      //Print() << "iComp " << iComp << std::endl;
	      //Print() << "iPt " << j << std::endl;
	      //Print() << soa.GetRealData(offset)[pindex] << std::endl;
	      //Print() << offset << std::endl;
	      fwrite(&(soa.GetRealData(offset)[pindex]),sizeof(Real),1,file);
	    }
	  }
	  //nStreamsCheck++; // sanity check to make sure we wrote number of streams anticipated
	}
      }
    }
  }
  
  /*
  ParallelDescriptor::ReduceIntMin(minId);
  ParallelDescriptor::ReduceIntMax(maxId);
  
  if (nStreams!=nStreamsCheck)
    std::cout << "(nStreams!=nStreamsCheck) : "
	      << nStreams << " != "  << nStreamsCheck << std::endl;
  */
  fclose(file);
  
  // Write header file with everything consistent across all processors
  ParallelDescriptor::ReduceIntSum(nLines);
  if (ParallelDescriptor::IOProcessor()) {
    fileName = outFile + "/Header";
    std::ofstream ofs(fileName.c_str());
    ofs << "Even odder-ball replacement for sampled streams" << std::endl;
    ofs << nProcs << std::endl;         // translates to number of files to read
    ofs << nLines/nPartsPerLoc << std::endl;       // number of locations
    ofs << nPartsPerLoc << std::endl;    //number of particles at each point
    ofs << nPts << std::endl; // number of points on each line
    ofs << pcomp << std::endl;      // number of variables
    for (int iComp=0; iComp<pcomp; ++iComp)
      ofs << outVarNames[iComp] << " ";
    ofs << std::endl;
    if (!faceData.empty()) {
      ofs << faceData.size() << std::endl;
      ofs.write((char*)faceData.dataPtr(),sizeof(int)*faceData.size());
      ofs << '\n';      
    }
    ofs.close();
  }
}
