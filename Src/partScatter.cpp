#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <ScatterPC.H>

using namespace amrex;
static
Vector<dim3>
GetSeedLocations (const ScatterParticleContainer& spc, Vector<int>& faceData)
{
  Vector<dim3> locs;

  ParmParse pp;
  int nc=pp.countval("oneSeedPerCell");
  int ni=pp.countval("isoFile");
  int ns=pp.countval("seedLoc");
  int nrL=pp.countval("seedRakeL");
  int nrR=pp.countval("seedRakeR");
  AMREX_ALWAYS_ASSERT((nc>0) ^ ((ni>0) ^ ((ns>0) ^ ((nrL>0) && nrR>0))));
  if (nc>0)
  {
    int finestLevel = spc.numLevels() - 1;
    std::vector< std::pair<int,Box> > isects;
    FArrayBox mask;
    for (int lev=0; lev<=finestLevel; ++lev)
    {
      const auto& geom = spc.Geom(lev);
      const auto& dx = geom.CellSize();
      const auto& plo = geom.ProbLo();

      BoxArray baf;
      if (lev < finestLevel) {
        baf = BoxArray(spc.ParticleBoxArray(lev+1)).coarsen(spc.GetParGDB()->refRatio(lev));
      }
      for (MFIter mfi = spc.MakeMFIter(lev); mfi.isValid(); ++mfi)
      {
        const Box& tile_box  = mfi.tilebox();
        if (AMREX_SPACEDIM<3 || tile_box.contains(IntVect(AMREX_D_DECL(0,50,107)))) {

          mask.resize(tile_box,1);
          mask.setVal(1);
          if (lev < finestLevel) {
            isects = baf.intersections(tile_box);
            for (const auto& p : isects) {
              mask.setVal(0,p.second,0,1);
            }
          }

          for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
          {
            if (mask(iv,0) > 0)
            {
              locs.push_back({AMREX_D_DECL(plo[0] + (iv[0] + 0.5)*dx[0],
                                           plo[1] + (iv[1] + 0.5)*dx[1],
                                           plo[2] + (iv[2] + 0.5)*dx[2])});
            }
          }
        }
      }
    }
  }
  else if (ni>0)
  {
    // Read in isosurface
    std::string isoFile; pp.get("isoFile",isoFile);
    if (ParallelDescriptor::IOProcessor())
      std::cerr << "Reading isoFile... " << isoFile << std::endl;

    std::ifstream ifs;
    ifs.open(isoFile.c_str());
    // AJA added dummy line read; sometimes time, sometimes `decimated'
    std::string topline;
    std::getline(ifs,topline);
    std::string line;
    std::getline(ifs,line);
    auto surfNames = Tokenize(line,std::string(", "));
    int nCompSeedNodes = surfNames.size();
    int nElts, nodesPerElt;
    ifs >> nElts;
    ifs >> nodesPerElt;

    FArrayBox tnodes;
    tnodes.readFrom(ifs);
    int nSeedNodes = tnodes.box().numPts();

    Real* ndat = tnodes.dataPtr();
    for (int i=0; i<nSeedNodes; ++i)
    {
      int o=i*nCompSeedNodes;
      locs.push_back({AMREX_D_DECL(ndat[o+0], ndat[o+1], ndat[o+2])});
    }
    tnodes.clear();

    faceData.resize(nElts*nodesPerElt);
    ifs.read((char*)faceData.dataPtr(),sizeof(int)*faceData.size());
    ifs.close();
  }
  else if (pp.countval("seedLoc")>0)
  {
    Vector<Real> loc(BL_SPACEDIM);
    pp.getarr("seedLoc",loc,0,BL_SPACEDIM);
    locs.push_back({AMREX_D_DECL(loc[0], loc[1], loc[2])});
  }
  else
  {
    int seedRakeNum;
    pp.get("seedRakeNum",seedRakeNum);
    AMREX_ALWAYS_ASSERT(seedRakeNum >= 2);
    Vector<Real> locL(BL_SPACEDIM), locR(BL_SPACEDIM);
    pp.getarr("seedRakeL",locL,0,BL_SPACEDIM);
    pp.getarr("seedRakeR",locR,0,BL_SPACEDIM);

    for (int i=0; i<seedRakeNum; ++i) {
      locs.push_back({AMREX_D_DECL(locL[0] + (i/double(seedRakeNum-1))*(locR[0] - locL[0]),
                                   locL[1] + (i/double(seedRakeNum-1))*(locR[1] - locL[1]),
                                   locL[2] + (i/double(seedRakeNum-1))*(locR[2] - locL[2]))});
    }
  }
  return locs;
}

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    ParmParse pp;

    std::string infile; pp.get("infile",infile);

    int nGrad = pp.countval("gradientField");
    AMREX_ALWAYS_ASSERT(nGrad == 0 || nGrad == AMREX_SPACEDIM);    
    Vector<std::string> gradVarNames(nGrad);
    pp.queryarr("gradientField",gradVarNames);
    const bool cartesian = (nGrad == 0);
    int plane = -1;
    pp.query("plane",plane);
    int nVars= pp.countval("vars");
    Vector<std::string> inVarNames(nVars);
    Vector<int> orders(nVars);
    pp.getarr("orders",orders);
    pp.getarr("vars",inVarNames);
    PlotFileData pf(infile);
    int finestLevel = pf.finestLevel();
    Vector<Geometry> geoms(finestLevel+1);
    Vector<BoxArray> grids(finestLevel+1);
    Vector<DistributionMapping> dms(finestLevel+1);
    Vector<int> ratios(finestLevel);
    const int nComp = nGrad + nVars;// + AMREX_SPACEDIM;
    
    // get base for output files
    std::string outfile = infile;
    pp.query("outfile",outfile);

    int writeParticles(0);
    pp.query("writeParticles",writeParticles);
    std::string particlefile = outfile+"_particles";
    pp.query("particlefile",particlefile);
    int writeStreams(0);
    pp.query("writeStreams",writeStreams);
    std::string streamfile = outfile+"_stream";
    pp.query("streamfile",streamfile);
    //
    int writeStreamBin(1);
    pp.query("writeStreamBin",writeStreamBin);
    std::string streamBinfile = outfile+"_streamBin";
    pp.query("streamBinfile",streamBinfile);
    
    Vector<std::string> outVarNames = {AMREX_D_DECL("X","Y","Z"),"R"};
    for (int n = 0; n < nVars; n++) {
      outVarNames.push_back(inVarNames[n]);
    }

    IntVect pp_is_per;
    pp.getarr("is_per",pp_is_per);
    Array<int,AMREX_SPACEDIM> is_per = {AMREX_D_DECL(pp_is_per[0],pp_is_per[1],pp_is_per[2])};
    RealBox rb(pf.probLo(),pf.probHi());

    int Nlev = finestLevel + 1;
    Vector<Vector<MultiFab>> pfdata(Nlev);
    for (int lev=0; lev<Nlev; ++lev) {
      geoms[lev].define(pf.probDomain(lev),rb,pf.coordSys(),is_per);
      grids[lev] = pf.boxArray(lev);
      dms[lev] = pf.DistributionMap(lev);
      if (lev < finestLevel) ratios[lev] = pf.refRatio(lev);

      pfdata[lev].resize(nComp);
      Print() << "Loading data on level " << lev << std::endl;
      for (int d=0; d<nGrad; ++d) {
        pfdata[lev][d] = pf.get(lev,gradVarNames[d]);
      }
      for (int n = 0; n < nComp; n++) {
	pfdata[lev][n+nGrad] = pf.get(lev,inVarNames[n]);
      }
    }
    //does this just need hardcoding?
    int nGrow = 3;
    pp.query("nGrow",nGrow);
      
    Real time=pf.time();
    PhysBCFunctNoOp f;
    PCInterp cbi;
    BCRec bc;
    AMREX_ALWAYS_ASSERT(nGrow>=1);
    Vector<MultiFab> indata(Nlev);
    for (int lev=0; lev<Nlev; ++lev) {
      indata[lev].define(grids[lev],dms[lev],nComp,nGrow);
      for (int d=0; d<nComp; ++d) {
        if (lev==0) {
          FillPatchSingleLevel(indata[lev],time,{&pfdata[lev][d]},{time},0,d,1,geoms[0],f,0);
        }
        else
        {
          FillPatchTwoLevels(indata[lev],time,{&pfdata[lev-1][d]},{time},{&pfdata[lev][d]},{time},0,d,1,
                             geoms[lev-1],geoms[lev],f,0,f,0,ratios[lev-1]*IntVect::Unit,&cbi,{bc},0);
        }
      }
      indata[lev].FillBoundary(geoms[lev].periodicity());
    }

    int Nsteps = 50;
    pp.query("Nsteps",Nsteps);
    int particlesPerLoc = 1; pp.query("particlesPerLoc",particlesPerLoc);
    ScatterParticleContainer spc(Nsteps,particlesPerLoc,geoms,dms,grids,ratios,nComp,pp_is_per,outVarNames,plane);
    
    Vector<int> faceData;
    auto locs = GetSeedLocations(spc,faceData);
    if (writeStreamBin == 1 && faceData.empty()) {
      Print() << "Writing stream binary without surface definition - not a problem if there is no surface given :)" << std::endl;
    }

    // Initialise particles
    Print() << "Initialising particles..." << std::endl;
    spc.InitParticles(locs);
    
    //Check if particles initialised fine
    if (spc.OK()) {
      Print() << "SPC is happy with initialisation :-)" << std::endl;
    }
       
    // Interpolate at start
    Print() << "Interpolation at the seed points..." << std::endl;
    spc.InterpDataAtLocation(0,indata,orders);
    
    Print() << "Computing streams and interpolating..." << std::endl;
    Real stepSize = 0.1; pp.query("stepSize",stepSize);
    AMREX_ALWAYS_ASSERT(stepSize>=0 && stepSize<=0.5);
    for (int step=0; step<Nsteps-1; ++step)
    {
      // find next location
      spc.ComputeNextLocation(step,stepSize);

      // interpolate all data
      spc.InterpDataAtLocation(step+1,indata,orders);

    }

    // check in again
    if (spc.OK()) {
      Print() << "SPC is happy afterwards :-)" << std::endl;
    }
        
    if (writeParticles) {
      Print() << "Writing particles in plotfile to " << particlefile << std::endl;
      spc.WritePlotFile(particlefile, "particles");
    }
    
    if (writeStreams) {
      Print() << "Writing streamlines in Tecplot ascii format to " << streamfile << std::endl;
      spc.WriteScatterAsTecplot(streamfile);
    }
    
    //
    // Write streamBin
    //
    if (writeStreamBin) {
      Print() << "Writing streamlines as binary " << streamBinfile << std::endl;
      spc.WriteScatterAsBinary(streamBinfile,faceData);
    }
    
  }
  Finalize();
  return 0;
}
