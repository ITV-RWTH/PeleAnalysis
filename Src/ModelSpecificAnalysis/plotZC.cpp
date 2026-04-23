#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <PelePhysics.H>

//#include <mechanism.H>

using namespace amrex;

static
void
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile infile=f1 [options] \n\tOptions:\n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile,std::string("/"));
  return tokens[tokens.size()-1];
}

int
main (int   argc,
      char* argv[])
{
  Initialize(argc,argv);
  {
    if (argc < 2)
      print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc,argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(false);

    std::string plotFileName; pp.get("infile",plotFileName);
    std::string fuelName = "H2"; pp.query("fuelName",fuelName);
    //Real s = 8.0; pp.query("stoichRatio",s);
    //Real Y_O_air = 0.233; pp.query("YO2Air",Y_O_air);
    //Real S = s/Y_O_air;
    //const Real Z_st = Y_O_air/(s+Y_O_air);
    std::string productName = "H2O"; pp.query("productName",productName);
    int clipProgress = 0; pp.query("clipProgress",clipProgress);
    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;
    Vector<std::string> spec_names;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
    auto eos = pele::physics::PhysicsType::eos();
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();    
    const std::string spName= "Y("+fuelName+")";
    //const std::string oxName= "Y(O2)";
    const std::string prodName= "Y("+productName+")";

    constexpr int nCompIn = NUM_SPECIES;
    constexpr int nCompOut = 3;
    Vector<std::string> outNames(nCompOut);
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    
    for (int i=0; i<NUM_SPECIES; ++i)
    {
      destFillComps[i] = i;
      inNames[i] =  "Y(" + spec_names[i] + ")";
    }
    //out
    constexpr int idZlocal = 0; // Z out here
    constexpr int idCFlocal = 1; // CF out here
    constexpr int idCPlocal = 2; // CP out here
    outNames[idZlocal] = "Z";
    outNames[idCFlocal] = "CF";
    outNames[idCPlocal] = "CP";       

    Vector<MultiFab> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);
    RealBox rb(&(amrData.ProbLo()[0]),&(amrData.ProbHi()[0]));
    constexpr int nGrow = 0;
    Vector<Real> YFminmax = {1e100,-1e100};
    Vector<Real> YPminmax = {1e100,-1e100};    
    Vector<Box> probDomain = amrData.ProbDomain();

    const int YFcomp = amrData.StateNumber("Y("+fuelName+")")-amrData.StateNumber("Y("+spec_names[0]+")");
    const int YPcomp = amrData.StateNumber("Y("+productName+")")-amrData.StateNumber("Y("+spec_names[0]+")");
    
    for (int lev=0; lev<Nlev;++lev) {
      Real min,max;	
      amrData.MinMax(probDomain[lev], "Y("+fuelName+")", lev,min,max);
      if (YFminmax[0] > min) {
	YFminmax[0] = min;
      }
      if (YFminmax[1] < max) {
	YFminmax[1] = max;
      }      
      amrData.MinMax(probDomain[lev], "Y("+productName+")",lev,min,max);
      if (YPminmax[0] > min) {
	YPminmax[0] = min;
      }
      if (YPminmax[1] < max) {
	YPminmax[1] = max;
      }
    }
    //TODO: add some verb and manual querying
    pp.queryarr("YFminmax",YFminmax);
    pp.queryarr("YPminmax",YPminmax);
    Real Zfu = -1.0;
    Real Zox = -1.0;

    Real YF[NUM_SPECIES], YO[NUM_SPECIES];
    for (int i = 0; i < NUM_SPECIES; ++i) {
      YF[i] = 0.0;
      YO[i] = 0.0;
      if (spec_names[i] == "O2") {
	YO[i] = 0.233;
      }
      if (spec_names[i] == "N2") {
	YO[i] = 0.767;
      }
      if (i == YFcomp) {
	YF[i] = 1.0;
      }
    }

    
    // Detailed chem - compute Bilger coefficients
    // Only interested in CHON -in that order. Compute Bilger weights
    Array<amrex::Real, 4> Beta_mix;
    Real atwCHON[4] = {0.0};
    pele::physics::eos::atomic_weightsCHON<pele::physics::PhysicsType::eos_type>(atwCHON);
    Beta_mix[0] = (atwCHON[0] != 0.0) ? 2.0 / atwCHON[0] : 0.0;
    Beta_mix[1] = (atwCHON[1] != 0.0) ? 1.0 / (2.0 * atwCHON[1]) : 0.0;
    Beta_mix[2] = (atwCHON[2] != 0.0) ? -1.0 / atwCHON[2] : 0.0;
    Beta_mix[3] = 0.0;
    
    // Compute each species weight for the Bilger formulation based on elemental
    // compo Only interested in CHON -in that order.
    amrex::Array<amrex::Real, NUM_SPECIES> spec_Bilger_fact;

    int ecompCHON[NUM_SPECIES * 4];
    pele::physics::eos::element_compositionCHON<pele::physics::PhysicsType::eos_type>(ecompCHON);
    amrex::Real mwt[NUM_SPECIES];
    eos.molecular_weight(mwt);
    Zfu = 0.0;
    Zox = 0.0;
    for (int i = 0; i < NUM_SPECIES; ++i) {
      spec_Bilger_fact[i] = 0.0;
      for (int k = 0; k < 4; k++) {
	spec_Bilger_fact[i] +=
	  Beta_mix[k] * (ecompCHON[i * 4 + k] * atwCHON[k] / mwt[i]);
      }
      Zfu += spec_Bilger_fact[i] * YF[i];
      Zox += spec_Bilger_fact[i] * YO[i];
    }
    const Real denom_inv = 1.0 / (Zfu - Zox);
    //amrex::GpuArray<amrex::Real, NUM_SPECIES> fact_Bilger;
    //for (int n = 0; n < NUM_SPECIES; ++n) {
    //  fact_Bilger[n] = a_pelelm->spec_Bilger_fact[n];
    //}

    
    for (int lev=0; lev<Nlev; ++lev) {
      const BoxArray ba = amrData.boxArray(lev);
      const Vector<Real>& delta = amrData.DxLevel()[lev];
      const DistributionMapping dm(ba);
      MultiFab indata(ba,dm,nCompIn,nGrow);
      outdata[lev].define(ba,dm,nCompOut,nGrow);
      
      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata,lev,inNames,destFillComps); //Problem
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,0,&(is_per[0]));
      Print() << "Data has been read for level " << lev << std::endl;
      auto in_ma = indata.const_arrays();
      auto out_ma = outdata[lev].arrays();
      amrex::ParallelFor(outdata[lev],[=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {	
	out_ma[box_no](i,j,k,idCFlocal) = 1.0 - in_ma[box_no](i,j,k,YFcomp)/YFminmax[1];
	out_ma[box_no](i,j,k,idCPlocal) = in_ma[box_no](i,j,k,YPcomp)/YPminmax[1];
	Real Zloc = 0.0;
	for (int n = 0; n < NUM_SPECIES; ++n) {
	  Zloc += in_ma[box_no](i,j,k,n) * spec_Bilger_fact[n];
	}
	out_ma[box_no](i,j,k,idZlocal) = (Zloc - Zox)*denom_inv;
      });
      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_ZC");
    Print() << "Writing new data to " << outfile << std::endl;
    const bool verb = false;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(outfile,Nlev,GetVecOfConstPtrs(outdata),outNames,geoms, 0.0, isteps, refRatios);
  }
  Finalize();
  return 0;
}
