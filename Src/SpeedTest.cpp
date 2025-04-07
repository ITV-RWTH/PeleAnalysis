#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

// Time measurement 
#include <chrono>

using namespace amrex;

int
main (int   argc,
      char* argv[])
{
  amrex::Initialize(argc,argv);
  {


    // ---------------------------------------------------------------------
    // Set defaults input values
    // ---------------------------------------------------------------------
    std::string gradVar       = "temp";
    std::string infile        = "";  
    int finestLevel           = 1000;
    int nAuxVar               = 0;

    // ---------------------------------------------------------------------
    // ParmParse
    // ---------------------------------------------------------------------
    ParmParse pp;

    pp.get("infile",infile);
    Print() << "infile = " << infile << std::endl; 
    pp.query("gradVar",gradVar);
    Print() << "This is gradVar: " << gradVar << std::endl;
    pp.query("finestLevel",finestLevel);

    // std::string path;
    // pp.get("path",path);
    // Print() << "output path = " << path << std::endl; 


    // Initialize DataService
    Print() << "Initialising DataService!" << std::endl;
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if( ! dataServices.AmrDataOk()) {
      Print() << "Data not okay!" << std::endl;
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }

    Print() << "Initialised amrData!" << std::endl;
    AmrData& amrData = dataServices.AmrDataRef();
    Print() << "Initialised DataService!" << std::endl;

    // Plotfile global infos
    Print() << "Initialising global infos!" << std::endl;
    finestLevel = std::min(finestLevel,amrData.FinestLevel());
    int Nlev = finestLevel + 1;
    Print() << "finestLevel = " << finestLevel << "!" << std::endl;
    const Vector<std::string>& plotVarNames = amrData.PlotVarNames();
    RealBox rb(&(amrData.ProbLo()[0]), 
               &(amrData.ProbHi()[0]));

    // Gradient variable
    Print() << "Looking for gradient variable!" << std::endl;
    int idC = -1;
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == gradVar) idC = i;
    }
    if (idC<0) {
      Print() << "Cannot find " << gradVar << " data in pltfile \n";
    }

    // Auxiliary variables
    Print() << "Initialising Auxiliary variables!" << std::endl;
    nAuxVar = pp.countval("Aux_Variables");
    Vector<std::string> AuxVar(nAuxVar);
    for(int ivar = 0; ivar < nAuxVar; ++ivar) { 
         pp.get("Aux_Variables", AuxVar[ivar],ivar);
    }

    // ---------------------------------------------------------------------
    // Variables index management
    // ---------------------------------------------------------------------
    Print() << "Index Management!" << std::endl;
    const int idCst = 0;
    int nCompIn = idCst + 1;
    Vector<std::string> inVarNames(nCompIn);
    inVarNames[idCst] = plotVarNames[idC];
    Print() << "invarNames = " << inVarNames[0] <<std::endl;

    if (nAuxVar>0)
    {
        inVarNames.resize(nCompIn+nAuxVar);
        for (int ivar=0; ivar<nAuxVar; ++ivar) {
            if ( amrData.StateNumber(AuxVar[ivar]) < 0 ) {
               amrex::Abort("Unknown auxiliary variable name: "+AuxVar[ivar]);
            }
            inVarNames[nCompIn] = AuxVar[ivar];
            nCompIn ++;
        } 
    }

    Vector<int> destFillComps(nCompIn);
    for (int i=0; i<nCompIn; ++i) {
      destFillComps[i] = i;
    }
    Print() << "Size of destFillComps = " << destFillComps.size() << std::endl;;

    const int idGr = nCompIn;
    const int nCompOut = idGr + AMREX_SPACEDIM +1 ; // 1 component stores the ||gradT||

    // Check symmetry/periodicity in given coordinate direction
    Vector<int> sym_dir(AMREX_SPACEDIM,0);
    pp.queryarr("sym_dir",sym_dir,0,AMREX_SPACEDIM);  

    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Print() << is_per[idim] << " ";
    }
    Print() << "\n";


    int coord = 0;

    // ---------------------------------------------------------------------
    // Let's start the real work
    // ---------------------------------------------------------------------
    Print() << "Start of gradient computation!" << std::endl;
    Vector<MultiFab> state(Nlev);
    Vector<Geometry> geoms(Nlev);
    Vector<BoxArray> grids(Nlev);
    Vector<DistributionMapping> dmap(Nlev);
    const int nGrow = 1;

    // Read data on all the levels
    Print() << "Reading data on all levels!" << std::endl;
    const Real strt_time_io = ParallelDescriptor::second();
    for (int lev=0; lev<Nlev; ++lev) {
      const BoxArray ba = amrData.boxArray(lev);
      grids[lev] = ba;
      dmap[lev] = DistributionMapping(ba);
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
      state[lev].define(grids[lev], dmap[lev], inVarNames.size(), nGrow);
      state[lev].setVal(0.0);  // Reduces risk of nan or inf caes for fully covered cells if EB is used

      Print() << "Reading data for level: " << lev << std::endl;
      amrData.FillVar(state[lev], lev, inVarNames, destFillComps);
    }
    const Real end_time_io = ParallelDescriptor::second();
    Real io_time = end_time_io - strt_time_io;
    Print() << "Duration: " << io_time << " ms" << std::endl;  
    
    // Alternative load of data 
    PlotFileData pf(infile);
    const Real pf_strt_time_io = ParallelDescriptor::second();
    Vector<MultiFab> pfdata(Nlev);
    Vector<Geometry> pf_geoms(Nlev);
    Vector<BoxArray> pf_grids(Nlev);
    Vector<DistributionMapping> pf_dmap(Nlev);

    for (int lev=0; lev<Nlev; ++lev) {
      // Load data into Vector<Vector<MultiFab>> 
      const BoxArray pf_ba = pf.boxArray(lev);
      pf_grids[lev] = pf_ba;
      pf_dmap[lev] = pf.DistributionMap(lev);
      pf_geoms[lev] = Geometry(pf.probDomain(lev),&rb,coord,&(is_per[0]));
      pfdata[lev].define(pf_grids[lev], pf_dmap[lev], nCompOut, nGrow);
      pfdata[lev].setVal(0.0);
      for (int n=0; n< inVarNames.size(); ++n) {
        const MultiFab& src = pf.get(lev, inVarNames[n]);
        MultiFab::Copy(pfdata[lev], src, 0, n, 1, 0);
        Print() << "inVarNames " << inVarNames[n] << std::endl;
      }
      Print() << "...done reading the plotfile data at level " << lev << "..." << std::endl;
    }
  
    const Real pf_end_time_io = ParallelDescriptor::second();
    Real pf_io_time = pf_end_time_io - pf_strt_time_io;
    Print() << "Duration: " << pf_io_time << " ms" << std::endl; 

  } // amrex::Initialize



  amrex::Finalize();
  return 0;
}
