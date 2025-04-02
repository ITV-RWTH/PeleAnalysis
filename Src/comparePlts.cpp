// takes two hopefully identical plotfiles and tells you if you fucked up 

#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

static void print_usage (int, char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infiles=<reference and change> outfile=<> vars=<>" ;
  exit(1);
}

int main (int argc, char* argv[])
{
  amrex::Initialize(argc,argv);

  if (argc < 2)
    print_usage(argc,argv);
  
  ParmParse pp;

  // get infile names and count
  int nfiles(pp.countval("infiles"));
  Vector<std::string> infiles; infiles.resize(nfiles);
  pp.getarr("infiles",infiles);

  if (nfiles != 2){
    amrex::Abort("Tool is only designed for two infiles, please adjust");
  }
  // get and count variables to copy
  int nvars(pp.countval("vars"));
  int id_comp_last_a  = 0;
  int id_comp_last_b  = 0;
  Vector<std::string> vars; vars.resize(nvars);
  pp.getarr("vars",vars);  
  Vector<std::string> new_vars = vars; 
  Vector<std::string> names;
  // outfile name
  std::string outfile;
  pp.get("outfile",outfile);
  
  DataServices::SetBatchMode();
  Amrvis::FileType fileType(Amrvis::NEWPLT);

//setting up for reading pltfiles 
  DataServices dataServices0(infiles[0], fileType); 
  DataServices dataServices1(infiles[1], fileType); 

  if (!dataServices0.AmrDataOk() || !dataServices1.AmrDataOk())
      DataServices::Dispatch(DataServices::ExitRequest, NULL);

  AmrData& amrData0 = dataServices0.AmrDataRef();
  AmrData& amrData1 = dataServices1.AmrDataRef();
  
  // getting the finest level (or whatever user sets)
  // aborting if finest level dont match
  int finestLevel = -1;
  finestLevel = amrData0.FinestLevel();
  pp.query("finestLevel",finestLevel);
  // user changed finest lev
  if (
     (finestLevel > amrData1.FinestLevel()) || 
     (finestLevel > amrData0.FinestLevel())
     )
      amrex::Abort("Warning: Declared finest level does not exist");
  // user did not change finest lev  
  if (
     (amrData0.FinestLevel() != amrData1.FinestLevel()) && 
     (amrData0.FinestLevel() != finestLevel)
     )
      amrex::Abort("Warning: The Finest Level of the infiles do not match");

  int Nlev = finestLevel + 1;

  // setting up periodicity
  Vector<int> is_per(AMREX_SPACEDIM,1);
  pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
  Print() << "Periodicity assumed for this case: ";
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    Print() << is_per[idim] << " ";
  }
  Print() << "\n";
  
  // setting up pltfile boxes
  RealBox rb(&(amrData0.ProbLo()[0]), 
	     &(amrData0.ProbHi()[0]));
  Vector<MultiFab*> fileData0(Nlev);
  Vector<MultiFab*> fileData1(Nlev);
  Vector<MultiFab*> allOnes(Nlev);
  Vector<Geometry> geoms(Nlev);
  int coord = 0;
  // create the geometry of MultiFabs, the "boxArray" Parent does'nt matter
  for (int lev=0; lev<Nlev; ++lev)
  {
      const DistributionMapping dm(amrData0.boxArray(lev));
      geoms[lev] = Geometry(amrData0.ProbDomain()[lev],&rb,coord,&(is_per[0]));
      fileData0[lev] = new MultiFab(amrData0.boxArray(lev),dm,nvars,0);
      fileData1[lev] = new MultiFab(amrData1.boxArray(lev),dm,nvars,0);
      allOnes[lev] = new MultiFab(amrData1.boxArray(lev),dm,nvars,0);
  }

  // data structure to create identical var tables
  const Vector<std::string>& plot0_VarNames = amrData0.PlotVarNames();
  const Vector<std::string>& plot1_VarNames = amrData1.PlotVarNames();

  Vector<int> idcomp0;
  Vector<int> idcomp1;

  // loop through all available vars
  for (int i = 0; i < nvars; i++){
    const std::string& var = vars[i];
    int idx0 = -1; 
    int idx1 = -1;

    // loop through first infile to find var number i
    for (int j = 0; j < plot0_VarNames.size(); j++){
      // if var was found, save the position in idx
      if (plot0_VarNames[j] == var){
        idx0 = j;
        break;
      } 
    }

    // loop through second infile to finde var number i
    for (int k = 0; k < plot1_VarNames.size(), k++){
      // if var was found, save the position in idx
      if (plot1_VarNames[k] == var){
        idx1 = k;
        break;
      }
    }

    // if first idx is unchanged: var was not found, abort
    if (idx0 < 0){
      Print() << "Error: Variable " << var << " could not be located in first infile \n";
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }

    // if second idx is unchanged: var was not found, abort
    if (idx1 < 0){
      Print() << "Error: Variable " << var << " could not be located in second infile \n";
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }

    // save the position for each file, in which the var lays
    idcomp0.push_back(idx0);
    idcomp1.push_back(idx1);
    // now we have two tables with indexes ponting to the right var for each infile
  } 
 
  // fill in the Multifabs with the infile Data
 for (int lev = 0; lev < Nlev; lev++){
  for (int i = 0; i < idcomp0.size(), i++){
    fileData0[lev]->ParallelCopy(amrData0.GetGrids(lev,idcomp0[0]),0,i,1);
    fileData1[lev]->ParallelCopy(amrData1.GetGrids(lev,idcomp1[0]),0,i,1);
    allOnes[lev]->setVal(1.0);
  }
 }

// creating an array containing the var names
 for (int i = 0; i < idcomp0.size(), i++){
  names.push_back(plot0_VarNames[idcomp0[i]] + "_rel_diff");
 }
  
  // calculate the difference: hopefully
 for (int lev = 0; lev < Nlev; ++lev) {
   for (int i = 0; i < nvars; i++){
     fileData0[lev]->Devide(fileData0[lev], fileData1[lev], i, i, 1);
     allOnes[lev]->minus(fileData0[lev], i, 1);
   }
 }
  
  // write pltfile
  Vector<int> isteps(Nlev, 0);
  Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
  amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(allOnes), names, geoms, 0.0, isteps, refRatios);
  
  amrex::Finalize();
  return 0;
}
