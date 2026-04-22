#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "Calculates a progress variable from a pltFile as C = (Sum(specNames) - unburntVal)/(burntVal - unburntVal)\n";
  std::cerr << "usage:\n";
  std::cerr << argv[0] << "inputs infile=<i> speciesNames=<s> unburntVal=<u> burntVal=<b> [options] \n\tOptions:\n";
  std::cerr << "\t     infile=<i> where <i> is a pltfile\n";
  std::cerr << "\t     speciesNames=<s> where <s> are all the species to be included in the progress variable\n";
  std::cerr << "\t     unburntVal=float This is the sum of the included species in the unburnt unburntVal[DEF->0]\n";
  std::cerr << "\t     burntVal=float This is the value of included species in the burnt comp[DEF->1]\n";
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

    // ---------------------------------------------------------------------
    // ParmParse
    // ---------------------------------------------------------------------

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc,argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(true);

    std::string plotFileName; pp.get("infile",plotFileName);

    Vector<std::string> species;
    int nSpec = pp.countval("speciesNames"); 
    species.resize(nSpec);
    pp.getarr("speciesNames",species,0,nSpec);

    // Default values
    Real unburnt = 0; pp.get("unburntVal", unburnt);
    Real burnt = 1; pp.get("burntVal", burnt);


    // DataServices (Reads infile into amrData)
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if( ! dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
      // ^^^ this calls ParallelDescriptor::EndParallel() and exit()
    }
    AmrData& amrData = dataServices.AmrDataRef();

    // Get Plt File Info like number of levels, periodicity, 
    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Print() << is_per[idim] << " ";
    }

    RealBox rb(&(amrData.ProbLo()[0]),&(amrData.ProbHi()[0]));

    int nCompIn = amrData.NComp();
   
    Vector<std::string> inNames = amrData.PlotVarNames();

    // Get ids of relevant species
    Vector<int> idY(nSpec, -1);
    for (int j=0; j<nSpec; ++j)
    {
      for (int i=0; i<inNames.size(); ++i)
      {
        if (inNames[i] == species[j]) idY[j] = i;
      }
      if (idY[j]<0) {
        amrex::Abort("Species " + species[j] + " not found in plotfile");
      }
    }

    Vector<MultiFab> outdata(Nlev); // same levels as indata
    Vector<Geometry> geoms(Nlev);
    int nGrow = 0;
    Vector<std::string> outNames = {"specSum", "progVar"};
 
    // Read all comps of indata TODO Can you do this better? 
    Vector<int> destFillComps(nCompIn);
    for (int i=0; i<nCompIn; ++i) destFillComps[i] = i;

    for (int lev = 0; lev < Nlev; ++lev) {
      
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);

      outdata[lev] = MultiFab(ba,dm,outNames.size(),nGrow);
      MultiFab indata(ba,dm,nCompIn,nGrow); //important: indata needs to have space for all variables from plt files! -> initialisation with nCompIn

      int coord = 0;
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb, coord, &(is_per[0]));      
 
      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata,lev,inNames,destFillComps);

      Print() << "Data has been read for level " << lev << std::endl;
      
      // Prepare efficient kernel work
      amrex::Gpu::DeviceVector<int> d_idY(idY.size());
      // 2. Copy host to device
      amrex::Gpu::copy(amrex::Gpu::hostToDevice,idY.begin(), idY.end(),d_idY.begin());

      // 3. Get raw device pointer
      int const* idY_d = d_idY.data();
      Real denom = burnt - unburnt;
      
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(indata,TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
	  const Box& bx = mfi.tilebox();
	  auto const& out_a = outdata[lev].array(mfi);
	  auto const& in_a = indata.array(mfi);
	  amrex::ParallelFor(bx, [=]
			     AMREX_GPU_DEVICE (int i, int j, int k) noexcept
				 {
				   // calculate the sum of species
				   Real sum = 0.0_rt;
				   for (int s=0; s<nSpec; ++s)
 				   {
				     sum += in_a(i,j,k,idY_d[s]);
				   }
				   // calculate the progress variable
				   out_a(i,j,k, 0) = sum;
				   out_a(i,j,k,1) = (out_a(i,j,k,0) - unburnt)/denom;
				 });
	}

    }
    
    std::string outfile(getFileRoot(plotFileName) + "_prog");
    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(outdata), outNames,
                                   geoms, 0.0, isteps, refRatios);

  }
  Finalize();
  return 0;
}
