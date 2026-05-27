#include <iostream>
#include <string>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
    std::cerr << "usage:\n";
    std::cerr << argv[0] << " infile=<filename>\n"
	      << "   vars= var1 var2 var3 : Variable list\n"
	      << "   useminmax%i= min max : Use min/max for var %i (count from 1)\n"
              << "   infile= plt1 plt2    : List of plot files" << std::endl;
    exit(1);
}

static
std::string
ProtectSlashes (const std::string str)
{
    std::string s = str;

    std::vector<int> where;

    const char* cstr = s.c_str();

    for (int i = 0; i < s.size(); i++)
        if (cstr[i] == '/')
            where.push_back(i);

    for (int i = where.size() - 1; i >= 0; i--)
        s.replace(where[i], 1, "_");

    return s;
}

int
main (int   argc,
      char* argv[])
{
    amrex::Initialize(argc,argv);

    if (argc < 2)
        print_usage(argc,argv);

    ParmParse pp;

    if (pp.contains("help"))
        print_usage(argc,argv);

    bool verbose = false;
    if (ParallelDescriptor::IOProcessor())
	verbose = true;

    if (pp.contains("verbose"))
        AmrData::SetVerbose(true);

    // Get plot file
    int nPlotFiles(pp.countval("infile"));
    if(nPlotFiles <= 0) {
	std::cerr << "Bad nPlotFiles:  " << nPlotFiles << std::endl;
	std::cerr << "Exiting." << std::endl;
	DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    if (verbose)
	std::cout << "Processing " << nPlotFiles << " plotfiles..." << std::endl;

    // Make an array of srings containing paths of input plot files
    Vector<std::string> plotFileNames(nPlotFiles);
    for(int iPlot = 0; iPlot < nPlotFiles; ++iPlot) {
	pp.get("infile", plotFileNames[iPlot], iPlot);
	if (verbose)
	    std::cout << "   " << plotFileNames[iPlot] << std::endl;
    }

    // Get finest level argument
    int inFinestLevel(-1);
    pp.query("finestLevel",inFinestLevel);

    // Go to pink & white past maximum
    int goPastMax(1);
    pp.query("goPastMax",goPastMax);

    // Variables to load
    int nVars(pp.countval("vars"));
    Vector<std::string> whichVar(nVars);
    if (verbose)
      std::cout << "Variable list:" << std::endl;
    for(int v=0; v<nVars; v++) {
	pp.get("vars", whichVar[v], v);
	if (verbose)
	  std::cout << "   " << whichVar[v] << std::endl;
    }
    // Copy the names of the variable for output filenames, replacing dodgy characters
    Vector<std::string> whichVarOut(nVars);
    for (int iVar=0; iVar<nVars; iVar++) {
	whichVarOut[iVar] = ProtectSlashes(whichVar[iVar]);
	if (verbose)
	  std::cout << whichVar[iVar] << " -> " << whichVarOut[iVar] << std::endl;
    }
    
    // More random initialisation stuff
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    // DataPtrArray for each plot
    Vector<AmrData *> amrDataPtrArray(nPlotFiles);
    Vector<Real>      time(nPlotFiles);

    // Loop over files
    for (int iPlot=0; iPlot<nPlotFiles; iPlot++) {

	// Open file and get an amrData pointer
	std::string infile = plotFileNames[iPlot];
	if (verbose)
	  std::cout << "\nOpening " << infile << "..." << std::endl;
	DataServices dataServices(infile, fileType);
	if (!dataServices.AmrDataOk())
	    DataServices::Dispatch(DataServices::ExitRequest, NULL);
	
	amrDataPtrArray[iPlot] = &(dataServices.AmrDataRef());

	if (verbose)
	    std::cout << "   ...done." << std::endl;
    
	// Check the names of the variables are present in the plotfile
	Vector<int> destFills(nVars);
	for (int v=0; v<nVars; v++) {
	    destFills[v] = v;
	    if (amrDataPtrArray[iPlot]->StateNumber(whichVar[v])<0) {
		std::string message="Bad variable name ("+whichVar[v]+")";
		amrex::Error(message.c_str());
	    }
	}
    
	// Get finest level
	int finestLevel = amrDataPtrArray[iPlot]->FinestLevel();    
	if (inFinestLevel>-1 && inFinestLevel<finestLevel) {
	    finestLevel = inFinestLevel;
	    if (verbose)
	      std::cout << "Finest level: " << finestLevel << std::endl;
	}
	int nLevels = finestLevel+1;

	// Get probDomain
	Vector<Box> probDomain = amrDataPtrArray[iPlot]->ProbDomain();

	// Get min/max for each component
	Vector<Real> vMin(nVars);
	Vector<Real> vMax(nVars);
	for (int iVar=0; iVar<nVars; iVar++) {
	    vMin[iVar]=1e100;
	    vMax[iVar]=-vMin[iVar];
	    for (int iLevel=0; iLevel<nLevels; iLevel++) {
		Real min,max;
		amrDataPtrArray[iPlot]->MinMax(probDomain[iLevel], whichVar[iVar], iLevel, min, max);
		if (vMin[iVar]>min) vMin[iVar]=min;
		if (vMax[iVar]<max) vMax[iVar]=max;
	    }
	    std::cout << "Min/Max from file var " << iVar << ": " << vMin[iVar] << " / " << vMax[iVar] << std::endl;
	}

	for (int iVar=0; iVar<nVars; iVar++) {
	    char argName[12];
	    sprintf(argName,"useminmax%i",iVar+1);
	    int nMinMax = pp.countval(argName);
	    if (nMinMax>0) {
		if (nMinMax != 2) {
		    amrex::Abort("Need to specify 2 values for useMinMax");
		} else {
		    pp.get(argName, vMin[iVar], 0);
		    pp.get(argName, vMax[iVar], 1);
		    if (verbose)
			std::cout << "Min/Max from command line var" << iVar+1 << " (" << whichVar[iVar] << ") using min/max: " << vMin[iVar] << " / " << vMax[iVar] << std::endl;
		}
	    }
	}
	
	// Make the fab
	if (verbose)
	    std::cout << "Making fab..." << std::endl;
	Box tempBox(probDomain[finestLevel]);
#if (AMREX_SPACEDIM==3)
	int xslice(-1);
	pp.query("xslice",xslice);
	int yslice(-1);
	pp.query("yslice",yslice);
	int zslice(-1);
	pp.query("zslice",zslice);
	int islice(-1);
	int sliceDir(-1);
	std::string sliceStr;
	if (xslice>-1) {
	  islice=xslice;
	  sliceDir=Amrvis::XDIR;
	  sliceStr="X";
	}
	if (yslice>-1) {
	  if (islice==-1) {
	    islice=yslice;
	    sliceDir=Amrvis::YDIR;
	    sliceStr="Y";
	  } else {
	    amrex::Error("Can only specify one slice");
	  }
	}
	if (zslice>-1) {
	  if (islice==-1) {
	    islice=zslice;
	    sliceDir=Amrvis::ZDIR;
	    sliceStr="Z";
	  } else {
	    amrex::Error("Can only specify one slice");
	  }
	}
	// default to yslice=0
	if (sliceDir==-1) {
	  sliceDir = 1;
	  islice   = 0;
	  sliceStr="Y";
	}
	char isliceStr[8];
	sprintf(isliceStr,"%i",islice);
	if (verbose)
	  std::cout << "   " << sliceStr << "slice" << isliceStr << std::endl;	
	tempBox.setSmall(sliceDir,islice);
	tempBox.setBig(sliceDir,islice);
#endif
	BoxArray domainBoxArray(1);
	domainBoxArray.set(0,tempBox);
	MultiFab mf;
	DistributionMapping dm(domainBoxArray);
	mf.define(domainBoxArray, dm, nVars, 0);

	// Load data
	if (verbose)
	    std::cout << "Loading data..." << std::endl;
	
	amrDataPtrArray[iPlot]->FillVar(mf, finestLevel, whichVar, destFills);

	for (int n=0; n<nVars; n++)
	    amrDataPtrArray[iPlot]->FlushGrids(amrDataPtrArray[iPlot]->StateNumber(whichVar[n]));

	// Construct image buffer
	int ix, kx;
#if (AMREX_SPACEDIM==2)
	ix = probDomain[finestLevel].length(Amrvis::XDIR);
	kx = probDomain[finestLevel].length(Amrvis::YDIR);
#endif
#if (AMREX_SPACEDIM==3)
	switch (sliceDir) {
	case 0:
	  ix = probDomain[finestLevel].length(Amrvis::YDIR);
	  kx = probDomain[finestLevel].length(Amrvis::ZDIR);
	  break;
	case 1:
	  ix = probDomain[finestLevel].length(Amrvis::XDIR);
	  kx = probDomain[finestLevel].length(Amrvis::ZDIR);
	  break;
	case 2:
	  ix = probDomain[finestLevel].length(Amrvis::XDIR);
	  kx = probDomain[finestLevel].length(Amrvis::YDIR);
	  break;
	}
#endif
	int *buff = (int*)malloc(3*ix*kx*sizeof(int));
	if (verbose)
	    std::cout << "Constructing image buffer..." << std::endl;
	for(MFIter ntmfi(mf); ntmfi.isValid(); ++ntmfi) {
	  const FArrayBox &myFab = mf[ntmfi];
	  for (int iVar=0; iVar<nVars; iVar++) {
	    const Real *dat = myFab.dataPtr(iVar);
	    for (int k=0; k<kx; k++) {
	      for (int i=0; i<ix; i++) {
		int cell = k*ix+i;
		int bc   = ((kx-k-1)*ix+i)*3;
		Real val = dat[cell];
		Real colour = fmax(0.,fmin(1.5,(dat[cell]-vMin[iVar])/(vMax[iVar]-vMin[iVar])));
		if (colour<0.125) {
		  buff[bc]   = 0;                         
		  buff[bc+1] = 0;
		  buff[bc+2] = (int)((colour+0.125)*1020.);
		} else if (colour<0.375)  {
		  buff[bc]   = 0;
		  buff[bc+1] = (int)((colour-0.125)*1020.);
		  buff[bc+2] = 255;
		} else if (colour<0.625)  {
		  buff[bc]   = (int)((colour-0.375)*1020.);
		  buff[bc+1] = 255;
		  buff[bc+2] = (int)((0.625-colour)*1020.);
		} else if (colour<0.875)  {
		  buff[bc]   = 255;
		  buff[bc+1] = (int)((0.875-colour)*1020.);
		  buff[bc+2] = 0;
		} else if (colour<1.000)  {
		  buff[bc]   = (int)((1.125-colour)*1020.);
		  buff[bc+1] = 0;
		  buff[bc+2] = 0;
		} else if (goPastMax==1) {
		  if (colour<1.125)  {
		    buff[bc]   = (int)((colour-0.875)*1020.);
		    buff[bc+1] = 0;
		    buff[bc+2] = (int)((colour-1.000)*1020.);
		  } else if (colour<1.250) {
		    buff[bc]   = 255;
		    buff[bc+1] = 0;
		    buff[bc+2] = (int)((colour-1.000)*1020.);
		  } else if (colour<1.500)  {
		    buff[bc]   = 255;
		    buff[bc+1] = (int)((colour-1.250)*1020.);
		    buff[bc+2] = 255;
		  } else { // default if above 1.5 with goPastMax==1
		    buff[bc]   = 255;
		    buff[bc+1] = 255;
		    buff[bc+2] = 255;
		  }
		} else { // default if above 1 with goPastMax==0
		  buff[bc]   = 128;
		  buff[bc+1] = 0;
		  buff[bc+2] = 0;	      
		}
	      }
	    }
	    ParallelDescriptor::ReduceIntMax(buff,ParallelDescriptor::IOProcessorNumber());
	    if (ParallelDescriptor::IOProcessor()) {
	      unsigned char *buff_char=(unsigned char*)malloc(3*ix*kx*sizeof(char));
	      //convert to unsigned char
	      for (int idx = 0; idx < 3*ix*kx; idx++) {
		buff_char[idx] = (unsigned char)buff[idx];
	      } 
	      // Output file
	      if (verbose)
		std::cout << "Outputting file..." << std::endl;
	      
	      //std::string filename=plotFileNames[iPlot]+"/"+whichVar[iVar] // origional
	      std::string filename=plotFileNames[iPlot]+whichVar[iVar] // mod
#if (AMREX_SPACEDIM==3)
		+sliceStr+"slice"+isliceStr
#endif
		+".ppm";
	      FILE *file = fopen(filename.c_str(),"w");
	      fprintf(file,"P6\n%i %i\n255\n",ix,kx);
	      fwrite(buff_char,ix*kx*3,sizeof(unsigned char),file);
	      fclose(file);
	      
	      if (verbose)
		std::cout << "   ...done." << std::endl;
	    }
	  } // iVar
	}
	
    } // iPlot
    amrex::Finalize();
    return 0;
}








