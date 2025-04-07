#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
// #include <AMReX_BLFort.H>
// #include <AMReX_MLMG.H>
#include <AMReX_MLPoisson.H>
// #include <AMReX_MLABecLaplacian.H>

// //EB
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_MLEBABecLap.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_MLMG.H>

using namespace amrex;

static
void 
print_usage (int,
             char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile=<plotfilename> \n\tOptions:\n\tis_per=<L M N> gradVar=<name>\n";
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
  amrex::Initialize(argc,argv);
  {
    if (argc < 2) {
      print_usage(argc,argv);
    }

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

    if (pp.contains("help")) {
      print_usage(argc,argv);
    }

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
    for (int lev=0; lev<Nlev; ++lev) {
      const BoxArray ba = amrData.boxArray(lev);
      grids[lev] = ba;
      dmap[lev] = DistributionMapping(ba);
      geoms[lev] = Geometry(amrData.ProbDomain()[lev],&rb,coord,&(is_per[0]));
      state[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      state[lev].setVal(0.0);  // Reduces risk of nan or inf caes for fully covered cells if EB is used

      Print() << "Reading data for level: " << lev << std::endl;
      amrData.FillVar(state[lev], lev, inVarNames, destFillComps);
      Print() << "Fill Boundary!" << std::endl;
      state[lev].FillBoundary(idCst,1,geoms[lev].periodicity());
    }



//------------------------------------------------------------------------------------------
// EB SECTION (Insert your own geometry + necessary parsing )
//------------------------------------------------------------------------------------------
#ifdef AMREX_USE_EB

Print() << "Start building EB!" << std::endl;
BL_PROFILE("PeleLMeX::makeEBGeometry()");
// Building EB (Needs to be changed to individual EB)

  //ParmParse pp("prob");
  Real prechamber_side_x, prechamber_side_y,cylinderRadius, cylinderLength;

  // Parse Information to EB
  pp.get("cylinder_radius",cylinderRadius);
  pp.get("cylinder_length",cylinderLength);
  pp.get("prechamber_side_x",prechamber_side_x);
  pp.get("prechamber_side_y",prechamber_side_y);

  //pp.get("required_coarsening_level",required_coarsening_level)
  //pp.get("max_coarsening_level",max_coarsening_level)

  int max_coarsening_level = 100;
  int required_coarsening_level = (geoms.size()) - 1;
  int eb_ref_level = finestLevel; 

  int ngrow = 4;
  bool build_coarse_level_by_coarsening = true;
  bool extend_domain_face = EB2::ExtendDomainFace();
  int num_coarsen_opt;
  pp.get("eb2.num_coarsen_opt",num_coarsen_opt);

  Print() << "geoms-size():"<< geoms.size() << std::endl;


  RealArray prechamberLo,prechamberHi,cylinderCentre,chamberLo,chamberHi;

  // Calculate Positions of shapes
  AMREX_D_TERM(prechamberLo[0] = geoms[eb_ref_level].ProbLo(0);,prechamberLo[1] = -prechamber_side_y/2.0;,prechamberLo[2] = -prechamber_side_y/2.0;); 
  Print() << prechamberLo << std::endl; // Just debugging or function?
  AMREX_D_TERM(prechamberHi[0] = geoms[eb_ref_level].ProbLo(0)+prechamber_side_x;,prechamberHi[1] = prechamber_side_y/2.0;,prechamberHi[2] = prechamber_side_y/2.0;);
  Print() << prechamberHi << std::endl;
  AMREX_D_TERM(cylinderCentre[0] = prechamberHi[0]+cylinderLength/2.0;,cylinderCentre[1] = 0.0;,cylinderCentre[2] = 0.0;);
  Print() << cylinderCentre << std::endl;
  AMREX_D_TERM(chamberLo[0] = geoms[eb_ref_level].ProbLo(0) + prechamber_side_x + cylinderLength;,chamberLo[1] = geoms[eb_ref_level].ProbLo(1);,chamberLo[2] = geoms[eb_ref_level].ProbLo(2););
  AMREX_D_TERM(chamberHi[0] = geoms[eb_ref_level].ProbHi(0);,chamberHi[1] = geoms[eb_ref_level].ProbHi(1);,chamberHi[2] = geoms[eb_ref_level].ProbHi(2););


  // Initialisation of the different shapes see https://amrex-codes.github.io/amrex/doxygen/index.html
  EB2::BoxIF prechamber(prechamberLo, prechamberHi,false);
  EB2::CylinderIF cylinder(cylinderRadius,cylinderLength,0,cylinderCentre,false);
  EB2::BoxIF chamber(chamberLo,chamberHi,false);

  auto domaininv = EB2::makeUnion(prechamber,cylinder,chamber); // Combines the shapes to one domain (automatic detects typ of variable)
  auto domain = EB2::makeComplement(domaininv); // Complement of an object. E.g. a sphere with fluid on outside becomes a sphere with fluid inside.
  
  // Build your geometry shop using EB2::makeShop
  auto gshop = EB2::makeShop(domain); 

  // Build geom using EB2::Build
  EB2::Build(gshop, geoms[eb_ref_level], required_coarsening_level, max_coarsening_level);
  Print() << "Building process finished!" << std::endl;



//-------------------------------------------------------------------------------------------------

  Print() << "Setting up EBFactory!" << std::endl;
  Vector<std::unique_ptr<EBFArrayBoxFactory>> eb_factory(Nlev);
  for (int lev = 0; lev < Nlev; ++lev) {
      const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
      const EB2::Level& eb_level = eb_is.getLevel(geoms[lev]);
      //eb_factory[lev] = makeEBFabFactory(geoms[lev], grids[lev], dmap[lev], {nGrow}, EBSupport::full);
      eb_factory[lev] = std::make_unique<EBFArrayBoxFactory>
          (eb_level, geoms[lev], grids[lev], dmap[lev], Vector<int>{2,2,2}, EBSupport::full);
  }

#endif


// Solver Section
#ifdef AMREX_USE_EB
    Print() << "Settings for Solver!" << std::endl;
    LPInfo info_apply;
    info_apply.setMaxCoarseningLevel(0);  
    MLEBABecLap poisson_eb(geoms, grids, dmap, info_apply, amrex::GetVecOfConstPtrs(eb_factory)); 
    poisson_eb.setVerbose(4);  

    poisson_eb.setMaxOrder(4);
    
    //Poisson like solver able to handle EB'S
    poisson_eb.setScalars(0.0, 1.0);   
    for (int lev = 0; lev <= finestLevel; ++lev) {
      poisson_eb.setBCoeffs(lev, -1.0);
    }

    Print() << "Setting Boundary Conditions" << std::endl;
    std::array<LinOpBCType, AMREX_SPACEDIM> lo_bc;
    std::array<LinOpBCType, AMREX_SPACEDIM> hi_bc;
    for (int idim = 0; idim< AMREX_SPACEDIM; idim++){
       if (is_per[idim] == 1) {
          lo_bc[idim] = hi_bc[idim] = LinOpBCType::Periodic;
       } else {
          if (sym_dir[idim] == 1) {
             lo_bc[idim] = hi_bc[idim] = LinOpBCType::reflect_odd;
          } else {
             lo_bc[idim] = hi_bc[idim] = LinOpBCType::Neumann;
          }
       }
    }
    
    poisson_eb.setDomainBC(lo_bc, hi_bc);



    // Need to apply the operator to ensure CF consistency with composite solve
    int nGrowGrad = 0;                   // No need for ghost face on gradient
    Vector<Array<MultiFab,AMREX_SPACEDIM> > grad(Nlev);
    Vector<std::unique_ptr<MultiFab>> phi;
    Vector<MultiFab> laps;
    for (int lev = 0; lev < Nlev; ++lev) {
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
         const auto& ba = grids[lev];
            grad[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                                dmap[lev], 1, nGrowGrad, MFInfo(), *eb_factory[lev]);

      }    
      phi.push_back(std::make_unique<MultiFab> (state[lev],amrex::make_alias,idCst,1));
      poisson_eb.setLevelBC(lev, phi[lev].get());
      laps.emplace_back(grids[lev], dmap[lev], 1, 1);
    }

    Print() << "Getting Fluxes!" << std::endl;
    MLMG mlmg_eb(poisson_eb);
    mlmg_eb.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
    mlmg_eb.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCentroid);
    //mlmg_eb.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);
#else

//------------------------------------------------------------------------------------------------
// Alternative if USE_EB = FALSE
//------------------------------------------------------------------------------------------------

    // Get face-centered gradients from MLMG 
    Print() << "Old Initialisation of MLMG!" << std::endl;
    LPInfo info;
    // info.setAgglomeration(1);
    // info.setConsolidation(1);
    info.setAgglomeration(true);
    info.setConsolidation(true);
    info.setMaxCoarseningLevel(0);
    info.setMetricTerm(false);
    info.setMaxCoarseningLevel(0);
    MLPoisson poisson({geoms}, {grids}, {dmap}, info);


    Print() << "Setting Boundary Conditions" << std::endl;
    std::array<LinOpBCType, AMREX_SPACEDIM> lo_bc;
    std::array<LinOpBCType, AMREX_SPACEDIM> hi_bc;
    for (int idim = 0; idim< AMREX_SPACEDIM; idim++){
       if (is_per[idim] == 1) {
          lo_bc[idim] = hi_bc[idim] = LinOpBCType::Periodic;
       } else {
          if (sym_dir[idim] == 1) {
             lo_bc[idim] = hi_bc[idim] = LinOpBCType::reflect_odd;
          } else {
             lo_bc[idim] = hi_bc[idim] = LinOpBCType::Neumann;
          }
       }
    }
    poisson.setDomainBC(lo_bc, hi_bc);
    poisson.setLevelBC(0, nullptr);
    Print() << "BC Set!" << std::endl;

    // Need to apply the operator to ensure CF consistency with composite solve
    int nGrowGrad = 0;                   // No need for ghost face on gradient
    Vector<Array<MultiFab,AMREX_SPACEDIM> > grad(Nlev);
    Vector<std::unique_ptr<MultiFab>> phi;
    Vector<MultiFab> laps;
    for (int lev = 0; lev < Nlev; ++lev) {
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
         const auto& ba = grids[lev];
         grad[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                                dmap[lev], 1, nGrowGrad);
      }    
      phi.push_back(std::make_unique<MultiFab> (state[lev],amrex::make_alias,idCst,1));
      poisson.setLevelBC(lev, phi[lev].get());
      laps.emplace_back(grids[lev], dmap[lev], 1, 1);


    }

    Print() << "Getting Fluxes!" << std::endl;
    MLMG mlmg(poisson);
    mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
    mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);


#endif

  Print() << "Starting gradalias etc" << std::endl;
    for (int lev = 0; lev < Nlev; ++lev) {
        // Convert to cell avg gradient
        MultiFab gradAlias(state[lev], amrex::make_alias, idGr, AMREX_SPACEDIM);
        EB_average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));  // if grad[lev] isn't build with ebfactory the function is equal to average_face_to_cellcenter(...);

        #ifdef AMREX_USE_EB
          // Not needed when using MLEBABecLap
        #else
          gradAlias.mult(-1.0);
        #endif
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(state[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {    
           const Box& bx = mfi.tilebox();
           auto const& grad_a   = gradAlias.const_array(mfi);
           auto const& gradMag  = state[lev].array(mfi,idGr+AMREX_SPACEDIM);
           amrex::ParallelFor(bx, [=]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {    
              gradMag(i,j,k) = std::sqrt(AMREX_D_TERM(  grad_a(i,j,k,0) * grad_a(i,j,k,0),
                                                      + grad_a(i,j,k,1) * grad_a(i,j,k,1),
                                                      + grad_a(i,j,k,2) * grad_a(i,j,k,2)));
           });  
        } 
    }

    // ---------------------------------------------------------------------
    // Write the results
    // ---------------------------------------------------------------------
    Vector<std::string> nnames(nCompOut);
    for (int i=0; i<nCompIn; ++i) {
      nnames[i] = inVarNames[i];
    }
    nnames[idGr+0] = gradVar + "_gx";
    nnames[idGr+1] = gradVar + "_gy";
#if AMREX_SPACEDIM==3
    nnames[idGr+2] = gradVar + "_gz";
#endif  
    nnames[idGr+AMREX_SPACEDIM] = "||grad"+ gradVar+ "||";
    std::string outfile(getFileRoot(infile) + "_gt"); pp.query("outfile",outfile);
    //std::string outfile(path + "_gt"); pp.query("outfile",outfile);

    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    //amrex::WriteMultiLevelPlotfile("volFrac_output", Nlev, GetVecOfConstPtrs(volFracMF), varnames, geoms, 0.0, isteps, refRatios);
    VisMF::SetNOutFiles(4);
    Real time = amrData.Time();
    amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(state), nnames,
                                   geoms, time, isteps, refRatios);
                 
  }


  amrex::Finalize();
  return 0;
}
