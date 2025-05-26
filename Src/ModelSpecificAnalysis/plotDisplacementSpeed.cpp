#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_BCRec.H>
#include <AMReX_Interpolater.H>

#include <AMReX_MLMG.H>

#include <mechanism.H>
#include <PelePhysics.H>

#ifdef AMREX_USE_EB
#include <AMReX_MLEBABecLap.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H> 
#include <pelelmex_prob_parm.H>
#include <PeleLMeX_EBUserDefined.H>
#else
#include <AMReX_MLPoisson.H>
#include <AMReX_MLABecLaplacian.H>
#endif

using namespace amrex;

pele::physics::PeleParams<pele::physics::transport::TransParm<
  pele::physics::PhysicsType::eos_type,
  pele::physics::PhysicsType::transport_type>>
  trans_parms;

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

    // ---------------------------------------------------------------------
    // Set defaults input values
    // ---------------------------------------------------------------------
    std::string fuelName  = "H2";
    int verbose           = 0; 
    Real sL               = 1.0; // laminar flame speed
    int do_soret          = 0;   
    int do_wbar           = 1;
    int n_files           = 4; 



    // ---------------------------------------------------------------------
    // ParmParse
    // ---------------------------------------------------------------------

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc,argv);

    pp.query("verbose",verbose);
    
    pp.query("amr.n_files",n_files);  // Changes how many files the written pltfile contains

    std::string plotFileName; pp.get("infile",plotFileName);
    PlotFileData pf(plotFileName);

    trans_parms.initialize();

    int finestLevel = pf.finestLevel();
    pp.query("finestLevel",finestLevel);
    int Nlev = finestLevel + 1;

    pp.get("fuelName",fuelName);

    std::string progressVar = "Y("+fuelName+")";
    Print () << "ProgressVar: " << progressVar << std::endl;
    pp.query("sL",sL);
    ParmParse pptrans("transport");
    pptrans.query("use_soret",do_soret);
    pp.query("use_wbar",do_wbar);
    int idYin = -1;
    int idTin = -1;
    int idRin = -1;
    int idCin = -1;
    int idEin = -1;

    Vector<std::string> spec_names;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);

    int fuelID = -1;
    for (int i = 0; i < spec_names.size(); ++i) {
        if (spec_names[i] == fuelName) {
            fuelID = i;
            break;
        }
    }

    if (fuelID < 0) {
        amrex::Abort("Species '" + fuelName + "' not found in mechanism!");
    }


    const Vector<std::string>& plotVarNames = pf.varNames();
    const std::string spName= "Y(" + spec_names[0] + ")";
    const std::string TName = "temp";
    const std::string RName = "density";
    const std::string CName = "I_R(" + fuelName + ")";
    const std::string EName = "rhoh";
    
    for (int i=0; i<plotVarNames.size(); ++i)
    {
      if (plotVarNames[i] == spName) idYin = i;
      if (plotVarNames[i] == TName)  idTin = i;
      if (plotVarNames[i] == RName)  idRin = i;
      if (plotVarNames[i] == CName)  idCin = i;
      if (plotVarNames[i] == EName)  idEin = i;

    }
    if (idYin<0 || idTin<0 || idRin<0 || idCin<0)
      Abort("Cannot find required data in pltfile");

    const int nCompIn  = NUM_SPECIES+4; //all species, temperature, density, chemical source term for fuel
    const int nCompOut = 6;
    Vector<std::string> outNames(nCompOut);  
    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    const int idYlocal = 0; // Xs start here
    const int idClocal = NUM_SPECIES;
    const int idTlocal = NUM_SPECIES+1;   // T starts here
    const int idRlocal = NUM_SPECIES+2; // R starts here
    const int idElocal = NUM_SPECIES+3;
    for (int i=0; i<NUM_SPECIES; ++i) {
      destFillComps[i] = idYlocal + i;
      inNames[i] =  "Y(" + spec_names[i] + ")";
    }
    destFillComps[idClocal] = idClocal;
    destFillComps[idTlocal] = idTlocal;
    destFillComps[idRlocal] = idRlocal;
    destFillComps[idElocal] = idElocal; 
    inNames[idTlocal] = TName;
    inNames[idRlocal] = RName;
    inNames[idClocal] = CName;    // "I_R(" + fuelName + ")";
    inNames[idElocal] = EName;
    outNames[0] = "Sd("+fuelName+")_dY";
    outNames[1] = "Sd("+fuelName+")_dT";
    outNames[2] = "Sd("+fuelName+")_dW";
    outNames[3] = "Sd("+fuelName+")_C";
    outNames[4] = "Sd("+fuelName+")";
    outNames[5] = "Sd("+fuelName+")rhoh";

    int idProglocal=-1;
    for (int i = 0; i < nCompIn; i++) {
      if (progressVar == inNames[i]) {
	idProglocal=i;
	break;
      }
    }
    if (idProglocal == -1) {
      Abort("Progress variable not being loaded");	    
    }

    Vector<int> sym_dir(AMREX_SPACEDIM,0);
    pp.queryarr("sym_dir",sym_dir,0,AMREX_SPACEDIM);

    Vector<int> is_per(AMREX_SPACEDIM,1);
    pp.queryarr("is_per",is_per,0,AMREX_SPACEDIM);
    Print() << "Periodicity assumed for this case: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Print() << is_per[idim] << " ";
    }
    Print() << "\n";
    BCRec gradVarBC;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        gradVarBC.setLo(idim,BCType::foextrap);
        gradVarBC.setHi(idim,BCType::foextrap);
        if ( is_per[idim] ) {
            gradVarBC.setLo(idim, BCType::int_dir);
            gradVarBC.setHi(idim, BCType::int_dir);
        }
    }

    int coord = 0;

    amrex::RealBox real_box({AMREX_D_DECL(pf.probLo()[0], pf.probLo()[1], pf.probLo()[2])},
                            {AMREX_D_DECL(pf.probHi()[0], pf.probHi()[1], pf.probHi()[2])});
                      


    Vector<Geometry> geoms(Nlev);
    Vector<BoxArray> grids(Nlev);
    Vector<DistributionMapping> dmap(Nlev);
                                                                                               
    for (int lev=0; lev<Nlev; ++lev) {
      const BoxArray ba = pf.boxArray(lev);
      grids[lev] = ba;
      dmap[lev] = pf.DistributionMap(lev);
      geoms[lev] = Geometry(pf.probDomain(lev),&real_box,coord,&(is_per[0]));
    
    }

//------------------------------------------------------------------------------------------------------------------------
// Initialise MultiFabs and read data
//------------------------------------------------------------------------------------------------------------------------

    Vector<MultiFab> coeffMF(Nlev); //holds transport coefficients
    Vector<MultiFab> gradMF(Nlev); //holds all the gradients, need gradY, gradW and gradT
    Vector<MultiFab> grad2MF(Nlev); //holds all the second derivatives
    Vector<MultiFab> deriveMF(Nlev); //just holds mean molecular weight   
    Vector<MultiFab> outdata(Nlev);
    Vector<MultiFab> indata(Nlev);
    const int nGrow = 1;
    const int nCoeffs = 3; //standard, soret and molecular weight coeffs for H2
    const int nGrads = 3*AMREX_SPACEDIM; //gradFuel gradWbar gradT
    const int n2Grads = 3*AMREX_SPACEDIM*AMREX_SPACEDIM; //each diffusive flux has D dims
    const int nDerive = 1 + 3*AMREX_SPACEDIM; //Wbar and the 3 different diffusive fluxes 

    // Read data on all the levels                                                                                                   
    for (int lev=0; lev<Nlev; ++lev) {

      indata[lev].define(grids[lev], dmap[lev], nCompIn, nGrow);
      coeffMF[lev].define(grids[lev], dmap[lev], nCoeffs, nGrow);
      deriveMF[lev].define(grids[lev], dmap[lev], nDerive, nGrow);
      gradMF[lev].define(grids[lev],dmap[lev], nGrads, nGrow);
      grad2MF[lev].define(grids[lev],dmap[lev], n2Grads, nGrow);
      outdata[lev].define(grids[lev], dmap[lev], nCompOut, nGrow);
      
      // Get input state data
      if (verbose) Print() << "Reading data for level " << lev << "\n";
      indata[lev].setVal(0.0);
      for (int n=0; n< inNames.size(); ++n) {
          const MultiFab& src = pf.get(lev, inNames[n]);
          MultiFab::Copy(indata[lev], src, 0, n, 1, 0);
        }
      indata[lev].FillBoundary(0,nCompIn,geoms[lev].periodicity());


      // Setting remaining initial values to 0 
      coeffMF[lev].setVal(0.0);
      deriveMF[lev].setVal(0.0);
      gradMF[lev].setVal(0.0);
      grad2MF[lev].setVal(0.0);
      outdata[lev].setVal(0.0);

    } // lev-loop

//---------------------------------------------------------------------------
// Building EB
//---------------------------------------------------------------------------

#ifdef AMREX_USE_EB

    if (verbose) Print() << "Start building EB!" << std::endl;
    BL_PROFILE("PeleLMeX::makeEBGeometry()");

    int max_coarsening_level = 100;
    int req_coarsening_level = static_cast<int>(geoms.size()) - 1;

    // Read the geometry type and act accordingly
    ParmParse ppeb2("eb2");
    std::string geom_type;
    ppeb2.get("geom_type", geom_type);

    // At what level should the EB be generated ?
    // Default : max_level
    int max_lvl_eb = finestLevel;
    ppeb2.query("max_level_generation", max_lvl_eb);

    // Generate the EB data at max_lvl_eb
    if (geom_type == "UserDefined") {
        EBUserDefined(
        geoms[max_lvl_eb], req_coarsening_level, max_coarsening_level);
    } else {
        // If geom_type is not an AMReX recognized type, it'll crash.
        EB2::Build(
        geoms[max_lvl_eb], req_coarsening_level, max_coarsening_level);
    }

    // Setting up an eb_factory for the solver to use
    if (verbose) Print() << "Setting up EBFactory!" << std::endl;
    Vector<std::unique_ptr<EBFArrayBoxFactory>> eb_factory(Nlev);
    for (int lev = 0; lev < Nlev; ++lev) {
        const EB2::IndexSpace& eb_is = EB2::IndexSpace::top();
        const EB2::Level& eb_level = eb_is.getLevel(geoms[lev]);
        eb_factory[lev] = std::make_unique<EBFArrayBoxFactory>
            (eb_level, geoms[lev], grids[lev], dmap[lev], Vector<int>{2,2,2}, EBSupport::full);
    }

#endif // AMREX_USE_EB

//------------------------------------------------------------------------------------------------------------------------
// Compute Coefficients
//------------------------------------------------------------------------------------------------------------------------

    for (int lev=0; lev<Nlev; ++lev)
    {
      //First compute transport properties and specific enthalpy
      // Get the transport data pointer
      auto const* ltransparm = trans_parms.device_parm();
      
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(indata[lev], amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const Box& bx = mfi.tilebox();
        Array4<Real const> const& Y_a = indata[lev].const_array(mfi,idYlocal);
        Array4<Real const> const& T_a = indata[lev].const_array(mfi,idTlocal);
        Array4<Real const> const& rho_a = indata[lev].const_array(mfi,idRlocal);
        Array4<Real> const& I_R = indata[lev].array(mfi,idClocal); // I_R(fuelName)
	
        Array4<Real> const& w_a = deriveMF[lev].array(mfi,0);
        Array4<Real> const& spec_coeff_a = coeffMF[lev].array(mfi,0);
        Array4<Real> const& soret_coeff_a = coeffMF[lev].array(mfi,1);
        Array4<Real> const& molar_coeff_a = coeffMF[lev].array(mfi,2);

        // For EB
      #ifdef AMREX_USE_EB
        Array4<const Real> const&  volFracBox = eb_factory[lev]->getVolFrac()[mfi].array();
      #endif

        // Debugging
        Array4<Real> const& sd_a = outdata[lev].array(mfi,0);

	amrex::ParallelFor(bx, [=]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

        #ifdef AMREX_USE_EB
          Real volFrac = volFracBox(i,j,k,0);
        #else
          Real volFrac = 1; 
        #endif

          if(volFrac > 0) {

	  
            Real Yloc[NUM_SPECIES] = {0.0};
            Real Dcoeff[NUM_SPECIES] = {0.0};
            Real chi[NUM_SPECIES] = {0.0};
            Real imw[NUM_SPECIES] = {0.0};

          
            Real ci[NUM_SPECIES] = {0.0};   // species concentration (molecular)
            Real wdot[NUM_SPECIES] = {0.0}; // production rate
            
            for (int n=0; n<NUM_SPECIES; ++n) {
              Yloc[n] = Y_a(i,j,k,n);
                  }

            //Dummies unused due to falses
            Real lam,mu_dummy,xi_dummy,mmw;
            //get D, lam, and chi
            const bool get_xi = false;
            const bool get_mu = false;
            const bool get_lam = false;
            const bool get_Ddiag = true;
            const bool get_chi = (do_soret == 1);
            pele::physics::transport::SimpleTransport::transport(get_xi,get_mu,get_lam,get_Ddiag,get_chi,T_a(i,j,k),rho_a(i,j,k),Yloc,Dcoeff,chi,mu_dummy,xi_dummy, lam,ltransparm);

            //inverse molecular weights
            get_imw(imw);
            //mean molecular weight
            CKMMWY(Yloc,mmw);
            mmw *= 1.0e-3; //cgs->mks
            w_a(i,j,k) = mmw;
            
            for (int n=0; n< NUM_SPECIES; n++) {
              imw[n] *= 1.0e3; //cgs->mks
              Dcoeff[n] *= 1.0e-1; //cgs->mks   // Why 1e-1? Which unit is Dcoeff? kg/(m*s)?
              chi[n] *= Dcoeff[n]; //chi->theta
            }


            // species concentration
            CKYTCR(rho_a(i,j,k), T_a(i,j,k), Yloc, ci);
            for (int n=0; n< NUM_SPECIES; n++) {
              //ci[n] *= 1e3; // correction for units
              ci[n] *= 1e3*1e-6; // correction + mks -> cgs
            }

            // production rate
            //productionRate(wdot, ci, T_a(i,j,k))
            //CKWC(T_a(i,j,k),ci*1e6,wdot); // ci in mks -> ci*1e6 for cgs
            //productionRate(wdot, ci, T_a(i,j,k)); // is in mks
            //Real wdot_local = wdot[fuelID]; //1.0e-6;  // cgs -> mks

            CKWC(T_a(i,j,k),ci,wdot); // returns molar production rate
            Real wdot_local = wdot[fuelID]*(1/imw[fuelID])*1e6;  // convert to mks and to mass production rate as kg/(m^3*s)


            //I_R(i,j,k) = wdot_local;  // overwrite the reaction rate
            
            //Enthalpy flux due to diff diff cofficient
                  //for (int n = 0; n < NUM_SPECIES; ++n) {

            spec_coeff_a(i,j,k) = Dcoeff[fuelID] * mmw * imw[fuelID];
            soret_coeff_a(i,j,k) = 0.664*chi[fuelID]/T_a(i,j,k); //divide by T here for convenience
            molar_coeff_a(i,j,k) = (do_wbar == 1) ? Dcoeff[fuelID] * Yloc[fuelID] * imw[fuelID] : 0.0;
            //}

          } else {
            // Setting values outside of Domain to 0
            w_a(i,j,k) = 0; 
            spec_coeff_a(i,j,k) = 0;
            soret_coeff_a(i,j,k) = 0;
            molar_coeff_a(i,j,k) = 0;
          }

	});
      }

      if(verbose) Print() << "Coefficients, enthalpy and mmw derived for level " << lev << std::endl;
    }


//------------------------------------------------------------------------------------------------------------------------
// Compute Gradients 
//------------------------------------------------------------------------------------------------------------------------


#ifdef AMREX_USE_EB


    if(verbose) Print() << "Setting up Solver!" << std::endl;
    LPInfo info_apply;
    info_apply.setMaxCoarseningLevel(0);  
    MLEBABecLap poisson(geoms, grids, dmap, info_apply, amrex::GetVecOfConstPtrs(eb_factory));
    poisson.setVerbose(4);  
    poisson.setMaxOrder(4);
    
    //Poisson like solver for EB
    poisson.setScalars(0.0, 1.0);   
    for (int lev = 0; lev <= finestLevel; ++lev) {
      poisson.setBCoeffs(lev, -1.0);
    }

    if(verbose) Print() << "Setting Boundary Conditions" << std::endl;
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
    if(verbose) Print() << "BC's set and solver setup finished" << std::endl;

#else

    if(verbose) Print() << "Setting up Solver!" << std::endl;
    // Get face-centered gradients from MLMG                                                                                          
    LPInfo info;
    info.setAgglomeration(1);
    info.setConsolidation(1);
    info.setMetricTerm(false);
    info.setMaxCoarseningLevel(0);
    MLPoisson poisson({geoms}, {grids}, {dmap}, info);
    poisson.setMaxOrder(4);


    if(verbose) Print() << "Setting Boundary Conditions" << std::endl;
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
    if(verbose) Print() << "BC's set and solver setup finished" << std::endl;
#endif   


    // Need to apply the operator to ensure CF consistency with composite solve                                                       
    int nGrowGrad = 0;                   // No need for ghost face on gradient   

    Vector<Array<MultiFab,AMREX_SPACEDIM> > grad(Nlev);
    Vector<std::unique_ptr<MultiFab>> phi;
    Vector<MultiFab> laps;
    
    //for (int n =0; n<NUM_SPECIES; n++) {
    if(verbose) Print() << "Calculating grad"+inNames[idProglocal] << std::endl;
    for (int lev = 0; lev < Nlev; ++lev) {
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
        const auto& ba = grids[lev];
        #ifdef AMREX_USE_EB
          grad[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                    dmap[lev], 1, nGrowGrad, MFInfo(), *eb_factory[lev]);
        #else
          grad[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                    dmap[lev], 1, nGrowGrad);
        #endif
      }

      phi.push_back(std::make_unique<MultiFab>(indata[lev],amrex::make_alias,idProglocal,1));
      poisson.setLevelBC(lev, phi[lev].get());
      laps.emplace_back(grids[lev], dmap[lev], 1, 1);

    }

    MLMG mlmg(poisson);
    mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
    #ifdef AMREX_USE_EB
      mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCentroid);
    #else
      mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);
    #endif
    phi.clear();

    for (int lev = 0; lev < Nlev; ++lev) {
      MultiFab gradAlias(gradMF[lev], amrex::make_alias, 0, AMREX_SPACEDIM); //put the gradient in here
      
      #ifdef AMREX_USE_EB
        EB_average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
        // Multiplication with (-1) not needed
      #else
        average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
        gradAlias.mult(-1.0);
      #endif
    }
    //}


    if(verbose) Print() << "Calculating grad"+inNames[idTlocal] << std::endl;
    for (int lev = 0; lev < Nlev; ++lev) {
      phi.push_back(std::make_unique<MultiFab>(indata[lev],amrex::make_alias,idTlocal,1));
      poisson.setLevelBC(lev, phi[lev].get());
    }
    MLMG mlmg2(poisson); 
    mlmg2.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
    #ifdef AMREX_USE_EB
      mlmg2.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCentroid);
    #else
      mlmg2.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);
    #endif
    phi.clear();
    for (int lev = 0; lev < Nlev; ++lev) {
      MultiFab gradAlias(gradMF[lev], amrex::make_alias, AMREX_SPACEDIM, AMREX_SPACEDIM); //put the gradient in here
      #ifdef AMREX_USE_EB
        EB_average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
        // Multiplication with (-1) not needed
      #else
        average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
        gradAlias.mult(-1.0);
      #endif
    }
      
    if(verbose) Print() << "Calculating W gradient..." << std::endl;
    //for (int n = 0; n < 2; n++) {
    for (int lev = 0; lev < Nlev; ++lev) {
      phi.push_back(std::make_unique<MultiFab>(deriveMF[lev],amrex::make_alias,0,1));
      poisson.setLevelBC(lev, phi[lev].get());
    }
    MLMG mlmg3(poisson);
    mlmg3.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
    #ifdef AMREX_USE_EB
      mlmg3.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCentroid);
    #else
      mlmg3.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);
    #endif
    phi.clear();
    for (int lev = 0; lev < Nlev; ++lev) {
      MultiFab gradAlias(gradMF[lev], amrex::make_alias, 2*AMREX_SPACEDIM, AMREX_SPACEDIM); // Alias on gradMF[lev] including the components starting from 2*AMREX_SPACEDIM and containing AMREX_SPACEDIM components
      #ifdef AMREX_USE_EB
        EB_average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
        // Multiplication with (-1) not needed
      #else
        average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
        gradAlias.mult(-1.0);
      #endif
    }
    

//------------------------------------------------------------------------------------------------------------------------
// Combining Coefficients and Gradients
//------------------------------------------------------------------------------------------------------------------------

    //So we have the coefficients and gradients in all directions, lets combine into the fluxes for all three directions.
    //Store in derive MF
    for (int lev=0; lev<Nlev; ++lev)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(indata[lev], amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const Box& bx = mfi.tilebox();
	Array4<Real> const& spec_coeff_a = coeffMF[lev].array(mfi,0);
	Array4<Real> const& soret_coeff_a = coeffMF[lev].array(mfi,1);
	Array4<Real> const& molar_coeff_a = coeffMF[lev].array(mfi,2);
	Array4<Real> const& gradY = gradMF[lev].array(mfi,0);
	Array4<Real> const& gradT = gradMF[lev].array(mfi,AMREX_SPACEDIM);
	Array4<Real> const& gradW = gradMF[lev].array(mfi,2*AMREX_SPACEDIM);

	Array4<Real> const& flux_a = deriveMF[lev].array(mfi,1);

  // For EB
#ifdef AMREX_USE_EB
  Array4<const Real> const&  volFracBox = eb_factory[lev]->getVolFrac()[mfi].array();
#endif

  Array4<Real> const& sd_a = outdata[lev].array(mfi,0);

	amrex::ParallelFor(bx, [=]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

        #ifdef AMREX_USE_EB
          Real volFrac = volFracBox(i,j,k,0);
        #else
          Real volFrac = 1; 
        #endif

          if(volFrac > 0 ) {

	  for (int d = 0; d < AMREX_SPACEDIM; d++) {
	    flux_a(i,j,k,d) = spec_coeff_a(i,j,k)*gradY(i,j,k,d);           
	    flux_a(i,j,k,d+AMREX_SPACEDIM) = soret_coeff_a(i,j,k)*gradT(i,j,k,d);
	    flux_a(i,j,k,d+2*AMREX_SPACEDIM) = molar_coeff_a(i,j,k)*gradW(i,j,k,d);

	  }

          // } else {
          //   // Setting values to 0 outside of valid domain -> Better way would be to initialise flux_a with 0 when loading data
          //   for (int d = 0; d < AMREX_SPACEDIM; d++) {
          //     flux_a(i,j,k,d) = 0;
          //     flux_a(i,j,k,d+AMREX_SPACEDIM) = 0;
          //     flux_a(i,j,k,d+2*AMREX_SPACEDIM) = 0;
          //   }
          }

	});

			   
    
      }
    }

//------------------------------------------------------------------------------------------------------------------------
// Divergence of Fluxes
//------------------------------------------------------------------------------------------------------------------------

    //we now have the diffusive fluxes, gradients and the chemical source term. Let's compute divergence of diffusive flux...
    for (int n =0; n<AMREX_SPACEDIM; n++) {


      if(verbose) Print() << "Calculating gradFspecies"+std::to_string(n) << std::endl;
      for (int lev = 0; lev < Nlev; ++lev) {
	phi.push_back(std::make_unique<MultiFab>(deriveMF[lev],amrex::make_alias,n+1,1)); //plus 1 because 1 is W
	poisson.setLevelBC(lev, phi[lev].get());
	laps.emplace_back(grids[lev], dmap[lev], 1, 1);
      }
      MLMG mlmg(poisson);
      mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
    #ifdef AMREX_USE_EB
      mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCentroid);
    #else
      mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);
    #endif
      phi.clear();
      for (int lev = 0; lev < Nlev; ++lev) {
	MultiFab gradAlias(grad2MF[lev], amrex::make_alias, n*AMREX_SPACEDIM, AMREX_SPACEDIM); //put the gradient in here
	    #ifdef AMREX_USE_EB
        EB_average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
        // Multiplication with (-1) not needed
      #else
        average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
        gradAlias.mult(-1.0);
      #endif
      }
    }

//------------------------------------------------------------------------------

    for (int n =0; n<AMREX_SPACEDIM; n++) {
      if (verbose) Print() << "Calculating gradFsoret"+std::to_string(n) << std::endl;
      for (int lev = 0; lev < Nlev; ++lev) {
	phi.push_back(std::make_unique<MultiFab>(deriveMF[lev],amrex::make_alias,n+1+AMREX_SPACEDIM,1)); //plus 1 because 1 is W
	poisson.setLevelBC(lev, phi[lev].get());
	laps.emplace_back(grids[lev], dmap[lev], 1, 1);
      }
      MLMG mlmg(poisson);
      mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
    #ifdef AMREX_USE_EB
      mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCentroid);
    #else
      mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);
    #endif
      phi.clear();
      for (int lev = 0; lev < Nlev; ++lev) {
	MultiFab gradAlias(grad2MF[lev], amrex::make_alias, n*AMREX_SPACEDIM + AMREX_SPACEDIM*AMREX_SPACEDIM, AMREX_SPACEDIM); //put the gradient in here
	    #ifdef AMREX_USE_EB
        EB_average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
        // Multiplication with (-1) not needed
      #else
        average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
        gradAlias.mult(-1.0);
      #endif
      }
    }

//------------------------------------------------------------------------------

    for (int n =0; n<AMREX_SPACEDIM; n++) {
      if(verbose) Print() << "Calculating gradFmolar"+std::to_string(n) << std::endl;
      for (int lev = 0; lev < Nlev; ++lev) {
	phi.push_back(std::make_unique<MultiFab>(deriveMF[lev],amrex::make_alias,n+1+2*AMREX_SPACEDIM,1)); //plus 1 because 1 is W
	poisson.setLevelBC(lev, phi[lev].get());
	laps.emplace_back(grids[lev], dmap[lev], 1, 1);
      }
      MLMG mlmg(poisson);
      mlmg.apply(GetVecOfPtrs(laps), GetVecOfPtrs(phi));
    #ifdef AMREX_USE_EB
      mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCentroid);
    #else
      mlmg.getFluxes(GetVecOfArrOfPtrs(grad), GetVecOfPtrs(phi), MLMG::Location::FaceCenter);
    #endif
      phi.clear();
      for (int lev = 0; lev < Nlev; ++lev) {
	MultiFab gradAlias(grad2MF[lev], amrex::make_alias, n*AMREX_SPACEDIM + 2*AMREX_SPACEDIM*AMREX_SPACEDIM, AMREX_SPACEDIM); //put the gradient in here
	    #ifdef AMREX_USE_EB
        EB_average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
        // Multiplication with (-1) not needed
      #else
        average_face_to_cellcenter(gradAlias, 0, GetArrOfConstPtrs(grad[lev]));
        gradAlias.mult(-1.0);
      #endif
      }
    }

//------------------------------------------------------------------------------------------------------------------------
// Compute Flame Displacement Speed
//------------------------------------------------------------------------------------------------------------------------

    //now we have all the second derivatives, lets combine everything for Sd
    for (int lev=0; lev<Nlev; ++lev)
    {


#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(indata[lev], amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

              const Box& bx = mfi.tilebox();
	Array4<Real> const& sd_a = outdata[lev].array(mfi,0);

	Array4<Real> const& gradY_a = gradMF[lev].array(mfi,0);
	Array4<Real> const& diff_spec_a = grad2MF[lev].array(mfi,0);
	Array4<Real> const& diff_soret_a = grad2MF[lev].array(mfi,AMREX_SPACEDIM*AMREX_SPACEDIM);
	Array4<Real> const& diff_molar_a = grad2MF[lev].array(mfi,2*AMREX_SPACEDIM*AMREX_SPACEDIM);
	
	Array4<Real> const& c_a = indata[lev].array(mfi,idClocal);
	Array4<Real> const& rho_a = indata[lev].array(mfi,idRlocal);
	Array4<Real> const& rhoh_a = indata[lev].array(mfi,idElocal);

#ifdef AMREX_USE_EB
  Array4<const Real> const&  volFracBox = eb_factory[lev]->getVolFrac()[mfi].array();
#endif

	amrex::ParallelFor(bx, [=]
        AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {

      #ifdef AMREX_USE_EB
          Real volFrac = volFracBox(i,j,k,0);
      #else
          Real volFrac = 1; 
      #endif
    
    Real maggrad = 0;
	  Real lapdiff_spec = 0;
	  Real lapdiff_soret = 0;
	  Real lapdiff_molar = 0;  

    if(volFrac > 0 ) {
    maggrad = std::sqrt(AMREX_D_TERM(gradY_a(i,j,k,0)*gradY_a(i,j,k,0),+gradY_a(i,j,k,1)*gradY_a(i,j,k,1),+gradY_a(i,j,k,2)*gradY_a(i,j,k,2)));
	  lapdiff_spec = AMREX_D_TERM(diff_spec_a(i,j,k,0), + diff_spec_a(i,j,k,AMREX_SPACEDIM + 1), + diff_spec_a(i,j,k,2*AMREX_SPACEDIM+2));
	  lapdiff_soret = AMREX_D_TERM(diff_soret_a(i,j,k,0), + diff_soret_a(i,j,k,AMREX_SPACEDIM + 1), + diff_soret_a(i,j,k,2*AMREX_SPACEDIM+2));
	  lapdiff_molar = AMREX_D_TERM(diff_molar_a(i,j,k,0), + diff_molar_a(i,j,k,AMREX_SPACEDIM + 1), + diff_molar_a(i,j,k,2*AMREX_SPACEDIM+2));

    }

    // was maggrad >1e-2
	  if (maggrad > 1e-2 && volFrac > 0) {
	    sd_a(i,j,k,0) = -lapdiff_spec/(rho_a(i,j,k)*maggrad);
	    sd_a(i,j,k,1) = -lapdiff_soret/(rho_a(i,j,k)*maggrad); 
	    sd_a(i,j,k,2) = -lapdiff_molar/(rho_a(i,j,k)*maggrad);
	    sd_a(i,j,k,3) = -c_a(i,j,k)/(rho_a(i,j,k)*maggrad);
	    sd_a(i,j,k,4) = sd_a(i,j,k,0) + sd_a(i,j,k,1) + sd_a(i,j,k,2) + sd_a(i,j,k,3);
	    sd_a(i,j,k,5) = sd_a(i,j,k,4)*rhoh_a(i,j,k);
      
	  } //else {
	  //   sd_a(i,j,k,0) = 0.0;
	  //   sd_a(i,j,k,1) = 0.0;
	  //   sd_a(i,j,k,2) = 0.0;
	  //   sd_a(i,j,k,3) = 0.0;
	  //   sd_a(i,j,k,4) = 0.0;
	  //   sd_a(i,j,k,5) = 0.0;
	  // }

      // sd_a(i,j,k,0) = 0.0;
	    // sd_a(i,j,k,1) = 0.0;
	    // sd_a(i,j,k,2) = 0.0;
	    // sd_a(i,j,k,3) = 0.0;
	    // sd_a(i,j,k,4) = 0.0;
	    // sd_a(i,j,k,5) = 0.0;

    // for(int ii=0; ii<6;ii++) {
    //     if( std::isnan(sd_a(i,j,k,ii))) {
    //       Print() << "sd_a is |inf| for level = " << lev << "and ii:" << ii << " , and volFrac= " << volFrac
    //     << " and x,y,z = " << xlo << "," << ylo << "," << zlo << " , and density = " 
    //     << rho_a(i,j,k) << " , and maggrad = " << maggrad <<std::endl; 
    //     Print() << "lapdiff_spec = " << lapdiff_spec << ", lapdiff_soret = " << lapdiff_soret 
    //     << ", and lapdiff_molar = " << lapdiff_molar << std::endl; 
    //     Print() << "gradY_a = " << gradY_a(i,j,k,0) << std::endl;
    //     }  
    // }

	}); // ParallelFor


      } // MFIter
    } //lev-loop
    
    std::string outfile(getFileRoot(plotFileName) + "_sd");
    pp.query("outfile",outfile);
    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev-1,{AMREX_D_DECL(2, 2, 2)});
    Real time = pf.time();
    VisMF::SetNOutFiles(n_files);
    amrex::WriteMultiLevelPlotfile(outfile, Nlev, GetVecOfConstPtrs(outdata), outNames,
                                   geoms, time, isteps, refRatios);
  }

  ParallelDescriptor::Barrier();
  Finalize();
  return 0;
}
