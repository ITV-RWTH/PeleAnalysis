#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <PltFileManager.H>

using namespace amrex;

static void
print_usage(int, char* argv[])
{
  std::cerr << "Usage:\n"
            << "  " << argv[0] << " infiles=FILE1 FILE2 ... [OPTIONS]\n\n"

            << "Required arguments:\n"
            << "  infiles=LIST       List of AMReX plotfiles to average\n\n"

            << "Options:\n"
            << "  -h, --help         Show this help message\n\n"

            << "For optional arguments and examples, refer to the documentation\n"
            << "or see PeleAnalysis/Src/InputSamples/\n";

  std::exit(1);
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
     // ParmParse
     // ---------------------------------------------------------------------
     ParmParse pp;

     if (pp.contains("help")) {
       print_usage(argc,argv);
     }

     // Arbitrary number of input files will be combined and averaged
     int nf = pp.countval("infiles");
     AMREX_ALWAYS_ASSERT(nf>0);
     Vector<std::string> plotFileNames; pp.getarr("infiles",plotFileNames,0,nf);
     // into a single output file
     std::string outfile("plt_averaged");
     pp.query("outfile",outfile);

     // Vairables to keep - will keep all if not specified
     int nvar = pp.countval("variables");
     Vector<std::string> variableNames;
     pp.queryarr("variables",variableNames,0,nvar);
     bool all_vars = nvar > 0 ? false : true;
     Vector<Vector<int>> var_idxs;

     // Maximum number of levels to keep
     int output_max_level = 1000;
     pp.query("output_max_level", output_max_level);
     output_max_level +=1; // account for base level

     // Max grid size in output data
     int output_max_grid_size = 32;
     pp.query("output_max_grid_size", output_max_grid_size);

     // Type of interpolation to do
     int interp_type = 1;
     pp.query("interp_type",interp_type);

     // Control over individual operations
     int do_divide = 1;
     pp.query("do_divide",do_divide);

     int do_average = 1;
     pp.query("do_average",do_average);

     int do_variance = 1;
     pp.query("do_variance",do_variance);

     // Validation: must compute at least average or variance
     if (do_average == 0 && do_variance == 0) {
       amrex::Abort("At least one of do_average or do_variance must be 1");
     }

     Print() << "Options: do_average=" << do_average << " do_variance=" << do_variance << " do_divide=" << do_divide << std::endl;

     // ---------------------------------------------------------------------
     // Execute
     // ---------------------------------------------------------------------

     // First load the metadata of each input plt file
     Print() << "Loading plt file metadata..." << std::endl;
     Vector<pele::physics::pltfilemanager::PltFileManager*> plt_file_data(nf);
     int nlevels = 0;
     for (int i = 0; i < plt_file_data.size(); ++i) {
       plt_file_data[i] = new pele::physics::pltfilemanager::PltFileManager(plotFileNames[i]);
       nlevels = max(nlevels, plt_file_data[i]->getNlev());
       // Verify we have the right vairables
       if (all_vars) {
         if (i==0) {
           variableNames = plt_file_data[i]->getVariableList();
           nvar = variableNames.size();
         } else {
           Vector<std::string> variableNamesTest = plt_file_data[i]->getVariableList();
           if (variableNamesTest.size() != nvar) {
             amrex::Abort("All plt files must have same number of variables unless variable list is specified. File: " + plotFileNames[i]);
           }
           for (int var = 0; var < nvar; ++var) {
             if (variableNames[var] != variableNamesTest[var]) {
               amrex::Abort("All plt files must have same variables unless variable list is specified. File: " + plotFileNames[i]);
             }
           }
         }
       } else {
         Vector<int> var_idx_loc;
         Vector<std::string> variableNamesPlt = plt_file_data[i]->getVariableList();
         for (int var = 0; var < nvar; ++var) {
           int pvar;
           for (pvar = 0; pvar < variableNamesPlt.size(); ++pvar) {
             if (variableNames[var] == variableNamesPlt[pvar]) {
               var_idx_loc.push_back(pvar);
               break;
             }
           }
           if (pvar == variableNamesPlt.size()) {
             amrex::Abort("Variable '" + variableNames[var] + "' not found in file: " + plotFileNames[i]);
           }
         }
         var_idxs.push_back(var_idx_loc);
       }
     }
     nlevels = min(nlevels, output_max_level);
     Print() << " -> Combining " << nf << " files across " << nlevels << " levels" << std::endl;
     // Find density index:
     Vector<std::string> variableNamesPlt = plt_file_data[0]->getVariableList();
     int idRho = -1;
     for (int var = 0; var < variableNamesPlt.size(); ++var) {
       if (variableNamesPlt[var] == "density") idRho = var;
     }
     if (idRho == -1) amrex::Abort("Density not found in file with index 0: " + plotFileNames[0]);

     // Loop over each input file to get union of boxes on each level
     // On any level with all same BoxArray, we use that without modification
     Print() << "Finding the combined grids..." << std::endl;
     Vector<BoxArray> combined_boxes;
     Vector<Geometry> level_geometries;
     Vector<int> boxarray_all_same(nlevels, 1);
     for (int i = 0; i < plt_file_data.size(); ++i) {
       int nlevels_file = min(nlevels, plt_file_data[i]->getNlev());
       for (int lev = 0; lev < nlevels_file; ++lev) {
         // Verify we have the same geometry
         if (level_geometries.size() <= lev) {
           level_geometries.push_back(plt_file_data[i]->getGeom(lev));
         } else {
           bool same_domain = AlmostEqual(plt_file_data[i]->getGeom(lev).ProbDomain(), level_geometries[lev].ProbDomain());
           bool same_box = plt_file_data[i]->getGeom(lev).Domain() == level_geometries[lev].Domain();
           if (!(same_domain and same_box)) {
             amrex::Abort("All plt files must have the same geometry");
           }
         }
         // Now combine the boxes - if not the same
         if (combined_boxes.size() <= lev) {
           combined_boxes.push_back(BoxArray(plt_file_data[i]->getGrid(lev)));
         } else {
           BoxList boxlist_file{plt_file_data[i]->getGrid(lev)};
           BoxList boxlist_combined{combined_boxes[lev]};
           if (boxlist_combined != boxlist_file) {
             boxlist_combined.catenate(boxlist_file);
             combined_boxes[lev] = BoxArray(boxlist_combined);
             combined_boxes[lev].removeOverlap();
             boxarray_all_same[lev] = 0;
           }
         }
       }
     }

     // Calculate number of output components based on what we're computing
     int ncomp_output = 1; // rho_mean
     if (do_average) ncomp_output += nvar;
     if (do_variance) ncomp_output += nvar;

     // Create the data structures to read in the data and keep running sums
     Vector<MultiFab> running_data(nlevels);
     Vector<MultiFab> running_data2(nlevels);
     Vector<MultiFab> running_rho(nlevels);
     Vector<MultiFab> tmp_data(nlevels);
     Vector<MultiFab> tmp_data2(nlevels);
     Vector<MultiFab> tmp_rho(nlevels);
     Vector<MultiFab> variance(nlevels);
     Vector<MultiFab> combined_data(nlevels);
     Vector<IntVect> refRatios(nlevels-1);
     for (int lev = 0; lev < nlevels; ++lev) {
       if (!boxarray_all_same[lev]) {
         combined_boxes[lev].maxSize(output_max_grid_size);
       }
       DistributionMapping dmap = DistributionMapping(combined_boxes[lev]);
       running_data[lev].define(combined_boxes[lev], dmap, nvar, 0);
       running_data2[lev].define(combined_boxes[lev], dmap, nvar, 0);
       running_rho[lev].define(combined_boxes[lev], dmap, 1, 0);
       tmp_data[lev].define(combined_boxes[lev], dmap, nvar, 0);
       tmp_data2[lev].define(combined_boxes[lev], dmap, nvar, 0);
       tmp_rho[lev].define(combined_boxes[lev], dmap, 1, 0);
       variance[lev].define(combined_boxes[lev], dmap, nvar, 0);
       combined_data[lev].define(combined_boxes[lev], dmap, ncomp_output, 0);
       running_data[lev].setVal(0.0);
       running_data2[lev].setVal(0.0);
       running_rho[lev].setVal(0.0);
       if (lev > 0) {
         int rr = int(level_geometries[lev-1].CellSize(0) / level_geometries[lev].CellSize(0));
         refRatios[lev-1] = {AMREX_D_DECL(rr,rr,rr)};
       }
     }
     // Fillpatch tmp_data from each pltfile and add to running data
     Print() << "Fillpatching and combining..." << std::endl;
     for (int i = 0; i < plt_file_data.size(); ++i) {
       Print() << "   working on file " << plotFileNames[i] << " (" << i+1 << "/" << plt_file_data.size() << ")" << std::endl;
       for (int lev = 0; lev < nlevels; ++lev) {
         plt_file_data[i]->fillPatchFromPlt(lev, level_geometries[lev], idRho, 0, 1, tmp_rho[lev], interp_type);
         if (all_vars) {
           plt_file_data[i]->fillPatchFromPlt(lev, level_geometries[lev], 0, 0, nvar, tmp_data[lev], interp_type);
         } else {
           for (int var = 0; var < nvar; ++var) {
             plt_file_data[i]->fillPatchFromPlt(lev, level_geometries[lev], var_idxs[i][var], var, 1, tmp_data[lev], interp_type);
           }
         }
         MultiFab::Copy(tmp_data2[lev], tmp_data[lev], 0, 0, nvar, 0);
         for (int var = 0; var < nvar; ++var) {
           MultiFab::Multiply(tmp_data[lev], tmp_rho[lev], 0, var, 1, 0);
           MultiFab::Multiply(tmp_data2[lev], tmp_data2[lev], var, var, 1, 0);
           MultiFab::Multiply(tmp_data2[lev], tmp_rho[lev], 0, var, 1, 0);
         }
         MultiFab::Add(running_data[lev], tmp_data[lev], 0, 0, nvar, 0);
         MultiFab::Add(running_data2[lev], tmp_data2[lev], 0, 0, nvar, 0);
         MultiFab::Add(running_rho[lev], tmp_rho[lev], 0, 0, 1, 0);
       }
       delete plt_file_data[i];
     }

     // Divide by number of files to get average
     Real factor = 1.0 / Real(nf);
     for (int lev = 0; lev < nlevels; ++lev) {
       running_data[lev].mult(factor);
       running_data2[lev].mult(factor);
       running_rho[lev].mult(factor);

       // Compute variance BEFORE dividing (need both rho-weighted quantities)
       // var(p) = E[p^2] - E[p]^2
       if (do_variance == 1) {
         MultiFab::Copy(variance[lev], running_data[lev], 0, 0, nvar, 0);
         MultiFab::Multiply(variance[lev], variance[lev], 0, 0, nvar, 0);
         variance[lev].mult(-1.0);
         MultiFab::Add(variance[lev], running_data2[lev], 0, 0, nvar, 0);
       }

       // Divide by rho_mean if requested (do this AFTER variance computation)
       if (do_divide == 1) {
         for (int var = 0; var < nvar; ++var) {
           MultiFab::Divide(running_data[lev], running_rho[lev], 0, var, 1, 0);
           MultiFab::Divide(running_data2[lev], running_rho[lev], 0, var, 1, 0);
           MultiFab::Divide(variance[lev], running_rho[lev], 0, var, 1, 0);
         }
       }
     }

     // Build output variable names based on what is being computed
     Vector<std::string> outputVariableNames;
     outputVariableNames.push_back("rho_mean");

     if (do_average) {
       for (int var = 0; var < nvar; ++var) {
         std::string variableName = variableNames[var];
         if (do_divide == 1) {
           outputVariableNames.push_back(variableName + "_favre_mean");
         } else {
           outputVariableNames.push_back(variableName + "_mean");
         }
       }
     }

     if (do_variance) {
       for (int var = 0; var < nvar; ++var) {
         std::string variableName = variableNames[var];
         if (do_divide == 1) {
           outputVariableNames.push_back(variableName + "_favre_variance");
         } else {
           outputVariableNames.push_back(variableName + "_variance");
         }
       }
     }

     // Populate combined_data with the computed fields
     for (int lev = 0; lev < nlevels; ++lev) {
       int icomp = 0;

       // Add rho_mean
       MultiFab::Copy(combined_data[lev], running_rho[lev], 0, icomp, 1, 0);
       icomp += 1;

       // Add averages if computed
       if (do_average) {
         MultiFab::Copy(combined_data[lev], running_data[lev], 0, icomp, nvar, 0);
         icomp += nvar;
       }

       // Add variances if computed
       if (do_variance) {
         MultiFab::Copy(combined_data[lev], variance[lev], 0, icomp, nvar, 0);
         icomp += nvar;
       }
     }

     // Save the final plt file
     Print() << "Saving final plt file..." << std::endl;
     Vector<int> stepidx(nlevels,0);
     WriteMultiLevelPlotfile(
       outfile,
       nlevels,
       GetVecOfConstPtrs(combined_data),
       outputVariableNames,
       level_geometries,
       0.0,
       stepidx,
       refRatios
     );
     Print() << "Done." << std::endl;
   }
   amrex::Finalize();
   return 0;
}
