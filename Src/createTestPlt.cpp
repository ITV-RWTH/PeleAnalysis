#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        // Read parametes
        ParmParse pp;
        
        // Geometry parms
        Vector<int> n_cell(AMREX_SPACEDIM);
        pp.getarr("amr.n_cell", n_cell, 0, AMREX_SPACEDIM);
        
        // Domain parms
        RealBox real_box;
        Vector<Real> pp_prob_x(AMREX_SPACEDIM, 0.0);
        pp.getarr("geometry.prob_lo", pp_prob_x, 0, AMREX_SPACEDIM);
        Array<Real,AMREX_SPACEDIM> prob_lo = {AMREX_D_DECL(pp_prob_x[0],pp_prob_x[1],pp_prob_x[2])};
        pp.getarr("geometry.prob_hi", pp_prob_x, 0, AMREX_SPACEDIM);
        Array<Real,AMREX_SPACEDIM> prob_hi = {AMREX_D_DECL(pp_prob_x[0],pp_prob_x[1],pp_prob_x[2])};
        
        for (int i = 0; i < AMREX_SPACEDIM; i++) {
            real_box.setLo(i, prob_lo[i]);
            real_box.setHi(i, prob_hi[i]);
        }
        
        // Periodicity
        IntVect pp_is_per;
        pp.getarr("geometry.is_periodic",pp_is_per);
        Array<int,AMREX_SPACEDIM> is_per = {AMREX_D_DECL(pp_is_per[0],pp_is_per[1],pp_is_per[2])};
        
        // Coordinate system
        int coord = 0;
        pp.query("geometry.coord_sys", coord);

        // Only cartesioan supported
        AMREX_ALWAYS_ASSERT(coord==0)
        
        // Create geometry
        IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
        IntVect domain_hi(AMREX_D_DECL(n_cell[0]-1, 
                                        n_cell[1]-1, 
                                        n_cell[2]-1));
        Box domain(domain_lo, domain_hi);
        
        Geometry geom(domain, real_box, coord, is_per);
        
        // Create BoxArray and DistributionMapping 
        BoxArray ba(domain);
        
        // Dist domain
        int max_grid_size = 32;
        pp.query("amr.max_grid_size", max_grid_size);
        ba.maxSize(max_grid_size);
        int ngrow = 0; // Ghost cells
        pp.query("amr.ngrow", ngrow);
        
        DistributionMapping dm(ba);
        
        // Get components
        int ncomp = pp.countval("field.comps");
        Vector<std::string> comps(ncomp);
        if (ncomp>0) {
            pp.getarr("field.comps",comps);
            Print() << "Components:" << std::endl;
            for (int i=0; i<ncomp; ++i) {
                Print() << "   " << i << ": " << comps[i] << std::endl;
            }
        } else {
            Abort("No test fields specified...");
        }
        
        // Create MultiFAB
        MultiFab mf(ba, dm, ncomp, ngrow);
        
        // MultiFAB mit Testdaten füllen
        for (MFIter mfi(mf); mfi.isValid(); ++mfi) {
            const Box& box = mfi.validbox();
            auto const& fab = mf.array(mfi);
            
            amrex::ParallelFor(box, ncomp,
            [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                // Beispiel: Sinusförmige Testdaten
                fab(i,j,k,n) = std::sin(2.0*M_PI*i/n_cell[0]) * 
                               std::cos(2.0*M_PI*j/n_cell[1]);
            });
        }
        
        // var names for plotfile
        Vector<std::string> varnames(ncomp);
        for (int n = 0; n < ncomp; ++n) {
            varnames[n] = "var_" + std::to_string(n);
        }
        
        // write plot file
        std::string plotfile_name = "pltTestFile";
        pp.query("plotfile_name", plotfile_name);
        
        WriteSingleLevelPlotfile(plotfile_name, mf, varnames, geom, 0.0, 0);
        
        amrex::Print() << "Test-Plotfile erstellt: " << plotfile_name << "\n";
    }
    amrex::Finalize();
    return 0;
}
