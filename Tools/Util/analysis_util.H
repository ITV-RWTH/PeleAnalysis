#ifndef analysis_util_H
#define analysis_util_H

#include <istream>
#include <string>
#include <vector>

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

namespace analysis_util {

// --- Global variables ---
extern int nlev;
extern amrex::Vector<amrex::BoxArray> grids;
extern amrex::Vector<amrex::DistributionMapping> dmap;
extern amrex::Vector<int> ref_ratio;
extern amrex::Vector<amrex::Geometry> geoms;
extern amrex::IntVect periodicity;
extern bool initialized;

// --- Initialization ---
void init(
  const amrex::Vector<amrex::MultiFab>& a_mf,
  const amrex::Vector<amrex::Geometry>& a_geoms);

// --- AMR utilities ---
[[nodiscard]]
amrex::Vector<std::unique_ptr<amrex::iMultiFab>> get_covered_mf();

// --- Integration tools ---
[[nodiscard]]
std::unique_ptr<amrex::Gpu::ManagedVector<amrex::Real>> integrate(
  const amrex::Vector<amrex::MultiFab>& a_mf,
  int scomp,
  int ncomp,
  const amrex::Vector<int>& axes_to_integrate);

// --- Gradient tools ---
void gradient(
  const amrex::Vector<amrex::MultiFab> a_mf,
  const int scomp,
  const int ncomp,
  amrex::Vector<amrex::MultiFab> grad_mf);

// --- String / file utilities ---

// Strip directory prefix: "path/to/plt00001" -> "plt00001"
std::string get_file_root(const std::string& infile);

// Read a single title/label line from a stream
std::string parse_title(std::istream& is);

// Read and tokenize a variable names line from a stream
std::vector<std::string> parse_var_names(std::istream& is);

// Find component index of a named variable; aborts if not found and
// abort_if_not_found=true
int find_var_index(
  const amrex::Vector<std::string>& var_names,
  const std::string& name,
  bool abort_if_not_found = true);

// --- I/O tools ---

struct PlotfileData
{
  amrex::Vector<amrex::MultiFab> mf;
  amrex::Vector<amrex::Geometry> geoms;
  amrex::Vector<int> ref_ratios;
  amrex::Real time;
  int n_lev;
  amrex::Vector<std::string> var_names; // all variables present in the file
};

// Load selected variables from an AMReX plotfile into MultiFabs + Geometries.
// The returned PlotfileData can be passed directly to init().
[[nodiscard]]
PlotfileData read_plotfile(
  const std::string& infile,
  const amrex::Vector<std::string>& var_names,
  int finest_level = 1000,
  int n_grow = 0,
  const amrex::Vector<int>& is_per = {});

// Write multi-level MultiFabs to an AMReX plotfile.
// ref_ratios defaults to all-2 if empty.
void write_plotfile(
  const std::string& outfile,
  const amrex::Vector<amrex::MultiFab>& mf,
  const amrex::Vector<std::string>& var_names,
  const amrex::Vector<amrex::Geometry>& geoms,
  amrex::Real time = 0.0,
  const amrex::Vector<int>& ref_ratios = {});

// --- MEF file tools ---

struct MEFData
{
  std::string title;
  std::vector<std::string> var_names;
  amrex::FArrayBox nodes;          // component-major layout [nNodes x nComp]
  amrex::Vector<int> connectivity; // flat array [n_elts * nodes_per_elt]
  int n_elts;
  int nodes_per_elt;
};

[[nodiscard]]
MEFData read_mef(const std::string& infile);

void write_mef(const std::string& outfile, const MEFData& data);

} // namespace analysis_util

#endif // analysis_util_H
