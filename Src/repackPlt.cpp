#include <filesystem>
#include <fstream>
#include <algorithm>
#include <cstring>
#include <regex>
#include <set>

#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_MultiFab.H"
#include "AMReX_VisMF.H"
#include "AMReX_Utility.H"

using namespace amrex;
namespace fs = std::filesystem;

// Rewrite existing AMReX plotfiles with a user-chosen number of data files
// per MultiFab (n_files). Each input is written to <outdir>/<same name>.  Only
// the level data (Level_*/<mf>_H and Level_*/<mf>_D_*) is re-packed; every
// other file or directory (Header, job_info, particles, ...) is copied
// verbatim, so the output is identical to the input apart from the data-file
// layout.

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "Usage:\n"
    << "  " << argv[0] << " infiles=DIR [DIR ...] [OPTIONS]\n\n"

    << "Required arguments:\n"
    << "  infiles=DIR ...    One or more AMReX plotfiles to re-pack\n\n"

    << "Options:\n"
    << "  outdir=DIR         Directory for the re-packed plotfiles, which\n"
    << "                     keep their names (DEF: repack_nf<n_files>)\n"
    << "  n_files=INT        Number of data files per MultiFab (DEF: 1).\n"
    << "                     Capped at the number of MPI ranks.\n"
    << "  float_type=32|64   Precision of the written data (DEF: same as "
       "input)\n"
    << "  verbose=1          Print progress information\n"
    << "  -h, --help         Show this help message\n\n"
    << "Visit PeleAnalysis/Src/InputSamples for examples or refer to "
    << "the documentation.\n";

  std::exit(1);
}

// Bytes per real of the data described by a MultiFab header <mfpath>_H.
// Version_v1 headers keep the RealDescriptor in each FAB header inside the
// data files; the NoFabHeader* versions store it as the last line of the
// header file.  Returns 0 if it cannot be determined.
static int
detectRealBytes(const std::string& mfpath, int vers)
{
  // RealDescriptor: ((<n>, (<bits> <exp bits> <mantissa bits> ...)),(<bytes>,
  // (<byte order>))); the first format entry is the number of bits per real.
  static const std::regex rd_re(R"(\(\(\s*\d+\s*,\s*\(\s*(\d+))");
  std::smatch m;
  if (vers == VisMF::Header::Version_v1) {
    // First FabOnDisk entry names a data file; its FAB header is at offset.
    std::ifstream hdr(mfpath + "_H");
    std::string line;
    while (std::getline(hdr, line)) {
      std::istringstream ls(line);
      std::string tag, fname;
      Long offset;
      if ((ls >> tag >> fname >> offset) && tag == "FabOnDisk:") {
        fs::path dfile = fs::path(mfpath).parent_path() / fname;
        std::ifstream df(dfile, std::ios::binary);
        df.seekg(offset);
        std::string fabline;
        std::getline(df, fabline);
        if (std::regex_search(fabline, m, rd_re)) {
          return std::stoi(m[1]) / 8;
        }
        return 0;
      }
    }
  } else {
    std::ifstream hdr(mfpath + "_H");
    std::string line, last;
    while (std::getline(hdr, line)) {
      if (line.rfind("((", 0) == 0) {
        last = line;
      }
    }
    if (std::regex_search(last, m, rd_re)) {
      return std::stoi(m[1]) / 8;
    }
  }
  return 0;
}

static int
readHeaderVersion(const std::string& mfpath)
{
  std::ifstream hdr(mfpath + "_H");
  int vers = VisMF::Header::Undefined_v1;
  hdr >> vers;
  return vers;
}

// Re-pack one plotfile <infile> into <outfile> (which must not exist yet or
// is moved aside by UtilCreateCleanDirectory).
static void
repackOne(
  const std::string& infile,
  const std::string& outfile,
  int float_type,
  int verbose)
{
  // ------------------------------------------------------------------
  // Classify the content of infile: level dirs hold MultiFabs to repack,
  // everything else is copied as is.
  // ------------------------------------------------------------------
  static const std::regex level_re(R"(Level_\d+)");
  std::vector<std::string> levelDirs;
  std::vector<fs::path> copyEntries; // relative to infile
  std::vector<std::string> mfNames;  // e.g. Level_0/Cell
  for (const auto& e : fs::directory_iterator(infile)) {
    const std::string name = e.path().filename().string();
    if (e.is_directory() && std::regex_match(name, level_re)) {
      levelDirs.push_back(name);
      std::set<std::string> mfs;
      for (const auto& le : fs::directory_iterator(e.path())) {
        const std::string lname = le.path().filename().string();
        if (
          le.is_regular_file() && lname.size() > 2 &&
          lname.compare(lname.size() - 2, 2, "_H") == 0) {
          mfs.insert(lname.substr(0, lname.size() - 2));
        }
      }
      // Anything in the level dir that does not belong to a MultiFab
      for (const auto& le : fs::directory_iterator(e.path())) {
        const std::string lname = le.path().filename().string();
        bool isMF = false;
        for (const auto& mf : mfs) {
          if (lname == mf + "_H" || lname.rfind(mf + "_D_", 0) == 0) {
            isMF = true;
            break;
          }
        }
        if (!isMF) {
          copyEntries.push_back(fs::path(name) / lname);
        }
      }
      for (const auto& mf : mfs) {
        mfNames.push_back(name + "/" + mf);
      }
    } else {
      copyEntries.push_back(fs::path(name));
    }
  }
  std::sort(mfNames.begin(), mfNames.end());
  if (mfNames.empty()) {
    Abort("No MultiFab data found in " + infile);
  }

  // ------------------------------------------------------------------
  // Create output directory tree and copy non-MultiFab content
  // ------------------------------------------------------------------
  UtilCreateCleanDirectory(outfile, true);
  if (ParallelDescriptor::IOProcessor()) {
    for (const auto& d : levelDirs) {
      fs::create_directories(fs::path(outfile) / d);
    }
    for (const auto& rel : copyEntries) {
      if (verbose) {
        Print() << "Copying " << rel.string() << "\n";
      }
      fs::copy(
        fs::path(infile) / rel, fs::path(outfile) / rel,
        fs::copy_options::recursive | fs::copy_options::copy_symlinks);
    }
  }
  ParallelDescriptor::Barrier();

  // ------------------------------------------------------------------
  // Re-pack every MultiFab, one at a time
  // ------------------------------------------------------------------
  const auto origVersion = VisMF::GetHeaderVersion();
  for (const auto& mfname : mfNames) {
    const std::string inpath = infile + "/" + mfname;
    const std::string outpath = outfile + "/" + mfname;

    int vers = readHeaderVersion(inpath);
    if (
      vers < VisMF::Header::Version_v1 ||
      vers > VisMF::Header::NoFabHeaderFAMinMax_v1) {
      Abort("Unknown VisMF header version in " + inpath + "_H");
    }
    int nbytes = float_type / 8;
    if (nbytes == 0) {
      if (ParallelDescriptor::IOProcessor()) {
        nbytes = detectRealBytes(inpath, vers);
      }
      ParallelDescriptor::Bcast(
        &nbytes, 1, ParallelDescriptor::IOProcessorNumber());
      if (nbytes != 4 && nbytes != 8) {
        Abort(
          "Could not detect precision of " + inpath +
          "; set float_type=32 or 64");
      }
    }

    MultiFab mf;
    VisMF::Read(mf, inpath);

    VisMF::SetHeaderVersion(static_cast<VisMF::Header::Version>(vers));
    FArrayBox::setFormat(
      nbytes == 4 ? FABio::FAB_NATIVE_32 : FABio::FAB_NATIVE);
    VisMF::Write(mf, outpath);

    if (verbose) {
      Print() << "Re-packed " << mfname << ": " << mf.boxArray().size()
              << " boxes, " << mf.nComp() << " comps, " << 8 * nbytes
              << "-bit, header v" << vers << "\n";
    }
  }
  VisMF::SetHeaderVersion(origVersion);
}

int
main(int argc, char* argv[])
{
  Initialize(argc, argv);
  {
    if (argc < 2) {
      print_usage(argc, argv);
    } else if (
      (std::strcmp(argv[1], "-h") == 0) ||
      (std::strcmp(argv[1], "--help") == 0)) {
      print_usage(argc, argv);
    }

    ParmParse pp;

    int verbose = 0;
    pp.query("verbose", verbose);

    const int nin = pp.countval("infiles");
    if (nin < 1) {
      Abort("No input plotfiles given; set infiles = plt1 [plt2 ...]");
    }
    Vector<std::string> infiles(nin);
    pp.getarr("infiles", infiles, 0, nin);
    for (auto& f : infiles) {
      while (f.size() > 1 && f.back() == '/') {
        f.pop_back();
      }
    }

    int n_files = 1;
    pp.query("n_files", n_files);
    if (n_files < 1) {
      Abort("n_files must be >= 1");
    }

    std::string outdir = "repack_nf" + std::to_string(n_files);
    pp.query("outdir", outdir);
    while (outdir.size() > 1 && outdir.back() == '/') {
      outdir.pop_back();
    }

    int float_type = 0; // 0: same as input
    pp.query("float_type", float_type);
    if (float_type != 0 && float_type != 32 && float_type != 64) {
      Abort("float_type must be 32 or 64");
    }

    // ------------------------------------------------------------------
    // Validate all inputs before writing anything
    // ------------------------------------------------------------------
    Vector<std::string> outfiles(nin);
    std::set<std::string> basenames;
    for (int i = 0; i < nin; ++i) {
      const std::string& in = infiles[i];
      if (
        !fs::is_directory(in) || !fs::exists(fs::path(in) / "Header") ||
        !fs::is_directory(fs::path(in) / "Level_0")) {
        Abort("infile " + in + " is not an AMReX plotfile directory");
      }
      const std::string base = fs::path(in).filename().string();
      if (!basenames.insert(base).second) {
        Abort(
          "Two input plotfiles share the name " + base +
          "; they would collide in outdir " + outdir);
      }
      outfiles[i] = outdir + "/" + base;
      if (
        fs::exists(outfiles[i]) &&
        fs::equivalent(fs::path(in), fs::path(outfiles[i]))) {
        Abort(
          "outdir " + outdir + " contains the input " + in +
          "; choose a different outdir");
      }
    }

    VisMF::SetNOutFiles(n_files);
    if (ParallelDescriptor::IOProcessor() && VisMF::GetNOutFiles() < n_files) {
      Print() << "WARNING: n_files=" << n_files << " is larger than the number"
              << " of MPI ranks; writing " << VisMF::GetNOutFiles()
              << " file(s) per MultiFab instead\n";
    }

    if (ParallelDescriptor::IOProcessor()) {
      fs::create_directories(outdir);
    }
    ParallelDescriptor::Barrier();

    for (int i = 0; i < nin; ++i) {
      Print() << "[" << i + 1 << "/" << nin << "] " << infiles[i] << " -> "
              << outfiles[i] << "\n";
      repackOne(infiles[i], outfiles[i], float_type, verbose);
    }

    Print() << "Wrote " << nin << " plotfile(s) to " << outdir
            << " with at most " << VisMF::GetNOutFiles()
            << " data file(s) per MultiFab\n";
  }
  Finalize();
  return 0;
}
