#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

using namespace amrex;

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "Usage:\n"
    << "  " << argv[0] << " infiles=\"f1.mef\" \"f2.mef\" ... xvar=NAME yvar=NAME nBins=N [OPTIONS]\n\n"
    << "Required arguments:\n"
    << "  infiles=FILE ...     MEF input files (space-separated list)\n"
    << "  xvar=NAME            Variable name for x-axis (as in MEF file)\n"
    << "  yvar=NAME            Variable name for y-axis (as in MEF file)\n"
    << "  nBins=N              Number of histogram bins per axis\n\n"
    << "Options:\n"
    << "  xvar_out=NAME        Output label for x-variable (default: xvar)\n"
    << "  yvar_out=NAME        Output label for y-variable (default: yvar)\n"
    << "  xmin=VAL xmax=VAL    x-axis bounds (auto-detected if omitted)\n"
    << "  ymin=VAL ymax=VAL    y-axis bounds (auto-detected if omitted)\n"
    << "  -h                   Show this help message\n\n"
    << "Output: MEF_JPDFAverage_{xvar_out}_{yvar_out}/\n"
    << "  Pdf_{xvar_out}_x.dat          x bin centers\n"
    << "  Pdf_{yvar_out}_x.dat          y bin centers\n"
    << "  Pdf_{xvar_out}_{yvar_out}.dat 2D joint PDF matrix\n\n"
    << "Visit PeleAnalysis/Src/InputSamples for examples.\n";
  std::exit(1);
}

static std::vector<std::string>
parseVarNames(std::istream& is)
{
  std::string line;
  std::getline(is, line);
  return Tokenize(line, std::string(", "));
}

static std::string
parseTitle(std::istream& is)
{
  std::string line;
  std::getline(is, line);
  return line;
}

static void
read_iso(
  const std::string& infile,
  FArrayBox& nodes,
  Vector<int>& faceData,
  int& nElts,
  std::vector<std::string>& names,
  std::string& label)
{
  std::ifstream ifs;
  ifs.open(infile.c_str(), std::ios::in | std::ios::binary);
  if (!ifs.good())
    amrex::Abort("Cannot open MEF file: " + infile);

  label = parseTitle(ifs);
  names = parseVarNames(ifs);
  const int nCompSurf = names.size();

  int nodesPerElt;
  ifs >> nElts;
  ifs >> nodesPerElt;
  ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

  std::string header;
  std::getline(ifs, header);

  auto p1 = header.rfind("((");
  auto p2 = header.find("))", p1);
  std::string box_str = header.substr(p1, p2 - p1 + 2);
  auto first_close = box_str.find(')');
  auto second_open = box_str.find('(', first_close);
  auto second_close = box_str.find(')', second_open);
  std::string stop_tuple = box_str.substr(second_open + 1, second_close - second_open - 1);
  int stop0 = std::stoi(stop_tuple.substr(0, stop_tuple.find(',')));
  const int nNodes = stop0 + 1;

  std::vector<Real> raw(nNodes * nCompSurf);
  ifs.read((char*)raw.data(), sizeof(Real) * raw.size());

  nodes.resize(Box(IntVect::TheZeroVector(),
                   IntVect(AMREX_D_DECL(nNodes - 1, 0, 0))),
               nCompSurf);
  for (int j = 0; j < nCompSurf; ++j) {
    Real* dst = nodes.dataPtr(j);
    for (int i = 0; i < nNodes; ++i)
      dst[i] = raw[i * nCompSurf + j];
  }

  faceData.resize(nElts * nodesPerElt, 0);
  ifs.read((char*)faceData.dataPtr(), sizeof(int) * faceData.size());
}

static int
find_comp(const std::vector<std::string>& names, const std::string& comp)
{
  for (int i = 0; i < (int)names.size(); ++i) {
    if (names[i] == comp)
      return i;
  }
  Print() << "comp: " << comp << std::endl;
  amrex::Abort("Component not found in names.");
  return -1;
}

static Real
elt_center_val(
  const FArrayBox& nodes,
  const Vector<int>& faceData,
  int k,
  int nodesPerElt,
  int comp)
{
  const Real* dat = nodes.dataPtr(comp);
  Real val = 0;
  for (int i = 0; i < nodesPerElt; ++i)
    val += dat[faceData[k * nodesPerElt + i] - 1];
  return val / nodesPerElt;
}

static Real
triangle_area(
  const FArrayBox& nodes, const Vector<int>& faceData, int k, int nodesPerElt)
{
  const Real* xdat = nodes.dataPtr(0);
  const Real* ydat = nodes.dataPtr(1);
#if AMREX_SPACEDIM == 3
  const Real* zdat = nodes.dataPtr(2);
  int n0 = faceData[k * nodesPerElt + 0] - 1;
  int n1 = faceData[k * nodesPerElt + 1] - 1;
  int n2 = faceData[k * nodesPerElt + 2] - 1;
  Real ax = xdat[n1] - xdat[n0], ay = ydat[n1] - ydat[n0],
       az = zdat[n1] - zdat[n0];
  Real bx = xdat[n2] - xdat[n0], by = ydat[n2] - ydat[n0],
       bz = zdat[n2] - zdat[n0];
  return 0.5 * std::sqrt(
                 (ay * bz - az * by) * (ay * bz - az * by) +
                 (az * bx - ax * bz) * (az * bx - ax * bz) +
                 (ax * by - ay * bx) * (ax * by - ay * bx));
#else
  int n0 = faceData[k * nodesPerElt + 0] - 1;
  int n1 = faceData[k * nodesPerElt + 1] - 1;
  return std::sqrt(
    (xdat[n1] - xdat[n0]) * (xdat[n1] - xdat[n0]) +
    (ydat[n1] - ydat[n0]) * (ydat[n1] - ydat[n0]));
#endif
}

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  {
    if (argc < 2)
      print_usage(argc, argv);

    ParmParse pp;

    if (pp.contains("help") || pp.contains("h"))
      print_usage(argc, argv);

    Vector<std::string> infiles;
    pp.getarr("infiles", infiles);
    if (infiles.empty())
      amrex::Abort("No input files specified.");

    std::string xvar, yvar;
    pp.get("xvar", xvar);
    pp.get("yvar", yvar);

    std::string xvar_out = xvar, yvar_out = yvar;
    pp.query("xvar_out", xvar_out);
    pp.query("yvar_out", yvar_out);

    int nBins;
    pp.get("nBins", nBins);

    bool have_xmin = pp.contains("xmin");
    bool have_xmax = pp.contains("xmax");
    bool have_ymin = pp.contains("ymin");
    bool have_ymax = pp.contains("ymax");

    Real xmin = std::numeric_limits<Real>::max();
    Real xmax = -std::numeric_limits<Real>::max();
    Real ymin = std::numeric_limits<Real>::max();
    Real ymax = -std::numeric_limits<Real>::max();

    if (have_xmin) pp.get("xmin", xmin);
    if (have_xmax) pp.get("xmax", xmax);
    if (have_ymin) pp.get("ymin", ymin);
    if (have_ymax) pp.get("ymax", ymax);

    const bool need_pass1 = !have_xmin || !have_xmax || !have_ymin || !have_ymax;

    if (need_pass1) {
      std::cout << "Pass 1: scanning " << infiles.size()
                << " file(s) for variable bounds..." << std::endl;
      for (const auto& f : infiles) {
        FArrayBox nodes;
        Vector<int> faceData;
        int nElts;
        std::vector<std::string> names;
        std::string label;
        read_iso(f, nodes, faceData, nElts, names, label);
        const int nodesPerElt = faceData.size() / nElts;
        const int xcomp = find_comp(names, xvar);
        const int ycomp = find_comp(names, yvar);
        for (int k = 0; k < nElts; ++k) {
          Real xc = elt_center_val(nodes, faceData, k, nodesPerElt, xcomp);
          Real yc = elt_center_val(nodes, faceData, k, nodesPerElt, ycomp);
          if (!have_xmin) xmin = std::min(xmin, xc);
          if (!have_xmax) xmax = std::max(xmax, xc);
          if (!have_ymin) ymin = std::min(ymin, yc);
          if (!have_ymax) ymax = std::max(ymax, yc);
        }
      }
      std::cout << "  " << xvar << " in [" << xmin << ", " << xmax << "]" << std::endl;
      std::cout << "  " << yvar << " in [" << ymin << ", " << ymax << "]" << std::endl;
    }

    const Real dx = (xmax - xmin) / nBins;
    const Real dy = (ymax - ymin) / nBins;

    Vector<Real> hist(nBins * nBins, 0.0);

    std::cout << "Pass 2: accumulating histogram..." << std::endl;
    for (const auto& f : infiles) {
      std::cout << "  [" << f << "]" << std::endl;
      FArrayBox nodes;
      Vector<int> faceData;
      int nElts;
      std::vector<std::string> names;
      std::string label;
      read_iso(f, nodes, faceData, nElts, names, label);
      const int nodesPerElt = faceData.size() / nElts;
      const int xcomp = find_comp(names, xvar);
      const int ycomp = find_comp(names, yvar);
      for (int k = 0; k < nElts; ++k) {
        Real area = triangle_area(nodes, faceData, k, nodesPerElt);
        Real xc = elt_center_val(nodes, faceData, k, nodesPerElt, xcomp);
        Real yc = elt_center_val(nodes, faceData, k, nodesPerElt, ycomp);
        int ix = (int)std::floor((xc - xmin) / dx);
        int iy = (int)std::floor((yc - ymin) / dy);
        ix = std::max(0, std::min(nBins - 1, ix));
        iy = std::max(0, std::min(nBins - 1, iy));
        hist[ix * nBins + iy] += area;
      }
    }

    Real total = 0;
    for (int i = 0; i < nBins * nBins; ++i)
      total += hist[i];
    if (total <= 0)
      amrex::Abort("Histogram is empty — check variable names and file contents.");
    for (int i = 0; i < nBins * nBins; ++i)
      hist[i] /= total;

    const std::string outdir =
      "MEF_JPDFAverage_" + xvar_out + "_" + yvar_out;
    if (!amrex::UtilCreateDirectory(outdir, 0755))
      amrex::Abort("Could not create output directory: " + outdir);

    {
      const std::string fname = outdir + "/Pdf_" + xvar_out + "_x.dat";
      std::cout << "Saving " << xvar_out << " bin centers to " << fname << std::endl;
      FILE* fp = fopen(fname.c_str(), "w");
      for (int i = 0; i < nBins; ++i)
        fprintf(fp, "%e\n", xmin + dx * (0.5 + (Real)i));
      fclose(fp);
    }

    {
      const std::string fname = outdir + "/Pdf_" + yvar_out + "_x.dat";
      std::cout << "Saving " << yvar_out << " bin centers to " << fname << std::endl;
      FILE* fp = fopen(fname.c_str(), "w");
      for (int i = 0; i < nBins; ++i)
        fprintf(fp, "%e\n", ymin + dy * (0.5 + (Real)i));
      fclose(fp);
    }

    {
      const std::string fname =
        outdir + "/Pdf_" + xvar_out + "_" + yvar_out + ".dat";
      std::cout << "Saving jPDF to " << fname << std::endl;
      FILE* fp = fopen(fname.c_str(), "w");
      for (int ix = 0; ix < nBins; ++ix) {
        for (int iy = 0; iy < nBins; ++iy)
          fprintf(fp, "%e ", hist[ix * nBins + iy]);
        fprintf(fp, "\n");
      }
      fclose(fp);
    }

    std::cout << "Done." << std::endl;
  }
  amrex::Finalize();
  return 0;
}
