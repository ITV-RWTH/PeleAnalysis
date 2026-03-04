#include <string>
#include <iostream>
#include <set>
#include <map>
#include <vector>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

using namespace amrex;

static void
print_usage(int, char* argv[])
{
  std::cerr << "Usage:\n"
            << "  " << argv[0] << " infile=FILE outfile=FILE [OPTIONS]\n\n"

            << "Required arguments:\n"
            << "  infile=FILE        MEF file\n"
            << "  outfile=FILE        MEF file\n\n"

            << "Options:\n"
            << "  -h, --help         Show this help message\n\n"
            << "Visit PeleAnalysis/Src/InputSamples for examples or refer to "
            << "the documentation.\n";

  std::exit(1);
}

void read_iso(
  const std::string& infile,
  FArrayBox& nodes,
  Vector<int>& faceData,
  int& nElts,
  std::vector<std::string>& names,
  std::string& label);

void write_iso(
  const std::string& outfile,
  const FArrayBox& nodes,
  const Vector<int>& faceData,
  int nElts,
  const std::vector<std::string>& names,
  const std::string& label);

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

Real
surface_area(
  const FArrayBox& nodes, const Vector<int>& faceData, int nodesPerElt)
{
  Vector<const Real*> dat(AMREX_SPACEDIM);
  for (int i = 0; i < dat.size(); ++i)
    dat[i] = nodes.dataPtr(i);

  AMREX_ASSERT(nodesPerElt == 3); // not general enough for anything else
  int nElts = faceData.size() / nodesPerElt;
  Real Area = 0;
  Vector<Real> p0(AMREX_SPACEDIM), p1(AMREX_SPACEDIM), p2(AMREX_SPACEDIM);
  for (int k = 0; k < nElts; ++k) {
    // set points of
    for (int i = 0; i < 3; ++i) {
      p0[i] = dat[i][faceData[k * nodesPerElt + 0] - 1];
      p1[i] = dat[i][faceData[k * nodesPerElt + 1] - 1];
      p2[i] = dat[i][faceData[k * nodesPerElt + 2] - 1];
    }

    Area +=
      0.5 *
      sqrt(
        pow(
          (p1[1] - p0[1]) * (p2[2] - p0[2]) - (p1[2] - p0[2]) * (p2[1] - p0[1]),
          2)

        +
        pow(
          (p1[2] - p0[2]) * (p2[0] - p0[0]) - (p1[0] - p0[0]) * (p2[2] - p0[2]),
          2)

        +
        pow(
          (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0]),
          2));
  }
  return Area;
}

int
find_comp(const std::vector<std::string>& names, const std::string& comp)
{
  for (int i = 0; i < names.size(); ++i) {
    if (names[i] == comp)
      return i;
  }
  Print() << "comp: " << comp << std::endl;
  amrex::Abort("Component not found in names.");
  return -1; // unreachable, but suppresses compiler warning
}

void
trim_surface(
  const Vector<std::string>& comps,
  const Vector<std::string>& signs,
  const Vector<Real>& vals,
  FArrayBox& nodes,
  Vector<int>& faceData,
  int nodesPerElt,
  std::vector<std::string> names)
{
  const Box& nbox = nodes.box();
  const Real* xdat = nodes.dataPtr(0);
  const Real* ydat = nodes.dataPtr(1);
#if AMREX_SPACEDIM == 3
  const Real* zdat = nodes.dataPtr(2);
#endif

  // Include the radius to the comp logic as comp = -1
  const int nc = comps.size();
  Vector<const Real*> dat(nc);
  for (int i = 0; i < nc; ++i) {
    if (comps[i] == "RXY" || comps[i] == "RXZ" || comps[i] == "RXZ") {
      // RXY, RXZ, and RXZ do not have a corresponding comp.
      // Inserting a nullptr here and catching it later.
      dat[i] = nullptr;
    } else {
      dat[i] = nodes.dataPtr(find_comp(names, comps[i]));
    }
  }
    
  int cnt = 0;
  int cnt_new = 0;
  std::vector<int> nodeIdx;
  for (IntVect iv = nbox.smallEnd(); iv <= nbox.bigEnd();
       nbox.next(iv), cnt++) {
    bool remove_this_node = false;
    for (int i = 0; i < nc; ++i) {
      Real data;
      if (comps[i] == "RXY") {
        Real x = xdat[cnt];
        Real y = ydat[cnt];
        data = std::sqrt(x * x + y * y);
#if AMREX_SPACEDIM == 3
      } else if (comps[i] == "RXZ") {
        Real x = xdat[cnt];
        Real z = zdat[cnt];
        data = std::sqrt(x * x + z * z);
      } else if (comps[i] == "RYZ") {
        Real y = ydat[cnt];
        Real z = zdat[cnt];
        data = std::sqrt(y * y + z * z);
#endif
      } else {
        data = dat[i][cnt];
      }
      if (signs[i] == "lt") {
        remove_this_node |= (data < vals[i]);
      } else if (signs[i] == "le") {
        remove_this_node |= (data <= vals[i]);
      } else if (signs[i] == "gt") {
        remove_this_node |= (data > vals[i]);
      } else if (signs[i] == "ge") {
        remove_this_node |= (data >= vals[i]);
      } else if (signs[i] == "eq") {
        remove_this_node |= (data == vals[i]);
      } else {
        Print() << "i,signs[i] " << i << "," << signs[i] << std::endl;
        amrex::Abort("Bad signs data. Use one of [lt,le,gt,ge,eq]");
      }
    }

    nodeIdx.push_back(remove_this_node ? -1 : cnt_new++);
  }
  int nNodesNEW = cnt_new;

  // Build new nodes structure with bad points removed
  Box newBox(IntVect::TheZeroVector(), (nNodesNEW - 1) * amrex::BASISV(0));
  int nNodesOLD = nbox.numPts();
  int nComp = nodes.nComp();
  FArrayBox newNodes(newBox, nComp);
  Vector<Real*> odp(nComp);
  Vector<Real*> ndp(nComp);
  for (int j = 0; j < nComp; ++j) {
    odp[j] = nodes.dataPtr(j);
    ndp[j] = newNodes.dataPtr(j);
  }
  cnt_new = 0;
  for (int i = 0; i < nNodesOLD; ++i) {
    int newNode = nodeIdx[i];
    if (newNode >= 0)
      for (int j = 0; j < nComp; ++j)
        ndp[j][newNode] = odp[j][i];
  }
  nodes.resize(newBox, nComp);
  nodes.copy(newNodes);
  newNodes.clear(); // Make some space...not really necessary, but what the heck

  const int nElts = faceData.size() / nodesPerElt;
  AMREX_ASSERT(nElts * nodesPerElt == faceData.size()); // Idiot check

  // Remove elements that refer to removed nodes
  std::vector<int> newFaceData;
  cnt_new = 0;
  for (int i = 0; i < nElts; ++i) {
    int offset = nodesPerElt * i;
    bool eltGood = true;
    for (int j = 0; j < nodesPerElt; ++j)
      eltGood &=
        (nodeIdx[faceData[offset + j] - 1] >=
         0); // Remember that faceData is 1-based

    if (eltGood) {
      int new_offset = nodesPerElt * cnt_new;
      for (int j = 0; j < nodesPerElt; ++j)
        newFaceData.push_back(
          nodeIdx[faceData[offset + j] - 1] +
          1); // Make sure new faceData is 1-based
      cnt_new++;
    }
  }

  faceData.resize(cnt_new * nodesPerElt);
  for (int i = 0; i < faceData.size(); ++i)
    faceData[i] = newFaceData[i];
}



void
remove_unused_nodes(FArrayBox& nodes, Vector<int>& faceData, int nodesPerElt)
{
  const Box& nbox = nodes.box();
  int nNodesOLD = nbox.numPts();

  std::vector<bool> removeNode(nNodesOLD, true);
  for (int i = 0; i < faceData.size(); ++i) {
    int node = faceData[i]; // These are 1-based node numbers
    AMREX_ASSERT(node <= nNodesOLD && node > 0);
    removeNode[node - 1] = false; // Indexing here is zero-based
  }

  Vector<int> nodeIdx(nNodesOLD);
  int nNodesNEW = 0;
  for (int i = 0; i < nNodesOLD; ++i) {
    nodeIdx[i] = (removeNode[i] ? -1 : nNodesNEW++); // zero-based
  }

  if (nNodesNEW != nNodesOLD) {
    // Build new nodes structure with bad points removed
    Box newBox(IntVect::TheZeroVector(), (nNodesNEW - 1) * amrex::BASISV(0));
    int nComp = nodes.nComp();
    FArrayBox newNodes(newBox, nComp);
    Vector<Real*> odp(nComp);
    Vector<Real*> ndp(nComp);
    for (int j = 0; j < nComp; ++j) {
      odp[j] = nodes.dataPtr(j);
      ndp[j] = newNodes.dataPtr(j);
    }
    int cnt_new = 0;
    for (int i = 0; i < nNodesOLD; ++i) {
      int newNode = nodeIdx[i];
      if (newNode >= 0)
        for (int j = 0; j < nComp; ++j)
          ndp[j][newNode] = odp[j][i];
    }
    nodes.resize(newBox, nComp);
    nodes.copy(newNodes);
    newNodes
      .clear(); // Make some space...not really necessary, but what the heck

    const int nElts = faceData.size() / nodesPerElt;
    AMREX_ASSERT(nElts * nodesPerElt == faceData.size()); // Idiot check

    // Redefine elements with new numbering
    std::vector<int> newFaceData;
    cnt_new = 0;
    for (int i = 0; i < nElts; ++i) {
      int offset = nodesPerElt * i;
      int new_offset = nodesPerElt * cnt_new;
      for (int j = 0; j < nodesPerElt; ++j) {
        if (
          nodeIdx[faceData[offset + j] - 1] < 0 ||
          nodeIdx[faceData[offset + j] - 1] >= nNodesOLD)
          Print() << "messed up" << std::endl;

        newFaceData.push_back(
          nodeIdx[faceData[offset + j] - 1] +
          1); // Make sure new faceData is 1-based
      }
      cnt_new++;
    }

    faceData.resize(cnt_new * nodesPerElt);
    for (int i = 0; i < faceData.size(); ++i)
      faceData[i] = newFaceData[i];

    if (ParallelDescriptor::IOProcessor())
      Print() << "Removed " << nNodesOLD - nNodesNEW << " unusued nodes"
                << std::endl;
  }
}

void
triangle_area(
  const Real* x,
  const Real* y,
  const Real* z,
  const Vector<int>& faceData,
  int nElts,
  int nNodes,
  Real* area)
{
  Real R1[3], R2[3], R3[3];
  const Real* loc[3];
  loc[0] = x;
  loc[1] = y;
  loc[2] = z;

  for (int i = 0; i < nElts; ++i) {
    int offset = i * 3;
    int pa = faceData[offset + 0] - 1;
    int pb = faceData[offset + 1] - 1;
    int pc = faceData[offset + 2] - 1;
    for (int j = 0; j < 3; ++j) {
      R1[j] = loc[j][pb] - loc[j][pa];
      R2[j] = loc[j][pc] - loc[j][pa];
    }
    R3[0] = R1[1] * R2[2] - R2[1] * R1[2];
    R3[1] = R1[2] * R2[0] - R2[2] * R1[0];
    R3[2] = R1[0] * R2[1] - R2[0] * R1[1];

    Real res = 0;
    for (int j = 0; j < 3; ++j)
      res += R3[j] * R3[j];
    area[i] = 0.5 * std::sqrt(res);
  }
}

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);

  if (argc < 2) {
    print_usage(argc, argv);
  } else if (
    (std::strcmp(argv[1], "-h") == 0) ||
    (std::strcmp(argv[1], "--help") == 0)) {
    print_usage(argc, argv);
  }

  ParmParse pp;

  int nElts;
  std::string infile;
  pp.get("infile", infile);
  std::string outfile;
  pp.get("outfile", outfile);

  FArrayBox nodes;
  Vector<int> faceData;
  std::vector<std::string> names;
  std::string label;
  read_iso(infile, nodes, faceData, nElts, names, label);
  int nodesPerElt = faceData.size() / nElts;
  AMREX_ASSERT(nodesPerElt * nElts == faceData.size());

  Print() << "Surface area before: " << surface_area(nodes, faceData, nodesPerElt)
       << '\n';

  int nc = pp.countval("comps");
  
  Vector<std::string> comps(nc);
  Vector<std::string> signs(nc);
  Vector<Real> vals(nc);
  if (nc > 0) {
    comps.resize(nc);
    pp.getarr("comps", comps, 0, nc);

    int ns = pp.countval("signs");
    AMREX_ALWAYS_ASSERT(ns == nc);
    pp.getarr("signs", signs, 0, nc);

    int nv = pp.countval("vals");
    AMREX_ALWAYS_ASSERT(nv == nc);
    pp.getarr("vals", vals, 0, nc);
  } else {
    amrex::Abort("No triming tarfet in inputs");
  }

  trim_surface(comps, signs, vals, nodes, faceData, nodesPerElt, names);
  
  Print() << "Surface area after: " << surface_area(nodes, faceData, nodesPerElt)
       << '\n';

  nElts = faceData.size() / nodesPerElt;
  AMREX_ASSERT(nElts * nodesPerElt == faceData.size());

  Vector<int> remComps;
  int nrc = pp.countval("remComps");
  if (nrc > 0) {
    remComps.resize(nrc);
    pp.getarr("remComps", remComps, 0, nrc);

    int nCompsNew = names.size() - nrc;
    FArrayBox newNodes(nodes.box(), nCompsNew);
    int cnt = 0;
    Vector<std::string> newNames(nCompsNew);
    for (int i = 0; i < names.size(); ++i) {
      bool keepThisComp = true;
      for (int j = 0; j < remComps.size(); ++j) {
        if (remComps[j] == i) {
          keepThisComp = false;
        }
      }
      if (keepThisComp) {
        newNames[cnt] = names[i];
        newNodes.copy(nodes, i, cnt, 1);
        cnt++;
        AMREX_ASSERT(cnt < nCompsNew);
      }
    }

    nodes.clear();
    names.clear();
    names = newNames;
    nodes.resize(newNodes.box(), newNodes.nComp());
    nodes.copy(newNodes);
  }

  bool do_area_stats = false;
  pp.query("do_area_stats", do_area_stats);
  if (do_area_stats) {
    Vector<Real> area(nElts);

    int idX = 0;
    int idY = 1;
    int idZ = 2;

    int nNodes = nodes.box().numPts();
    triangle_area(
      nodes.dataPtr(idX), nodes.dataPtr(idY), nodes.dataPtr(idZ), faceData,
      nElts, nNodes, area.dataPtr());

    // Find min, max
    Real area_min = area[0];
    Real area_max = area[0];
    for (int i = 0; i < area.size(); ++i) {
      area_min = std::min(area_min, area[i]);
      area_max = std::max(area_max, area[i]);
    }

    Print() << "  Triangle area min, max: " << area_min << " , " << area_max
         << '\n';
  }

  remove_unused_nodes(nodes, faceData, nodesPerElt);

  write_iso(outfile, nodes, faceData, nElts, names, label);

  amrex::Finalize();
  return 0;
}

void
write_iso(
  const std::string& outfile,
  const FArrayBox& nodes,
  const Vector<int>& faceData,
  int nElts,
  const std::vector<std::string>& names,
  const std::string& label)
{
  // Rotate data to vary quickest on component
  int nCompSurf = nodes.nComp();
  int nNodes = nodes.box().numPts();

  FArrayBox tnodes(nodes.box(), nCompSurf);
  const Real** np = new const Real*[nCompSurf];
  for (int j = 0; j < nCompSurf; ++j)
    np[j] = nodes.dataPtr(j);
  Real* ndat = tnodes.dataPtr();
  for (int i = 0; i < nNodes; ++i) {
    for (int j = 0; j < nCompSurf; ++j) {
      ndat[j] = np[j][i];
    }
    ndat += nCompSurf;
  }
  delete[] np;

  std::ofstream ofs;
  ofs.open(outfile.c_str(), std::ios::out | std::ios::trunc | std::ios::binary);
  ofs << label << std::endl;
  for (int i = 0; i < nCompSurf; ++i) {
    ofs << names[i];
    if (i < nCompSurf - 1)
      ofs << " ";
    else
      ofs << std::endl;
  }
  int nodesPerElt = faceData.size() / nElts;
  ofs << nElts << " " << nodesPerElt << std::endl;
  tnodes.writeOn(ofs);
  ofs.write((char*)faceData.dataPtr(), sizeof(int) * faceData.size());
  ofs.close();
}

void
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
  label = parseTitle(ifs);
  names = parseVarNames(ifs);
  const int nCompSurf = names.size();

  int nodesPerElt;
  ifs >> nElts;
  ifs >> nodesPerElt;

  FArrayBox tnodes;
  tnodes.readFrom(ifs);
  const int nNodes = tnodes.box().numPts();

  // "rotate" the data so that the components are 'in the right spot for fab
  // data'
  nodes.resize(tnodes.box(), nCompSurf);
  Real** np = new Real*[nCompSurf];
  for (int j = 0; j < nCompSurf; ++j)
    np[j] = nodes.dataPtr(j);

  Real* ndat = tnodes.dataPtr();
  for (int i = 0; i < nNodes; ++i) {
    for (int j = 0; j < nCompSurf; ++j) {
      np[j][i] = ndat[j];
    }
    ndat += nCompSurf;
  }
  delete[] np;
  tnodes.clear();

  faceData.resize(nElts * nodesPerElt, 0);
  ifs.read((char*)faceData.dataPtr(), sizeof(int) * faceData.size());
}
