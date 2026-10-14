#include <algorithm>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_Reduce.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

using namespace amrex;

// Maximum number of table parameters (dimensions) supported
constexpr int MAX_TABLE_DIM = 4;

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "Looks up tabulated quantities from an OpenFOAM-format chemistry table "
       "(as\nwritten by the chemtable_generator) at the state given by "
       "plotfile\nvariables, by multilinear interpolation.\n\n"
    << "Usage:\n"
    << "  " << argv[0]
    << " infile=FILE tableDir=DIR tables=NAMES coordVars=NAMES [OPTIONS]\n\n"

    << "Required arguments:\n"
    << "  infile=FILE              AMReX plotfile\n"
    << "  tableDir=DIR             directory holding tableProperties and the\n"
    << "                           <name><tableSuffix> files\n"
    << "  tables=NAMES             tables to look up, e.g. DC Dz_filt\n"
    << "  coordVars=NAMES          one plotfile variable per table parameter,\n"
    << "                           in the order of 'params' in "
       "tableProperties\n\n"

    << "Options:\n"
    << "  outfile=FILE             output plotfile [DEF: <infile>_table]\n"
    << "  tableSuffix=STR          table file suffix [DEF: _prem_table]\n"
    << "  outPrefix=STR            prefix of the output names [DEF: none]\n"
    << "  keepVars=NAMES           variables copied unchanged to the output\n"
    << "  finestLevel=N            finest AMR level to process [DEF: "
       "finest]\n"
    << "  n_files=N                max number of plotfile data files\n"
    << "  -h, --help               show this help message\n\n"
    << "States outside the table range are clamped to its boundary; the\n"
    << "fraction of clamped cells is reported per parameter.\n\n"

    << "Visit PeleAnalysis/Src/InputsSamples for examples or refer to "
    << "the documentation.\n";

  std::exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile, std::string("/"));
  return tokens[tokens.size() - 1];
}

// Split an OpenFOAM dictionary into tokens: comments are dropped and the
// characters ( ) { } ; are tokens of their own
static std::vector<std::string>
tokenizeFoam(const std::string& fileName)
{
  std::ifstream ifs(fileName);
  if (!ifs.good()) {
    Abort("Cannot open table file " + fileName);
  }
  std::stringstream buffer;
  buffer << ifs.rdbuf();
  const std::string text = buffer.str();

  std::vector<std::string> tokens;
  std::string cur;
  auto flush = [&]() {
    if (!cur.empty()) {
      tokens.push_back(cur);
      cur.clear();
    }
  };
  for (size_t i = 0; i < text.size(); ++i) {
    const char c = text[i];
    if (c == '/' && i + 1 < text.size() && text[i + 1] == '/') {
      flush();
      while (i < text.size() && text[i] != '\n') {
        ++i;
      }
    } else if (c == '/' && i + 1 < text.size() && text[i + 1] == '*') {
      flush();
      const size_t end = text.find("*/", i + 2);
      i = (end == std::string::npos) ? text.size() : end + 1;
    } else if (std::isspace(static_cast<unsigned char>(c))) {
      flush();
    } else if (c == '(' || c == ')' || c == '{' || c == '}' || c == ';') {
      flush();
      tokens.emplace_back(1, c);
    } else {
      cur += c;
    }
  }
  flush();
  return tokens;
}

static bool
isInteger(const std::string& s)
{
  return !s.empty() && std::all_of(s.begin(), s.end(), ::isdigit);
}

// Position of the first token equal to key that starts a list (is followed by
// an integer count)
static size_t
findList(
  const std::vector<std::string>& tokens,
  const std::string& key,
  const std::string& fileName)
{
  for (size_t i = 0; i + 1 < tokens.size(); ++i) {
    if (tokens[i] == key && isInteger(tokens[i + 1])) {
      return i + 1;
    }
  }
  Abort("List " + key + " not found in " + fileName);
  return 0;
}

static void
expectToken(
  const std::vector<std::string>& tokens,
  size_t& pos,
  const std::string& tok,
  const std::string& fileName)
{
  if (pos >= tokens.size() || tokens[pos] != tok) {
    Abort(
      "Malformed table file " + fileName + ": expected '" + tok + "', found '" +
      (pos < tokens.size() ? tokens[pos] : "EOF") + "'");
  }
  ++pos;
}

// Read a (nested) list "N ( ... )" starting at pos into values, checking the
// count of every nesting level against shape
static void
readNestedList(
  const std::vector<std::string>& tokens,
  size_t& pos,
  const std::vector<int>& shape,
  int depth,
  std::vector<Real>& values,
  const std::string& fileName)
{
  if (pos >= tokens.size() || !isInteger(tokens[pos])) {
    Abort("Malformed table file " + fileName + ": expected a list count");
  }
  const int n = std::stoi(tokens[pos++]);
  if (n != shape[depth]) {
    Abort(
      "Table " + fileName + " has " + std::to_string(n) +
      " entries along parameter " + std::to_string(depth) + ", expected " +
      std::to_string(shape[depth]));
  }
  expectToken(tokens, pos, "(", fileName);
  for (int i = 0; i < n; ++i) {
    if (depth == static_cast<int>(shape.size()) - 1) {
      values.push_back(std::stod(tokens[pos++]));
    } else {
      readNestedList(tokens, pos, shape, depth + 1, values, fileName);
    }
  }
  expectToken(tokens, pos, ")", fileName);
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

    if (pp.contains("verbose")) {
      AmrData::SetVerbose(true);
    }

    std::string infile;
    pp.get("infile", infile);
    std::string outfile(getFileRoot(infile) + "_table");
    pp.query("outfile", outfile);

    std::string tableDir;
    pp.get("tableDir", tableDir);
    std::string tableSuffix = "_prem_table";
    pp.query("tableSuffix", tableSuffix);
    std::string outPrefix;
    pp.query("outPrefix", outPrefix);

    const int nTables = pp.countval("tables");
    if (nTables < 1) {
      Abort("At least one table must be given with tables");
    }
    Vector<std::string> tables(nTables);
    pp.getarr("tables", tables, 0, nTables);

    const int nCoord = pp.countval("coordVars");
    Vector<std::string> coordVars(nCoord);
    if (nCoord > 0) {
      pp.getarr("coordVars", coordVars, 0, nCoord);
    }

    const int nKeep = pp.countval("keepVars");
    Vector<std::string> keepVars(nKeep);
    if (nKeep > 0) {
      pp.getarr("keepVars", keepVars, 0, nKeep);
    }

    // ---------------------------------------------------------------------
    // Read the table on the I/O rank and broadcast it
    // ---------------------------------------------------------------------
    int nDim = 0;
    Vector<int> shape;
    Vector<std::string> params;
    Vector<Real> axes; // all parameter axes, one after the other
    Vector<Real> data; // all tables, one after the other
    if (ParallelDescriptor::IOProcessor()) {
      const std::string propFile = tableDir + "/tableProperties";
      const auto propTokens = tokenizeFoam(propFile);
      size_t pos = findList(propTokens, "params", propFile);
      nDim = std::stoi(propTokens[pos++]);
      expectToken(propTokens, pos, "(", propFile);
      for (int d = 0; d < nDim; ++d) {
        std::string name = propTokens[pos++];
        name.erase(std::remove(name.begin(), name.end(), '"'), name.end());
        params.push_back(name);
      }
      expectToken(propTokens, pos, ")", propFile);
      for (int d = 0; d < nDim; ++d) {
        pos = findList(propTokens, params[d], propFile);
        const int n = std::stoi(propTokens[pos]);
        std::vector<Real> axis;
        readNestedList(propTokens, pos, {n}, 0, axis, propFile);
        for (size_t i = 1; i < axis.size(); ++i) {
          if (axis[i] <= axis[i - 1]) {
            Abort("Axis " + params[d] + " is not strictly increasing");
          }
        }
        shape.push_back(n);
        axes.insert(axes.end(), axis.begin(), axis.end());
      }

      for (const auto& table : tables) {
        const std::string fileName = tableDir + "/" + table + tableSuffix;
        const auto tokens = tokenizeFoam(fileName);
        pos = findList(tokens, table + tableSuffix, fileName);
        std::vector<Real> values;
        readNestedList(
          tokens, pos, std::vector<int>(shape.begin(), shape.end()), 0, values,
          fileName);
        data.insert(data.end(), values.begin(), values.end());
        Print() << "Read table " << fileName << "\n";
      }
    }
    ParallelDescriptor::Bcast(
      &nDim, 1, ParallelDescriptor::IOProcessorNumber());
    shape.resize(nDim);
    ParallelDescriptor::Bcast(
      shape.data(), nDim, ParallelDescriptor::IOProcessorNumber());
    BroadcastStringArray(
      params, ParallelDescriptor::MyProc(),
      ParallelDescriptor::IOProcessorNumber(),
      ParallelDescriptor::Communicator());
    Long tableSize = 1;
    int nAxes = 0;
    for (int d = 0; d < nDim; ++d) {
      tableSize *= shape[d];
      nAxes += shape[d];
    }
    axes.resize(nAxes);
    data.resize(tableSize * nTables);
    ParallelDescriptor::Bcast(
      axes.data(), nAxes, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::Bcast(
      data.data(), data.size(), ParallelDescriptor::IOProcessorNumber());

    if (nDim > MAX_TABLE_DIM) {
      Abort(
        "Tables with more than " + std::to_string(MAX_TABLE_DIM) +
        " parameters are not supported");
    }
    if (nCoord != nDim) {
      std::string msg = "coordVars must name one variable per table "
                        "parameter:";
      for (const auto& p : params) {
        msg += " " + p;
      }
      Abort(msg);
    }
    Print() << "Table parameters:\n";
    {
      int off = 0;
      for (int d = 0; d < nDim; ++d) {
        Print() << "  " << params[d] << " <- " << coordVars[d] << ": "
                << shape[d] << " points in [" << axes[off] << ", "
                << axes[off + shape[d] - 1] << "]\n";
        off += shape[d];
      }
    }

    // Copy the table to the device
    Gpu::DeviceVector<Real> axes_d(nAxes);
    Gpu::DeviceVector<Real> data_d(data.size());
    Gpu::copy(Gpu::hostToDevice, axes.begin(), axes.end(), axes_d.begin());
    Gpu::copy(Gpu::hostToDevice, data.begin(), data.end(), data_d.begin());
    GpuArray<int, MAX_TABLE_DIM> nPts{};
    GpuArray<int, MAX_TABLE_DIM> axisOff{};
    GpuArray<Long, MAX_TABLE_DIM> stride{};
    {
      int off = 0;
      Long s = tableSize;
      for (int d = 0; d < nDim; ++d) {
        nPts[d] = shape[d];
        axisOff[d] = off;
        off += shape[d];
        s /= shape[d];
        stride[d] = s;
      }
    }

    // ---------------------------------------------------------------------
    // Plotfile
    // ---------------------------------------------------------------------
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices dataServices(infile, fileType);
    if (!dataServices.AmrDataOk()) {
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel", finestLevel);
    finestLevel = std::min(finestLevel, amrData.FinestLevel());
    const int Nlev = finestLevel + 1;

    // Input: coordinates, then the carried variables
    Vector<std::string> inNames = coordVars;
    for (const auto& var : keepVars) {
      inNames.push_back(var);
    }
    for (const auto& var : inNames) {
      if (amrData.StateNumber(var) < 0) {
        Abort("Variable " + var + " not found in file " + infile);
      }
    }
    const int nCompIn = inNames.size();
    Vector<int> destFillComps(nCompIn);
    for (int i = 0; i < nCompIn; ++i) {
      destFillComps[i] = i;
    }

    // Output: carried variables, then the tables
    const int nCompOut = nKeep + nTables;
    Vector<std::string> outNames = keepVars;
    for (const auto& table : tables) {
      outNames.push_back(outPrefix + table);
    }

    Vector<MultiFab> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);
    RealBox rb(&(amrData.ProbLo()[0]), &(amrData.ProbHi()[0]));
    int coord = 0;
    Vector<int> is_per(AMREX_SPACEDIM, 0);

    const Real* axes_p = axes_d.data();
    const Real* data_p = data_d.data();

    for (int lev = 0; lev < Nlev; ++lev) {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      geoms[lev] =
        Geometry(amrData.ProbDomain()[lev], &rb, coord, &(is_per[0]));

      MultiFab indata(ba, dm, nCompIn, 0);
      outdata[lev].define(ba, dm, nCompOut, 0);

      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata, lev, inNames, destFillComps);

      for (int n = 0; n < nKeep; ++n) {
        MultiFab::Copy(outdata[lev], indata, nDim + n, n, 1, 0);
      }

      // Report the states outside the table
      for (int d = 0; d < nDim; ++d) {
        const Real lo = axes[axisOff[d]];
        const Real hi = axes[axisOff[d] + nPts[d] - 1];
        ReduceOps<ReduceOpSum, ReduceOpSum> reduce_op;
        ReduceData<Long, Long> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        for (MFIter mfi(indata, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.tilebox();
          auto const& in = indata.const_array(mfi);
          reduce_op.eval(
            bx, reduce_data,
            [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple {
              const Real v = in(i, j, k, d);
              return {static_cast<Long>(v < lo), static_cast<Long>(v > hi)};
            });
        }
        auto counts = reduce_data.value();
        Long below = amrex::get<0>(counts);
        Long above = amrex::get<1>(counts);
        ParallelDescriptor::ReduceLongSum(below);
        ParallelDescriptor::ReduceLongSum(above);
        const Real nCells = static_cast<Real>(ba.numPts());
        Print() << "  " << coordVars[d] << ": " << 100.0 * below / nCells
                << " % of the cells below, " << 100.0 * above / nCells
                << " % above the table range\n";
      }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(outdata[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto const& in = indata.const_array(mfi);
        auto const& out = outdata[lev].array(mfi);
        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          // Lower bracketing index and weight of the upper point per
          // parameter, clamped to the table
          GpuArray<int, MAX_TABLE_DIM> idx{};
          GpuArray<Real, MAX_TABLE_DIM> w{};
          for (int d = 0; d < nDim; ++d) {
            const Real* ax = axes_p + axisOff[d];
            const int n = nPts[d];
            const Real v = in(i, j, k, d);
            if (n == 1 || v <= ax[0]) {
              idx[d] = 0;
              w[d] = 0.0;
            } else if (v >= ax[n - 1]) {
              idx[d] = n - 2;
              w[d] = 1.0;
            } else {
              int lo = 0;
              int hi = n - 1;
              while (hi - lo > 1) {
                const int mid = (lo + hi) / 2;
                if (ax[mid] <= v) {
                  lo = mid;
                } else {
                  hi = mid;
                }
              }
              idx[d] = lo;
              w[d] = (v - ax[lo]) / (ax[lo + 1] - ax[lo]);
            }
          }
          for (int t = 0; t < nTables; ++t) {
            const Real* tab = data_p + t * tableSize;
            Real val = 0.0;
            for (int corner = 0; corner < (1 << nDim); ++corner) {
              Real weight = 1.0;
              Long offset = 0;
              for (int d = 0; d < nDim; ++d) {
                const int up = (corner >> d) & 1;
                weight *= up ? w[d] : (1.0 - w[d]);
                offset += (idx[d] + (nPts[d] > 1 ? up : 0)) * stride[d];
              }
              if (weight != 0.0) {
                val += weight * tab[offset];
              }
            }
            out(i, j, k, nKeep + t) = val;
          }
        });
      }
      Print() << "Table lookup finished for level " << lev << std::endl;
    }

    // Cap the number of plotfile data files via the n_files option (AMReX)
    int n_files = amrex::VisMF::GetNOutFiles();
    pp.query("n_files", n_files);
    amrex::VisMF::SetNOutFiles(n_files);

    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev - 1);
    for (int lev = 0; lev < Nlev - 1; ++lev) {
      int r = amrData.RefRatio()[lev];
      refRatios[lev] = IntVect(AMREX_D_DECL(r, r, r));
    }
    amrex::WriteMultiLevelPlotfile(
      outfile, Nlev, GetVecOfConstPtrs(outdata), outNames, geoms,
      amrData.Time(), isteps, refRatios);
  }
  Finalize();
  return 0;
}
