#include <string>
#include <iostream>
#include <map>
#include <memory>
#include <regex>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_Parser.H>

using namespace amrex;

// Maximum number of distinct variables a single expression may reference
constexpr int MAX_EXPR_VARS = 16;

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "Evaluates user-defined expressions of plotfile variables and writes "
       "the\nresults to a new plotfile. Expressions are evaluated in the given "
       "order\nand may use the results of earlier expressions.\n\n"
    << "Usage:\n"
    << "  " << argv[0]
    << " infiles=FILE [FILE ...] outVars=NAME [NAME ...] "
       "expr.NAME=EXPR [OPTIONS]\n\n"

    << "Required arguments:\n"
    << "  infiles=FILE [FILE ...]         AMReX plotfile(s) on the same "
       "domain\n"
    << "  outVars=NAME [NAME ...]         names of the variables to compute\n"
    << "  expr.NAME=EXPR                  expression for each outVars NAME\n\n"

    << "Options:\n"
    << "  outfile=FILE                    output plotfile [DEF: "
       "<infiles[0]>_eval]\n"
    << "  alias=SYM:VAR [SYM:VAR ...]     expression symbol SYM for plotfile\n"
    << "                                  variable VAR, needed for names "
       "that\n"
    << "                                  are not identifiers, e.g. "
       "YH2:Y(H2)\n"
    << "  constants=SYM:VAL [...]         named constants (pi is predefined)\n"
    << "  keepVars=NAME [NAME ...]        input variables copied to the "
       "output\n"
    << "  writeVars=NAME [NAME ...]       subset of outVars to write [DEF: "
       "all]\n"
    << "  The cell-center coordinates are available as x, y, z unless a\n"
    << "  variable of that name exists.\n"
    << "  finestLevel=N                   finest AMR level to process [DEF: "
       "finest]\n"
    << "  n_files=N                       max number of plotfile data files\n"
    << "  -h, --help                      show this help message\n\n"

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

static bool
isIdentifier(const std::string& s)
{
  static const std::regex re("[A-Za-z_][A-Za-z0-9_]*");
  return std::regex_match(s, re);
}

static std::pair<std::string, std::string>
splitPair(const std::string& entry, const std::string& key)
{
  const auto pos = entry.find(':');
  if (pos == std::string::npos || pos == 0 || pos == entry.size() - 1) {
    Abort(key + " entry '" + entry + "' must have the form SYM:VALUE");
  }
  return {entry.substr(0, pos), entry.substr(pos + 1)};
}

// Where the value of an expression symbol is read from
struct SymbolSource
{
  enum Kind { Input, Output, Coord };
  Kind kind = Input; // plotfile variable, earlier result or coordinate
  int comp = -1;     // component in the input, output or coordinate data
};

template <int N>
void
evalExpression(
  Parser const& parser,
  MultiFab& out,
  int outComp,
  const Vector<const MultiFab*>& srcMF,
  const Vector<int>& srcComp)
{
  auto const exe = parser.compile<N>();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(out, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    auto const& res = out.array(mfi);
    GpuArray<Array4<Real const>, N> args;
    GpuArray<int, N> comps;
    for (int n = 0; n < N; ++n) {
      args[n] = srcMF[n]->const_array(mfi);
      comps[n] = srcComp[n];
    }
    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      GpuArray<double, N> v;
      for (int n = 0; n < N; ++n) {
        v[n] = args[n](i, j, k, comps[n]);
      }
      res(i, j, k, outComp) = exe(v);
    });
  }
}

template <int N>
void
dispatchExpression(
  int nVars,
  Parser const& parser,
  MultiFab& out,
  int outComp,
  const Vector<const MultiFab*>& srcMF,
  const Vector<int>& srcComp)
{
  if constexpr (N > MAX_EXPR_VARS) {
    Abort(
      "An expression may reference at most " + std::to_string(MAX_EXPR_VARS) +
      " variables");
  } else {
    if (nVars == N) {
      evalExpression<N>(parser, out, outComp, srcMF, srcComp);
    } else {
      dispatchExpression<N + 1>(nVars, parser, out, outComp, srcMF, srcComp);
    }
  }
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

    // ---------------------------------------------------------------------
    // Input files
    // ---------------------------------------------------------------------
    const int nFiles = pp.countval("infiles");
    if (nFiles < 1) {
      Abort("At least one file must be given with infiles");
    }
    Vector<std::string> infiles(nFiles);
    pp.getarr("infiles", infiles, 0, nFiles);

    std::string outfile(getFileRoot(infiles[0]) + "_eval");
    pp.query("outfile", outfile);

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    Vector<std::unique_ptr<DataServices>> dataServices(nFiles);
    for (int f = 0; f < nFiles; ++f) {
      dataServices[f] = std::make_unique<DataServices>(infiles[f], fileType);
      if (!dataServices[f]->AmrDataOk()) {
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
      }
    }
    AmrData& amrData0 = dataServices[0]->AmrDataRef();

    int finestLevel = amrData0.FinestLevel();
    pp.query("finestLevel", finestLevel);
    for (int f = 1; f < nFiles; ++f) {
      AmrData& amrData = dataServices[f]->AmrDataRef();
      finestLevel = std::min(finestLevel, amrData.FinestLevel());
      for (int lev = 0; lev <= finestLevel; ++lev) {
        if (amrData.ProbDomain()[lev] != amrData0.ProbDomain()[lev]) {
          Abort(
            "Domain of " + infiles[f] + " differs from " + infiles[0] +
            " at level " + std::to_string(lev));
        }
        if (amrData.boxArray(lev) != amrData0.boxArray(lev)) {
          Print() << "WARNING: grids of " << infiles[f] << " differ from "
                  << infiles[0] << " at level " << lev
                  << "; its data is interpolated onto the grids of "
                  << infiles[0] << "\n";
        }
      }
    }
    const int Nlev = finestLevel + 1;

    // ---------------------------------------------------------------------
    // Symbols: plotfile variables, aliases, constants
    // ---------------------------------------------------------------------
    // Symbol -> (file index, plotfile variable name), first file wins
    std::map<std::string, std::pair<int, std::string>> inputSymbols;
    for (int f = 0; f < nFiles; ++f) {
      for (const auto& name : dataServices[f]->AmrDataRef().PlotVarNames()) {
        if (isIdentifier(name) && inputSymbols.count(name) == 0) {
          inputSymbols[name] = {f, name};
        }
      }
    }
    auto findVariable = [&](const std::string& var) -> int {
      for (int f = 0; f < nFiles; ++f) {
        const auto& names = dataServices[f]->AmrDataRef().PlotVarNames();
        if (std::find(names.begin(), names.end(), var) != names.end()) {
          return f;
        }
      }
      return -1;
    };

    const int nAlias = pp.countval("alias");
    for (int n = 0; n < nAlias; ++n) {
      std::string entry;
      pp.get("alias", entry, n);
      const auto [sym, var] = splitPair(entry, "alias");
      if (!isIdentifier(sym)) {
        Abort("alias symbol '" + sym + "' is not a valid identifier");
      }
      const int f = findVariable(var);
      if (f < 0) {
        Abort("alias: variable " + var + " not found in the input files");
      }
      inputSymbols[sym] = {f, var};
    }

    std::map<std::string, Real> constants{{"pi", M_PI}};
    const int nConst = pp.countval("constants");
    for (int n = 0; n < nConst; ++n) {
      std::string entry;
      pp.get("constants", entry, n);
      const auto [sym, val] = splitPair(entry, "constants");
      if (!isIdentifier(sym)) {
        Abort("constant name '" + sym + "' is not a valid identifier");
      }
      constants[sym] = std::stod(val);
    }

    // ---------------------------------------------------------------------
    // Expressions
    // ---------------------------------------------------------------------
    const int nOut = pp.countval("outVars");
    if (nOut < 1) {
      Abort("At least one variable must be given with outVars");
    }
    Vector<std::string> outVars(nOut);
    pp.getarr("outVars", outVars, 0, nOut);
    for (const auto& name : outVars) {
      if (!isIdentifier(name)) {
        Abort("outVars name '" + name + "' is not a valid identifier");
      }
    }

    const int nKeep = pp.countval("keepVars");
    Vector<std::string> keepVars(nKeep);
    if (nKeep > 0) {
      pp.getarr("keepVars", keepVars, 0, nKeep);
    }

    // Input variables to load: (file, name), keepVars first
    Vector<std::pair<int, std::string>> loadVars;
    auto loadIndex = [&](int f, const std::string& var) -> int {
      for (int n = 0; n < loadVars.size(); ++n) {
        if (loadVars[n].first == f && loadVars[n].second == var) {
          return n;
        }
      }
      loadVars.push_back({f, var});
      return static_cast<int>(loadVars.size()) - 1;
    };
    Vector<int> keepIdx;
    for (const auto& var : keepVars) {
      const int f = findVariable(var);
      if (f < 0) {
        Abort("keepVars: variable " + var + " not found in the input files");
      }
      keepIdx.push_back(loadIndex(f, var));
    }

    Vector<Parser> parsers(nOut);
    Vector<Vector<SymbolSource>> exprSources(nOut);
    std::map<std::string, int> outputSymbols;
    const Vector<std::string> coordSymbols{AMREX_D_DECL("x", "y", "z")};
    bool useCoords = false;
    for (int n = 0; n < nOut; ++n) {
      // An unquoted expression with blanks arrives as several tokens
      const std::string key = "expr." + outVars[n];
      const int nTok = pp.countval(key.c_str());
      if (nTok < 1) {
        Abort("No expression given for " + outVars[n] + " (" + key + ")");
      }
      std::string expr;
      for (int t = 0; t < nTok; ++t) {
        std::string tok;
        pp.get(key.c_str(), tok, t);
        expr += (t > 0 ? " " : "") + tok;
      }
      try {
        parsers[n].define(expr);
      } catch (const std::runtime_error& e) {
        Abort(std::string(e.what()) + " (expr." + outVars[n] + ")");
      }
      for (const auto& [sym, val] : constants) {
        parsers[n].setConstant(sym, val);
      }
      Vector<std::string> regVars;
      for (const auto& sym : parsers[n].symbols()) {
        // Earlier outputs shadow input variables of the same name
        // and coordinates are used only for names not otherwise taken
        const auto coordIt =
          std::find(coordSymbols.begin(), coordSymbols.end(), sym);
        if (outputSymbols.count(sym) > 0) {
          regVars.push_back(sym);
          exprSources[n].push_back({SymbolSource::Output, outputSymbols[sym]});
        } else if (inputSymbols.count(sym) > 0) {
          regVars.push_back(sym);
          const auto& [f, var] = inputSymbols[sym];
          exprSources[n].push_back({SymbolSource::Input, loadIndex(f, var)});
        } else if (coordIt != coordSymbols.end()) {
          regVars.push_back(sym);
          exprSources[n].push_back(
            {SymbolSource::Coord,
             static_cast<int>(coordIt - coordSymbols.begin())});
          useCoords = true;
        }
      }
      parsers[n].registerVariables(regVars);
      outputSymbols[outVars[n]] = n;
      Print() << "  " << outVars[n] << " = " << expr << "\n";
    }

    Vector<std::string> writeVars = outVars;
    if (pp.countval("writeVars") > 0) {
      writeVars.resize(pp.countval("writeVars"));
      pp.getarr("writeVars", writeVars, 0, writeVars.size());
      for (const auto& name : writeVars) {
        if (outputSymbols.count(name) == 0) {
          Abort("writeVars: " + name + " is not one of outVars");
        }
      }
    }

    // ---------------------------------------------------------------------
    // Evaluate level by level
    // ---------------------------------------------------------------------
    const int nLoad = loadVars.size();
    const int nWrite = nKeep + writeVars.size();
    Vector<MultiFab> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);
    RealBox rb(&(amrData0.ProbLo()[0]), &(amrData0.ProbHi()[0]));
    int coord = 0;
    Vector<int> is_per(AMREX_SPACEDIM, 0);

    for (int lev = 0; lev < Nlev; ++lev) {
      const BoxArray ba = amrData0.boxArray(lev);
      const DistributionMapping dm(ba);
      geoms[lev] =
        Geometry(amrData0.ProbDomain()[lev], &rb, coord, &(is_per[0]));

      MultiFab indata(ba, dm, std::max(nLoad, 1), 0);
      MultiFab results(ba, dm, nOut, 0);

      // Cell-center coordinates, if an expression uses x, y or z
      MultiFab coords;
      if (useCoords) {
        coords.define(ba, dm, AMREX_SPACEDIM, 0);
        const auto geomData = geoms[lev].data();
        for (MFIter mfi(coords, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& bx = mfi.tilebox();
          auto const& xyz = coords.array(mfi);
          ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            const IntVect iv(AMREX_D_DECL(i, j, k));
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
              xyz(i, j, k, d) =
                geomData.ProbLo(d) + (iv[d] + 0.5) * geomData.CellSize(d);
            }
          });
        }
      }

      Print() << "Reading data for level " << lev << std::endl;
      for (int f = 0; f < nFiles; ++f) {
        Vector<std::string> names;
        Vector<int> comps;
        for (int n = 0; n < nLoad; ++n) {
          if (loadVars[n].first == f) {
            names.push_back(loadVars[n].second);
            comps.push_back(n);
          }
        }
        if (!names.empty()) {
          dataServices[f]->AmrDataRef().FillVar(indata, lev, names, comps);
        }
      }

      for (int n = 0; n < nOut; ++n) {
        Vector<const MultiFab*> srcMF;
        Vector<int> srcComp;
        for (const auto& src : exprSources[n]) {
          srcMF.push_back(
            src.kind == SymbolSource::Output  ? &results
            : src.kind == SymbolSource::Coord ? &coords
                                              : &indata);
          srcComp.push_back(src.comp);
        }
        try {
          if (srcMF.empty()) {
            results.setVal(parsers[n].compileHost<0>()(), n, 1);
          } else {
            dispatchExpression<1>(
              static_cast<int>(srcMF.size()), parsers[n], results, n, srcMF,
              srcComp);
          }
        } catch (const std::runtime_error& e) {
          Abort(std::string(e.what()) + " (expr." + outVars[n] + ")");
        }
      }

      outdata[lev].define(ba, dm, nWrite, 0);
      for (int n = 0; n < nKeep; ++n) {
        MultiFab::Copy(outdata[lev], indata, keepIdx[n], n, 1, 0);
      }
      for (int n = 0; n < writeVars.size(); ++n) {
        MultiFab::Copy(
          outdata[lev], results, outputSymbols[writeVars[n]], nKeep + n, 1, 0);
      }
      Print() << "Expressions evaluated for level " << lev << std::endl;
    }

    // ---------------------------------------------------------------------
    // Write the results
    // ---------------------------------------------------------------------
    Vector<std::string> outNames;
    for (const auto& var : keepVars) {
      outNames.push_back(var);
    }
    for (const auto& var : writeVars) {
      outNames.push_back(var);
    }

    // Cap the number of plotfile data files via the n_files option (AMReX)
    int n_files = amrex::VisMF::GetNOutFiles();
    pp.query("n_files", n_files);
    amrex::VisMF::SetNOutFiles(n_files);

    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev - 1);
    for (int lev = 0; lev < Nlev - 1; ++lev) {
      int r = amrData0.RefRatio()[lev];
      refRatios[lev] = IntVect(AMREX_D_DECL(r, r, r));
    }
    amrex::WriteMultiLevelPlotfile(
      outfile, Nlev, GetVecOfConstPtrs(outdata), outNames, geoms,
      amrData0.Time(), isteps, refRatios);
  }
  Finalize();
  return 0;
}
