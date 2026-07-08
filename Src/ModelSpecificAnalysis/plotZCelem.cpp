#include <string>
#include <iostream>
#include <set>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DataServices.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <PelePhysics.H>

using namespace amrex;

static void
print_usage(int, char* argv[])
{
  std::cerr << "usage:\n";
  std::cerr << argv[0] << " infile infile=f1 [options] \n\tOptions:\n";
  std::cerr << "\tfuelNames=H2 CH4       (list of fuel species)\n";
  std::cerr << "\tfuelMoleFracs=0.4 0.6  (mole fractions of each fuel species, "
               "must sum to 1)\n";
  exit(1);
}

std::string
getFileRoot(const std::string& infile)
{
  std::vector<std::string> tokens = Tokenize(infile, std::string("/"));
  return tokens[tokens.size() - 1];
}

int
main(int argc, char* argv[])
{
  Initialize(argc, argv);
  {
    if (argc < 2)
      print_usage(argc, argv);

    ParmParse pp;

    if (pp.contains("help"))
      print_usage(argc, argv);

    if (pp.contains("verbose"))
      AmrData::SetVerbose(false);

    std::string plotFileName;
    pp.get("infile", plotFileName);

    // ---- Multi-species fuel input ----
    // Default: pure H2
    Vector<std::string> fuelNames;
    Vector<Real> fuelMoleFracs;
    if (pp.countval("fuelNames") > 0) {
      pp.getarr("fuelNames", fuelNames);
      pp.getarr("fuelMoleFracs", fuelMoleFracs);
      if (fuelNames.size() != fuelMoleFracs.size())
        amrex::Abort("fuelNames and fuelMoleFracs must have the same length");
      // Normalise just in case
      Real sumX = 0.0;
      for (Real x : fuelMoleFracs)
        sumX += x;
      for (Real& x : fuelMoleFracs)
        x /= sumX;
    } else {
      // Legacy single-fuel path
      std::string singleFuel = "H2";
      pp.query("fuelName", singleFuel);
      fuelNames = {singleFuel};
      fuelMoleFracs = {1.0};
    }

    std::string productName = "H2O";
    pp.query("productName", productName);
    int clipProgress = 0;
    pp.query("clipProgress", clipProgress);
    Vector<int> is_per(AMREX_SPACEDIM, 1);
    pp.queryarr("is_per", is_per, 0, AMREX_SPACEDIM);
    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);

    DataServices dataServices(plotFileName, fileType);
    if (!dataServices.AmrDataOk())
      DataServices::Dispatch(DataServices::ExitRequest, NULL);
    AmrData& amrData = dataServices.AmrDataRef();

    int finestLevel = amrData.FinestLevel();
    pp.query("finestLevel", finestLevel);
    int Nlev = finestLevel + 1;

    Vector<std::string> spec_names;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
      spec_names);
    auto eos = pele::physics::PhysicsType::eos();

    // ---- Output components ----
    // Z, CF, CP  +  Z_H, Z_C, Z_O, Z_N
    constexpr int nCompIn = NUM_SPECIES;
    constexpr int nCompOut = 7; // Z  CF  CP  Z_H  Z_C  Z_O  Z_N
    constexpr int idZlocal = 0;
    constexpr int idCFlocal = 1;
    constexpr int idCPlocal = 2;
    constexpr int idZH = 3;
    constexpr int idZC = 4;
    constexpr int idZO = 5;
    constexpr int idZN = 6;

    Vector<std::string> outNames(nCompOut);
    outNames[idZlocal] = "Z";
    outNames[idCFlocal] = "CF";
    outNames[idCPlocal] = "CP";
    outNames[idZH] = "Z_H";
    outNames[idZC] = "Z_C";
    outNames[idZO] = "Z_O";
    outNames[idZN] = "Z_N";

    Vector<std::string> inNames(nCompIn);
    Vector<int> destFillComps(nCompIn);
    for (int i = 0; i < NUM_SPECIES; ++i) {
      destFillComps[i] = i;
      inNames[i] = "Y(" + spec_names[i] + ")";
    }

    // ---- EOS data ----
    amrex::Real mwt[NUM_SPECIES];
    eos.molecular_weight(mwt);

    int ecompCHON[NUM_SPECIES * 4];
    pele::physics::eos::element_compositionCHON<
      pele::physics::PhysicsType::eos_type>(ecompCHON);

    Real atwCHON[4] = {0.0};
    pele::physics::eos::atomic_weightsCHON<
      pele::physics::PhysicsType::eos_type>(atwCHON);

    // Bilger weights  (C, H, O, N  in CHON order: index 0=C,1=H,2=O,3=N)
    Array<amrex::Real, 4> Beta_mix;
    Beta_mix[0] = (atwCHON[0] != 0.0) ? 2.0 / atwCHON[0] : 0.0;         // C
    Beta_mix[1] = (atwCHON[1] != 0.0) ? 1.0 / (2.0 * atwCHON[1]) : 0.0; // H
    Beta_mix[2] = (atwCHON[2] != 0.0) ? -1.0 / atwCHON[2] : 0.0;        // O
    Beta_mix[3] = 0.0;                                                  // N

    amrex::Array<amrex::Real, NUM_SPECIES> spec_Bilger_fact;
    for (int i = 0; i < NUM_SPECIES; ++i) {
      spec_Bilger_fact[i] = 0.0;
      for (int k = 0; k < 4; k++)
        spec_Bilger_fact[i] +=
          Beta_mix[k] * (ecompCHON[i * 4 + k] * atwCHON[k] / mwt[i]);
    }

    // ---- Per-element mass-fraction weights for each species ----
    // w_elem[e][i] = (atoms_e_in_i * atw_e) / mwt_i   => mass fraction of
    // element e in species i CHON order: 0=C, 1=H, 2=O, 3=N
    amrex::Array2D<amrex::Real, 0, 3, 0, NUM_SPECIES - 1> w_elem;
    for (int e = 0; e < 4; ++e)
      for (int i = 0; i < NUM_SPECIES; ++i)
        w_elem(e, i) =
          (atwCHON[e] > 0.0) ? ecompCHON[i * 4 + e] * atwCHON[e] / mwt[i] : 0.0;

    // ---- Build fuel / oxidiser Y arrays from mole fractions ----
    // Convert mole fractions X_k to mass fractions Y_k:
    //   Y_k = X_k * W_k / sum_j(X_j * W_j)
    Real YF[NUM_SPECIES], YO[NUM_SPECIES];
    for (int i = 0; i < NUM_SPECIES; ++i) {
      YF[i] = 0.0;
      YO[i] = 0.0;
    }

    // Oxidiser: standard air
    for (int i = 0; i < NUM_SPECIES; ++i) {
      if (spec_names[i] == "O2")
        YO[i] = 0.233;
      if (spec_names[i] == "N2")
        YO[i] = 0.767;
    }

    // Fuel: multi-species mixture
    {
      // Build a map from species name to index for quick lookup
      std::map<std::string, int> specIdx;
      for (int i = 0; i < NUM_SPECIES; ++i)
        specIdx[spec_names[i]] = i;

      // Numerator  W_k * X_k for each fuel species
      Real denom = 0.0;
      Vector<int> fuelIdx(fuelNames.size(), -1);
      for (int f = 0; f < (int)fuelNames.size(); ++f) {
        auto it = specIdx.find(fuelNames[f]);
        if (it == specIdx.end())
          amrex::Abort("Fuel species not found in mechanism: " + fuelNames[f]);
        fuelIdx[f] = it->second;
        denom += fuelMoleFracs[f] * mwt[fuelIdx[f]];
      }
      for (int f = 0; f < (int)fuelNames.size(); ++f)
        YF[fuelIdx[f]] = fuelMoleFracs[f] * mwt[fuelIdx[f]] / denom;

      Print() << "Fuel composition (mass fractions):\n";
      for (int f = 0; f < (int)fuelNames.size(); ++f)
        Print() << "  Y(" << fuelNames[f] << ") = " << YF[fuelIdx[f]] << "\n";
    }

    // ---- Bilger mixture-fraction at pure fuel / pure oxidiser ----
    Real Zfu = 0.0, Zox = 0.0;
    for (int i = 0; i < NUM_SPECIES; ++i) {
      Zfu += spec_Bilger_fact[i] * YF[i];
      Zox += spec_Bilger_fact[i] * YO[i];
    }
    const Real denom_inv = 1.0 / (Zfu - Zox);
    Print() << "Zfu = " << Zfu << "  Zox = " << Zox << "\n";

    // ---- Find fuel-species index for CF (use first / dominant fuel species)
    // ---- For multi-fuel we track the species with the highest mass fraction
    // in the fuel stream.
    int YFcomp = -1;
    {
      Real maxYF = -1.0;
      for (int i = 0; i < NUM_SPECIES; ++i)
        if (YF[i] > maxYF) {
          maxYF = YF[i];
          YFcomp = i;
        }
    }
    const int YPcomp = amrData.StateNumber("Y(" + productName + ")") -
                       amrData.StateNumber("Y(" + spec_names[0] + ")");

    // ---- Global min/max for normalisation ----
    Vector<Box> probDomain = amrData.ProbDomain();
    Vector<Real> YFminmax = {1e100, -1e100};
    Vector<Real> YPminmax = {1e100, -1e100};
    for (int lev = 0; lev < Nlev; ++lev) {
      Real mn, mx;
      amrData.MinMax(
        probDomain[lev], "Y(" + spec_names[YFcomp] + ")", lev, mn, mx);
      YFminmax[0] = std::min(YFminmax[0], mn);
      YFminmax[1] = std::max(YFminmax[1], mx);
      amrData.MinMax(probDomain[lev], "Y(" + productName + ")", lev, mn, mx);
      YPminmax[0] = std::min(YPminmax[0], mn);
      YPminmax[1] = std::max(YPminmax[1], mx);
    }
    pp.queryarr("YFminmax", YFminmax);
    pp.queryarr("YPminmax", YPminmax);
    Print() << "YF (dominant fuel species " << spec_names[YFcomp]
            << ") range: [" << YFminmax[0] << ", " << YFminmax[1] << "]\n";

    // ---- Level loop ----
    RealBox rb(&(amrData.ProbLo()[0]), &(amrData.ProbHi()[0]));
    constexpr int nGrow = 0;
    Vector<MultiFab> outdata(Nlev);
    Vector<Geometry> geoms(Nlev);

    for (int lev = 0; lev < Nlev; ++lev) {
      const BoxArray ba = amrData.boxArray(lev);
      const DistributionMapping dm(ba);
      MultiFab indata(ba, dm, nCompIn, nGrow);
      outdata[lev].define(ba, dm, nCompOut, nGrow);

      Print() << "Reading data for level " << lev << std::endl;
      amrData.FillVar(indata, lev, inNames, destFillComps);
      geoms[lev] = Geometry(amrData.ProbDomain()[lev], &rb, 0, &(is_per[0]));
      Print() << "Data has been read for level " << lev << std::endl;

      auto in_ma = indata.const_arrays();
      auto out_ma = outdata[lev].arrays();

      amrex::ParallelFor(
        outdata[lev],
        [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
          // Progress variables
          out_ma[box_no](i, j, k, idCFlocal) =
            1.0 - in_ma[box_no](i, j, k, YFcomp) / YFminmax[1];
          out_ma[box_no](i, j, k, idCPlocal) =
            in_ma[box_no](i, j, k, YPcomp) / YPminmax[1];

          // Bilger mixture fraction Z
          Real Zloc = 0.0;
          for (int n = 0; n < NUM_SPECIES; ++n)
            Zloc += in_ma[box_no](i, j, k, n) * spec_Bilger_fact[n];
          out_ma[box_no](i, j, k, idZlocal) = (Zloc - Zox) * denom_inv;

          // Elemental mass fractions  Z_H, Z_C, Z_O, Z_N
          // CHON order: 0=C, 1=H, 2=O, 3=N
          Real ZC_loc = 0.0, ZH_loc = 0.0, ZO_loc = 0.0, ZN_loc = 0.0;
          for (int n = 0; n < NUM_SPECIES; ++n) {
            const Real Yn = in_ma[box_no](i, j, k, n);
            ZC_loc += w_elem(0, n) * Yn;
            ZH_loc += w_elem(1, n) * Yn;
            ZO_loc += w_elem(2, n) * Yn;
            ZN_loc += w_elem(3, n) * Yn;
          }
          out_ma[box_no](i, j, k, idZH) = ZH_loc;
          out_ma[box_no](i, j, k, idZC) = ZC_loc;
          out_ma[box_no](i, j, k, idZO) = ZO_loc;
          out_ma[box_no](i, j, k, idZN) = ZN_loc;
        });

      Print() << "Derive finished for level " << lev << std::endl;
    }

    std::string outfile(getFileRoot(plotFileName) + "_ZC");
    Print() << "Writing new data to " << outfile << std::endl;
    Vector<int> isteps(Nlev, 0);
    Vector<IntVect> refRatios(Nlev - 1, {AMREX_D_DECL(2, 2, 2)});
    amrex::WriteMultiLevelPlotfile(
      outfile, Nlev, GetVecOfConstPtrs(outdata), outNames, geoms, 0.0, isteps,
      refRatios);
  }
  Finalize();
  return 0;
}
