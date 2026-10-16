#ifndef RAYLEIGH_CROSS_SECTIONS_H
#define RAYLEIGH_CROSS_SECTIONS_H

// Rayleigh scattering cross-sections [m^2] per molecule, tabulated at 12
// wavelengths (157-694 nm). Generated from
// Rayleigh_cross-section_sheet.xls (33 species).
//
// Requires mechanism.H (NUM_SPECIES and <SPECIES>_ID macros). A
// NUM_SPECIES x n_lambda table is zero-initialised and filled at compile
// time for every tabulated species whose _ID is defined; all other
// mechanism species have sigma = 0.
//
// Interpolation is linear in log(sigma) vs log(lambda) (piecewise power
// law; exact for sigma ~ lambda^-4). Outside the table the nearest end
// point is scaled by lambda^-4. Species with zero data return 0.
// NB: some species have resonances within the table (O3, CN, CH, CH2,
// N2H, NO2 in the UV) where interpolation between points is crude.

//#include "mechanism.H"
#include <AMReX_REAL.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_Extension.H>
#include <cmath>
#include <PelePhysics.H>

namespace rayleigh {

constexpr int n_lambda = 12;

struct Table
{
  amrex::Real lam[n_lambda]; // [nm], ascending
  amrex::Real sig[NUM_SPECIES][n_lambda]; // [m^2], zero if no data
};

constexpr void
set_row(Table& t, int k, const amrex::Real (&s)[n_lambda])
{
  for (int i = 0; i < n_lambda; ++i) { t.sig[k][i] = s[i]; }
}

constexpr Table
make_table()
{
  Table t{}; // zero-initialised
  constexpr amrex::Real lam[n_lambda] = {157, 193, 248, 266, 338, 355, 444, 488, 514, 532, 633, 694};
  for (int i = 0; i < n_lambda; ++i) { t.lam[i] = lam[i]; }

#ifdef O_ID
  static_assert(O_ID < NUM_SPECIES, "O_ID out of range");
  set_row(t, O_ID, {3.069547e-30, 9.653826e-31, 2.990804e-31, 2.193067e-31, 7.854324e-32, 6.389379e-32, 2.524688e-32, 1.713031e-32, 1.385264e-32, 1.203794e-32, 5.936977e-33, 4.090246e-33}); // O
#endif
#ifdef O2_ID
  static_assert(O2_ID < NUM_SPECIES, "O2_ID out of range");
  set_row(t, O2_ID, {2.767962e-29, 4.945305e-30, 1.323971e-30, 9.563494e-31, 3.331307e-31, 2.700611e-31, 1.055948e-31, 7.144483e-32, 5.769877e-32, 5.010301e-32, 2.463440e-32, 1.695171e-32}); // O2
#endif
#ifdef O3_ID
  static_assert(O3_ID < NUM_SPECIES, "O3_ID out of range");
  set_row(t, O3_ID, {2.851928e-29, 9.071971e-30, 1.404860e-30, 7.744885e-31, 5.672591e-30, 2.796937e-30, 5.498772e-31, 3.361936e-31, 2.604359e-31, 2.211347e-31, 9.970259e-32, 6.638065e-32}); // O3
#endif
#ifdef OH_ID
  static_assert(OH_ID < NUM_SPECIES, "OH_ID out of range");
  set_row(t, OH_ID, {1.050156e-29, 2.738520e-30, 6.618108e-31, 4.661140e-31, 1.786107e-31, 1.429494e-31, 5.481660e-32, 3.697295e-32, 2.982040e-32, 2.587640e-32, 1.268888e-32, 8.723350e-33}); // OH
#endif
#ifdef C3H6_ID
  static_assert(C3H6_ID < NUM_SPECIES, "C3H6_ID out of range");
  set_row(t, C3H6_ID, {1.591699e-28, 1.847796e-28, 2.667967e-29, 1.825364e-29, 5.763581e-30, 4.614761e-30, 1.736203e-30, 1.162307e-30, 9.339996e-31, 8.087361e-31, 3.929563e-31, 2.691612e-31}); // C3H6propene
#endif
#ifdef AR_ID
  static_assert(AR_ID < NUM_SPECIES, "AR_ID out of range");
  set_row(t, AR_ID, {1.188313e-29, 4.194657e-30, 1.342180e-30, 9.881510e-31, 3.565265e-31, 2.902862e-31, 1.149996e-31, 7.807998e-32, 6.315898e-32, 5.489492e-32, 2.709128e-32, 1.866926e-32}); // AR
#endif
#ifdef C2H5_ID
  static_assert(C2H5_ID < NUM_SPECIES, "C2H5_ID out of range");
  set_row(t, C2H5_ID, {2.546817e-28, 2.964018e-29, 3.499817e-29, 7.899251e-30, 2.814809e-30, 2.238905e-30, 8.363310e-31, 5.599670e-31, 4.500996e-31, 3.898146e-31, 1.896384e-31, 1.299751e-31}); // C2H5
#endif
#ifdef CH_ID
  static_assert(CH_ID < NUM_SPECIES, "CH_ID out of range");
  set_row(t, CH_ID, {2.222675e-29, 9.870559e-30, 3.102621e-30, 1.942868e-30, 3.373356e-31, 3.789060e-30, 3.759088e-31, 1.953233e-31, 1.480010e-31, 1.248333e-31, 5.614723e-32, 3.764020e-32}); // CH
#endif
#ifdef CH2_ID
  static_assert(CH2_ID < NUM_SPECIES, "CH2_ID out of range");
  set_row(t, CH2_ID, {1.971209e-29, 6.070333e-29, 4.804409e-30, 2.831781e-30, 8.046953e-31, 6.400045e-31, 2.370010e-31, 1.581298e-31, 1.268840e-31, 1.097785e-31, 5.317545e-32, 3.638266e-32}); // CH2
#endif
#ifdef CH3_ID
  static_assert(CH3_ID < NUM_SPECIES, "CH3_ID out of range");
  set_row(t, CH3_ID, {3.973411e-29, 8.897818e-30, 1.223115e-29, 3.836045e-30, 9.338568e-31, 7.404601e-31, 2.735556e-31, 1.826365e-31, 1.466130e-31, 1.268844e-31, 6.154955e-32, 4.213962e-32}); // CH3
#endif
#ifdef CH4_ID
  static_assert(CH4_ID < NUM_SPECIES, "CH4_ID out of range");
  set_row(t, CH4_ID, {4.165867e-29, 1.195149e-29, 3.486919e-30, 2.527366e-30, 8.805707e-31, 7.133449e-31, 2.779241e-31, 1.877884e-31, 1.515497e-31, 1.315443e-31, 6.455595e-32, 4.438769e-32}); // CH4
#endif
#ifdef CN_ID
  static_assert(CN_ID < NUM_SPECIES, "CN_ID out of range");
  set_row(t, CN_ID, {2.967412e-29, 7.978716e-30, 2.164519e-30, 1.456428e-30, 2.896999e-31, 2.374075e-31, 1.153177e-30, 4.200091e-31, 2.857358e-31, 2.288768e-31, 8.589124e-32, 5.329075e-32}); // CN
#endif
#ifdef CO_ID
  static_assert(CO_ID < NUM_SPECIES, "CO_ID out of range");
  set_row(t, CO_ID, {5.456109e-29, 8.028222e-30, 2.168388e-30, 1.559682e-30, 5.363980e-31, 4.338992e-31, 1.683611e-31, 1.136444e-31, 9.167389e-32, 7.955186e-32, 3.900250e-32, 2.680766e-32}); // CO
#endif
#ifdef CO2_ID
  static_assert(CO2_ID < NUM_SPECIES, "CO2_ID out of range");
  set_row(t, CO2_ID, {3.458430e-29, 1.162644e-29, 3.611666e-30, 2.644543e-30, 9.422272e-31, 7.657356e-31, 3.014841e-31, 2.043234e-31, 1.651340e-31, 1.434540e-31, 7.064763e-32, 4.864372e-32}); // CO2
#endif
#ifdef H_ID
  static_assert(H_ID < NUM_SPECIES, "H_ID out of range");
  set_row(t, H_ID, {7.791985e-30, 1.634410e-30, 3.941870e-31, 2.760150e-31, 8.884940e-32, 7.114156e-32, 2.666653e-32, 1.781541e-32, 1.429998e-32, 1.237368e-32, 5.993069e-33, 4.099324e-33}); // H
#endif
#ifdef H2_ID
  static_assert(H2_ID < NUM_SPECIES, "H2_ID out of range");
  set_row(t, H2_ID, {4.029758e-30, 1.264072e-30, 3.741866e-31, 2.714100e-31, 9.454114e-32, 7.656870e-32, 2.979759e-32, 2.012467e-32, 1.623761e-32, 1.409219e-32, 6.911631e-33, 4.751134e-33}); // H2
#endif
#ifdef H2O_ID
  static_assert(H2O_ID < NUM_SPECIES, "H2O_ID out of range");
  set_row(t, H2O_ID, {9.589597e-30, 7.412702e-30, 1.245024e-30, 8.809252e-31, 2.966860e-31, 2.395648e-31, 9.255612e-32, 6.242140e-32, 5.033530e-32, 4.367125e-32, 2.139602e-32, 1.470293e-32}); // H2O
#endif
#ifdef H2O2_ID
  static_assert(H2O2_ID < NUM_SPECIES, "H2O2_ID out of range");
  set_row(t, H2O2_ID, {2.507603e-29, 9.281921e-30, 3.026340e-30, 2.143285e-30, 7.333443e-31, 5.935242e-31, 2.309365e-31, 1.560432e-31, 1.259402e-31, 1.093214e-31, 5.366968e-32, 3.690980e-32}); // H2O2
#endif
#ifdef HCN_ID
  static_assert(HCN_ID < NUM_SPECIES, "HCN_ID out of range");
  set_row(t, HCN_ID, {4.703627e-29, 1.296601e-29, 3.718811e-30, 2.688065e-30, 9.310444e-31, 7.536141e-31, 2.928434e-31, 1.977193e-31, 1.595090e-31, 1.384247e-31, 6.787573e-32, 4.665511e-32}); // HCN
#endif
#ifdef HO2_ID
  static_assert(HO2_ID < NUM_SPECIES, "HO2_ID out of range");
  set_row(t, HO2_ID, {7.497297e-30, 1.471313e-30, 3.067853e-30, 2.008034e-30, 6.163532e-31, 4.935745e-31, 1.865678e-31, 1.251846e-31, 1.007176e-31, 8.727325e-32, 4.254357e-32, 2.917998e-32}); // HO2
#endif
#ifdef N2_ID
  static_assert(N2_ID < NUM_SPECIES, "N2_ID out of range");
  set_row(t, N2_ID, {1.444835e-29, 5.070496e-30, 1.608047e-30, 1.181521e-30, 4.241851e-31, 3.451089e-31, 1.363576e-31, 9.250639e-32, 7.480023e-32, 6.499819e-32, 3.204720e-32, 2.207582e-32}); // N2
#endif
#ifdef C2H6_ID
  static_assert(C2H6_ID < NUM_SPECIES, "C2H6_ID out of range");
  set_row(t, C2H6_ID, {1.329190e-28, 3.565373e-29, 1.027923e-29, 7.440505e-30, 2.586022e-30, 2.094295e-30, 8.152482e-31, 5.507202e-31, 4.444027e-31, 3.857147e-31, 1.892473e-31, 1.301122e-31}); // C2H6
#endif
#ifdef CH3NO2_ID
  static_assert(CH3NO2_ID < NUM_SPECIES, "CH3NO2_ID out of range");
  set_row(t, CH3NO2_ID, {8.766899e-29, 9.928763e-29, 1.504067e-29, 1.050760e-29, 3.465635e-30, 2.789995e-30, 1.067937e-30, 7.183974e-31, 5.786061e-31, 5.016597e-31, 2.450831e-31, 1.682317e-31}); // CH3NO2
#endif
#ifdef CH3OH_ID
  static_assert(CH3OH_ID < NUM_SPECIES, "CH3OH_ID out of range");
  set_row(t, CH3OH_ID, {4.073711e-29, 2.188640e-29, 5.745657e-30, 4.136386e-30, 1.428917e-30, 1.156788e-30, 4.500783e-31, 3.040484e-31, 2.453621e-31, 2.129657e-31, 1.045093e-31, 7.185987e-32}); // CH3OH
#endif
#ifdef N_ID
  static_assert(N_ID < NUM_SPECIES, "N_ID out of range");
  set_row(t, N_ID, {1.032491e-29, 2.567745e-30, 7.162871e-31, 5.160747e-31, 1.777302e-31, 1.437653e-31, 5.576158e-32, 3.763166e-32, 3.035312e-32, 2.633794e-32, 1.290893e-32, 8.871649e-33}); // N
#endif
#ifdef NO_ID
  static_assert(NO_ID < NUM_SPECIES, "NO_ID out of range");
  set_row(t, NO_ID, {9.749492e-30, 6.786651e-30, 1.783305e-30, 1.256587e-30, 4.235694e-31, 3.422519e-31, 1.325947e-31, 8.950585e-32, 7.220816e-32, 6.266503e-32, 3.073652e-32, 2.113155e-32}); // NO
#endif
#ifdef NO2_ID
  static_assert(NO2_ID < NUM_SPECIES, "NO2_ID out of range");
  set_row(t, NO2_ID, {4.521953e-27, 1.835112e-29, 4.344292e-30, 2.921345e-30, 8.410600e-31, 5.326045e-31, 4.655148e-31, 2.854674e-31, 2.241007e-31, 1.920103e-31, 9.052259e-32, 6.152590e-32}); // NO2
#endif
#ifdef C2H2_ID
  static_assert(C2H2_ID < NUM_SPECIES, "C2H2_ID out of range");
  set_row(t, C2H2_ID, {1.105669e-28, 3.745709e-29, 8.194124e-30, 5.751965e-30, 1.882003e-30, 1.512047e-30, 5.740412e-31, 3.850884e-31, 3.097299e-31, 2.683247e-31, 1.306301e-31, 8.953917e-32}); // C2H2
#endif
#ifdef C2H4_ID
  static_assert(C2H4_ID < NUM_SPECIES, "C2H4_ID out of range");
  set_row(t, C2H4_ID, {8.416556e-29, 1.004164e-28, 1.264423e-29, 8.657099e-30, 2.724577e-30, 2.179874e-30, 8.179335e-31, 5.471176e-31, 4.394750e-31, 3.804476e-31, 1.846752e-31, 1.264470e-31}); // C2H4
#endif
#ifdef HNO_ID
  static_assert(HNO_ID < NUM_SPECIES, "HNO_ID out of range");
  set_row(t, HNO_ID, {6.679462e-30, 8.457736e-30, 4.712378e-30, 2.760990e-30, 7.863549e-31, 6.271060e-31, 2.339060e-31, 1.562898e-31, 1.254883e-31, 1.085526e-31, 5.236647e-32, 3.542235e-32}); // HNO
#endif
#ifdef N2H_ID
  static_assert(N2H_ID < NUM_SPECIES, "N2H_ID out of range");
  set_row(t, N2H_ID, {1.121937e-29, 7.348559e-28, 2.485517e-30, 1.816378e-30, 1.234106e-30, 9.435827e-31, 3.257502e-31, 2.148513e-31, 1.715980e-31, 1.479782e-31, 7.538398e-32, 4.914061e-32}); // N2H
#endif
#ifdef N2O_ID
  static_assert(N2O_ID < NUM_SPECIES, "N2O_ID out of range");
  set_row(t, N2O_ID, {7.039434e-29, 1.859748e-29, 5.261097e-30, 3.801781e-30, 1.313596e-30, 1.063448e-30, 4.131449e-31, 2.789061e-31, 2.250281e-31, 1.952549e-31, 9.574681e-32, 6.581146e-32}); // N2O
#endif
#ifdef NH_ID
  static_assert(NH_ID < NUM_SPECIES, "NH_ID out of range");
  set_row(t, NH_ID, {1.244195e-29, 4.568588e-30, 1.087793e-30, 7.372716e-31, 5.727363e-31, 3.509498e-31, 1.081489e-31, 7.122079e-32, 5.692737e-32, 4.914290e-32, 2.367878e-32, 1.618089e-32}); // NH
#endif
  return t;
}

// Cross-section [m^2] for species k (mechanism index) at wavelength lambda_nm [nm].
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
cross_section_wavelength(int k, amrex::Real lambda_nm) noexcept
{
  constexpr Table tab = make_table();
  const amrex::Real* lam = tab.lam;
  const amrex::Real* s = tab.sig[k];

  if (s[0] == amrex::Real(0.0)) { return amrex::Real(0.0); }
  if (lambda_nm <= lam[0]) {
    const amrex::Real r = lam[0] / lambda_nm;
    return s[0] * r * r * r * r;
  }
  if (lambda_nm >= lam[n_lambda - 1]) {
    const amrex::Real r = lam[n_lambda - 1] / lambda_nm;
    return s[n_lambda - 1] * r * r * r * r;
  }
  int i = 0;
  while (lambda_nm > lam[i + 1]) { ++i; }
  const amrex::Real t = std::log(lambda_nm / lam[i]) / std::log(lam[i + 1] / lam[i]);
  return s[i] * std::pow(s[i + 1] / s[i], t);
}

// Same, taking wavenumber nu [cm^-1] (lambda_nm = 1e7 / nu).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
cross_section_wavenumber(int k, amrex::Real nu_cm1) noexcept
{
  return cross_section_wavelength(k, amrex::Real(1.0e7) / nu_cm1);
}

// Fill sigma[0..NUM_SPECIES) at one wavelength [nm] (e.g. once, before a loop).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
cross_sections(amrex::Real lambda_nm, amrex::Real sigma[NUM_SPECIES]) noexcept
{
  for (int k = 0; k < NUM_SPECIES; ++k) {
    sigma[k] = cross_section_wavelength(k, lambda_nm);
  }
}

} // namespace rayleigh

#endif
