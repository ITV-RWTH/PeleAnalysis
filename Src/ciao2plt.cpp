#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <iostream>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_Print.H>

#include "H5Cpp.h"

using namespace amrex;
using namespace H5;

static void
print_usage(int, char* argv[])
{
  std::cerr
    << "usage:\n"
    << "  " << argv[0] << " infile=<file.h5> [OPTIONS]\n\n"
    << "Options:\n"
    << "  outfile=<name>          output plt name [default: plt_<stem>]\n"
    << "  layout=native|flat      HDF5 layout [default: native]\n"
    << "  vars=\"v1 v2 ...\"        fields to read [default: all]; flat layout:\n"
    << "                          dataset paths, e.g. \"SC/T flow/rho\"\n"
    << "  names=\"n1 n2 ...\"       flat layout: output names of vars\n"
    << "                          [default: path with '/' -> '_']\n"
    << "  coords=\"cx cy cz\"       flat layout: cell-centre coordinate\n"
    << "                          datasets [default: coords/xm coords/ym "
       "coords/zm]\n"
    << "  time_dset=<path>        flat layout: time dataset [default: "
       "flow/time1]\n"
    << "  interp_dx=<h>           flat layout: uniform spacing for a "
       "non-uniform\n"
    << "                          x or y direction [default: dx or dy of "
       "the\n"
    << "                          first uniform direction]\n"
    << "  interp_range=\"lo hi\"    flat layout: extent of the interpolated\n"
    << "                          direction [default: full]\n"
    << "  interp_align=0|1        flat layout: align the interpolated cells\n"
    << "                          with the finest source cells [default: 1]\n"
    << "  uniform_tol=<eps>       flat layout: relative tolerance of the\n"
    << "                          uniformity check [default: 1e-3]\n"
    << "  max_grid_size=<N>       AMReX box size [default: 32]\n"
    << "  per=\"0 0 0\"             periodicity (fallback if not in file)\n";
  std::exit(1);
}

// Strip directory and extension: /path/to/data.h5 -> data
static std::string
stem(const std::string& path)
{
  auto sep = path.rfind('/');
  std::string base = (sep == std::string::npos) ? path : path.substr(sep + 1);
  auto dot = base.rfind('.');
  return (dot == std::string::npos) ? base : base.substr(0, dot);
}

// H5Lexists emits diagnostics in HDF5 >=1.14 when an intermediate path
// component is absent. Suppress the error stack for the duration of the check.
static bool
h5_exists(hid_t loc, const char* name)
{
  H5E_auto_t old_func;
  void* old_data;
  H5Eget_auto(H5E_DEFAULT, &old_func, &old_data);
  H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
  htri_t ret = H5Lexists(loc, name, H5P_DEFAULT);
  H5Eset_auto(H5E_DEFAULT, old_func, old_data);
  return ret > 0;
}

// Find the single mesh group (skip IO-information).
static std::string
discover_mesh_group(const H5File& file)
{
  Group root = file.openGroup("/");
  std::string mesh;
  hsize_t n = root.getNumObjs();
  for (hsize_t i = 0; i < n; ++i) {
    std::string name = root.getObjnameByIdx(i);
    if (root.getObjTypeByIdx(i) == H5G_GROUP && name != "IO-information") {
      if (!mesh.empty())
        Abort("Multiple mesh groups found in HDF5 file; only one is supported");
      mesh = name;
    }
  }
  if (mesh.empty())
    Abort("No mesh group found in HDF5 file");
  Print() << "Mesh group: " << mesh << "\n";
  return mesh;
}

// List all field names: cv_data_real datasets + scalars/SC attribute names.
static Vector<std::string>
list_fields(const H5File& file, const std::string& mesh)
{
  Vector<std::string> names;

  Group cvg = file.openGroup("/" + mesh + "/data/cv_data_real");
  hsize_t n = cvg.getNumObjs();
  for (hsize_t i = 0; i < n; ++i)
    if (cvg.getObjTypeByIdx(i) == H5G_DATASET)
      names.push_back(cvg.getObjnameByIdx(i));

  std::string scpath = "/" + mesh + "/data/scalars";
  if (h5_exists(file.getId(), scpath.c_str())) {
    Group scalg = file.openGroup(scpath);
    if (h5_exists(scalg.getId(), "SC")) {
      DataSet sc = scalg.openDataSet("SC");
      hsize_t dims[4] = {};
      sc.getSpace().getSimpleExtentDims(dims);
      for (int i = 1; i <= (int)dims[0]; ++i) {
        Attribute attr = sc.openAttribute("Index " + std::to_string(i));
        std::string val;
        attr.read(attr.getStrType(), val);
        names.push_back(val);
      }
    }
  }
  return names;
}

// Read simulation time from globals_r0/time (required).
static Real
read_time(const H5File& file, const std::string& mesh)
{
  std::string path = "/" + mesh + "/data/globals_r0/time";
  if (!h5_exists(file.getId(), path.c_str()))
    Abort("Time dataset not found at " + path);
  double t = 0.0;
  file.openDataSet(path).read(&t, PredType::NATIVE_DOUBLE);
  Print() << "Time = " << t << "\n";
  return (Real)t;
}

// Read periodicity from sd_info if present; fall back to per= input parameter.
static Vector<int>
read_periodicity(const H5File& file, const std::string& mesh)
{
  Vector<int> per(3, 0);
  std::string sdpath = "/" + mesh + "/geometry/sd_info";
  if (h5_exists(file.getId(), sdpath.c_str())) {
    Group sd = file.openGroup(sdpath);
    if (h5_exists(sd.getId(), "xper")) {
      auto ri = [&](const std::string& name) {
        int v = 0;
        sd.openDataSet(name).read(&v, PredType::NATIVE_INT);
        return v;
      };
      per[0] = ri("xper");
      per[1] = ri("yper");
      per[2] = ri("zper");
      Print() << "Periodicity (from file): " << per[0] << " " << per[1] << " "
              << per[2] << "\n";
      return per;
    }
  }
  ParmParse pp;
  pp.queryarr("per", per, 0, 3);
  return per;
}

// Read coordinate system type from sd_info/icyl if present; fall back to
// coord_sys= input.
static int
read_coord_sys(const H5File& file, const std::string& mesh)
{
  std::string sdpath = "/" + mesh + "/geometry/sd_info";
  if (h5_exists(file.getId(), sdpath.c_str())) {
    Group sd = file.openGroup(sdpath);
    if (h5_exists(sd.getId(), "icyl")) {
      int v = 0;
      sd.openDataSet("icyl").read(&v, PredType::NATIVE_INT);
      Print() << "Coord sys (from file): " << v << "\n";
      return v;
    }
  }
  int coord = 0;
  ParmParse pp;
  pp.query("coord_sys", coord);
  return coord;
}

// Build BoxArray and Geometry from grid node coordinates.
// Grid coords are node-based: N+1 nodes define N cells.
static void
setup_grid(
  const H5File& file,
  const std::string& mesh,
  const Vector<int>& per,
  int coord_sys,
  int max_grid_size,
  BoxArray& ba,
  Geometry& geom)
{
  std::string gp = "/" + mesh + "/geometry/grid/";

  auto read_coord = [&](const std::string& name) {
    DataSet ds = file.openDataSet(gp + name);
    hsize_t dim[1] = {};
    ds.getSpace().getSimpleExtentDims(dim);
    Vector<Real> v(dim[0]);
    ds.read(v.data(), PredType::NATIVE_DOUBLE);
    return v;
  };

  Vector<Real> x = read_coord("x");
  Vector<Real> y = read_coord("y");
  Vector<Real> z = read_coord("z");

  int nx = (int)x.size() - 1;
  int ny = (int)y.size() - 1;
  int nz = (int)z.size() - 1;

  Print() << "Grid: " << nx << " x " << ny << " x " << nz << " cells\n";

#if AMREX_SPACEDIM == 2
  if (nz != 1)
    Abort("DIM=2 but file has nz=" + std::to_string(nz) + " (expected 1)");
#endif

  Box domain(IntVect::Zero, IntVect(AMREX_D_DECL(nx - 1, ny - 1, nz - 1)));
  Real plo[AMREX_SPACEDIM] = {AMREX_D_DECL(x[0], y[0], z[0])};
  Real phi[AMREX_SPACEDIM] = {AMREX_D_DECL(x[nx], y[ny], z[nz])};
  RealBox rb(plo, phi);
#if AMREX_SPACEDIM != 2
  if (coord_sys != 0)
    Abort("coord_sys != 0 (cylindrical/RZ) requires a DIM=2 build");
#endif
  geom = Geometry(domain, &rb, coord_sys, per.data());

  ba = BoxArray(domain);
  ba.maxSize(max_grid_size);
}

// Read all requested fields into data.
// Fields come from cv_data_real (3D datasets) or scalars/SC (4D stacked array).
static void
read_fields(
  const H5File& file,
  const std::string& mesh,
  const Vector<std::string>& vars,
  MultiFab& data)
{
  std::string cvbase = "/" + mesh + "/data/cv_data_real/";
  std::string scpath = "/" + mesh + "/data/scalars/SC";
  bool has_sc = h5_exists(file.getId(), scpath.c_str());

  // Build scalar name -> 0-based index map.
  std::map<std::string, int> sc_map;
  if (has_sc) {
    DataSet sc = file.openDataSet(scpath);
    hsize_t dims[4] = {};
    sc.getSpace().getSimpleExtentDims(dims);
    for (int i = 1; i <= (int)dims[0]; ++i) {
      Attribute attr = sc.openAttribute("Index " + std::to_string(i));
      std::string name;
      attr.read(attr.getStrType(), name);
      sc_map[name] = i - 1;
    }
  }

  // Validate all vars exist before starting I/O.
  for (const auto& v : vars) {
    bool in_cv = h5_exists(file.getId(), (cvbase + v).c_str());
    bool in_sc = sc_map.count(v) > 0;
    if (!in_cv && !in_sc) {
      auto avail = list_fields(file, mesh);
      std::string msg = "Field '" + v + "' not found. Available:";
      for (const auto& a : avail)
        msg += " " + a;
      Abort(msg);
    }
  }

  for (int n = 0; n < (int)vars.size(); ++n) {
    const std::string& vname = vars[n];
    Print() << "Reading " << vname << "\n";

    bool is_cv = h5_exists(file.getId(), (cvbase + vname).c_str());
    DataSet ds = file.openDataSet(is_cv ? cvbase + vname : scpath);

    for (MFIter mfi(data); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const int* lo = box.loVect();
      const int* hi = box.hiVect();

      int nx_b = hi[0] - lo[0] + 1;
      int ny_b = hi[1] - lo[1] + 1;
#if AMREX_SPACEDIM == 3
      int nz_b = hi[2] - lo[2] + 1;
      int lo2 = lo[2];
#else
      int nz_b = 1;
      int lo2 = 0;
#endif
      Vector<Real> buf(nx_b * ny_b * nz_b);
      DataSpace fsp = ds.getSpace();

      if (is_cv) {
        hsize_t cnt[3] = {(hsize_t)nz_b, (hsize_t)ny_b, (hsize_t)nx_b};
        hsize_t off[3] = {(hsize_t)lo2, (hsize_t)lo[1], (hsize_t)lo[0]};
        fsp.selectHyperslab(H5S_SELECT_SET, cnt, off);
        DataSpace msp(3, cnt);
        ds.read(buf.data(), PredType::NATIVE_DOUBLE, msp, fsp);
      } else {
        int si = sc_map.at(vname);
        hsize_t cnt[4] = {1, (hsize_t)nz_b, (hsize_t)ny_b, (hsize_t)nx_b};
        hsize_t off[4] = {
          (hsize_t)si, (hsize_t)lo2, (hsize_t)lo[1], (hsize_t)lo[0]};
        fsp.selectHyperslab(H5S_SELECT_SET, cnt, off);
        hsize_t mcnt[3] = {(hsize_t)nz_b, (hsize_t)ny_b, (hsize_t)nx_b};
        DataSpace msp(3, mcnt);
        ds.read(buf.data(), PredType::NATIVE_DOUBLE, msp, fsp);
      }

      Array4<Real> const& arr = data.array(mfi);
      int lo0 = lo[0], lo1 = lo[1], nxb = nx_b, nyb = ny_b;
      amrex::ParallelFor(
        box, [arr, &buf, n, lo0, lo1, lo2, nxb,
              nyb] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          arr(i, j, k, n) =
            buf[(k - lo2) * nyb * nxb + (j - lo1) * nxb + (i - lo0)];
        });
    }
  }
}

// ---------------------------------------------------------------------------
// Flat layout (e.g. CIAO post-processing "*_plot.h5" files):
//   /<group>/<name>      3-D datasets (nz, ny, nx), C order, x fastest
//   coords/xm, ym, zm    cell-centre coordinates (N values for N cells)
//   flow/time1           time, 1-element dataset
// x and z (or x and y) must be uniform; at most one of x, y may be stretched,
// which is linearly interpolated onto a uniform spacing. The file is read in
// contiguous z-slabs (one per rank, efficient for large files), interpolated
// within the slab and then redistributed into max_grid_size boxes.
// ---------------------------------------------------------------------------

// 1-D linear interpolation map from source cell centres onto the target cell
// centres lo + (i + 1/2) h, i = 0..n-1 (constant extrapolation at the ends)
struct InterpMap
{
  bool identity = true;
  Vector<int> i0;  // left source index
  Vector<Real> w;  // weight of source index i0 + 1
};

static Vector<Real>
read_1d(const H5File& file, const std::string& path)
{
  if (!h5_exists(file.getId(), path.c_str()))
    Abort("Dataset not found: " + path);
  DataSet ds = file.openDataSet(path);
  hsize_t dim[1] = {};
  if (ds.getSpace().getSimpleExtentNdims() != 1)
    Abort("Expected a 1-D dataset: " + path);
  ds.getSpace().getSimpleExtentDims(dim);
  Vector<Real> v(dim[0]);
  ds.read(v.data(), PredType::NATIVE_DOUBLE);
  return v;
}

static bool
is_uniform(const Vector<Real>& c, Real tol)
{
  const int n = (int)c.size();
  if (n < 3)
    return true;
  const Real h = (c[n - 1] - c[0]) / (n - 1);
  for (int i = 1; i < n; ++i)
    if (std::abs((c[i] - c[i - 1]) - h) > tol * std::abs(h))
      return false;
  return true;
}

static InterpMap
build_interp(const Vector<Real>& c, Real lo, Real h, int n)
{
  InterpMap m;
  m.identity = false;
  m.i0.resize(n);
  m.w.resize(n);
  const int ns = (int)c.size();
  for (int i = 0; i < n; ++i) {
    const Real xc = lo + (i + 0.5) * h;
    if (xc <= c[0]) {
      m.i0[i] = 0;
      m.w[i] = 0.0;
    } else if (xc >= c[ns - 1]) {
      m.i0[i] = ns - 2;
      m.w[i] = 1.0;
    } else {
      int k = (int)(std::upper_bound(c.begin(), c.end(), xc) - c.begin()) - 1;
      k = std::min(std::max(k, 0), ns - 2);
      m.i0[i] = k;
      m.w[i] = (xc - c[k]) / (c[k + 1] - c[k]);
    }
  }
  return m;
}

static void
convert_flat(const H5File& file, ParmParse& pp, const std::string& outfile)
{
  // ---- Fields and output names
  Vector<std::string> vars, names;
  if (pp.countval("vars") == 0)
    Abort("layout=flat requires vars=\"group/dataset ...\"");
  pp.getarr("vars", vars);
  if (pp.countval("names") > 0) {
    pp.getarr("names", names);
    if (names.size() != vars.size())
      Abort("names must have as many entries as vars");
  } else {
    for (auto v : vars) {
      std::replace(v.begin(), v.end(), '/', '_');
      names.push_back(v);
    }
  }

  // ---- Coordinates
  Vector<std::string> cpath = {"coords/xm", "coords/ym", "coords/zm"};
  if (pp.countval("coords") > 0)
    pp.getarr("coords", cpath, 0, 3);
  Real tol = 1.0e-3;
  pp.query("uniform_tol", tol);

  Vector<Real> c[3];
  for (int d = 0; d < 3; ++d)
    c[d] = read_1d(file, "/" + cpath[d]);
  int nsrc[3];
  bool uni[3];
  for (int d = 0; d < 3; ++d) {
    nsrc[d] = (int)c[d].size();
    uni[d] = is_uniform(c[d], tol);
    if (nsrc[d] < 2 && !(AMREX_SPACEDIM == 2 && d == 2))
      Abort("Need at least 2 cells in direction " + std::to_string(d));
  }
#if AMREX_SPACEDIM == 2
  if (nsrc[2] != 1)
    Abort("DIM=2 but file has nz=" + std::to_string(nsrc[2]));
  uni[2] = true;
#endif
  if (!uni[2])
    Abort("z must be uniform (only x or y can be interpolated)");
  if (!uni[0] && !uni[1])
    Abort("At most one of x and y may be non-uniform");

  // ---- Target grid: uniform directions keep their cells, a non-uniform
  // direction is resampled onto interp_dx over interp_range
  Real plo[3], phi[3], h[3];
  int ntgt[3];
  InterpMap map[3];
  int dint = -1;
  for (int d = 0; d < 3; ++d) {
    if (uni[d]) {
      h[d] = (nsrc[d] > 1) ? (c[d][nsrc[d] - 1] - c[d][0]) / (nsrc[d] - 1)
                           : 1.0;
      plo[d] = c[d][0] - 0.5 * h[d];
      phi[d] = c[d][nsrc[d] - 1] + 0.5 * h[d];
      ntgt[d] = nsrc[d];
    } else {
      dint = d;
    }
  }
  if (dint >= 0) {
    const auto& cc = c[dint];
    const int ns = nsrc[dint];
    Real hi_ = h[dint == 0 ? 1 : 0];
    pp.query("interp_dx", hi_);
    Real lo = cc[0] - 0.5 * (cc[1] - cc[0]);
    Real hi = cc[ns - 1] + 0.5 * (cc[ns - 1] - cc[ns - 2]);
    if (pp.countval("interp_range") > 0) {
      Vector<Real> r;
      pp.getarr("interp_range", r, 0, 2);
      lo = r[0];
      hi = r[1];
    }
    // Align the target cells with the finest source cells: where the source
    // spacing equals interp_dx (e.g. a uniform core) the data is copied
    // exactly instead of being interpolated between shifted cells
    int align = 1;
    pp.query("interp_align", align);
    if (align) {
      int kmin = 0;
      for (int k = 1; k < ns - 1; ++k)
        if (cc[k + 1] - cc[k] < cc[kmin + 1] - cc[kmin])
          kmin = k;
      const Real face = cc[kmin] - 0.5 * hi_;
      lo = face - std::floor((face - lo) / hi_) * hi_;
    }
    const int n = std::max(1, (int)std::floor((hi - lo) / hi_ + 1.0e-6));
    h[dint] = hi_;
    plo[dint] = lo;
    phi[dint] = lo + n * hi_;
    ntgt[dint] = n;
    map[dint] = build_interp(cc, lo, hi_, n);
    Real hmin = cc[1] - cc[0], hmax = hmin;
    for (int i = 2; i < ns; ++i) {
      hmin = std::min(hmin, cc[i] - cc[i - 1]);
      hmax = std::max(hmax, cc[i] - cc[i - 1]);
    }
    Print() << "Direction " << dint << " is non-uniform (spacing " << hmin
            << " ... " << hmax << "): interpolated onto " << n
            << " cells of " << hi_ << " in [" << plo[dint] << ", "
            << phi[dint] << "]\n";
  }
  Print() << "Grid: " << ntgt[0] << " x " << ntgt[1] << " x " << ntgt[2]
          << " cells, dx = " << h[0] << " " << h[1] << " " << h[2] << "\n";

  // ---- Time
  std::string tpath = "flow/time1";
  pp.query("time_dset", tpath);
  Real time = 0.0;
  {
    if (!h5_exists(file.getId(), ("/" + tpath).c_str()))
      Abort("Time dataset not found: " + tpath);
    DataSet ds = file.openDataSet("/" + tpath);
    if (ds.getSpace().getSimpleExtentNpoints() != 1)
      Abort("Time dataset must hold one value: " + tpath);
    double t = 0.0;
    ds.read(&t, PredType::NATIVE_DOUBLE);
    time = (Real)t;
  }
  Print() << "Time = " << time << " (" << tpath << ")\n";

  // ---- Geometry and box arrays
  Vector<int> per(3, 0);
  pp.queryarr("per", per, 0, 3);
  int max_grid_size = 32;
  pp.query("max_grid_size", max_grid_size);

  Box domain(
    IntVect::Zero, IntVect(AMREX_D_DECL(ntgt[0] - 1, ntgt[1] - 1, ntgt[2] - 1)));
  RealBox rb(
    {AMREX_D_DECL(plo[0], plo[1], plo[2])},
    {AMREX_D_DECL(phi[0], phi[1], phi[2])});
  Geometry geom(domain, &rb, 0, per.data());

  BoxArray ba(domain);
  ba.maxSize(max_grid_size);
  DistributionMapping dm(ba);
  MultiFab data(ba, dm, (int)vars.size(), 0);

  // Slabs along z that cover full x-y planes: one contiguous hyperslab read
  // per slab and field (2-D: a single slab)
  const int nprocs = ParallelDescriptor::NProcs();
  int ks = 1;
  BoxArray ba_slab(domain);
#if AMREX_SPACEDIM == 3
  {
    ks = std::max(1, (ntgt[2] + nprocs - 1) / nprocs);
    IntVect msz = domain.length();
    msz[2] = ks;
    ba_slab.maxSize(msz);
  }
#endif
  Vector<int> pmap(ba_slab.size());
  for (int i = 0; i < (int)pmap.size(); ++i)
    pmap[i] = i % nprocs;
  DistributionMapping dm_slab(pmap);
  MultiFab slab(ba_slab, dm_slab, 1, 0);
  Print() << ba.size() << " boxes, " << ba_slab.size() << " read slabs of "
          << ks << " planes on " << nprocs << " ranks\n";

  for (int n = 0; n < (int)vars.size(); ++n) {
    const std::string path = "/" + vars[n];
    if (!h5_exists(file.getId(), path.c_str()))
      Abort("Dataset not found: " + vars[n]);
    DataSet ds = file.openDataSet(path);
    DataSpace fsp = ds.getSpace();
    const int rank = fsp.getSimpleExtentNdims();
    if (rank != 3)
      Abort("Expected a 3-D dataset: " + vars[n]);
    hsize_t dims[3] = {};
    fsp.getSimpleExtentDims(dims);
    if (
      (int)dims[0] != nsrc[2] || (int)dims[1] != nsrc[1] ||
      (int)dims[2] != nsrc[0])
      Abort(
        "Dataset " + vars[n] + " does not match the coordinates (nz, ny, nx)");
    Print() << "Reading " << vars[n] << " -> " << names[n] << "\n";

    for (MFIter mfi(slab); mfi.isValid(); ++mfi) {
      const Box& box = mfi.validbox();
      const auto lo = lbound(box);
      const auto hi = ubound(box);
      const int klo = lo.z;
      const int nk = hi.z - lo.z + 1;
      const hsize_t sx = nsrc[0], sy = nsrc[1];

      Vector<double> buf((size_t)nk * sy * sx);
      hsize_t cnt[3] = {(hsize_t)nk, sy, sx};
      hsize_t off[3] = {(hsize_t)klo, 0, 0};
      DataSpace fs = ds.getSpace();
      fs.selectHyperslab(H5S_SELECT_SET, cnt, off);
      DataSpace msp(3, cnt);
      ds.read(buf.data(), PredType::NATIVE_DOUBLE, msp, fs);

      Array4<Real> const& arr = slab.array(mfi);
      for (int k = lo.z; k <= hi.z; ++k) {
        const size_t koff = (size_t)(k - klo) * sy * sx;
        for (int j = lo.y; j <= hi.y; ++j) {
          int j0 = j, j1 = j;
          Real wy = 0.0;
          if (!map[1].identity) {
            j0 = map[1].i0[j];
            j1 = j0 + 1;
            wy = map[1].w[j];
          }
          const double* r0 = &buf[koff + (size_t)j0 * sx];
          const double* r1 = &buf[koff + (size_t)j1 * sx];
          for (int i = lo.x; i <= hi.x; ++i) {
            int i0 = i, i1 = i;
            Real wx = 0.0;
            if (!map[0].identity) {
              i0 = map[0].i0[i];
              i1 = i0 + 1;
              wx = map[0].w[i];
            }
            const Real a = (1.0 - wx) * r0[i0] + wx * r0[i1];
            const Real b = (1.0 - wx) * r1[i0] + wx * r1[i1];
            arr(i, j, k) = (1.0 - wy) * a + wy * b;
          }
        }
      }
    }
    data.ParallelCopy(slab, 0, n, 1);
  }
  slab.clear();

  int n_files = amrex::VisMF::GetNOutFiles();
  pp.query("n_files", n_files);
  amrex::VisMF::SetNOutFiles(n_files);

  Print() << "Writing " << outfile << "...\n";
  WriteSingleLevelPlotfile(outfile, data, names, geom, time, 0);
}

// Native NGA2/CIAO layout (see the documentation)
static void
convert_native(
  const H5File& file, ParmParse& pp, const std::string& outfile,
  int max_grid_size)
{
  std::string mesh = discover_mesh_group(file);

  Vector<std::string> vars;
  if (pp.countval("vars") > 0) {
    pp.getarr("vars", vars);
  } else {
    vars = list_fields(file, mesh);
    Print() << vars.size() << " fields:";
    for (const auto& v : vars)
      Print() << " " << v;
    Print() << "\n";
  }

  Real time = read_time(file, mesh);
  Vector<int> per = read_periodicity(file, mesh);
  int coord_sys = read_coord_sys(file, mesh);

  BoxArray ba;
  Geometry geom;
  setup_grid(file, mesh, per, coord_sys, max_grid_size, ba, geom);

  MultiFab data(ba, DistributionMapping(ba), (int)vars.size(), 0);
  read_fields(file, mesh, vars, data);

  // Cap the number of plotfile data files via the n_files option (AMReX)
  int n_files = amrex::VisMF::GetNOutFiles();
  pp.query("n_files", n_files);
  amrex::VisMF::SetNOutFiles(n_files);

  Print() << "Writing " << outfile << "...\n";
  WriteSingleLevelPlotfile(outfile, data, vars, geom, time, 0);
}

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  {
    const Real strt_time = amrex::second();

    if (argc < 2)
      print_usage(argc, argv);

    ParmParse pp;
    if (pp.contains("help"))
      print_usage(argc, argv);

    std::string infile;
    pp.get("infile", infile);

    std::string outfile = "plt_" + stem(infile);
    pp.query("outfile", outfile);

    int max_grid_size = 32;
    pp.query("max_grid_size", max_grid_size);

#ifdef AMREX_USE_MPI
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(
      fapl_id, ParallelDescriptor::Communicator(), MPI_INFO_NULL);
    H5File file(
      infile.c_str(), H5F_ACC_RDONLY, H5::FileCreatPropList::DEFAULT,
      H5::FileAccPropList(fapl_id));
    H5Pclose(fapl_id);
#else
    H5File file(infile.c_str(), H5F_ACC_RDONLY);
#endif

    std::string layout = "native";
    pp.query("layout", layout);
    if (layout == "flat") {
      convert_flat(file, pp, outfile);
    } else if (layout == "native") {
      convert_native(file, pp, outfile, max_grid_size);
    } else {
      Abort("layout must be native or flat");
    }

    Real elapsed = amrex::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(elapsed);
    Print() << "Done. Run time: " << elapsed << " s\n";
  }
  amrex::Finalize();
  return 0;
}
