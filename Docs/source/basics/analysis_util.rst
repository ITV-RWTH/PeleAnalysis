.. highlight:: c++

analysis_util API Reference
***************************

``analysis_util`` is a C++ utility library that wraps the low-level AMReX
plotfile and MEF I/O routines behind a concise, MPI-safe interface.
All symbols live in the ``analysis_util`` namespace.
Include the header with:

::

   #include <analysis_util.H>

Enable the library in your ``GNUmakefile`` by setting ``USE_UTILS = TRUE``
(done automatically when ``EBASE = template``).


Data Structures
---------------

PlotfileData
~~~~~~~~~~~~

Returned by ``read_plotfile``.  Holds everything needed to work with
a multi-level dataset and to call ``write_plotfile``.

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Field
     - Type
     - Description
   * - ``mf``
     - ``Vector<MultiFab>``
     - One ``MultiFab`` per AMR level; components correspond to the requested variables.
   * - ``geoms``
     - ``Vector<Geometry>``
     - Domain geometry for each level.
   * - ``ref_ratios``
     - ``Vector<int>``
     - Refinement ratio between consecutive levels (length ``n_lev - 1``).
   * - ``time``
     - ``Real``
     - Simulation time stored in the plotfile.
   * - ``n_lev``
     - ``int``
     - Number of levels loaded.
   * - ``var_names``
     - ``Vector<string>``
     - All variable names present in the file (not just the ones requested).

MEFData
~~~~~~~

Returned by ``read_mef`` / passed to ``write_mef``.  Represents an
unstructured surface or polyline in Marc's Element Format.

.. list-table::
   :header-rows: 1
   :widths: 25 20 55

   * - Field
     - Type
     - Description
   * - ``title``
     - ``string``
     - Label stored on the first line of the file.
   * - ``var_names``
     - ``vector<string>``
     - Per-component variable names.
   * - ``nodes``
     - ``FArrayBox``
     - Node data: ``Box`` is ``(0:nNodes-1, 0:0, 0:0)`` with ``nComp`` components.
   * - ``connectivity``
     - ``Vector<int>``
     - Flat array of length ``n_elts * nodes_per_elt``; 1-based node indices.
   * - ``n_elts``
     - ``int``
     - Number of elements.
   * - ``nodes_per_elt``
     - ``int``
     - Nodes per element (e.g. 3 for triangles, 2 for line segments).


Functions
---------

init
~~~~

**What it does**

Populates the global state used by ``get_covered_mf``, ``integrate``,
and ``gradient``: the number of levels, box arrays, distribution maps,
geometries, periodicity, and refinement ratios.
Must be called once before using any of those three functions.
Calling it a second time with the same data is a no-op.

**Inputs**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``a_mf``
     - ``Vector<MultiFab>`` — one MultiFab per AMR level (typically ``data.mf``).
   * - ``a_geoms``
     - ``Vector<Geometry>`` — one Geometry per level (typically ``data.geoms``).

**Output**

None (modifies global state).

**Example**

::

   auto data = analysis_util::read_plotfile(infile, {"density"}, finestLevel);
   analysis_util::init(data.mf, data.geoms);


----

read_plotfile
~~~~~~~~~~~~~

**What it does**

Opens an AMReX plotfile, reads the requested variables into one
``MultiFab`` per AMR level, and returns a ``PlotfileData`` struct.
MPI-safe: redistributes boxes when the domain has fewer boxes than
MPI ranks to avoid ``FillVar`` deadlocks.

**Inputs**

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Description
   * - ``infile``
     - required
     - Path to the plotfile directory.
   * - ``var_names``
     - required
     - Variables to load; each name must exist in the file (aborts otherwise).
   * - ``finest_level``
     - ``1000``
     - Cap on the highest AMR level to read; clamped to the file's actual finest level.
   * - ``n_grow``
     - ``0``
     - Number of ghost cells to allocate in the returned ``MultiFab``.
   * - ``is_per``
     - ``{}`` (all 0)
     - Periodicity flags, one per spatial dimension.

**Output**

``PlotfileData`` — see structure description above.

**Example**

::

   Vector<std::string> vars = {"density", "velocity_x"};
   Vector<int> is_per(AMREX_SPACEDIM, 1);   // fully periodic domain
   auto data = analysis_util::read_plotfile("plt00100", vars, /*finestLevel=*/2,
                                            /*n_grow=*/1, is_per);
   Print() << "Loaded " << data.n_lev << " level(s), t = " << data.time << "\n";


----

write_plotfile
~~~~~~~~~~~~~~

**What it does**

Writes a multi-level set of ``MultiFab`` objects to an AMReX plotfile.
MPI-safe: adds ``Barrier`` calls before and after the collective write,
and redistributes boxes when the domain has fewer boxes than MPI ranks.

**Inputs**

.. list-table::
   :header-rows: 1
   :widths: 30 20 50

   * - Parameter
     - Default
     - Description
   * - ``outfile``
     - required
     - Output directory name (created by AMReX).
   * - ``mf``
     - required
     - ``Vector<MultiFab>`` to write, one per level.
   * - ``var_names``
     - required
     - Variable name for each component; length must equal ``mf[0].nComp()``.
   * - ``geoms``
     - required
     - Domain geometry for each level.
   * - ``time``
     - ``0.0``
     - Simulation time to embed in the plotfile header.
   * - ``ref_ratios``
     - ``{}`` (all 2)
     - Refinement ratio between consecutive levels (length ``n_lev - 1``).

**Output**

None (writes to disk).

**Example**

::

   const std::string outfile = analysis_util::get_file_root(infile) + "_result";
   analysis_util::write_plotfile(outfile, data.mf, {"density", "velocity_x"},
                                 data.geoms, data.time, data.ref_ratios);


----

read_mef
~~~~~~~~

**What it does**

Reads a binary MEF file from disk and returns a ``MEFData`` struct
containing the title, variable names, node data, and element connectivity.

**Inputs**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``infile``
     - Path to the MEF file.

**Output**

``MEFData`` — see structure description above.

**Example**

::

   auto mef = analysis_util::read_mef("surface.mef");
   Print() << "Title: " << mef.title << "\n";
   Print() << mef.n_elts << " elements, " << mef.nodes_per_elt << " nodes/elt\n";


----

write_mef
~~~~~~~~~

**What it does**

Writes a ``MEFData`` struct to disk in binary MEF format.
Only the IO rank writes; all ranks synchronize via a ``Barrier``
before returning so that a subsequent ``read_mef`` on any rank sees
the complete file.

**Inputs**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``outfile``
     - Output file path.
   * - ``data``
     - ``MEFData`` to write.

**Output**

None (writes to disk).

**Example**

::

   analysis_util::MEFData mef;
   mef.title        = "my surface";
   mef.var_names    = {"x", "y", "z", "temperature"};
   mef.n_elts       = nTri;
   mef.nodes_per_elt = 3;
   // ... fill mef.nodes and mef.connectivity ...
   analysis_util::write_mef("surface_out.mef", mef);


----

get_covered_mf
~~~~~~~~~~~~~~

**What it does**

Builds an integer mask ``MultiFab`` at each AMR level.
A cell is marked ``1`` (uncovered) if it is not covered by a finer
level, and ``0`` (covered) otherwise.
Used internally by ``integrate`` to avoid double-counting.

**Inputs**

None (uses global state set by ``init``).

**Output**

``Vector<unique_ptr<iMultiFab>>`` — one mask MultiFab per level.

**Example**

::

   analysis_util::init(data.mf, data.geoms);
   auto mask = analysis_util::get_covered_mf();
   // mask[lev] == 1 where lev is the finest data at that location


----

integrate
~~~~~~~~~

**What it does**

Computes the volume integral of one or more components across all AMR
levels, correctly excluding fine-covered cells.
Partial integrals along selected axes are supported: e.g. integrating
over x and z leaves a 1-D array indexed by y.
Result is reduced across all MPI ranks before returning.

**Inputs**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``a_mf``
     - ``Vector<MultiFab>`` — data to integrate (one per level).
   * - ``scomp``
     - First component index to integrate.
   * - ``ncomp``
     - Number of components to integrate, starting at ``scomp``.
   * - ``axes_to_integrate``
     - ``Vector<int>`` of axis indices to sum over: ``0``=x, ``1``=y, ``2``=z. Pass all axes for a scalar result.

**Output**

``unique_ptr<Gpu::ManagedVector<Real>>`` of length
``ncomp * (Nx if x not integrated, else 1) * (Ny if y not integrated, else 1) * (Nz if z not integrated, else 1)``,
where Nx/Ny/Nz are the cell counts of the finest level.

**Example — scalar volume integral**

::

   analysis_util::init(data.mf, data.geoms);
   // Integrate component 0 over all axes -> one value per component
   auto result = analysis_util::integrate(data.mf, /*scomp=*/0, /*ncomp=*/1,
                                          {0, 1, 2});
   Print() << "Integral = " << (*result)[0] << "\n";

**Example — line-average along y**

::

   // Integrate over x and z; result is an array indexed by y
   auto profile = analysis_util::integrate(data.mf, 0, 1, {0, 2});
   // (*profile)[j] is the x-z integral at y-index j


----

gradient
~~~~~~~~

**What it does**

Computes the cell-centered gradient of ``ncomp`` scalar fields using
AMReX's ``MLMG`` / ``MLPoisson`` solver with 4th-order accuracy.
Boundary conditions are set automatically from the periodicity and
an optional ``sym_dir`` ParmParse flag for symmetric (reflect-odd)
directions.

**Inputs**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Parameter
     - Description
   * - ``a_mf``
     - ``Vector<MultiFab>`` — source data; **must have at least 1 ghost cell**.
   * - ``scomp``
     - First component to differentiate.
   * - ``ncomp``
     - Number of scalar fields to differentiate.
   * - ``grad_mf``
     - ``Vector<MultiFab>`` — output; must have ``AMREX_SPACEDIM * ncomp`` components.

**Output**

None (fills ``grad_mf`` in-place).
Components are ordered as ``(df0/dx, df0/dy, df0/dz, df1/dx, ...)``.

**Example**

::

   analysis_util::init(data.mf, data.geoms);
   // Read with n_grow=1 so gradient has ghost cells
   auto data = analysis_util::read_plotfile(infile, {"temperature"},
                                            finestLevel, /*n_grow=*/1, is_per);
   analysis_util::init(data.mf, data.geoms);

   Vector<MultiFab> grad_mf(data.n_lev);
   for (int lev = 0; lev < data.n_lev; ++lev)
     grad_mf[lev].define(data.mf[lev].boxArray(),
                         data.mf[lev].DistributionMap(),
                         AMREX_SPACEDIM, 0);
   analysis_util::gradient(data.mf, /*scomp=*/0, /*ncomp=*/1, grad_mf);
   // grad_mf[lev] now holds (dT/dx, dT/dy [, dT/dz])


----

get_file_root
~~~~~~~~~~~~~

**What it does**

Strips the directory prefix from a file path, returning only the
final component.

**Input**

Path string (e.g. ``"runs/case1/plt00100"``).

**Output**

Root filename string (e.g. ``"plt00100"``).

**Example**

::

   std::string root = analysis_util::get_file_root("runs/case1/plt00100");
   // root == "plt00100"
   const std::string outfile = root + "_processed";


----

find_var_index
~~~~~~~~~~~~~~

**What it does**

Searches a list of variable names for a given name and returns its
index.  Aborts with a helpful message if the name is not found
(unless ``abort_if_not_found = false``, in which case it returns ``-1``).

**Inputs**

.. list-table::
   :header-rows: 1
   :widths: 35 15 50

   * - Parameter
     - Default
     - Description
   * - ``var_names``
     - required
     - List of names to search (e.g. ``data.var_names``).
   * - ``name``
     - required
     - Variable name to find.
   * - ``abort_if_not_found``
     - ``true``
     - If ``true``, aborts when the name is absent; if ``false``, returns ``-1``.

**Output**

``int`` — 0-based component index, or ``-1`` if not found and
``abort_if_not_found = false``.

**Example**

::

   int iTemp = analysis_util::find_var_index(data.var_names, "temperature");
   // Use iTemp as scomp in integrate() or as an array index into data.mf


----

parse_title / parse_var_names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**What they do**

Low-level stream parsers for reading the header section of MEF-style
binary files.

- ``parse_title(is)`` — reads and returns one line from ``is`` as the title string.
- ``parse_var_names(is)`` — reads the next line from ``is`` and splits it on
  spaces and commas, returning the tokens as a ``vector<string>``.

These are called internally by ``read_mef``; you only need them when
writing a custom file reader that follows the same header convention.

**Example**

::

   std::ifstream ifs("custom.mef");
   std::string title        = analysis_util::parse_title(ifs);
   auto        var_names    = analysis_util::parse_var_names(ifs);


----

Complete minimal example
------------------------

The following is a minimal tool that reads a plotfile, doubles the
``density`` field in-place, and writes the result as a new plotfile.
It mirrors the pattern used in ``Src/template.cpp``.

::

   #include <AMReX_AmrData.H>
   #include <AMReX_MultiFab.H>
   #include <AMReX_ParmParse.H>
   #include <analysis_util.H>

   using namespace amrex;

   int main(int argc, char* argv[])
   {
     Initialize(argc, argv);
     {
       ParmParse pp;
       std::string infile; pp.get("infile", infile);

       Vector<std::string> vars = {"density"};
       Vector<int> is_per(AMREX_SPACEDIM, 1);

       // 1. Read
       auto data = analysis_util::read_plotfile(infile, vars,
                                                /*finestLevel=*/1000,
                                                /*n_grow=*/0, is_per);
       analysis_util::init(data.mf, data.geoms);

       // 2. Process (in-place kernel)
       for (int lev = 0; lev < data.n_lev; ++lev) {
   #ifdef AMREX_USE_OMP
   #pragma omp parallel if (Gpu::notInLaunchRegion())
   #endif
         for (MFIter mfi(data.mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
           const Box& bx = mfi.tilebox();
           auto const& a = data.mf[lev].array(mfi);
           amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
             a(i, j, k, 0) *= 2.0;   // double the density
           });
         }
       }

       // 3. Write
       const std::string outfile = analysis_util::get_file_root(infile) + "_doubled";
       analysis_util::write_plotfile(outfile, data.mf, vars, data.geoms, data.time);
       Print() << "Wrote " << outfile << "\n";
     }
     Finalize();
     return 0;
   }
