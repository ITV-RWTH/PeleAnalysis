.. highlight:: bash


partScatter
***********

Sample structure functions from a plotfile. From each of a set of "seed"
locations, ``partScatter`` fires ``particlesPerLoc`` rays in evenly distributed
directions and marches a particle along each one in fixed steps. At every step
it interpolates the requested variables to the particle position and stores the
increment relative to the seed point, raised to a user-chosen exponent:

.. math::

   \Delta_n(\mathbf{x}, \mathbf{e}, r) =
     \big[\, q(\mathbf{x} + r\,\mathbf{e}) - q(\mathbf{x}) \,\big]^{n}

for variable :math:`q`, seed point :math:`\mathbf{x}`, unit ray direction
:math:`\mathbf{e}`, separation :math:`r` and exponent :math:`n` (the ``orders``
input). Averaging the output over seeds and ray directions gives the
angle-averaged structure function :math:`S_n(r)`; the raw per-ray samples are
written out so that the averaging, and any conditioning on the seed location,
can be done in post-processing.

Unlike :doc:`partStream`, the paths are straight lines rather than integral
curves of a vector field, so no vector field is read and no integration scheme
is involved. The tool shares ``partStream``'s seeding options and output
formats.

This tool must be compiled with ``DIM = 3``: the ray-direction generators are
genuinely three-dimensional and have no meaningful two-dimensional counterpart.

Usage: ::

   ./partScatter3d.gnu.MPI.ex infile=FILE vars=LIST orders=LIST is_per=I J K [OPTIONS]

Help: ::

   ./partScatter3d.gnu.MPI.ex -h

Example: ::

   ./partScatter3d.gnu.MPI.ex ./InputSamples/partScatter.inp


Tool Options
############

Input control
~~~~~~~~~~~~~

::

   #------------------- INPUT CONTROL --------------------------------------------------------
   infile = plt00000                          # Plot file to sample
   vars = x_velocity y_velocity z_velocity    # Variables to sample along the rays
   orders = 2 2 2                             # Exponent applied to each entry of vars

`vars` names the plotfile variables to sample; `orders` gives the exponent
applied to each increment and must contain exactly one value per entry of
`vars`. Exponent 2 gives the second-order structure function, 3 the third, and
so on. Different exponents may be used for different variables in the same run.
All AMR levels present in `infile` are read.

Ray options
~~~~~~~~~~~

::

   #------------------- RAY OPTIONS ----------------------------------------------------------
   Nsteps = 50                                # DEF: 50; Number of points recorded along each ray
   stepSize = 0.1                             # DEF: 0.1; Step as a fraction of the finest cell size
   particlesPerLoc = 100                      # DEF: 1; Number of rays fired from each seed location
   plane = -1                                 # DEF: -1; Sphere of directions, or circle normal to axis

Each ray records `Nsteps` points, the first of which is the seed point itself.
The step length is `stepSize` times the cell size on the finest level of the
plotfile, and must lie in :math:`[0, 0.5]`, so the largest separation sampled is
``(Nsteps-1) * stepSize`` finest-level cells.

`particlesPerLoc` sets how many rays leave each seed. With `plane` negative the
directions are spread over the unit sphere using a Fibonacci lattice, which
gives a near-uniform covering for any count and is the right choice for an
angle-averaged structure function. Setting `plane` to 0, 1 or 2 instead spreads
the directions evenly around the circle normal to that axis, which is useful
when the statistics are wanted in a particular plane, for example to separate
directions parallel and perpendicular to a mean shear or a magnetic field.

Seed points
~~~~~~~~~~~

::

   #------------------- SEED POINTS ----------------------------------------------------------
   oneSeedPerCell = 0                         # [0, 1], DEF: 0; One seed in each uncovered cell
   isoFile = plt00000_surf.mef                # DEF: None; One seed per node of the surface
   seedLoc = 0.1 0.0 0.3                      # DEF: None; A single seed at this location
   seedRakeNum = 10                           # DEF: None; Evenly spaced seeds along a line segment
   seedRakeL = 0.0 0.0 0.0                    # Left endpoint for seedRake
   seedRakeR = 0.5 0.5 0.5                    # Right endpoint for seedRake

Exactly one of the four seeding modes must be given; they behave as in
:doc:`partStream`:

1. `oneSeedPerCell` places one seed in each cell not covered by a finer level.
   This gets very expensive in 3D, and note the MPI caveat below.
2. `isoFile` places one seed at each node of an MEF surface. The surface
   connectivity is copied into the binary output header, so the samples can be
   related back to the surface in post-processing.
3. `seedLoc` places a single seed at the given coordinates.
4. `seedRakeNum` places that many evenly spaced seeds along the segment from
   `seedRakeL` to `seedRakeR`; it must be at least 2.

Grid control
~~~~~~~~~~~~

::

   #------------------- GRID CONTROL ---------------------------------------------------------
   is_per = 1 1 0                             # Sets case periodicity
   nGrow = 3                                  # DEF: 3; Grow cells. Must be >= 1

`is_per` is required and sets the periodicity used to build the geometry, one
flag per direction; it is independent of whatever the original simulation used.
Rays crossing a periodic boundary wrap around, so separations larger than the
domain will alias onto it. In a non-periodic direction, particles that leave the
domain are dropped and simply contribute no further samples, which means the
sample count falls off with separation near such a boundary.

Output formats
~~~~~~~~~~~~~~

::

   #------------------- OUTPUT CONTROL -------------------------------------------------------
   outfile = plt00000                         # DEF: infile; Name base of output files
   writeParticles = 0                         # [0, 1], DEF: 0; Write particles as plt file
   particlefile = plt00000_particles          # DEF: outfile + "_particles"
   writeStreams = 0                           # [0, 1], DEF: 0; Write rays in Tecplot ascii format
   streamfile = plt00000_stream               # DEF: outfile + "_stream"
   writeStreamBin = 1                         # [0, 1], DEF: 1; Write rays as binary
   streamBinfile = plt00000_streamBin         # DEF: outfile + "_streamBin"

Note that `writeStreamBin` defaults to 1 here, so the binary output is written
unless it is explicitly disabled.

Each output records ``AMREX_SPACEDIM + 1 + nVars`` quantities per point along a
ray: the coordinates ``X``, ``Y``, ``Z`` of the point, the separation ``R`` from
the seed, and then one increment per entry of `vars`. At the seed point the
increment is zero by construction and ``R`` is zero.

`writeStreamBin` produces a directory containing one ``str_<rank>.bin`` per MPI
rank plus a ``Header``. Each binary file starts with the number of rays it
contains, written as an ``int``; then, for each ray, a zero-based ray id
followed by the samples for that ray ordered component-major, that is all
``Nsteps`` values of the first quantity, then all ``Nsteps`` values of the
second, and so on. The ray id encodes the seed and the direction: seed index is
``id / particlesPerLoc`` and direction index is ``id % particlesPerLoc``. The
``Header`` carries the number of files, the number of seed locations, the rays
per location, the points per ray, the number of quantities and their names, and
the surface connectivity when `isoFile` was used.

`writeStreams` writes the same data as one Tecplot ASCII file per rank, with one
``ZONE`` per ray; this is convenient for inspecting a handful of rays but does
not scale.

.. note::

   When ``writeParticles`` is enabled, the particle plotfile is written through
   the AMReX particle I/O layer, whose data-file count is controlled by the
   native ``particles.particles_nfiles`` option (read directly by AMReX) rather
   than the ``n_files`` option used by the grid-based tools.

.. warning::

   ``oneSeedPerCell`` currently only seeds every cell when run on a single MPI
   rank. The seed list is built per rank from that rank's own boxes and is then
   strided again by rank, so on ``N`` ranks only about ``1/N`` of the cells are
   seeded. The other three seeding modes build the same list on every rank and
   distribute correctly. :doc:`partStream` shares this behaviour.

Memory
~~~~~~

Every particle carries its whole ray, so the storage is
``particlesPerLoc * nSeeds * Nsteps * (AMREX_SPACEDIM + 1 + nVars)`` reals. This
grows quickly: raising `particlesPerLoc` to improve the angular average costs
memory linearly. Prefer many seeds with a moderate number of directions over a
few seeds with very many.
