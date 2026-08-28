# Changelog

All notable changes to FVQual are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added — non-uniform lateral refinement, and a better-controlled answer

Uniform refinement refines the wrong place. `planLateralCells(maxLanes,
minCellWidth)` chooses lanes per segment so no cell falls below a stated width:
DeGray's 79 m tributary arms keep one lane while the 5530 m water near the dam —
where the overturning cell under investigation lives — gets up to eight. The mesh
is non-uniform but still **conforming**: where the lane count changes, the shared
boundary carries the union of both sides' lane positions and the coarser cell is
a polygon walking all of them, so every edge has exactly two faces and no node
hangs. Up to 7 nodes per face on DeGray, 41 such faces at a 250 m minimum.

**The cost this buys.** At a 700 m minimum DeGray becomes 67 columns and needs
**7218** internal steps against **6770** at one lane — 6.6% more. Uniform two
lanes needs 3.1×, for half the columns.

**And it changes the answer.** Volume-weighted warming below 20 m at the dam,
days 64.5–79.5, every run through the same mesh path so the jitter cancels:

| | warming below 20 m | vs control |
|---|---|---|
| one lane (laterally averaged) | +1.724 °C | — |
| adaptive, 67 columns, ≤7 lanes in the wide water | +1.695 °C | **−1.7%** |
| uniform two lanes | +1.972 °C | +14.4% |

The mesh-path jitter — the same model through two builders that differ only in
face-node ordering — is about **0.08 °C**. The adaptive difference is 0.029 °C,
*below* that floor: refining the wide water where the overturning cell lives
produces **no detectable change** in how much heat reaches the hypolimnion. The
uniform result, 0.25 °C, is about three times the floor and is now the suspect
number: uniform two lanes cuts the narrow arms to 39 m cells, and that is where
its extra transport is likely to come from rather than from anything about
lateral averaging near the dam.

So the schematisation hypothesis remains unsupported, and on a better-controlled
test than the previous one. Still fifteen days of early spring against a bias
that develops over months.

### Measured — widening DeGray does not reduce its interior heat transport

The residual +3 to +4 °C interior warm bias has been attributed to the
laterally-averaged schematisation: one segment wide, so every return path is
forced through the vertical. A synthetic basin supported that — 39% less heat
below 20 m at eight cells wide. `toLateralMeshDefinition()` makes the same test
possible on DeGray's own bathymetry, inventing no survey data.

**It does not reproduce.** One lane against two, identical forcing, profiles taken
at the deepest lane of the dam segment so both sample the same 54.5 m column:

| | warming below 20 m, days 64.5–79.5 |
|---|---|
| one lane (laterally averaged) | +1.751 °C |
| two lanes | +1.972 °C |

Two lanes carries **12.7% more** heat below 20 m, not less. Day-to-day scatter in
the difference is roughly half the signal, so the sign is more trustworthy than
the magnitude — but it is the opposite sign from the hypothesis.

**What this does not settle.** Two lanes over fifteen days of early spring is a
weak test of a claim about eight cells over a season, and the bias it is meant to
explain develops over months. The honest statement is that the schematisation
hypothesis is unsupported so far, not refuted.

**Why the stronger test was not run.** DeGray's narrowest active segment is 79 m
wide, so four lanes makes 19.8 m cells against kilometre-long segments and the
barotropic Courant limit collapses; two lanes already needs 3.1× the internal
steps for the same period. Uniform lateral refinement refines the narrow upstream
segments hardest, which is exactly where it is least wanted — the overturning
cell in question lives in the wide water near the dam. Testing this properly
needs *non-uniform* lateral refinement, and that is the next piece of work rather
than more of this one.

### Added — F1.5 conveyance: momentum reconciled with the hypsographic depth

Hypsography corrected *storage* in F1.5; momentum was left behind. Manning's drag
goes as `h^(-4/3)`, and the depth it was handed is `m_h = eta − bed`, the depth at
the column's **deepest point**. Almost none of a natural section is that deep, so
using it under-estimates friction, and worst where the hypsography is steepest.

- `BarotropicSolver::conveyanceDepth()` — the single depth carrying the section's
  actual conveyance, derived from the existing storage curve with **no free
  parameter**: `h_c = [(1/A) ∫ (eta − z)^{5/3} dA]^{3/5}`, the compound-section
  conveyance integral. A prismatic column returns its full depth exactly, so
  models without a curve are bit-for-bit unchanged. The three depths order as
  `V/A < h_c < h_max`, and the correction moves momentum from the top of that
  range to the middle. For DeGray's shape `h_c ≈ 0.39 H`, so the drag coefficient
  rises about fourfold.

**On DeGray that fourfold change moves nothing, and that was worth measuring.**
The commit that introduced it warned that "any hypsographic result predating this
will move." Re-run to day 150 against the 24 measured casts, with the friction
wiring the only difference between the two runs:

| | RMSE | 5–10 m | 20–30 m | 30–60 m |
|---|---|---|---|---|
| friction on `h_max` (before) | 3.59 °C | −1.11 | +3.22 | +3.65 |
| friction on `h_c` (after) | 3.59 °C | −1.11 | +3.22 | +3.65 |

Every band identical to two decimals over 333 measurements. Whole-run maxima:
surface temperature **0.035 °C** (0.26%), bottom temperature 0.039 °C, volume
0.001%, heat content 0.018%. The daily forcing columns differ too, which is the
tell: the two runs diverge only by step-sequence jitter, not by physics.

The reason is that DeGray is an impoundment, not a river. Its barotropic
velocities are millimetres per second and Manning drag goes as `u²`, so the term
is negligible whatever depth it is given. The correction is real, and it will
matter in a riverine or shallow domain where friction actually controls the flow.
It does not matter here. The baseline run reproducing the historical 3.59 °C
table exactly is what makes that a comparison rather than a coincidence.
- Wired into the Manning term only. Continuity, the pressure gradient and the
  wave speeds still see `m_h`, because those depend on how much water is present
  and where its surface is, not on how well the section carries it.

**The bank pressure term is not needed, and that is a measured result.** The plan
named it, and in a porosity-form scheme it is mandatory: there the pressure flux
carries a porosity factor, its divergence does not vanish under a flat surface,
and `g h²/2 ∇φ` is what cancels it. FVQual is not written in that form —
hypsography enters through the volume-to-stage relation, so the flux terms never
see the curve. A lake at rest over an uneven bed, with differently shaped curves
*and* differing surface areas column to column, stays exact to **1e-13** with no
extra term. Adding one would introduce an imbalance rather than remove one.

The gate earns that conclusion by being shown to fail: scaling the pressure
correction by each column's storage ratio — putting the scheme into porosity form
without its companion term — produces a **0.587 m/s** spurious current, five
orders of magnitude past tolerance. It stands as the guard on anyone later moving
the flux terms to porosity form and forgetting what has to come with it.

**A correction to the first commit of this work.** Its header claimed the depth
handed to Manning was the storage average `V/A`, and concluded the code
*over*-estimated friction. `m_h` is the maximum depth, so the sign was backwards:
the code under-estimated it, and the fix increases drag rather than reducing it.
Found by reading how `m_h` is actually assigned instead of trusting the comment.

Eight gates. The two shapes with closed forms are checked against the Euler beta
function via `tgamma` — a different route from the quadrature under test. Three
ablations fail them: the 5/3 exponent replaced by 1, the bed-area term dropped,
and the friction wiring reverted. That last one is the informative one: with it
reverted, a shelved basin and a prismatic one produce **bit-identical** speeds,
which is the defect stated as a measurement.

### Investigated — the vertical Courant "self-reference" does not reproduce

Carried as an open defect: `diagnoseVerticalFlux` closes layer continuity with
`(V_now − V_old)/dt`, and `stableTimeStep()` reads the result to size the *next*
step, so a noisy volume difference divided by a vanishing `dt` could inflate the
diagnosed velocity and drive the step down further. Measured rather than assumed
(`verification/interfaces/vertical_cfl_selfreference.cpp`):

| | at rest | tilted surface, 600 s of seiche |
|---|---|---|
| dt probed | 100 s → 1e-4 s | cap 10 s → 0.25 s |
| max \|w\| | 1.2e-16 → 7.7e-39 m/s | 1.89e-6 → 1.78e-6 m/s |
| behaviour | *falls* with dt | flat to **1.06×** |

At rest the volume change is exactly zero — the stage does not move, so the two
volumes are bitwise equal — and there is nothing to amplify. With real motion the
ratio converges: a fortyfold change in step size moves it by six percent. The lag
is real but benign. **No code change**, because nothing is broken.

The measurement took two corrections worth recording. The first tilted run took a
single step from rest, where nothing has accelerated, so it compared round-off
against round-off, reported a 16000× "spread", and concluded the defect was real
— entirely an artefact. And the step cap has to be shown to bind, or four runs
take identical steps and the table compares one calculation with itself; the
achieved step is now printed beside each cap.

### Added — F8.2c–f, coupling semantics carried across the boundary

Three interfaces added upstream to HydroCouple, all optional and discovered by
`dynamic_cast` rather than added to existing interfaces: these headers cross a
plugin boundary, and a new virtual on an interface changes a vtable that
already-compiled components — `CEQUALW2Component` among them — were built
against.

- `IValueSemantics` / `ValueKind` (F8.2c) — intensive, extensive, flux, density.
- `ITemporalSemantics` / `TimeKind` (F8.2d) — instantaneous, interval mean, min,
  max, accumulated, with the interval length.
- `IVerticalCoordinate` / `VerticalCoordinateKind` (F8.2e) — sigma, z-level,
  z-star, hybrid, isopycnal, depth-below-surface. `ILayeredMeshComponentDataItem`
  moved upstream alongside it and now carries a vertical coordinate rather than
  three CF formula terms.

All three default to `Unknown`, and each header states that a consumer must not
guess from it. An `Unknown` that quietly meant "instantaneous" would recreate the
problem the enum exists to solve. The motivating defects are this project's own:
cloud cover in tenths read as a fraction, and a daily report drifting from noon
to five in the afternoon — both two sides disagreeing about what a number meant.

- **The vertical coordinate answers with elevations, not formula terms.**
  Publishing terms is right for a file, where a reader implements a handful of CF
  formulae once. At a coupling boundary it obliges every partner to redo the
  arithmetic, and a partner that gets the signs wrong produces a plausible
  profile and no error at all.
- `Component::SigmaVerticalCoordinate` implements it over `LayeredMesh`,
  delegating to the mesh's own `z(column, k)` rather than rederiving
  `eta + sigma·(bed − eta)`. That distinction is not cosmetic: with a minimum
  layer thickness enforced, `SigmaGridGenerator` lifts thin layers to the floor
  and takes the surplus from the thick ones, and on a shoaling arm the naive
  formula puts the top layer's base at −0.006 m where the mesh has it at
  −0.10 m — a sixteenfold error in the surface layer, which is the one a partner
  reads surface temperature from.
- Eight gates, all written through the `IVerticalCoordinate` base so they assert
  only what a partner can rely on. Each was confirmed to fail under ablation.
  The gate pinning the delegation initially used uniform layering, measured a
  disagreement of 3.6e-15, and failed rather than passing: under uniform sigma
  the thickness floor binds on every layer at once or none, so the mechanism it
  claimed to test was never engaged. It now uses stretched layering and the
  comment records why.
- `ICrossSection` / `CrossSectionKind` and `ILayeredNetworkComponentDataItem`
  (F8.2f) — the W2 / HEC-RAS / SWMM shape: a network resolved vertically at each
  station, with a channel section there. `ILayering` factored out of the mesh
  item and shared with the network item, so a consumer wanting values by layer
  with trustworthy elevations can `dynamic_cast` to it without learning whether
  it was handed a 3-D mesh or a laterally-averaged network — which is the
  FVQual-to-W2 exchange this hierarchy exists for.
- **The cross-section answers derived quantities, not survey points**, for a
  sharper reason than the vertical coordinate had: the three model families do
  not store the same thing. HEC-RAS and SWMM hold station-elevation pairs;
  CE-QUAL-W2 holds *width against elevation* and has no stations at all (see
  `W2Bathymetry::width(layer, segment)`); a closed conduit has a formula and
  neither. An interface built on points would have forced two of the three to
  invent data they do not have.
- `verification/interfaces/crosssection_agreement.cpp` implements the same
  trapezoid in both representations and compares them against the closed form at
  every stage from invert to top of bank. They agree to **3.6e-15**. Three
  ablations confirmed it can fail: rectangular instead of trapezoidal
  integration (64 failures), flow area collapsed onto storage area (1), and a
  perimeter ignoring the bank slope (32). Kept out of `src/` deliberately —
  FVQual is a 3-D mesh model and has no use for a cross-section class of its own.
- Flow area and storage area are separate and provably not collapsible: on the
  test section 48 of 80 m² conveys, 60%.
- `SDK::Spatial::LayeredMeshComponentDataItem<T>` (F8.2b) — the concrete item
  over the SDK's typed hyperslab plane, entity at dimension 0 and layer at
  dimension 1, so it stays readable through the polyhedral-surface base by a
  partner that knows nothing about layers. It borrows the vertical coordinate
  and hands it out live rather than snapshotting at construction: an item
  holding a copy of the elevations would answer confidently and wrongly for the
  rest of a run, and nothing about the values would look wrong. Five gates,
  three ablations. FVQual's `SigmaVerticalCoordinate` drops into it directly —
  verified by compiling the two together.
- **The specification had a hole, and implementing it twice found it.** At
  exactly the invert both implementations returned the bottom width for
  `topWidth` and `wettedPerimeter` while the closed form said zero. The limit
  from above is right — zero would hand a partner an infinity in `sqrt(gA/T)` at
  the wet front — and `ICrossSection` now says so rather than leaving each
  producer to decide.

### Added — F3.7 spatially varying forcing and F5.4 turbulence closure

- `Forcing::ForcingField` holds either one value or a vector, so meteorological
  forcing is uniform or per column chosen per variable, and the uniform case
  costs one number. `MeteorologyField::at(column)` materialises a `Meteorology`,
  so `surfaceHeatFlux` stays the pure function every F3.5 gate tests and an
  `operator=` from `Meteorology` keeps those gates textually unchanged.
- **Wind stress is derived in the engine** from speed and direction rather than
  supplied beside them, and reaches both modes. The DeGray driver had been
  setting `windSpeed` for the heat budget while `windStressX/Y` stayed at zero:
  wind evaporated water and moved none, and because `Ri = N²/max(shear, 1e-12)`
  has no shear without it, the whole mixing scheme sat at its 1e-7 background
  for a simulated year. Two consumers of one variable, one wired and one not.
- Wind direction convention taken from the reference rather than guessed. W2
  projects stress as `cos(PHI − PHI0)` (`w2_4_win.f90:543`) and `W2Bathymetry`
  lays the chain out along `(sin PHI0, cos PHI0)`, so a compass angle maps to
  `π/2 − angle` and FVQual's own projection reproduces W2's expression without
  this code deciding whether PHI is a from- or a towards-bearing.
- `Hydro::TurbulenceClosure`: k-epsilon with Galperin stability functions,
  behind `useTurbulenceClosure`, off by default. Ten gates against theory.
- `Testing::basinGrid` — an nx-by-ny basin. Every prior gate used a one-cell-wide
  channel strip, which is right for the 1-D benchmarks and wrong for anything
  that circulates.

### Fixed — the DeGray wind failure, and three defects found looking for it

- **The vertical Courant limit read a scaled velocity.** `m_verticalFlux` stores
  `volumeFlux/planArea`, a true velocity only where interface area equals plan
  area; under hypsography DeGray's near-bed interfaces are 0.8% of it, so the
  limit under-restricted by up to 125×, never bound, and near-bed momentum
  diverged to 1e70 m/s in six steps. This was the failure. Fixed at the point of
  use so the mass-consistency identity stayed out of the change.
- **A data race in the closure's tridiagonal workspace** — the second time this
  exact bug has been written here. F7.1b found it in the momentum and scalar
  solves; it survives unit testing both times because the gates are
  single-threaded.
- **Wind stress on the body-force path.** `τ/(ρH)` multiplied back by a depth
  that evolves through the subcycle became `τh/(ρH₀)`: a column gaining depth
  gained stress. Now on the depth-independent surface-stress path.
- **The barotropic Courant limit used the plan area at full pool** rather than
  the storage area at the current stage.

### Changed — what the measurements now say

- Wind stress fixes the epilimnion: the 5–10 m band goes from −3.91 to −1.11 °C,
  and without it the modelled April surface layer is 3 m against a measured 14 m.
- The residual interior warm bias is **advective, not mixing**: removing mixing
  entirely makes it worse (5.87 °C against 3.42). A one-segment-wide basin
  forces every wind-driven return path through the vertical, and the same basin
  eight cells wide transports 39% less heat below 20 m.
- k-epsilon is a modest improvement on Richardson damping, 3.42 against 3.59 °C.
  **An earlier entry claimed the opposite** — 8.95 against 8.35 — from runs that
  predated the Courant fix and were pinned at the time-step floor. A comparison
  between two schemes is only worth its numbers if the model underneath steps
  correctly.

### Added — Phase 0 (scaffolding)

- CMake 3.21+ build with vcpkg manifest, presets for Linux/macOS/Windows, and a
  dependency-light `no-deps` preset.
- MIT license, contribution guide, citation metadata, and
  `THIRD_PARTY_LICENSES.md` recording CE-QUAL-W2 provenance policy.
- GitHub Actions build-and-test workflow.

### Added — Phase 1 (mesh & data model)

- `HorizontalMesh`: tri-quad structure-of-arrays mesh with unique-edge
  extraction, face areas and centroids, edge normals/lengths, face-neighbour
  topology, and structural validation (F1.1).
- `LayeredMesh`: column x layer indexing (`cell(c,k) = c * layerCount + k`),
  terrain-following sigma interface elevations built through the SDK
  `SigmaGridGenerator`, per-cell thicknesses and volumes, and a Haney-number
  hydrostatic-consistency gate (F1.1).
- `SigmaSpec`: sigma layering configuration (uniform / Song-Haidvogel
  stretched / explicit), with validation (F1.2).
- `UgridReader`: reads UGRID-1.0 2-D tri-quad meshes with fill-value-padded
  mixed connectivity from NetCDF-4 into an SDK `MeshDefinition` (F1.2).
- `LayeredUgridWriter`: CF-1.11 / UGRID-1.0 writer with a vertical layer
  dimension and `ocean_sigma_coordinate` formula terms, implementing the SDK
  `IOutputWriter` contract so it rides the SDK `IOThread` (F1.4). Node
  elevations are carried through, so output files are complete model inputs
  and FVQual can re-read its own results.
- `StateStore`: allocation-stable SoA state arrays with a species registry;
  temperature and water age are reserved transported species (F1.3).
- `fvqual` command-line runner (mesh inspection and sigma-grid generation).

### Added — Phase 2 (external barotropic mode)

- `BarotropicSolver`: explicit cell-centred finite-volume shallow-water solver
  on the tri-quad mesh — the fast mode of the barotropic–baroclinic splitting
  (F2.1). HLLC face fluxes in the face-normal frame, Audusse hydrostatic
  reconstruction for exact well-balancedness, semi-implicit Manning friction,
  wind stress, positivity-preserving wet/dry, and reflective wall boundaries.
  No linear solver and no per-step allocation.
- Adaptive CFL time stepping with `advance()` substepping to land exactly on a
  target time — the contract the internal mode needs when it subcycles this
  solver (F2.2).
- Verification gates: lake-at-rest over irregular and emergent bathymetry
  (exact to 1e-13), mass closure below 1e-10 in a sloshing basin, seiche
  period against linear theory, wind setup against the analytic slope with an
  unforced control, and a Ritter dam-break refinement study confirming
  first-order convergence.

### Added — Phase 3 (3D hydrostatic internal mode)

- `EquationOfState`: freshwater density with the 4 °C maximum, UNESCO salinity
  terms, and suspended-solids excess density (F3.2).
- `IPressureClosure` seam with `HydrostaticPressureClosure` — the hook
  through which a non-hydrostatic correction can arrive without touching
  momentum assembly. Baroclinic pressure uses Green-Gauss with every face
  evaluated at the cell's own elevation, which is exactly zero for a
  horizontally uniform ocean.
- `InternalModeSolver`: mode-splitting driver with time-averaged flux
  accumulation, mode adjustment, flux adjustment, continuity-diagnosed
  vertical velocity, implicit vertical viscosity and diffusion (batched
  Thomas), conservative scalar transport, Richardson-damped vertical mixing
  with convective adjustment, and prescribed surface heat flux (F3.1, F3.3).
- `BarotropicSolver::setBodyAcceleration` — how the internal mode hands the
  external mode its share of the baroclinic forcing.

- Verification against theory: quiescent stratified basin over steep
  bathymetry stays bounded at 1.1e-2 m/s; two-layer internal seiche period
  within 19% of `T = 2L/sqrt(g' h1 h2/H)`; lock-exchange front within 19% of
  `0.5 sqrt(g' H)`.

### Fixed

- **The baroclinic pressure gradient integrated the full density instead of
  the anomaly `rho - rho0`.** Because
  `integral(z..eta) rho = rho0 (eta - z) + integral(z..eta) rho'`, the
  supposedly baroclinic term carried the entire barotropic surface-slope term,
  which the external mode then received as a duplicate of its own restoring
  force, lagged across a whole subcycling interval. A free surface that should
  have tilted by millimetres reached several metres. Invisible to every
  rest-state test, since a flat surface hides the term entirely.
- The external mode never received the depth-averaged baroclinic forcing, so
  the two modes were decoupled and the mode adjustment erased the barotropic
  response that should have developed.
- The internal time step estimated wave speed from each column's *vertical*
  density range, which is zero for a horizontal density front — a lock
  exchange received the maximum step and was unstable from the first step.
- The internal time step had no vertical Courant limit, though sigma layers
  are metres thick where cells are tens of metres wide.

### Added — Phase 5, partial (transport engines)

- `EulerianTransport`: conservative scalar advection with MUSCL reconstruction
  and Zalesak flux correction (F5.1). Inverse-distance least-squares cell
  gradients with Barth-Jespersen limiting horizontally, minmod slopes
  vertically, and a Zalesak limiter that admits only as much of the
  antidiffusive correction as the local bounds allow — so boundedness is a
  property of the algorithm rather than of a tolerance.
- `InternalModeOptions::advectionScheme` selects MUSCL+FCT (default) or
  first-order upwind, the latter retained as the baseline the accuracy tests
  measure against.
- Measured effect: the two-layer internal seiche period bias fell from +18.6%
  to **+7.9%**, and the L1 error of a translated square pulse more than halved.
  The lock-exchange front speed is unchanged, which is expected — a gravity
  current's speed is set by integral dynamics, and at the front itself the
  limiter correctly degrades to upwind. That front converges under mesh
  refinement instead: -11.1%, -6.2%, -4.3% at 100, 200 and 400 cells.

### Added — verification animations

- `fvqual-verify` runs the six benchmarks that have closed-form or theoretical
  solutions and emits animation frames; `verification/index.html` is a
  self-contained page that animates the model against each reference with live
  error metrics.

### Changed

- **ELM and RWPT transport engines deferred indefinitely.** ELM cannot
  conserve mass — interpolation at backtracked departure points leaks, and the
  remedies are approximate corrections to a structural defect — and closed
  constituent budgets are what a water-quality model is for. RWPT conserves
  mass exactly but delivers it as statistical noise that nonlinear kinetics
  amplify. Water age and source attribution will be transported species rather
  than particle diagnostics, which keeps every diagnostic on the same
  conservative footing as the constituents. The `TransportEngine` interface is
  not built: with one engine it would be an abstraction guessing at its own
  seam.

### Added — F3.4 / Phase 4 (inflows and withdrawals)

- `Forcing::SourceSet`: inflows and withdrawals as volumetric source terms.
  Vertical placement is the substantive part — **density-matched plunging
  inflow** (a tributary finds its own level, giving overflow, interflow, or
  underflow), **selective withdrawal** through an envelope whose half-height
  follows from the discharge and the column's buoyancy frequency, plus
  uniform, surface, bed, and fixed-elevation placements. `setDischarge`
  applies an operating schedule mid-run and refuses to reverse a source's
  physics.
- Wired through both modes: a per-column source on the barotropic continuity,
  with inflow arriving momentum-free (and so diluting) while outflow leaves at
  the local velocity; per-cell sources in the continuity integration; and the
  scalar source inside the transport's conservation update.
- Nine verification gates, including that the free-surface continuity residual
  stays below 1e-10 with sources active — the mass-consistency contract now
  covers them.

### Fixed

- **The scalar source was applied after transport instead of inside it.**
  Continuity had already carried the source into the vertical flux field, so
  the divergence changed the concentration regardless: a withdrawal removed
  water but not its constituent load, concentrating what remained, and the
  exported mass came out as zero.
- **The withdrawal envelope used the density gradient local to the outlet.** A
  deep port below a sharp thermocline sits in uniform cold water where that
  gradient vanishes, so the envelope was unbounded and the outlet appeared to
  draw from the surface. Now uses the column-scale buoyancy frequency.

### Added — PY0 (C API foundation)

- Public C API under `include/fvqual/capi/` — `fvqual_engine.h`,
  `fvqual_mesh.h`, `fvqual_sigma.h`, `fvqual_options.h`, `fvqual_sim.h` —
  exported from `fvqual_core`. Opaque incomplete-struct handles, integer error
  codes with thread-local detail messages, and a two-phase
  create/register/initialize lifecycle. This is the layer the Python bindings
  will bind to, and the ABI firewall that lets the C++ engine be
  refactored without breaking a compiled wheel.
- Simulations copy the mesh they are created from, so mesh handles carry no
  lifetime obligation.
- `tests/capi/test_capi.c`: 50 checks written in C and compiled by a C
  compiler, including a barotropic-seiche gate (volume drift 1.2e-15,
  continuity residual 5.4e-16), plus a strict `-std=c89 -pedantic-errors`
  header-compilation check so the ABI promise is enforced rather than stated.

### Added — Phase 4 (structures & operations)

- `Forcing::SourceSet`: inflows and withdrawals as volumetric source terms with
  density-matched (plunging) inflow, selective-withdrawal envelopes, and
  uniform/surface/bed/fixed-elevation placements, wired mass-consistently
  through both the barotropic and internal modes (F3.4, F4.1).
- `Forcing::StructureSet`: weirs with Villemonte submergence, orifices and
  gates with a linear opening control, pipes, and pumps with a start/stop
  deadband. Each structure drives a paired source, so a rating curve composes
  with selective withdrawal (F4.1).
- `BarotropicSolver::setColumnSourceUpdater`: a per-substep hook that re-rates
  level-driven sources against the free surface as the subcycle moves it. Rating a spillway once per internal step freezes its discharge at a
  head that may be minutes stale, which drains a reservoir below its own crest;
  re-rating each substep holds the level at the crest by physics rather than by
  a limiter. Transport receives the time-mean release over the step so the two
  modes move the same volume.
- `Forcing::OutletGroupSet`: multi-port outlet works operated to a target
  release temperature (F4.2). Ports are blended by the heat-conservation mixing
  law `f_w = (T* - T_c) / (T_w - T_c)` using the tightest bracketing pair;
  targets outside what the reservoir holds open the closest port and report the
  temperature actually achieved, and demand beyond the ports' combined capacity
  is reported as a shortfall rather than met by inventing water.
- `Structure::scheduled`: under a group schedule a structure's rating curve
  becomes its capacity limit rather than its discharge, which is what lets an
  outlet be operated to a target instead of merely observed.
- `SourceSet::withdrawalMean`: the envelope-weighted value a withdrawal would
  draw at a stated discharge. The discharge is an argument because the question
  has to be asked of ports that are currently shut, and a closed outlet has no
  envelope of its own.
- `Forcing::TransferSet`: internal transfers moving water between columns —
  pumped-storage returns, bypasses, destratification pumps (F4.3). The
  receiving inflow carries the composition its donor is actually drawing rather
  than a stated one, taken through the same cached vertical weights the
  withdrawal removes mass at, which makes constituent conservation an algebraic
  identity rather than an approximation.
- `SourceSet::weightedMean` (cached weights) alongside `withdrawalMean`
  (probed at a stated discharge); the two answer different questions and must
  not be conflated.
- Reservoir drawdown is verified against the closed-form weir law
  `h(t) = h0 / (1 + C L sqrt(h0) t / 2A)^2`, matching within 2-6% once the
  basin has had two gravity-wave transits to level.

### Changed

- **Selective withdrawal now follows CE-QUAL-W2's actual formulation**.
  The previous zone limits `(Q/N)^(1/3)` and `(Q/(W N))^(1/2)` dropped the
  published leading coefficients and evaluated the buoyancy frequency once at
  column scale; the reference evaluates it locally and implicitly, between the
  outlet layer and each candidate layer over their own separation, marching
  outward until the separation catches the computed thickness. The velocity
  profile is now parabolic in density anomaly relative to the outlet layer
  rather than in distance from it. Consequences: a port in a well-mixed
  hypolimnion draws uniformly across it and nothing across a sharp thermocline.
  No previous gate caught this, because all of them asserted the sign of
  selectivity rather than its magnitude.
- `Source::withdrawalAngle`: the effective horizontal angle of withdrawal,
  default pi for a port on a wide vertical dam face.
- The withdrawal velocity maximum now shifts off the port in an asymmetric zone
  (Bohan and Grace 1973, eq. 4: `y1 = H sin^2(1.57 z1/H)`), and the parabolic
  profile is normalised about that layer rather than the port layer. A
  symmetric zone returns the port exactly, so this generalises the previous
  behaviour.

### Added — Phase 2 (F2.3, open boundaries)

- `BarotropicSolver` open boundary conditions: prescribed-discharge `Inflow`,
  prescribed-elevation `Stage`, and zero-gradient `Transmissive`, assigned per
  boundary edge and settable mid-run. How many quantities may be imposed is
  counted by the characteristics rather than chosen: in subcritical flow
  exactly one leaves the domain, so one condition is imposed and the rest comes
  from the outgoing Riemann invariant `u + 2c`. Supercritical outflow takes
  nothing, supercritical inflow takes both.
- Subcritical inflow solves `-q/h + 2 sqrt(g h) = R` for the boundary depth by
  bisection; the left side is monotone in `h`, so the root is unique and the
  search cannot fail where Newton could.
- `boundaryFlux` reports volume flux per boundary, signed into the domain.

### Added — Phase 2 (F2.4, MacDonald benchmarks)

- MacDonald-type manufactured steady-flow benchmarks: a smooth depth profile
  with the bed integrated from the steady momentum equation so the profile is
  an exact solution. Verifies the bed-slope/friction/flux balance that Ritter
  cannot reach. The manufactured friction slope uses the wide-channel form
  `n^2 q^2 / h^(10/3)` to match the solver exactly, rather than the
  hydraulic-radius form, which would pose the problem for a different model.
- Measured: 0.57 mm mean depth error on a 2 m depth at 200 cells, reached both
  from the exact profile and from a flat start; refinement ratios 1.97 and 1.98
  across 50/100/200 cells, i.e. clean first-order convergence.

### Added — Phase 7 (F7.2, profiling baseline)

- `Hydro::StepProfile`: per-phase wall-clock accounting inside the internal
  step, with `profile()` / `resetProfile()`.
- `benchmarks/fvqual-bench`: a stratified-basin sweep reporting per-phase share,
  throughput, external substep count, and cost per column-substep.
- Measured baseline: the barotropic subcycle is 58-71% of runtime and scalar
  transport 20-28%; everything else is under 7%. ~108 external substeps per
  internal step at every mesh size, at a flat 73-74 ns per column-substep.

### Changed — Phase 7 (F7.1a)

- The barotropic edge sweep no longer scatters into adjacent cells. It is split
  into an edge-parallel flux pass and a cell-parallel gather over the mesh's
  existing `faceEdges` CSR adjacency, which is race-free by construction and is
  the form Kokkos and a GPU want. Costs ~7% serially for the extra edge traffic.
- OpenMP applied across the external mode, scalar transport, the pressure
  closure, and the internal-mode kernels, behind `FVQual_USE_OPENMP`. All gates
  pass threaded and unthreaded.
- Vertical-solve workspace is now per-thread (`layers x maxThreads()`); it was
  six shared member vectors, which is a silent-corruption race once the column
  loops are threaded. `maxRichardson` became a reduction for the same reason.
- Measured 1.87x overall at 102,400 cells on four threads (75.4 -> 40.2
  ms/step, 1.36 -> 2.55 Mcell-steps/s); pressure 3.47x, transport 2.70x,
  external 1.68x.

### Changed

- The barotropic gather assigns the derivative arrays outright rather than
  accumulating into them, so the per-substep `std::fill` of `m_dh`, `m_dhu`,
  and `m_dhv` was dead work left over from the scattering version. Removed;
  only the boundary-flux tally still needs clearing.

### Added — Phase 6 (F6.1a, kinetics numerics)

- `FVQual::Kinetics::Expression`: Shunting-yard parser and stack evaluator for
  reaction-rate expressions, with name binding and a bound fast path.
- `FVQual::Kinetics::integrate`: adaptive RK45 Cash-Karp integrator, the
  substrate for operator-split reaction kinetics.
- Both relicensed to MIT from the author's openswmm.engine (see
  THIRD_PARTY_LICENSES.md), and verified in this build rather than assumed from
  the origin: expressions against the same formulas evaluated by the compiler,
  the integrator against closed-form solutions for exponential decay, coupled
  sequential decay, and a stiff pair.
- **Parser hardening:** malformed expressions (`""`, `"2 +"`, `"* 3"`) parsed as
  valid, which for a user-authored rate law means a model that runs and is
  quietly wrong. `parse()` now rejects any postfix sequence that does not
  evaluate to exactly one value.

### Added — Phase 6 (F6.1b, reaction network)

- `FVQual::Kinetics::ReactionNetwork`: a rate law per species, written as an
  expression over the other species and named constants, resolved once at
  `compile()` and integrated per cell. The network is data, not code, so adding
  a constituent does not mean recompiling the engine.
- Temperature is not a special case: it is a transported species addressable by
  its registry id, so Arrhenius dependence is something the modeller writes.
- Wired into `InternalModeSolver::step()` as a Lie-split operator after
  transport and before the density update, with a `StepProfile::reactions`
  timer and `compileReactions()` binding the network to the solver's own
  registry — passing a different registry would resolve to the wrong species.
- A network with rate laws that has not been compiled refuses the step up
  front, rather than failing partway through and leaving the state transported
  but unreacted.
- `ReactionDiagnostics` reports the smallest concentration reached, the cells
  that went negative, and the cells the integrator could not converge.
  Negatives are **not** clamped: a clamp restores plausibility while destroying
  the mass balance, and hides the rate law that caused it.
- Reaction conservation is deliberately not enforced. Transport conservation is
  an engine identity; whether reactions conserve depends on whether the rate
  laws sum to zero across species, which is the modeller's arithmetic. The
  gates verify that a closed network conserves, that an open one is *allowed*
  to lose mass, and that a zero-order sink crosses zero and says so.
- Gates: first-order decay against its closed form both in isolation and
  through the full solver; `A -> B -> C` checked against its analytical
  solution and its invariant every step; and a linear decay in a stirred basin
  against the bulk law `M(t) = M₀ exp(-kt)`, which holds under any flow field
  and so tests that transport and reaction compose. Each was proved to have
  teeth by ablation — removing the operator from the step, removing the
  pre-step guard, removing the unresolved-name check, and sharing one
  workspace slice across threads each break the gates that claim them.
- Reaction workspace is re-sized if the host raises the thread count after
  `compile()`; the previous sizing would have had threads writing past the end.

### Added — Phase 6 (F6.2a, light field and derived variables)

- `FVQual::Kinetics::LightField`: Beer-Lambert irradiance down each column, with
  clear-water attenuation, per-species self-shading, and surface absorption
  reported rather than discarded. Because attenuation is piecewise constant in a
  layer and the exponential composes exactly, the interface irradiances equal
  `I0 exp(-Kd z)` to round-off at any layer count — so the gates are written at
  1e-12 against exact values rather than against a discretization tolerance.
- **Cell values are layer means, not top-of-layer values.** Over a layer of
  optical depth `tau` the mean is smaller by `(1 - e^-tau)/tau`, which in a 2 m
  layer at Kd = 0.6 is a factor of 1.4. Using the top value would bias
  production high everywhere, worst in the surface layer where most of it
  happens. `meanAttenuationFactor` is written to stay accurate as `tau` vanishes,
  which is an ordinary case: a sigma grid's top layer thins whenever the surface
  drops.
- A third category of quantity in a rate law, after species and constants:
  engine-computed `DerivedFields` addressed by reserved name — `__LIGHT__`,
  `__DEPTH__`, `__THICKNESS__`. Referencing one that was not supplied is
  refused; ablation showed that check is a bounds guard rather than the
  diagnostic it was written as, since the integration loop indexes the span
  directly.
- Wired into the step: light is marched from the just-transported field
  immediately before the reactions, so self-shading reflects the constituents
  the reactions are about to act on, and only when a rate law reads it.
- Gates include light-limited growth through the full solver where every layer
  grows at its own rate against its own closed form — which a single
  column-averaged irradiance would fail while still passing any "growth
  happened" check.
- Two of my own test errors, recorded because both were wrong expectations
  rather than defects: a linearised expectation held to a tolerance a thousand
  times below its own truncation error, and a `__DEPTH__ * __THICKNESS__` gate
  that multiplication's commutativity would have let pass with the two fields
  swapped.

### Added — Phase 6 (F6.2b, reaction scope)

- `RateScope { WaterColumn, Bed }`. A benthic rate law applies only in the cell
  touching the bed, and is written per unit **area** — g/m²/s, which is how a
  sediment oxygen demand is actually quoted, because it is a property of the
  sediment surface and not of the water above it. The engine divides by the bed
  cell's thickness. Asking the modeller to write `/ __THICKNESS__` instead would
  make the benthic demand silently proportional to the layer count, and it looks
  plausible either way.
- The duplicate-law key is now (species, scope) rather than species, which
  revisits the F6.1b refusal deliberately: a water-column law and a bed law are
  distinct processes, not an ambiguity. Two laws in the *same* scope are still
  refused.
- Where the two overlap — a bed cell is part of the water column too — the rates
  **add into one derivative** rather than being applied in sequence. Sequential
  application would introduce a splitting error between two terms of the same
  ODE, for nothing. The integrated state vector is therefore indexed by species
  rather than by rate law.
- The solver supplies the bed-cell mask, so the network never needs to know the
  mesh's indexing convention.
- **A duplicate mechanism found and removed.** The areal conversion factor was
  zeroed away from the bed, which enforced scope a second time; ablating the
  scope check itself then changed no result, so nothing was testing it. The
  factor is now computed unconditionally and the check is the only mechanism —
  ablating it fails four gates.
- Gates: a bed law touches only the benthic cell; an areal flux removes the same
  mass per unit area at 4, 8 and 16 layers (the layer-count invariance a missing
  thickness division breaks); water-column and bed laws sum against the closed
  form of `dC/dt = -kC - F/h`; and the solver puts the flux at the *deepest*
  cell, checked against depth rather than a layer index so a mask marking the
  surface fails.

### Added — Phase 6 (F6.2c, settling and the benthic handoff)

- `FVQual::Transport::SettlingSet`: per-species downward settling, applied after
  transport and before the reactions so a benthic rate law acts on material that
  arrived this step.
- **Solved implicitly, by a single forward substitution.** Downward-only
  upwinding makes the matrix lower triangular, so the implicit solve costs the
  same as an explicit update and buys three things: no CFL constraint (the
  binding case is a thin surface sigma layer, not a fast-settling species —
  0.1 m and 10 m/d caps an explicit step near 900 s), unconditional positivity,
  and exact conservation, because what leaves the column telescopes to
  `w dt C_bed`, which is what the bed receives. The cost is extra numerical
  diffusion, which is recorded rather than glossed.
- Settled material accumulates in a benthic species measured per unit **area**,
  matching the `RateScope::Bed` convention: a benthic store expressed as a
  concentration in the bottom cell would change value whenever the layering did.
  A transported receiver is refused, since the sediment would be advected away.
- With no receiver configured, the mass is tallied and reported. A model may
  legitimately not track sediment; what it may not do is lose mass silently,
  because an unreported loss is indistinguishable from a conservation bug.
- Gates: water-plus-sediment conserved to 1e-12 checked every step; positive and
  conservative at **100x** the explicit stability limit (the gate that earns the
  implicit solve — an explicit ablation fails it); the single-layer closed form
  of the scheme; L1 convergence at first order against the exact translated
  profile; and direction, pinned separately because a sign error would give a
  perfectly conservative, perfectly stable, upside-down model.

### Added — Phase 6 (F6.3a, surface scope)

- `RateScope::Surface`, the mirror of `Bed` at the other end of the column, for
  air-water exchange. It is areal without needing to be declared so: reaeration
  is `K_L (C_sat - C)`, and m/s times g/m³ is already g/m²/s, so it reuses the
  same areal-to-volumetric division.
- The bed-cell mask generalises into a `CellPosition` bit field, which lets a
  **one-layer column be both surface and bed**. That is not a degenerate case to
  exclude — it is a shallow pond, where reaeration and sediment demand really do
  act on the same water — and a gate pins it.
- Gates: a surface law and a bed law on one species act at opposite ends of the
  same column and sum in a one-layer column; surface fluxes are layer-count
  invariant like benthic ones; and the solver places exchange at the
  *shallowest* cell, checked against depth rather than a layer index. Ablations:
  treating Surface as WaterColumn fails three gates, and a solver mask that
  marks the bed as the surface fails the wiring gate.

### Added — Phase 6 (F6.4 tier 1, constituent templates)

- `FVQual::Kinetics` templates: ready-made constituent models built on the
  reaction network. `WaterAgeTemplate` and `FirstOrderDecayTemplate` first —
  water age, decaying tracers, and conservative constituents.
- **Builders, not a second engine.** Every template installs readable rate-law
  strings through the same public API a user would call, so the result is
  inspectable and editable afterwards and no template can express anything a
  hand-written network could not.
- **Parameters are per day, not per second.** Every reference a modeller
  consults quotes rates per day, and SI in the struct invites a silent factor of
  86400 — which produces output that looks entirely plausible. The conversion
  happens inside the template and is gated: ablating it fails two tests.
- Water age's surface condition is a **relaxation, not a reset**, and the
  documentation says so rather than claiming zero. A hard assignment inside a
  split operator would make the surface contribution depend on layer thickness
  and hand the integrator a discontinuity every step. The relaxation carries an
  explicit `__THICKNESS__` factor to cancel the areal division a Surface scope
  applies — without it the reset strength would depend on the sigma grid, which
  is what the layer-count gate checks.
- A zero decay rate installs **no rate law at all**, so a conservative
  constituent is carried by transport and never reaches the integrator.
- `InternalModeSolver::speciesIndex()` and `species()` expose the registry
  read-only, so templates can resolve names without being able to bind a network
  against the wrong registry.

### Added — Phase 6 (F6.4 tier 2, the dissolved-oxygen suite)

- `DissolvedOxygenTemplate`: oxygen, carbonaceous demand, ammonium, nitrate,
  phosphate and algae, in the form QUAL2K and CE-QUAL-W2 share — reaeration at
  the surface, CBOD decay and nitrification in the water column, sediment demand
  at the bed, and algae that produce oxygen growing and consume it respiring.
  Every part of Phase 6 is load-bearing here: surface and bed scopes, the light
  field with self-shading, and settling with a benthic handoff.
- **The stoichiometry is the point, and the budgets are how it is verified.**
  Growth, respiration, mortality and nitrification all move nitrogen between
  pools, and each must carry the same `g N per g algae` or the model quietly
  manufactures nitrogen. The gates check the budgets rather than the
  coefficients, because a budget that closes is evidence and a coefficient that
  looks right is not.
- **Nutrient uptake is split without a regulariser.** `NH4 / (NH4 + NO3)` is
  0/0 in water containing neither, so the split reuses the half-saturation
  denominator limitation already carries — `NH4 / (Ksn + NH4 + NO3)` — and the
  two uptake terms then sum to exactly the nitrogen growth requires. Ablating
  this to an epsilon-regularised form breaks the nitrogen budget, which is why
  it is written this way.
- Gates: Streeter-Phelps against its analytical deficit in a well-mixed column;
  `DO - CBOD` conserved without exchange, which is the statement that the
  ultimate deficit equals the initial demand; nitrification costing 4.57 g O₂
  per g N; nitrogen and phosphorus budgets closed every step over active growth;
  growth equal to the product of its limitation factors each computed
  independently; and no growth in the dark, which catches a light factor
  written upside down.
- Deliberately absent, and documented as such rather than defaulted quietly:
  temperature- and elevation-dependent oxygen saturation, separate organic
  nitrogen and phosphorus pools, and Steele photoinhibition.
- Reaeration is a **transfer velocity** (m/day), not a first-order rate (1/day).
  A rate builds in an assumed depth, and the model would then reaerate
  differently for no physical reason when the reservoir is drawn down.

### Added — Phase 6 (F6.3b, text-configured reaction networks)

- `FVQual::Kinetics::ReactionFile`: a `.rxn` text format with `[SPECIES]`,
  `[CONSTANTS]`, `[RATES]`, `[SETTLING]` and `[LIGHT]` sections, so changing a
  constituent model no longer means recompiling. This is the other half of
  making the network data rather than code.
- **Numeric fields are arithmetic expressions**, evaluated by the same parser
  the rate laws use: `0.23/86400` rather than `2.6620370370e-6`. It is the one
  place the file can undo the per-day trap without guessing which constants in
  an arbitrary rate law are rates.
- Failures report a line number and what was wrong, and parsing commits to the
  object only on success — a rejected file leaves nothing half-applied.
- Two format flaws that gates caught before the format shipped, both recorded in
  the header rather than quietly patched:
  - Numeric fields may not contain spaces, because a settling record has an
    optional field *after* the number and whitespace is the only delimiter.
    `TSS 1.0 / 86400 TSS_BED` parsed the receiver as part of the velocity.
  - The scope on a rate law is **mandatory**. An optional field in front of a
    free-form expression is unparseable: `DO   0 - k * CBOD` is
    indistinguishable from a law scoped to something named `0`. Requiring it
    also makes a model file state where each process acts.
- The central gate is equivalence rather than parsing: a file expressing
  first-order decay and the C++ template expressing the same model must produce
  bit-identical concentrations, and both must equal the closed form.

### Added — Phase 2 (F2.3b, 3-D transport across open boundaries)

- The internal mode no longer refuses open boundaries. `bindOpenBoundaries()`
  turns each (boundary, adjacent column) pair into a `Forcing::Source` whose
  discharge is refreshed each step from the volume the external mode actually
  passed — so boundary water rides the placement, composition and conservation
  machinery a tributary inflow already uses. This unblocks tributary inflows and
  so Phase 8.
- **A source, not an edge flux.** The external mode solves a depth-averaged
  problem, so all it knows about a boundary is the depth-integrated volume that
  crossed it; how that spreads over the column is a modelling choice either way.
  A column source makes the choice explicit and selectable rather than burying
  it in a flux reconstruction, and keeps mode consistency for free because the
  vertical velocity is already diagnosed from continuity *including* sources.
- **Inflow is density-placed by default**, so a cold tributary plunges to its own
  density rather than mixing across the surface — the behaviour that
  distinguishes a stratified reservoir model from a mixed one, and what W2 does.
  `Placement::Uniform` remains selectable per boundary.
- **Outflow always leaves uniformly.** A density-matched search has no meaning
  for water leaving, and a boundary has no outlet elevation or width, so the
  selective-withdrawal treatment a real port gets is not available to it. Model
  such an opening as a structure instead.
- **An unset boundary carries ambient values**, so it is transparent to
  constituents rather than a source of distilled water — with a gate proving a
  uniform field crosses unchanged. `setBoundaryValue()` prescribes a composition
  per (boundary, species).
- `Source::externallyAccounted` marks water the external mode has already moved,
  keeping it out of the column total handed back to the depth-averaged solver.
- Three defects the gates caught, all of which would have run and looked
  plausible: boundary water counted twice (once as flux, once as source); the
  density placement applied to outflow, leaving the mode-consistency residual at
  1.2e-2 against 1e-15; and the ambient composition read before the placement
  weights were current, evaluating to zero and importing distilled water.

### Added — Phase 3 (F3.5, meteorological surface heat budget)

- `FVQual::Hydro::surfaceHeatFlux()`: net shortwave with albedo, atmospheric
  longwave with cloud correction, back radiation, latent and sensible heat, from
  a `Meteorology` record settable each step.
- **The flux now depends on the water temperature it is setting**, which is the
  point. A prescribed flux has no feedback and will drive a column past boiling;
  a computed one is self-limiting, because back radiation and evaporation both
  rise with surface temperature. `equilibriumTemperature()` exposes where the
  terms cancel, solved by bisection since the balance is monotone.
- **Penetrating shortwave is a distributed source, not a surface flux.**
  Longwave, latent and sensible exchange happen in the top millimetres; sunlight
  does not, and the depth over which it is absorbed is what builds a thermocline.
  The penetrating part is marched with the same Beer-Lambert code the kinetics
  use, and the energy absorbed in a cell is the irradiance difference across it.
  Sunlight reaching the bed is absorbed there rather than lost.
- The solar field is kept **separate from the kinetics light field**: growth
  responds to the photosynthetically active band and is shaded by algae, while
  heating follows the whole spectrum through clear water. One field would force
  a single attenuation to stand for both.
- Inert until `enableMeteorology()`, so existing models keep their prescribed
  `surfaceHeatFlux` unchanged — with a gate for that.
- Gates are behavioural rather than arithmetic: the budget is self-limiting over
  -5 to 40 °C; a shallow column reaches its equilibrium temperature from *both*
  directions; energy entering matches the flux to 1%; sunlight warms the
  interior with a stably ordered profile. Ablations — all shortwave to the
  surface, frozen back radiation, flipped latent sign — each break the gates
  that claim them.

### Added — Phase 1 (F1.5, depth-varying plan area)

- `LayeredMesh::setHypsography()`: a per-column plan area that varies with
  elevation, so cell volumes and the volume-elevation curve reflect a real
  water body rather than a prism.
- **Why it matters, measured on DeGray's actual bathymetry:** the deepest
  section runs from 6418 m wide at the surface to 75 m at the bed, a ratio of
  85, and a prism on surface area overstates the reservoir's volume by
  **2.07x** — 3.11e9 m³ against 1.51e9. Almost all of that error is in the
  hypolimnion, which is the water whose volume governs stratification strength
  and oxygen depletion, so a validation built on prisms would measure the mesh
  rather than the engine.
- Cell volume is the **exact integral** of the piecewise-linear area curve, split
  at every sample inside the cell. A midpoint rule is exact only where the cell
  contains no kink, and a cell straddling a break in a surveyed curve is
  precisely where it is not — there is a gate for that, and it fails under a
  midpoint ablation.
- Curves are clamped outside their samples rather than extrapolated: a linear
  fit below the deepest sample runs the area negative.
- A curve that widens downward is refused, because it would make the storage
  curve non-monotone and a free surface would no longer follow uniquely from a
  volume.
- **Nothing changes without a curve.** A column with no hypsography produces
  bit-for-bit its previous volumes, gated with `EXPECT_DOUBLE_EQ`, and the whole
  existing suite is unmoved — the barotropic solver branches so the
  constant-area path is the original code, untouched.
- **The external mode evolves volume, not depth, on a hypsographic column.**
  Updating stage with the area at the start of the step leaves an error of order
  `(dEta)² dA/dz` against the exact integral the mesh computes: measured, the
  mode-consistency residual sat at **7.7e-7** rather than round-off. Advancing
  volume and inverting the storage curve — a closed-form quadratic root within
  each linear segment, not bisection — makes both modes compute the same
  integral, and the residual returns to 1e-15.
- Gates: mode consistency holds with the surface actively rising on a steeply
  tapered mesh; still water over sloping bathymetry stays still, so
  well-balancedness survives the change; and the stage after a known inflow
  matches the analytic inverse of the storage curve, with a check that a
  constant-area answer would have been visibly different.

### Added — Phase 8 (F8.1, CE-QUAL-W2 bathymetry ingest)

- `FVQual::IO::W2Bathymetry` reads a W2 bathymetry CSV — segment lengths,
  orientations, friction, initial water surface, and per-layer widths — and
  builds a FVQual mesh from it. Thirty years of reservoir schematisations live
  in this format, and reading one is what turns "passes its analytical
  benchmarks" into "reproduces a real reservoir".
- **The vertical datum is a required parameter, not a guess.** The layer rows
  carry a `K` column and the header names an `ELEV` column, but no example in
  the W2 distribution populates it: elevations come from the control file's
  location card as `EBOT`. For DeGray that is 66.35 m, which with 36 two-metre
  layers puts the initial surface at 123.80 m holding **7.29e8 m³ over
  5.09e7 m²** — matching the reservoir's published figures, which is how the
  interpretation was checked.
- **A segment becomes one column with a storage curve**, using F1.5. W2 is
  laterally averaged; FVQual is three-dimensional with a plan area that varies
  with elevation, so `width(z) x DLX` carries the storage the schematisation
  encodes without pretending the reservoir is a box.
- **The curve is sampled at layer midpoints**, and that is the load-bearing
  choice. W2's width is a step function of elevation; a hypsographic curve is
  piecewise linear; and the integral of a linear function over a layer equals
  its midpoint value times the thickness — so each layer's volume comes out at
  exactly `B x DLX x H`. Sampling at layer boundaries averages adjacent widths
  and loses storage wherever the reservoir narrows, which for DeGray is
  everywhere. Ablating midpoint to boundary sampling fails the volume gate.
- Verified against the real file rather than only a synthetic one: FVQual's mesh
  integrates DeGray's storage to **0.03%** of what W2's own cell arithmetic
  gives (7.29223e8 against 7.28973e8), on a sigma grid that deliberately does
  not align with W2's layers — so the agreement is a statement about the storage
  curve, not about indexing. The prism mesh the same geometry would have
  produced is **2.40x** too large.

### Added — Phase 3 (F3.6, computed solar radiation)

- `clearSkySolarRadiation()` and `solarPosition()`: solar declination, hour
  angle and elevation, extraterrestrial irradiance modulated by orbital
  eccentricity, an air-mass transmittance, and a cloud reduction.
- **Needed because CE-QUAL-W2 computes shortwave rather than reading it** when
  its `SROC` switch is off, which is how DeGray is configured — its meteorology
  file has `JDAY, TAIR, TDEW, WIND, PHI, CLOUD` and no shortwave column at all.
  A heat budget drivable only by a measured record cannot run that model.
- The evaporative wind function generalises from `a + bU` to `a + bU^c`, because
  calibrated wind functions are routinely quadratic — DeGray's is
  `9.2 + 0.46 W²` — and a linear form cannot follow one at both low and high
  wind. At 8 m/s the quadratic term is six times the constant.
- Gates use the exact reference values solar geometry has: declination at the
  solstices and equinoxes, twelve-hour days at the equator all year, polar day
  and night above the Arctic circle, symmetry about solar noon, and a clear
  midsummer daily total of 20-40 MJ/m² at DeGray's latitude.

### Added — Phase 8 (DeGray Reservoir, end to end)

- `fvqual-degray` runs CE-QUAL-W2's DeGray example through FVQual for March to
  December 1980: real bathymetry, meteorology with computed solar radiation,
  density-placed inflow, withdrawal at the dam, and the initial profile from
  `vpr.npt`. Results and a full account are in `verification/degray/`.
- The reservoir starts nearly isothermal at 7.4/6.5 °C, stratifies through
  spring, holds a cold hypolimnion all summer while the surface peaks at
  **30.9 °C** on day 197, deepens its thermocline from 5.2 m at midsummer to
  13.0 m in October, and turns over by December to an isothermal ~9.9 °C.
  Maximum stratification 23.8 °C. Every one of those is what DeGray does.
- The driver also writes the whole profile at the dam every day
  (`degray_profiles.csv`), which is what the depth-time contour is drawn from.
  Surface and bottom temperatures say a reservoir stratified; only the profile
  says whether it stratified the way a reservoir does.
- `verification/degray/plot_degray.py` draws five figures and
  `verify_figures.py` restates every claim they make as an inequality over the
  same CSV data, exiting non-zero if one is not true. Annotated numbers are
  computed from the data rather than typed, so they cannot drift; the checks
  exist for the words, which do not recompute themselves.
- **The stage-storage figure is an independent cross-check, not decoration.**
  The plot script re-derives DeGray's storage curve in Python straight from
  `bth1.csv`, with no reference to anything the C++ computed, and plots every
  day of the run against it. They agree to **0.22% at worst**, and `dV/dz`
  recovered from the run is 4.98e7 m² against the reservoir's 5.09e7 m² plan
  area at full pool — the hypsographic mesh and the barotropic mode agreeing
  through a curve neither of them was told.
- **It is not a code-to-code validation**, and the distinction matters: the W2
  distribution ships DeGray's inputs but no output, so there is nothing to
  compare against. It is judged on physics, which catches gross errors and would
  not catch a subtle formulation difference. Nor does it exercise
  three-dimensionality — a laterally-averaged schematisation is one segment
  wide, so the mesh is a chain of columns.

### Fixed — found by the DeGray run, not by the test suite

- **The evaporative wind function was roughly eight times too weak**: 10 W/m² at
  a 1000 Pa vapour deficit where the physical value is nearer 80. Evaporation is
  the dominant cooling term at high temperature, so equilibrium temperature ran
  tens of degrees high and peak surface temperature reached **124 °C**. The
  defaults are now CE-QUAL-W2's `9.2 + 0.46 W²` converted to SI.
  **247 unit tests passed over this**, because every heat-budget gate computed
  its expectation from the same function it was testing — the budget was still
  self-limiting, still conservative, still reached "its" equilibrium. A new gate
  compares the terms against magnitudes from outside the code.
- **Cloud cover units at the W2 boundary.** W2 records cloud in tenths (0–10);
  the budget wants a fraction. The cloud factor is quadratic, so a value of 10
  multiplied the atmospheric longwave by **eighteen**, and peak surface
  temperature was still 76 °C after the wind function was fixed. Converted at
  the boundary where the conventions meet, and `surfaceHeatFlux` now clamps —
  it had not, while `clearSkySolarRadiation` did, and two functions in the same
  file disagreeing about what a cloud fraction is was the underlying fault.
- **The daily report drifted five hours through the run.** The schedule was
  `nextWrite = day + 1.0`, rebasing on the time actually reached instead of
  advancing a fixed schedule, so each step's overshoot carried forward — up to
  five minutes a day under the 300 s cap. The report walked from 12:05 in March
  to 17:00 in December, sampling a steadily later hour of the diurnal cycle and
  aliasing it into the seasonal signal. The trajectory was never affected, the
  same 119,750 steps either way, but the reported peak surface temperature fell
  from 31.4 °C to 30.9 °C once the sample was consistent. Found by asserting
  that a figure caption reading "near solar noon" was true of the timestamps.

### Known limitations

- An open boundary cannot do selective withdrawal, and its inflow composition is
  uniform over the boundary rather than specified per layer.

- **Hypsography corrects storage but not conveyance.** With plan area varying
  inside a cell the depth-averaged equations become the sub-grid, or porous,
  shallow-water equations, in which storage and conveyance are separate
  functions of stage and momentum carries a pressure term from the sloping
  banks. Continuity is now exact; momentum still uses the constant-area form.
  That is defensible for a reservoir, where the dynamics are dominated by
  storage and buoyancy rather than fast shallow flow, and would not be for a
  flood wave down a valley. Horizontal face areas likewise still use the surface
  edge length at every layer.

- Meteorology is uniform over the water body, and the atmospheric longwave uses
  Swinbank's empirical emissivity rather than a measured incoming-longwave
  record. Wind-function coefficients are exposed because they are the least
  universal part of any heat budget and are routinely recalibrated per site.

- The expression parser binds unary minus tighter than exponentiation, so
  `-3 ^ 2` is 9, as in Excel, rather than -9 as in Python and mathematical
  notation. Inherited from the shared openswmm lineage and pinned by a test
  rather than changed, because changing it would silently give a different
  answer to an expression shared between the two projects.

- Hoisting a single OpenMP parallel region around the barotropic substep, to
  replace its five per-loop regions, was implemented and then **reverted**. An
  A/B of the two binaries measured it as neutral at one thread (ratio 0.998)
  with no usable multi-thread signal on the available hardware. It adds macros
  and a non-obvious constraint — the affected functions must then be called
  from inside a region — for no demonstrated benefit, so the simpler form
  stands. Worth retrying where scaling can actually be measured.
- Parallel scaling is measured only on a development container that exposes no
  core topology and where four threads still underperform three on some phases.
  The numbers above are indicative, not a characterisation of the code, and
  must be retaken on hardware whose topology is known.
- The external mode scales worst (1.68x) because its parallel regions are small
  and entered ~108 times per internal step. Hoisting a single parallel region
  around the subcycle is the next step.
- `computeLayerFluxes` and `diagnoseVerticalFlux` remain serial (1.0 ms/step).

- Open boundaries are barotropic-only. The internal mode **refuses to
  initialize** when any is present, because scalar transport across an opening
  is not yet defined — an inflow needs a per-layer composition and an outflow
  must export ambient values (F2.3b). Refusing beats silently transporting
  nonsense.

- Selective withdrawal implements the classic zone march plus the shifted
  velocity maximum. The V5.0 boundary-interference bisection is **deliberately
  not implemented**: the distributed technical note and the shipped v5
  source contradict each other on whether the withdrawal angle multiplies or
  divides the zone limit and on the boundary-interference coefficient (8/2
  versus 1/2), and the bisection's solved distance is defined from the
  interference boundary but consumed as a distance from the port. The classic
  coefficient-of-2 treatment stays until real W2 output can settle it.
- The DeGray comparison against W2 output is formulation-level only so far; the
  whole-reservoir run needs bathymetry ingest and meteorology and belongs with
  the Phase 8 validation release.

- The dry-bed front in a Ritter dam break trails theory by roughly 15 m
  independent of mesh refinement. The depth field converges cleanly at first
  order; the front speed does not, because the vanishing tip is a degenerate
  point for a first-order scheme. MUSCL reconstruction is the fix and arrives
  with the transport engines in Phase 5.
- Boundaries are reflective walls only. Inflows, withdrawals, and open
  boundaries arrive with the internal mode in Phase 3 (F3.4).
- Tiered local time stepping is deferred until the global-time-step
  solver can serve as its conservation reference.
- The barotropic solver itself is still first order: MUSCL reconstruction was
  added to *scalar transport*, not to the shallow-water depth and momentum, so
  the Ritter dry-front lag is unchanged. Extending reconstruction to the
  external mode is separate work.
- Spurious currents over steep bathymetry are bounded but not negligible
  (~1e-2 m/s for a 40 m step under a thermocline). The dominant residual is
  the below-bed extrapolation where adjacent columns have very different
  depths; hybrid z-sigma is the planned mitigation.
- The internal mode requires every column to be wet; three-dimensional
  wetting and drying, inflows, and withdrawals are deferred (F3.4).

[Unreleased]: https://github.com/hydrocouple/FVQual/commits/main
