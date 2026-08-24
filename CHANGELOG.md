# Changelog

All notable changes to FVQual are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

### Known limitations

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
