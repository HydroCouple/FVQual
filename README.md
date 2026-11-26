# FVQual

[![Build and Test](https://github.com/hydrocouple/FVQual/actions/workflows/build_and_test.yml/badge.svg)](https://github.com/hydrocouple/FVQual/actions/workflows/build_and_test.yml)
[![Documentation](https://github.com/hydrocouple/FVQual/actions/workflows/build_docs.yml/badge.svg)](https://github.com/hydrocouple/FVQual/actions/workflows/build_docs.yml)
[![API docs](https://img.shields.io/badge/API%20docs-Doxygen-blue)](https://hydrocouple.github.io/FVQual/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-00599C.svg?logo=c%2B%2B&logoColor=white)](https://en.cppreference.com/w/cpp/20)
[![CMake](https://img.shields.io/badge/build-CMake%203.21%2B-064F8C.svg?logo=cmake&logoColor=white)](CMakeLists.txt)
[![C API](https://img.shields.io/badge/C%20API-C89--clean-555555.svg)](include/fvqual/capi/)
[![Python bindings](https://img.shields.io/badge/Python-C%20ABI%20bindings-555555.svg)](python/)

**A fully 3D, unstructured, finite-volume hydrodynamic and water-quality model —
the HydroCouple-native successor to CE-QUAL-W2.**

FVQual solves 3D hydrostatic free-surface flow and multispecies water quality on
an unstructured **tri-quad** horizontal mesh with **terrain-following sigma**
layers, using explicit finite-volume kernels written to be portable to GPUs via
Kokkos. It is designed as a HydroCouple 2.0 model component so it composes with
OpenSWMM networks, groundwater, and stream-temperature components.

> **Status: pre-alpha.** Complete and verified: the mesh and data model, the
> explicit barotropic mode with characteristic open boundaries, the 3D
> hydrostatic internal mode, MUSCL + FCT scalar transport, and reservoir
> operations — inflows and withdrawals with selective withdrawal, rating-curve
> structures, multi-port outlet blending, and internal transfers. Shared-memory
> threading is in place behind `FVQual_USE_OPENMP`. Reaction kinetics, a
> HydroCouple ABI 4 component, conservative interval inflows, complete
> in-memory checkpoints and Python C ABI bindings are implemented. Native
> mixed SWMM meshes can be imported without splitting quads. The worked
> [SWMM reservoir example](examples/swmm_reservoir/) includes recorded forcing,
> a live one-way surface-boundary bridge, native 1D SWMM outfall coupling and
> uncertainty experiments. The native CPU coupling is qualified through
> Composer on synthetic heat/species cases and a six-hour Bellinge event;
> its receiving basin and Bellinge inlet temperature are explicitly synthetic.
> [Scheduled reservoir operations](examples/reservoir_operations/) provide
> contextual Basic/Advanced editors, shared component/C API configuration,
> process-resolved water/heat/species budgets, scalar ledger plots, recorded
> layer prisms and transects. A worked mixed-mesh case qualifies CLI/Composer
> parity and independent irregular-schedule budgets.
> GPU execution, hydraulic feedback, coordinated durable restart and full
> solver derivatives remain unqualified. `CHANGELOG.md` records scientific
> accuracy and known limitations; a short DeGray run is not field validation.

## Why

CE-QUAL-W2 is the workhorse reservoir water-quality model, but it is laterally
averaged, structured-grid, serial Fortran with global module state, and carries
a monolithic kinetics module. FVQual targets the same problems in full 3D with a
modern, parallel, composable architecture:

| CE-QUAL-W2 | FVQual |
|---|---|
| Laterally averaged (2D longitudinal-vertical) | Fully 3D unstructured tri-quad |
| Structured segment x layer grid | Terrain-following sigma over any tri-quad mesh |
| Serial Fortran, global state | C++20, instance-safe, thread-parallel, Kokkos-ready |
| Hard-coded kinetics | Declarative MSX-convention `.rxn` kinetics |
| Bespoke ASCII IO | CF-1.11 / UGRID-1.0 NetCDF + HDF5 |
| Standalone | HydroCouple 2.0 model component |

## Design decisions

The architecture rests on a small number of deliberate choices. Each is
documented with its rationale in the header that implements it:

- **Hydrostatic core** with a `PressureClosure` seam so a non-hydrostatic
  pressure-correction module can be added later without refactoring.
- **Explicit barotropic-baroclinic mode splitting** plus implicit vertical
  column solves — no linear-solver dependency in the base configuration.
- **Eulerian MUSCL + FCT transport** as the *only* water-quality engine. The Lagrangian alternatives were dropped on mass-conservation
  grounds: ELM cannot conserve mass, and RWPT conserves it as statistical noise
  that nonlinear kinetics amplify. Water age and source attribution are
  transported species instead.
- **Reaction kinetics reused from OpenSWMM's MSX-convention engine** — its
  integrator core is already mesh-agnostic.
- **UGRID NetCDF/HDF5 is the sole full-field results format**.
- **SI units internally**, converted only at IO boundaries.

## Building

Requires a C++20 compiler, CMake 3.21+, and
[HydroCoupleSDK](https://github.com/cbuahin/HydroCoupleSDK) 2.0 built with
`-DBUILD_TOOLS=ON`.

```bash
# With an installed SDK
cmake --preset linux -DCMAKE_PREFIX_PATH=/path/to/sdk-install
cmake --build build --parallel
ctest --test-dir build --output-on-failure

# With an SDK source checkout (the usual path today — the SDK is not yet
# published to a package registry)
cmake --preset no-deps -DFVQual_SDK_SOURCE_DIR=../HydroCoupleSDK
cmake --build build --parallel
```

Presets: `linux`, `macos`, `windows` (vcpkg, NetCDF enabled) and `no-deps`
(no vcpkg, NetCDF off — useful for a quick core build).

| Option | Default | Purpose |
|---|---|---|
| `FVQual_BUILD_TESTS` | ON | Unit and verification suite (GoogleTest) |
| `FVQual_BUILD_APPS` | ON | `fvqual` command-line runner |
| `FVQual_BUILD_VERIFICATION` | ON | `fvqual-verify`, which regenerates the animated page |
| `FVQual_USE_NETCDF` | OFF | NetCDF-4 UGRID mesh input and layered output |
| `FVQual_USE_OPENMP` | OFF | OpenMP threading |
| `FVQual_SDK_SOURCE_DIR` | "" | Build HydroCoupleSDK in-tree from a source checkout |

## Layout

```
include/fvqual/     Public headers (mesh, hydro, transport, forcing, state, io)
include/fvqual/capi/  Public C API — the ABI firewall the bindings target
src/                Implementation — core library libfvqual-core
apps/fvqual/        Command-line runner
tests/              GoogleTest unit and verification suite
verification/       Animated benchmarks against analytical solutions
cmake/              Build helpers, package config, version template
```

## Verification

Physics lands with the test that pins it. The external mode holds still water
over irregular and emergent bathymetry to 1e-13, closes mass below 1e-10 while
sloshing, reproduces the analytic seiche period and wind-setup slope, and
converges at first order against Ritter's dam-break solution. The internal mode
holds stratified and homogeneous rest states to 1e-12 m/s, closes free-surface
continuity to 1e-14, preserves a uniform tracer, and reproduces the two-layer
internal seiche period to within 8% of theory and the lock-exchange front speed
to within 11%, the latter converging to 4% under mesh refinement.

```bash
ctest --test-dir build --output-on-failure
```

Numbers are one thing; watching a front arrive late is another.
[`verification/`](verification/) holds an animated page putting each benchmark
side by side with its analytical solution — open `verification/index.html`, or
regenerate it from the solver with `fvqual-verify`.

## Roadmap

The CPU hydro, transport, heat, reservoir structures, reaction kinetics and
HydroCouple component pathways are implemented. Current qualification focuses
on imported mesh geometry and numerical convergence, richer component
configuration and outputs, SWMM wrappers and scheduling, DeGray thermal
validation, and installation packaging. Python ensembles and central finite
differences are available. Analytic JVP/VJP covers reference kernels and actual
passive upwind transport/vertical diffusion on recorded hydraulic trajectories;
the controlled-column path also supports thermal/mixing/wind feedback, depth/area
scaling, changing-depth observations and blocked adjoint replay. Native nonlinear
multispecies kinetics and conditional species derivatives now integrate with the
SDK differentiable workflow and its logarithmic replay. General reservoir
thermal/momentum, settling and geometry derivatives remain unqualified. GPU work and full coupled derivatives
need their own verification evidence.
Verification gates and known limitations are recorded in
[`CHANGELOG.md`](CHANGELOG.md).

## Related projects

FVQual is part of the HydroCouple family and follows the same build, testing,
and documentation conventions as its siblings:

- **[HydroCouple](https://github.com/HydroCouple/HydroCouple)** — the
  component-based modeling interface definitions FVQual targets
  ([API docs](https://hydrocouple.github.io/HydroCouple/)).
- **[HydroCoupleSDK](https://github.com/HydroCouple/HydroCoupleSDK)** — the SDK
  FVQual builds on for the data plane, component lifecycle, and mesh IO
  ([API docs](https://hydrocouple.github.io/HydroCoupleSDK/)).
- **[openswmm.engine](https://github.com/HydroCouple/openswmm.engine)** — the
  SWMM engine whose public C API supplies the optional live surface-boundary
  bridge and whose MSX-convention reaction kernels informed FVQual kinetics
  ([Python docs](https://hydrocouple.github.io/openswmm.engine/)).

## License

MIT — see [`LICENSE`](LICENSE). CE-QUAL-W2 is itself MIT-licensed; FVQual's
attribution policy and provenance records are in
[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).
