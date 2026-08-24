# FVQual

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
> threading is in place behind `FVQual_USE_OPENMP`. **Not yet implemented:**
> water-quality kinetics, GPU acceleration, and the HydroCouple component.
> `CHANGELOG.md` records measured accuracy and known limitations. A public
> C API (`include/fvqual/capi/`) is the foundation for Python bindings.

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

Phases 0-3, the Eulerian transport engine (F5.1), and reservoir sources
including selective withdrawal are implemented. Next: rating-curve structures
(weirs, gates, pumps), Phase 6 (kinetics, reusing
OpenSWMM's MSX engine), Phase 7 (Kokkos), Phase 8 (HydroCouple component and
the CE-QUAL-W2 validation release).
Verification gates and known limitations are recorded in
[`CHANGELOG.md`](CHANGELOG.md).

## License

MIT — see [`LICENSE`](LICENSE). CE-QUAL-W2 is itself MIT-licensed; FVQual's
attribution policy and provenance records are in
[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).
