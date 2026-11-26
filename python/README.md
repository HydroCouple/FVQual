# FVQual Python

The package uses the public C ABI through standard-library `ctypes`. It
requires Python 3.10 or later and a shared FVQual core library. It has no NumPy,
Torch or Cython dependency. Configure with `BUILD_SHARED_LIBS=ON`, build with
the current SDK, set `PYTHONPATH` to this directory and `FVQUAL_LIBRARY` to the
core `.dylib`, `.so` or DLL. Development builds search the repository's default
`build` directory when no library is specified. Binary wheels are not supplied.

```python
from fvqual import Library, Mesh, Simulation

library = Library()
with Mesh.channel(library, columns=4, length=100, width=20, bed=-5) as mesh:
    with Simulation(mesh, layers=4, temperature=10,
                    species={"tracer": {"units": "kg/m3", "initial": 0}}) as model:
        token = model.save()
        model.advance(60)
        print(model.time, model.volume, model.mean_temperature)
        model.restore(token)
        model.release(token)
```

Context managers own and close handles even after exceptions. Simulations
copy their meshes. Serialize calls on one instance; separate instances can run
independently. Checkpoints are in memory, belong to their creating instance,
and are invalid after release or instance destruction.

`Mesh.from_arrays` accepts CSR triangle and native quad connectivity.
`Mesh.read` accepts UGRID or a `.mesh.json` package. Packages require metre
units and explicit bed datum; terrain cannot silently become bathymetry.
Package roughness and initial depth initialize the model.

`Simulation.apply` consumes already debited donor exports: unique frame ID,
elapsed start/end seconds, water in m3, heat in J relative to 0 degrees Celsius,
and transported species in kg. Declare every incoming transported species,
including zero loads. Repeating an identical ID returns the accepted receipt
without advancing or adding material again. A changed payload is refused.
Numerical failure rolls back the complete receiver interval.

The [worked example](../examples/swmm_reservoir/README.md) provides these runners:

- `fvqual.recorded`: forward runs, receipts and water/heat/species budgets.
- `fvqual.coupled`: optional native SWMM 2D boundary exports and recorded replay.
- `fvqual.uncertainty`: seeded independent forcing ensembles and central finite-difference sweeps.
- `fvqual.kernel_uncertainty`: analytic JVP/VJP checks and local variance versus Monte Carlo.

The three reference kernels cover decay, fully mixed heat and conservative
mixing. Actual solver derivatives now also cover passive upwind advection and
implicit vertical diffusion over multiple executed steps. This is a restricted
scalar trace; the conditional reactive component contract is described below.
It is not a full hydrodynamic/thermal adjoint. General whole-run sensitivities still use finite
differences. Inflow ensembles remain conditional on supplied SWMM forcing.

```python
from fvqual import Library, Mesh, Simulation
from fvqual.scalar_uncertainty import FixedObservation, local_covariance

library = Library()
with Mesh.channel(library) as mesh:
    with Simulation(mesh, species={"tracer": {"units": "kg/m3", "initial": .001}},
                    options={"advectionScheme": 0}) as model:
        model.begin_scalar_trace("tracer", max_steps=1000)
        model.advance(60)
        trace = model.capture_scalar_trace()

# An owned immutable trace outlives the model and shares completed records.
with trace:
    observation = FixedObservation((trace.steps + 1, trace.cells),
                                   [[(trace.steps, 0, 1)]])
    tangent = trace.jvp(initial=[.0001] * trace.cells)
    pullback = trace.vjp(observation.vjp([1]))
    covariance = local_covariance(trace, [{"initial": [.0001] * trace.cells}],
                                  observation)
```

`ScalarTrace` records the actual split order and the production tridiagonal
coefficients. JVP solves the tangent system; VJP solves its true transpose and
reverses the steps. Inputs are initial cell concentrations, per-step source
concentrations and a multiplier of the **recorded** vertical diffusion
operator. `source_ids` identify repeated source/interval contributions; gradients
for a forcing coefficient must accumulate all its occurrences. Uniform
initial-condition gradients accumulate all initial cells. Derivatives are at
the recorded baseline. `evaluate` replays other scalar inputs on the frozen
hydraulic trajectory. It does not rerun hydraulics or meteorology.

Recording refuses temperature, salinity, age, MUSCL/FCT, horizontal diffusion,
reactions, settling, internal transfers and open boundaries. External state
changes and a reached step limit are refused before the next step. An interval
failure restores the trace along with the solver; ordinary `advance` can retain
previous successful steps. Snapshot/restore/release remains in memory. Trace
storage grows with steps and mesh size; the reported bytes estimate payload,
not process peak memory. Logarithmic replay and the general SDK derivative
interfaces remain future work.

`FixedObservation` handles fixed spatial/time weights and unit conversions,
including their transpose and repeated observation contributions. Moving
depth/geometry weights require additional derivatives. `local_covariance`
accepts columns of a covariance factor, so correlations are explicit. It
excludes measurement error, structural discrepancy, numerical and mesh error.
Changing wind, heat, bathymetry or SWMM parameters is outside this trace's
derivative coverage. In particular, changing stable mixing options in a
thermal run is not equivalent to scaling a recorded Kz trajectory.

The [phase 6a review guide](../verification/degray/README.md#phase-6a-passive-scalar-sensitivity-review)
provides the Bellinge demonstration, qualification evidence and remaining gates.

The combined thermal/geometry/replay increment adds `ThermalColumn` and
`SigmaDepthObservation`. The column uses the actual `InternalModeSolver` for
forward evaluation. Shared production density, heat, wind and mixing kernels
provide analytic local derivatives; each native step is independently checked
against its differentiated column operator. Density, Richardson mixing, wind
shear and temperature feedback are active. Eighteen parameters include wind,
meteorology, mixing, optical attenuation, heat capacity and depth/area scaling.
The output blocks are temperature, x shear and y shear (velocity departures
from the depth mean), all surface to bed. Initial velocities are projected to
zero depth mean. Initial temperature is parameterized by sigma-layer index.

```python
from fvqual import Library, ThermalColumn, SigmaDepthObservation

with ThermalColumn(Library(), {"layers": 8, "duration_s": 3600}) as column:
    tangent = column.jvp(parameters={"wind_speed_m_s": .3, "depth_scale": .025})
    observation = SigmaDepthObservation(
        (column.time_count, column.state_size), column.layers, 10,
        [(column.time_count - 1, 2.6)])
    response = observation.jvp(column.states, tangent, depth_direction=.025)
    mapped = observation.vjp(column.states, [1])
    gradient = column.vjp(mapped["states"])
    gradient["parameters"]["depth_scale"] += mapped["depth_scale"]
```

Depth interpolation includes changing weights and their transpose. It refuses
extrapolation and geometry derivatives at a changing interpolation stencil.
Mixing/density/clamp boundaries are refused by the local derivative path;
`evaluate` still runs the native solver for ensembles across regimes. Forward
evaluation and failed derivatives leave the baseline unchanged. The adjoint
rebuilds one step's tape at a time and replays blocks between stored states;
`checkpoint_stride=1` stores every step, while zero chooses a square-root block
size. Reported tape/checkpoint counts are payload estimates, not process peak
memory or a logarithmic checkpoint algorithm.

The qualified column is closed and always wet, with zero Manning/Coriolis,
Richardson mixing and upwind transport. It has no sources, horizontal
circulation, salinity, kinetics, settling or SWMM derivatives. It does not
advertise the SDK's differentiable-component capability. The combined review
also supplies Bellinge **full native ensembles and explicitly named finite-
difference references**, preserving horizontal/thermal feedback there without
claiming an analytic Bellinge thermal adjoint. See the combined review section
of `verification/degray/README.md` for commands, GUI evaluation and evidence.

## Conditional reactive transport and SDK differentiation

Owned reactive traces chain native upwind transport and implicit vertical
diffusion through nonlinear multispecies reaction laws. The local analytic
tape follows the native adaptive Cash-Karp controller, including smooth
dependence of accepted internal step lengths. It refuses expression-domain,
scaling/max-error and switching boundaries, and unconverged reaction solves.
Forward ensemble evaluations remain available when local derivatives refuse.

```python
from fvqual import Library, Mesh, Simulation

deck = """[SPECIES]
A kg/m3 .0002
B kg/m3 .00005
[CONSTANTS]
k .00002
cref .001
[RATES]
A WATER -k*A*(1+B/cref)
B WATER k*A*(1+B/cref)
"""
with Mesh.channel(Library(), columns=2) as mesh:
    with Simulation(mesh, layers=3, options={"advectionScheme": 0},
                    operations={"reaction_deck": deck}) as sim:
        sim.begin_reactive_trace()
        sim.advance(3600)
        trace = sim.capture_reactive_trace()
    with trace:  # The owned history outlives the simulation.
        response = trace.jvp(constants={"k": 3e-6})
        a_response = trace.field(response, "A")
        weights = [[0.] * trace.state_size for _ in trace.times]
        weights[-1][0] = 1
        gradient = trace.vjp(weights)
        print(gradient["constants"]["k"])
```

Rows have shape `[steps+1, species*cells]`, with each species in face-major
layer order. `metadata` lists species/units and constant order. `initial` may be
a flat vector or a species-to-vector dictionary. `constants` uses named deck
constants; `evaluate` runs the native reaction operator with the recorded
hydraulic/transport trajectory. `vjp` returns named initial-state and constant
cotangents. A new trace is required after changing constants on a simulation.

`fvqual.reactive_uncertainty.local_covariance` propagates correlated covariance
factor modes through `FixedObservation`. Optional `gaussian_diagnostics`
requires NumPy and reports whitened sensitivity singular values, confounding
directions and a local Gaussian observation update with explicitly supplied
noise. These are conditional calculations; measurement error and structural
discrepancy are separate from parameter-only bands. The Bellinge synthetic
observation update is retained even though its withheld error worsens.

The loadable `FVQualComponent` exposes the existing HydroCouple
`IDifferentiableModelComponent` contract when the explicit
`conditional_species_derivatives` argument is 1 and preparation qualifies the
configuration. This selects upwind transport. Its declared states/outputs are
the supported pollutant/MSX concentrations; its differentiable argument is
`reaction_constants`, in deck units. SDK `DifferentiableWorkflow` carries state
cotangents through time using every-step or logarithmic checkpoints. Parameter
overrides are validated and applied before each update. Default components
retain their original transport scheme and do not advertise differentiation.

Hydraulics, temperature/salinity, geometry, Kz, source discharge/composition
and placement remain prescribed. Density-changing rate laws, self-shading,
settling, horizontal diffusion, internal transfers and open boundaries refuse
the conditional contract. No thermal/momentum, general geometry or SWMM
parameter adjoint is advertised. Checkpoint counts describe stored component
states; latest-interval operator history is bounded separately and process
peak memory is not qualified. The portable
`examples/reactive_transport/composition.json` exercises the new arguments.

CTest registers integration tests for shared builds with an available Python
interpreter. Set `SWMM_ENGINE_LIBRARY` to include the optional live bridge test.
Artifacts remain under `tests/output/python_*`. Zero-copy arrays, the complete
runtime control surface and binary distribution are outside this initial binding.
