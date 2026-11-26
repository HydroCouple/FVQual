# SWMM exports into a mixed-mesh reservoir

This synthetic case preserves two triangles and one quad, extrudes four sigma
layers and transfers water, heat and tracer into a reservoir with an explicit
synthetic bed datum. Urban terrain is not a surveyed reservoir bed.

From the repository, with a NetCDF-enabled build and the current SDK:

```sh
build/apps/fvqual import-mesh examples/swmm_reservoir/reservoir.2dm \
  examples/swmm_reservoir/reservoir.mesh.json \
  --crs LOCAL --datum synthetic-bed --elevation-kind bed
build/apps/fvqual run examples/swmm_reservoir/composition.json
```

The converter writes a manifest and UGRID companion, preserving node/face IDs,
roughness and initial depth. It refuses to overwrite existing packages.
The composition resolves input `@uri` paths relative to its folder. Temperature
CSV and a reopenable manifest appear under `output/`; time is elapsed seconds.
`inflow_series.json` is the recorded-forcing argument of the loadable component.
Composer's mesh editor previews cell types, datum and bed/terrain semantics
before applying the argument.

Load the built component directory into Composer, or run its headless route:

```sh
HydroCoupleComposer --run examples/swmm_reservoir/composition.json \
  --components build/components
```

The component exposes an explicit `epoch_julian_day` argument, defaulting to
1970-01-01 UTC for synthetic runs. Composer records physical Julian days for
timed components; the FVQual CLI records elapsed seconds. Set the epoch to
the real event start before combining calendar-based forcing or observations.

With a shared core and `PYTHONPATH=python`:

```sh
python -m fvqual.recorded examples/swmm_reservoir/recorded_event.json \
  tests/output/my_recorded_run
python -m fvqual.uncertainty examples/swmm_reservoir/recorded_event.json \
  tests/output/my_ensemble --samples 128 --seed 20261010
python -m fvqual.kernel_uncertainty tests/output/my_kernel_uncertainty.json
```

Output directories must be new. Results retain receipts, history, budget
residuals, input/binary hashes and revisions. Ensembles record distributions,
seed, every success/failure, prediction bands and finite-difference step sweeps.

The optional live example requires exact exchange times, surface heat transport
and accepted boundary receipts in the native SWMM library:

```sh
python -m fvqual.coupled examples/swmm_reservoir/recorded_event.json \
  examples/swmm_reservoir/surface_donor.inp tests/output/my_live_run \
  --swmm-library "$SWMM_ENGINE_LIBRARY"
```

The driver prescribes outward flow of 0.01 m3/s at a **2D surface boundary**,
converts TSS native mass to tracer kg explicitly, checks combined inventories
and requires recorded replay to match the live result. It never withdraws an
export a second time. Pending exports must be accepted before advancing the
donor; receiver failure stops the run because the donor has no rollback contract.

This is one-way surface coupling. The native 1D component qualification below
has a separate ownership contract; bidirectional ownership and receiving-stage
feedback remain separate work. JVP/VJP support covers three smooth reference kernels. Full-run
uncertainty uses ensembles and finite differences; no full hydrodynamic adjoint,
coupled SWMM gradient, field calibration or posterior uncertainty is claimed.

## Bellinge mesh qualification review

The October 10 mesh phase uses the existing Bellinge export at
`openswmm.gui/examples/bellinge_2d/BellingeSWMM_v021_nopervious.2dm`.
It contains **78,285 nodes and 155,688 triangles, with no quads**. Its source
terrain and EPSG:25832 reference are preserved in the imported package.
The separate synthetic mixed case checks native quad support.

The `MeshQualification` tests generate review artifacts under
`tests/output/mesh_review/`. With `FVQUAL_BELLINGE_MESH` set to the Bellinge
export, they also generate `tests/output/bellinge_mesh_phase1/`:

| Artifact | Purpose |
| --- | --- |
| `bellinge.terrain.mesh.json` and `.nc` | Original terrain, identities and source attributes; not reservoir bathymetry |
| `geometry.json` | Full geometry inspection, including closure, topology and boundary mappings |
| `bellinge.synthetic.mesh.json` and `.nc` | Same footprint/connectivity; explicitly synthetic flat bed at −5 m and wall boundaries |
| `bellinge.synthetic.composition.json` | Two-layer, 10 °C, always-wet review; two 0.02 s intervals |
| `bellinge.synthetic_output/run.json` and `temperature.csv` | Actual component execution and recorded temperatures |
| `../mesh_review/mixed.mesh.json` and `.nc` | 81 nodes, 64 triangles and 32 retained quads |

The synthetic review preserves original terrain and native boundary sections as
provenance. It does **not** predict Bellinge flooding or run coupled Bellinge
network hydraulics. Real receiving-reservoir bathymetry and physical boundary
ownership still need to be defined for coupling.

The reservoir results phase below supersedes the first mesh-only GUI handoff.
Recorded fields now carry owned mesh geometry, horizontal and vertical
references, physical quantities and time-dependent layer interfaces. Opening
the matching UGRID package companion also adopts its horizontal CRS. The saved
resource editor, completed-result selection and map/3D view handoff have been
corrected.

Numerical evidence is in `tests/output/mesh_corpus.json`,
`mesh_advection_convergence.json`, `mesh_diffusion_convergence.json`,
`mesh_surface_seiche.json` and `mesh_internal_seiche.json`. The tests cover
geometry defects, intentional holes/disconnected domains, still water,
conservation, transport bounds, tri/quad/mixed refinement, wave periods,
open-boundary budgets and checkpoint replay. These qualify the tested
always-wet CPU paths; 3D wetting/drying and full-model derivatives remain
separate gates.

## Reservoir settings and recorded results review

The current review files are in `tests/output/reservoir_gui_phase2/`:

| File | Configuration |
| --- | --- |
| `bellinge.composition.json` | Full Bellinge footprint, two layers, sigma `[0, 0.2, 1]`, initial 22/12 °C profile, conservative tracer 0.002 kg/m³ |
| `mixed.composition.json` | 64 triangles and 32 quads, three nonuniform layers, temperature varying across faces and layers, horizontal diffusivity 0.01 m²/s |
| `bellinge_output/run.json`, `mixed_output/run.json` | Completed two-step runs, 13 recorded fields and a shared mesh sidecar per run |
| `host_parity.json` | Composer and FVQual CLI comparison: all 4,044 mixed-case rows agree exactly |
| `phase_review.json` | Build/artifact hashes, checks and scope of the visual review |

Both cases run for two seconds with one-second exchange intervals. They use
synthetic bathymetry, wall boundaries and always-wet columns. The Bellinge case
is a geometry and reservoir-results exercise, **not a coupled Bellinge flood
prediction**. Initial uniform profiles make Bellinge a simple layer-selection
check; the mixed case provides spatially varying temperature.

With the rebuilt Composer and FVQual `build/components` loaded:

1. Open `bellinge.composition.json`. Select the reservoir on the **Composition**
   canvas. In the argument panel choose **Argument → mesh**, then
   **Mode → Basic**, then **Source → File or address → Preview mesh**. The saved relative URI should
   be visible.
   Check 78,285 nodes, 155,688 triangles, zero quads, EPSG:25832 and the
   `synthetic-flat-bed` datum. Use **Cancel** to retain the configuration.
2. Use the **Argument** dropdown to check duration, exchange interval, layer count,
   sigma interfaces and initial temperature. Choose **Mode → Basic** for scalar
   or profile editing; use **Source → Enter values** when a source choice is offered; tables save with **Apply values**. Switching argument or
   mode does not change saved settings. **Home → Save**, reopen the composition, then
   select **Home → Run**. Expect **Finished** and two completed steps.
3. The recorded temperature map opens automatically. **Layer Properties →
   Information** should report 155,688 features, ETRS89 / UTM zone 32N and
   extent `(579466.312, 6132060.139)` to `(587389.128, 6136464.184)`.
   The map may retain its separate display CRS; the layer's source CRS is
   EPSG:25832 without manual assignment.
4. In **Layer Properties → Rendering**, set **Vertical layer (surface = 1)**
   to 1 or 2, then select **Apply** or **OK**. Expect approximately 22 °C at
   the surface and 12 °C at the bed. The **Attributes** and **Plot** panels use
   the chosen layer. To enter a number, click the spin field and type it.
5. Select one cell with the **Map → Select** tool or select a row in
   **Attributes**, then open **Profile**. Expect the full two-point water-column
   profile: about 22 °C at −0.5 m and 12 °C at −3 m. Elevations are midpoints
   of the recorded interfaces in the synthetic datum, not layer numbers.
6. Step between the two times. Expect `2026-01-01 00:00:01` and
   `2026-01-01 00:00:02`. **Map → From Components** adds seven spatial fields
   from the matching completed run. In **Runs**, expand the run and component,
   choose a spatial field and select **Show on Map**. The six scalar inventory
   and accepted-transfer records are in the CSV and run catalogue.
7. Open `mixed.composition.json` and repeat the run. Confirm the three-layer
   selector and spatially varying values. Its mesh preview should show 64 teal
   triangles and 32 purple quads. Reopen `bellinge_output/run.json` through
   **Runs → Open Run** to check recorded results independently of a new run.

The named optional settings preserve the legacy geometry arguments when
omitted. Surface elevation accepts one value or one per face. Temperature
accepts a uniform value, a surface-to-bed profile or full face-major cells;
velocities accept uniform or full cell fields. Sigma interfaces run strictly
from surface 0 to bed 1 and must agree with any explicit layer count.
Conservative species map identifiers to initial concentrations in kg/m³.
Invalid shapes, dry columns, unsafe species identifiers and nonfinite values
are rejected before execution.

The outputs include temperature (°C), stage/depth (m), horizontal velocities
(m/s), cell and total water volumes (m³), species concentrations, heat inventory
(J relative to 0 °C), species inventories (kg where concentration units are
recognized), and accepted water/heat/species transfers. Recognized reaction
concentration units carry their dimensional scale; undeclared deck units stay
explicitly labelled as undeclared. Inventories and accepted transfers are not
a complete balance of meteorological, reaction, outlet and boundary terms.

Verification: FVQual **356 passed / one conditional skip**, SDK **425 passed /
one conditional skip**, and all **65 Composer suites passed** after the eight
graphics/local-server suites were rerun outside the sandbox. The final path
field change also passed the four affected Composer suites. DeGray and actual
Bellinge tests were enabled. Logs are `build/test-reservoir-phase-full.log`,
SDK `build-macos-full/test-reservoir-phase-full.log`, and Composer
`build/darwin-abi4/test-reservoir-phase-{full,unrestricted,final}.log`.

Remaining gates include complete meteorology/source/outlet/structure authoring,
full budget residuals, durable coupled
restart and feedback, recorded 3D prisms/transects, broader scientific
validation and full-solver derivatives. The new recorded layer/profile path
does not complete all of those capabilities. Older incompatible component
collection plugins are refused; only the rebuilt FVQual plugin is qualified
by this review.

## Native SWMM outfall coupling and Bellinge review

The phase 4 operations increment adds dedicated SWMM mapping and reservoir
editors, complete process ledgers, recorded prisms/sections and scalar plots.
The current GUI review is
`tests/output/swmm_coupling_phase4/bellinge/review.composition.json`; its
completed results are `inline_output/run.json`. See the
[operations review instructions](../reservoir_operations/README.md) for this
increment. The phase 3 evidence below describes the earlier native coupling
gate; the authoring and budget work it lists as remaining is delivered by
phase 4. Scientific validation, feedback, durable restart and actual-solver
derivatives remain open.

`live/composition.json` connects the SDK-free native OpenSWMM ABI 4 component's
`outfall_transfers` to FVQual's `integrated_inflows`. Its setup is in
`live/configuration.json`. Build OpenSWMM with `OPENSWMM_WITH_HYDROCOUPLE=ON`,
the same HydroCouple interface checkout/toolchain as Composer, and nlohmann-json
(the optional vcpkg `hydrocouple` feature). Load both component directories.

```sh
python3 verification/swmm_coupling.py \
  --composer ../HydroCoupleComposer/build/darwin-abi4/HydroCoupleComposer.app/Contents/MacOS/HydroCoupleComposer \
  --swmm-components ../openswmm.engine/build/fvqual-phase3/components \
  --engine-library ../openswmm.engine/build/fvqual-phase3/src/engine/libopenswmm.engine.dylib \
  --bellinge ../openswmm.gui/examples/bellinge_2d/BellingeSWMM_v021_nopervious.inp
```

The paths above are the qualified macOS development builds; select the matching
executable/library paths on other platforms. The script needs only Python's
standard library and a NetCDF-enabled FVQual build. It generates the Bellinge
review, executes both native components through the actual Composer host,
compares transfers with an independently driven public SWMM API instance and
requires exact recorded replay. It also checks distinct instances from one
factory and nine failure cases: epoch, weights, mesh order, unknown outfall,
dual forcing, count units, missing receiver species, reversal and multiple
receivers. Donor failure audits retain committed native boundary amounts.

The native ledger books SI amounts where SWMM books outfall continuity, before
report filtering; it does not withdraw them again. MG/L and UG/L exports map to
kg. The heat fixture uses native heat transport, including its native heat
capacity. A fixed-temperature assumption uses FVQual's 1000 kg/m³ and
4186 J/(kg·K) basis. Count-based concentrations, MSX/non-LEGACY export accounting,
internal routed outfalls, adapted packet connections, reversed flows and
multiple receiving components are refused. The receiver applies each interval
atomically and acknowledges its transfer ID before another native advance,
including the final interval. Numeric CSVs omit protocol strings; packet and
acknowledgment JSON files retain the full control records.

Current evidence in `tests/output/swmm_coupling_phase3/qualification.json`:

| Case | Result |
| --- | --- |
| Synthetic native heat and two pollutants | Five irregular intervals, 3.358561 m³; 400 recorded numeric values match replay exactly |
| Bellinge native runoff and 1D network | Nine outfalls, twenty irregular intervals through six hours, 12,162.107480 m³; 1,160 numeric values match replay exactly |
| Bellinge inventory closure | Maximum water residual 8.19×10⁻¹⁰ m³; heat residual 0.490 J against a ~10¹² J inventory (relative 1.11×10⁻¹³) |
| Failure and factory checks | Nine expected refusals; two independent native instances created by one component factory |
| FVQual | 359 passed, one conditional NetCDF-disabled branch skipped; actual Bellinge and DeGray fixtures enabled |
| Native engine and Composer | Receipt, continuity and checkpoint suites pass; five affected Composer suites pass, with simulation manager rechecked after final cleanup change |

The review lives under `tests/output/swmm_coupling_phase3/bellinge/`:

| File | Purpose |
| --- | --- |
| `composition.json` | Live Bellinge network → synthetic mixed-mesh receiving basin |
| `configuration.json`, `mapping.csv` | Explicit object/column/face mapping, normalized 0.75/0.25 weights, common epoch and temperature assumption |
| `bellinge_network.geojson`, `bellinge_outfalls.geojson` | 1,044 source network links and nine source outfalls in EPSG:25832 |
| `synthetic_transfer_mapping.geojson` | Eighteen review connections with outfall, receiver face and weight attributes; lines are not surveyed pipe alignments |
| `receiver.mesh.json` and `.nc` | Synthetic 100×100 m basin, −5 m flat bed, two triangles and one quad, three sigma layers |
| `provenance.json` | Original source hash, ownership and physical assumptions |
| `live_output/run.json`, `fields.csv` | Completed coupled run and numeric results |
| `qualification.json` | Mapping, inventory and live/replay evidence; names the specific donor audit directory |
| `donor_output/<run>/` | Setup, native boundary audits, exported packets, acknowledgments, replay forcing and native report |

To evaluate in the rebuilt GUI:

1. Open this phase's `composition.json` with both rebuilt component directories
   loaded. **Home → Run** should finish both components after twenty exchanges.
   The receiving temperature map opens automatically.
2. Use **Layer Properties → Rendering → Vertical layer (surface = 1)** to
   inspect all three layers. Compare the first and final event times; select
   one face and inspect **Profile** and **Plot**. Initial temperatures are
   22/16/12 °C. At 2012-06-29 06:01, column 0's profile is approximately
   20.237/17.497/13.270 °C; stage is about 1.216 m above its initial level.
   **Map → From Components**, or **Runs → Show on Map**, adds
   stage/depth and other completed spatial fields.
3. Use **View → Add Vector Data…** to add the network, outfall and mapping
   GeoJSONs. **Zoom to layer** on the network frames the full source model;
   zooming to the receiver temperature layer reveals the mixed basin. Mapping
   attributes identify the nine outfalls, destination faces and split weights.
4. Inspect `mapping.csv` and the setup file before accepting a physical mapping.
   Select either component on **Composition** to use the **Argument** dropdown
   and **Mode → Basic / Advanced**. SWMM Basic selects its setup file; detailed
   mapping edits currently use that JSON file. A dedicated mapping editor
   remains part of the next authoring phase.

Bellinge's original civil model start, 2012-06-29 00:01, is retained as the shared
Julian-day epoch; the source does not identify a timezone. Native 2D is disabled
before initialization so network exports belong to the separate receiver.
The six-hour coupling window ends before the original 48-hour native input;
the native report may therefore describe the early stop as interrupted.
The basin location, bathymetry, wall boundaries, always-wet state, mapping and
18 °C inlet temperature are synthetic. Bellinge has no pollutant or thermal
forcing in this model; those exports are exercised separately by the native
heat/species fixture. This phase qualifies software exchange, not a reservoir
forecast or field validation. Receiving-stage feedback, joint durable restart,
full network authoring, complete budgets and full-solver JVP/VJP remain open.
