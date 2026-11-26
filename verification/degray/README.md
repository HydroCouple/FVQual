# DeGray Reservoir, March–December 1980

The **phase 5c circulation and geometry review** at the end of this document is
the latest workflow. It preserves the earlier registered phase 5 evidence. The
earlier results and diagnosis sections record historical experiments; their
temperature statistics are not qualified results for the corrected solver.
The phase 5 audit found that thickness-based vertical diffusion and surface
heating did not conserve heat when storage area varied with depth. Their
finite-volume replacements use cell volumes and interface areas, with separate
analytical balance checks. Superseded runs are retained under
`tests/output/degray_phase5/pre_conservation_fix/` and excluded from acceptance.

An end-to-end run of CE-QUAL-W2's DeGray example through FVQual: real
bathymetry, real meteorology, real inflows and outflows, real initial
conditions.

    fvqual-degray <w2-input-directory> <output-directory>
    python3 plot_degray.py [output-directory] [bth1.csv]
    python3 verify_figures.py [output-directory]

The run writes `degray_daily.csv` (one row a day), `degray_profiles.csv` (the
whole profile at the dam, every day) and `degray_final_profile.csv`. The plot
script draws five figures from them; the verify script restates every claim the
figures make as an inequality over the same data and exits non-zero if one of
them is not true.

## What this is

**Not a code-to-code validation.** The W2 distribution ships DeGray's inputs
but no output — `el_stats.opt` contains only `NaN` — so there is nothing here
to compare against. This is an exercise judged on physics, which catches gross
errors and would not catch a subtle formulation difference.

**Not a test of three-dimensionality.** W2 is laterally averaged, so its
schematisation is one segment wide and the mesh built from it is a chain of 30
columns. The vertical physics is exercised fully; the horizontal
three-dimensionality FVQual exists for is not exercised at all.

## Configuration

| | |
|---|---|
| Mesh | 30 columns × 30 sigma layers, from `bth1.csv` |
| Vertical datum | `EBOT` = 66.35 m, from the control file's location card |
| Initial volume | 7.292e8 m³ over 5.092e7 m² at 123.80 m |
| Period | Julian day 64.5 to 358.7, 1980 |
| Meteorology | `met.npt`; solar radiation **computed** (W2 sets `SROC = OFF`) |
| Cloud cover | W2 records tenths; converted to a fraction at the boundary |
| Inflow | `qin_br1.npt` / `tin_br1.npt`, density-placed at the head |
| Outflow | `qot_br1.npt`, at the dam |

## Result

| Quantity | Value | Expected for DeGray |
|---|---|---|
| Initial surface / bottom | 7.4 / 6.5 °C | matches `vpr.npt` |
| Peak surface temperature | 30.9 °C, day 197 | ~30 °C in midsummer |
| Hypolimnion through summer | 6.5 → 8.4 °C | retains winter temperature |
| Maximum stratification | 23.8 °C, day 192 | ~22 °C |
| Thermocline, midsummer | 5.2 m mean | 5–15 m |
| Thermocline, October | 13.0 m mean | deepening |
| Final profile | isothermal ~9.9 °C | isothermal ~9–10 °C |
| Volume excursion | −11.9% then +1.4% | summer drawdown, autumn refill |

The reservoir starts nearly isothermal, stratifies through spring, holds a cold
hypolimnion all summer while the surface reaches 31 °C, then deepens its
thermocline through autumn and turns over by December. Surface temperature
tracks below the instantaneous equilibrium temperature throughout, which is what
thermal inertia does.

## Figures

| File | What it shows |
|---|---|
| `degray_isotherms.png` | Temperature against depth and time at the dam — the whole run in one frame |
| `degray_temperatures.png` | Surface and bottom against the equilibrium temperature, and the stratification cycle |
| `degray_profiles.png` | Five profiles through the year, as a thermal survey would report them |
| `degray_budgets.png` | Whole-reservoir heat content and volume |
| `degray_stage_storage.png` | The model's stage-storage against the raw bathymetry file's |

The stage-storage figure is the one that can fail. `plot_degray.py` re-derives
DeGray's storage curve in Python straight from `bth1.csv`, independently of
anything the C++ computed, and plots every day of the run against it. The two
agree to **0.22% at worst**, and the recovered `dV/dz` is 4.98e7 m² against the
reservoir's 5.09e7 m² plan area at full pool. The residual is systematically
positive and correlates with drawdown at r = −0.76, which is the 30-layer sigma
mesh sampling a storage curve the file resolves in 36 layers — a discretisation
difference, not a leak.

The steepest-gradient line on the isotherm figure is masked below 2 °C of
stratification. Spread 2 °C over 55 m and the tracker is choosing between layer
differences the size of the transport's own round-off. It still oscillates
through November, and that part is real: during turnover the profile has two
competing maxima — a sharp shallow one near 2 m and the eroding metalimnion base
near 20–40 m — and they trade places from day to day.

## Wind stress, and what measuring against the casts showed

Matched runs, 30 layers, against the measured casts in `prf_dam.npt`:

| | overall RMSE | 5–10 m | 20–30 m | 30–60 m |
|---|---|---|---|---|
| no wind stress | **2.29 °C** | −3.91 | −0.83 | −0.55 |
| Richardson + wind | 3.59 °C | −1.11 | +3.22 | +3.65 |
| k-epsilon + wind | **3.42 °C** | −1.23 | +2.95 | +3.39 |
| no mixing at all + wind | 5.87 °C | +0.84 | +5.78 | +6.21 |

**Wind stress fixes the epilimnion.** Without it the modelled April surface layer
is 3 m against a measured 14 m, and the 5–10 m band sits at −3.91 °C.

**The interior bias it costs is advective, not mixing.** Removing mixing entirely
makes the interior *worse*, so heat is being carried down by the wind-driven
overturning cell rather than mixed down; mixing partly counteracts it by holding
heat near the surface where the budget can shed it.

**And the cell is inflated by the schematisation.** A laterally-averaged grid is
one segment wide, so every return path is forced through the vertical. The same
synthetic basin under the same wind transports 39% less heat below 20 m when it
is eight cells wide instead of one. Whether that accounts for *all* of the
+3 to +4 °C cannot be decided on this dataset, which is one segment wide by
construction — the caveat at the top of this file, arriving as a number.

Note that no-wind still scores best overall. It gets there with an epilimnion a
fifth of the measured thickness, which is worth stating plainly: a better summary
statistic from a model that cannot stir its surface layer is not the better model.

## What this run found

Two defects that 247 passing unit tests did not, both in the heat budget:

- **The wind function was about eight times too weak** — 10 W/m² of evaporative
  cooling at a 1000 Pa vapour deficit where the physical value is nearer 80.
  Every heat-budget gate computed its expectation from the same function it was
  testing, so the budget was still self-limiting, still conservative, still
  reached "its" equilibrium. Only a comparison against outside numbers found it.
  Peak surface temperature was 124 °C.
- **Cloud cover units.** W2 records cloud in tenths, on a 0–10 scale; the
  budget wants a fraction. The cloud factor is quadratic, so a value of 10
  multiplied the atmospheric longwave by eighteen. Peak surface temperature was
  still 76 °C after the wind function was fixed.

Both now have gates: one comparing the terms against magnitudes from outside the
code, one ensuring an out-of-range cloud value cannot run away.

A third turned up when the figures were checked against their own data. The
daily report was scheduled as `nextWrite = day + 1.0`, rebasing on the time
actually reached rather than advancing a fixed schedule, so each step's overshoot
carried forward: the report drifted from 12:05 in March to 17:00 in December. The
run was sampling a steadily later hour of the diurnal cycle and aliasing it into
the seasonal signal. The trajectory was never affected — the same 119,750 steps
run either way — but the reported peak surface temperature fell from 31.4 °C to
30.9 °C once the sample was consistent, and the daily series lost a wobble that
had nothing to do with the reservoir. Caught by asserting that a caption reading
"near solar noon" was true of the timestamps, which is not where one expects to
find a timing bug.

## How small a difference is worth believing

DeGray amplifies round-off to a saturated few hundredths of a degree within a
single simulated day. Any single-run difference below about **0.1 °C** is not
interpretable, and several of the questions worth asking about this model produce
differences smaller than that.

The measurement. The mesh can be built by two routes — `toMeshDefinition` and
`toPlannedMeshDefinition` at one lane — that produce *geometrically identical*
meshes: node coordinates equal to the last bit, face areas agreeing to 5e-14
relative, the same 30 faces and 91 edges. They differ only in the order faces
list their nodes, which changes the order edges are numbered and therefore the
order fluxes are summed. Running both:

| | surface temperature difference |
|---|---|
| after one simulated day | 0.0039 °C |
| mean over 15 days | 0.021 °C |
| maximum over 15 days | 0.069 °C |

**This is not a defect and cannot be tuned away.** The model is deterministic —
two runs of the same binary on the same mesh are bit-identical. What differs is a
perturbation at the 1e-16 level, and the barotropic mode amplifies it to a few
hundredths of a degree within a day, after which it saturates rather than
continuing to grow. Compensated summation or a canonical edge order would change
the seed, not the amplification; the floor would come back.

Three things this qualifies, all measured this session:

| | difference | against the floor |
|---|---|---|
| Manning friction on the conveyance depth (F1.5b) | 0.035 °C | at the floor — no interpretable effect |
| adaptive lateral refinement, 67 columns | 0.029 °C | below the floor — no interpretable effect |
| uniform two lanes | 0.25 °C | 3–12x the floor — a real difference |

So the conclusions stand, and the first two stand for a sharper reason than "the
numbers looked close": the effect is smaller than the model's own reproducibility.

What to do instead of chasing it: compare quantities that average the noise out —
volume-weighted means over a season rather than a point at an instant, or an
ensemble of runs seeded by trivial perturbations, reporting the spread alongside
the difference. A single pair of runs cannot resolve a tenth of a degree here.

## Registered seasonal evaluation (phase 5)

The October 11 review completed eight full-season experiments. All passed the
registered numerical checks: the largest relative water residual was
3.92e-13 and the largest relative heat residual was 3.97e-12. The main build
passes 378 tests with one expected NetCDF-disabled-branch skip; the actual
DeGray and Bellinge fixtures are enabled. Five affected Composer suites pass,
including the full-season saved-data reader. Native visual inspection is
pending because the Mac is locked.

These are the errors on **970 common supported measurements** across all eight
cases; individual case statistics also retain their own support in `metrics.json`.

| Case | Common-support RMSE (°C) |
| --- | ---: |
| Baseline, 30 layers, 300 s ceiling | 7.5185 |
| 150 s ceiling | 7.5280 |
| 75 s ceiling | 7.5394 |
| 15 layers | 8.0483 |
| 60 layers | 7.4391 |
| Inferred two-lane geometry | 5.7520 |
| Wind-stress diagnostic | 2.6106 |
| Point-outlet diagnostic | 6.9038 |

No case is accepted as a field model. Timestep differences are much smaller
than the observed mismatch. Vertical changes, geometry and wind stress affect
the results, but neither these finite-grid comparisons nor the retrospective
autumn split establish asymptotic convergence or independent validation.
The baseline's own 991 supported points give RMSE 7.4450 °C and mean bias
+5.3251 °C. The wind diagnostic's own-support development/autumn RMSEs are
2.1915/3.2695 °C, and its mixed layer remains too shallow.

The thermal matrix used the archived executable/core in
`seasonal/review/executed_binary/` before the later chemistry balance additions.
Chemistry is inactive in those experiments. A separate one-day comparison of
the final main build matches all 1,864 daily, profile, field and budget rows
exactly; `final_build_thermal_parity/parity.json` records its hashes. The archive
preserves the actual seasonal binary and local dependency fingerprints.

`seasonal_validation.py` runs the current driver on the supplied historical
variable-flow records, registers their hashes and the 24 measured dam casts,
and writes per-cast/depth/split errors, independent water/source-heat audits,
process inventories and run provenance. Three casts precede the March 4 start;
14 in-season casts are development observations and seven from day 273.5
onward are reserved from future tuning. Earlier work already inspected all
observations: this is a retrospective reservation, **not blind validation**.
Only one measured station is supplied, so no station holdout is possible.

The driver writes the true initial state and lands exactly on daily reports,
forcing knots and the requested endpoint. Weather and flows remain linearly
interpolated W2 Julian-day inputs, sampled at the opening of each step. Cloud
tenths become fractions; the wind and evaporation reference heights remain
10 m and 2 m. Required forcing must cover the requested interval.

The comparison brackets time within daily reports and interpolates depth
between cell centres only. It refuses temporal/depth extrapolation and records
excluded observations. Surface measurements above the shallowest centre do
not silently become layer-centre temperatures. Cross-case error comparisons
use the intersection of supported observation depths.

```sh
python3 verification/degray/seasonal_validation.py \
  --driver build/apps/fvqual-degray \
  --inputs /path/to/v5/Examples/DeGray/InputFiles \
  --observed /path/to/v5/Examples/DeGray/prf_dam.npt \
  --output tests/output/degray_phase5/seasonal --jobs 4
```

Cases use 300/150/75 s timestep ceilings, 15/30/60 sigma layers, an inferred
two-lane geometry, and a diagnostic without momentum wind stress (weather
heat exchange remains active). The two-lane case is a geometry sensitivity,
not formal horizontal convergence: lateral depths are inferred from W2 width
curves, and changing lane count changes the deepest sampling column's bed.
W2's EBOT of 66.35 m is a file datum; the single-lane active dam bed is 68.35 m.

The supplied V5 control workbook describes a testing scenario with constant
12.6 m³/s input/output and a small gate release. It specifies a point outlet
at 115 m, width zero and withdrawal angle 3.14 radians. This differs from the
historical-flow, uniform-withdrawal driver. Keep those scenarios distinct. A
separate diagnostic tests that outlet with the historical flow records:

```sh
python3 verification/degray/seasonal_validation.py \
  --driver build/apps/fvqual-degray \
  --inputs /path/to/v5/Examples/DeGray/InputFiles \
  --observed /path/to/v5/Examples/DeGray/prf_dam.npt \
  --output tests/output/degray_phase5/outlet_diagnostic \
  --cases outlet115 --jobs 1
```

This is not reproduction of the workbook's complete V5 configuration. The
selective-withdrawal heat audit uses the shared envelope formula; the uniform
withdrawal audit uses an independently computed thickness-weighted mean.
The driver reports that distinction in `driver.json`.

After the seasonal runs finish, use a Python environment with numpy and
matplotlib to generate the scientific plots and native recorded views:

```sh
python3 verification/degray/review_artifacts.py \
  tests/output/degray_phase5/seasonal \
  --outlet-case tests/output/degray_phase5/outlet_diagnostic/outlet115
```

Review `review/profiles.png`, `seasonal_errors.png`, `refinement.png`,
`deep_temperature.png`, `budgets.png` and `common_support.json`. Numerical
qualification requires the full interval, finite output, no barotropic floor
violation and water/heat inventory residuals below 1e-8 relative to inventory.
These checks do not establish field acceptance; large thermal bias and lack
of asymptotic refinement must remain visible. No parameters are fitted.

### Review in Composer

Rebuild Composer, then launch its app with the Bellinge composition and
repeatable `--open-run` arguments for the saved DeGray runs. No producing
component library is required to reopen these DeGray results.

```sh
HydroCoupleComposer /path/to/bellinge/review.composition.json \
  --components /path/to/FVQual/build/components \
  --components /path/to/openswmm.engine/build/fvqual-phase3/components \
  --open-run /path/to/seasonal/baseline/gui_run.json \
  --open-run /path/to/seasonal/lanes2/gui_run.json \
  --open-run /path/to/seasonal/no_wind_stress/gui_run.json \
  --open-run /path/to/outlet_diagnostic/outlet115/gui_run.json
```

1. Expand a DeGray run and its component in **Runs**. Select
   `bottom_temperature_C` and choose **Plot** to inspect the saved seasonal
   history, including the initial state.
2. Select `temperature` and choose **Show**. Use **Full Extent**, open layer
   properties and choose **Graduated Colour**, coloring by value. Select a
   vertical layer, then move the time control through spring, summer and
   autumn. The clock represents the W2 1980 dates.
3. In the map, activate **Section**, then drag and release across the chain.
   Inspect the full saved column depth. Use **3D** for layer geometry.
   Footprints and lanes are W2 schematization; prismatic drawings do not
   depict hypsographic cell volumes or surveyed DeGray bathymetry.
4. Compare baseline, wind-stress diagnostic and point-outlet histories against
   the measured-profile figures. Check whether improvements persist on the
   reserved autumn dates and exceed numerical refinement differences.
5. Open the Bellinge saved run to review coupling ledgers. Its receiving basin
   remains synthetic; these DeGray results do not validate that receiver.
   Hide the previous mapped field before switching between DeGray and Bellinge;
   their coordinates and calendar ranges differ.

Next scientific work requires diagnosing the remaining thermal bias, testing
withdrawal/storage/pressure/mixing assumptions with the development dates,
and evaluating agreed acceptance tolerances against independent data.
Record velocity, mixing and heat-transfer-by-depth diagnostics to distinguish
wind-driven transport from pressure and conveyance assumptions before tuning.
Hypsographic particle deposition includes projected-bank interception in the
aggregate plan-area sediment inventory. Scoped water rates use the actual
surface or lowest bottom interface; reactions along intermediate sloping
banks and separate bank sediment compartments are not represented. Independent
two-box tests qualify these balances, rather than inferring them from the
temperature-only DeGray case.
Actual-solver JVP/VJP validation and CPU installation qualification follow;
finite-difference scenarios alone do not advertise a differentiated solver.

## Spring heat-transport audit

The registered 1980 day 64.5–100.5 experiments under
`tests/output/degray_thermal_audit/` diagnose the early deep warming. They
preserve the corrected full-season evidence under `degray_phase5/`.
All four experiments finish with finite output, no barotropic floor binding,
and independently closed cell/interface heat ledgers. Maximum relative cell
closure error is 2.28×10⁻¹¹; interface divergence closure is below 5.26×10⁻¹¹.
The normal-wind audit matches all 34,484 previously saved spring rows exactly.

| Spring case | Dam bottom at day 100.5 (°C) | Purpose |
| --- | ---: | --- |
| Normal wind | 11.87184 | Existing model with observation-only audit |
| Zero wind stress | 6.56146 | Remove momentum drag; retain weather heat exchange |
| Reduced shear mixing | 11.24929 | Suppress neutral shear viscosity/diffusivity; retain convection and background momentum viscosity |
| Closed basin with wind | 11.77247 | Remove prescribed inflow/withdrawal; retain wind and weather heat exchange |

The normal-wind bottom cell changes from 6.5 to 11.87184 °C. Its temperature
change is accounted for by advection plus sources (+7.06585 °C), diffusion
(−1.72463 °C), and changing storage (+0.03063 °C); direct bottom solar heating
is negligible. The closed basin warms nearly as much with zero prescribed
source loads. These results point toward the circulation/advective transport
path for the next investigation. They do not distinguish resolved circulation
from numerical advection mixing, or prove a pressure/conveyance defect.
The reduced-mixing switch also changes momentum shear viscosity, so it is not
an isolated scalar-diffusivity intervention. Zero wind stress is not a proposed
calibration. Spring errors use only the 96 supported observations in this
short interval and cannot replace the full-season comparison.

`degray_heat_transport.csv` holds cumulative cell heat contributions in J and
upward advective/diffusive interface amounts. Positive cell amounts add heat;
positive interface amounts move heat upward. All use the consistent reference
`rho_ref * Cp * temperature_C * volume`. Cell identities follow the moving sigma
grid, not fixed physical-depth bands. Temperature attribution explicitly
includes changing volume. Large opposing horizontal, vertical and withdrawal
amounts cancel; interpret their sum rather than an individual amount as a
temperature prediction. The saved `u_m_s`, `v_m_s`, `upward_relative_top_m_s`,
`kz_*` and `baroclinic_*` fields are last-step snapshots, not daily means.
The upward velocity is relative to the moving grid and uses the actual wet
interface area, not the plan-polygon normalization of the internal flux array.

Reproduce the diagnostic runs using the same input/observation paths as above:

```sh
python3 verification/degray/seasonal_validation.py \
  --driver build/apps/fvqual-degray --inputs /path/to/DeGray/InputFiles \
  --observed /path/to/DeGray/prf_dam.npt \
  --output tests/output/new_spring_audit \
  --cases baseline no_wind_stress reduced_shear_mixing \
  --end-day 100.5 --heat-audit --jobs 3
python3 verification/degray/seasonal_validation.py \
  --driver build/apps/fvqual-degray --inputs /path/to/DeGray/InputFiles \
  --observed /path/to/DeGray/prf_dam.npt \
  --output tests/output/new_closed_audit --cases no_sources \
  --end-day 100.5 --heat-audit --jobs 1
python3 verification/degray/heat_transport_review.py \
  tests/output/new_spring_audit --closed-basin tests/output/new_closed_audit
```

Plotting requires Matplotlib; execution and audit qualification use the standard
library. Choose fresh output directories. `registration.json` fingerprints the
model, inputs and protocol; each case saves its command, numerical qualification,
and CSV fingerprints. The executed driver/core/source archive is under
`degray_thermal_audit/executed/`; its replay instructions require the external
SDK/runtime and are not a portable installation. The archived driver used the
legacy `turbulent_mixing` flag for the reduction switch. The current driver
names this accurately as `shear_mixing_reduced` and explicitly records retained
convection; this metadata correction leaves numerical output unchanged.

The review figure is `spring/review/spring_diagnosis.png`. Per-case
`dam_bottom_attribution.csv` and `heat_audit_review.json` independently check the
attribution and interface divergence. `spring/review/status.json` summarizes all
four cases. Composer manifests contain the 13 new layered diagnostics and
physical dimensions for scalar/field temperature, volume, heat and momentum.

The rebuilt GUI opens the current Bellinge review composition and all four
spring cases. To evaluate it:

1. In **Runs**, expand a spring case → **degray**. Select
   **bottom_temperature_C** → **Plot series**. Repeat for each case and compare
   with the four-case review figure. The scalar plot is shown per selected run.
2. Select **temperature** → **Show on Map** → **Map** → **Full Extent**. In
   **Layer Properties**, choose **Symbology → Graduated**, then
   **Rendering → Vertical layer**. Layer 1 is the surface and layer 30 the bed.
   Use the first/last-time buttons to inspect March 4 and April 9, 1980 noon.
3. Scroll through the run's results to **u_m_s**, **kz_top_m2_s**,
   **baroclinic_x_m_s2**, **vertical_advection_J** and
   **vertical_diffusion_J**. Show one on the map, color by its value and inspect
   the same layer/time. Hide previous layers when switching fields.
4. Use **Section Line** and drag/release across the mapped chain to inspect a
   full saved column section; use **3D** to inspect layer geometry. The W2
   footprint drawing does not depict surveyed bathymetry or hypsographic storage.
5. Select the Bellinge run to inspect SWMM exported-water/heat and receiver
   budgets. Hide DeGray layers first because the coordinates and calendar differ.
   The receiver remains synthetic; the coupling check does not qualify a physical
   Bellinge receiving water body.

The phase 5c investigation below implements the next scientific gate described
at this earlier checkpoint:
closed tapered-basin benchmarks with independently specified storage and face
conveyance; momentum and pressure-rest balance checks; transport refinement and
an intervention that changes scalar mixing without changing viscosity. Diagnose
the effect before selecting any physics correction. Rerun a newly registered
seasonal matrix after a justified correction. Actual-solver differentiation,
uncertainty qualification and ABI 4 CPU installation remain subsequent gates.

## Phase 5c: circulation, momentum and transport geometry

The October 13 review is under `tests/output/phase5c/`. The isolated analytical
tests identified and corrected thickness-based vertical viscosity, surface wind
loading and bottom stress on tapered storage. The corrected operator uses
interface contact area and actual cell volume. Quadratic bed drag freezes the
same vector speed for both velocity components. A steep-step test also identified
pressure extrapolation beneath a neighbouring bed; that portion of a face is
now a solid wall, and one-sided boundary density reconstruction preserves a
horizontally uniform linear profile.

`tests/test_momentumgeometry.cpp` supplies the independent two-box, stress,
rotation, replay, coefficient-isolation and pressure-rest gates. Stable
Richardson diffusivity and viscosity have separate diagnostic multipliers.
Unstable convective adjustment is retained; non-default multipliers are refused
with the k-epsilon closure. The coefficients are independent at the same state;
subsequent feedback through changed velocity, stratification and heat exchange
is part of the experiment. These interventions are not calibrated models.

The optional per-column momentum audit closes the implicit vertical operator's
impulse and backward-Euler kinetic-energy balance, including viscous, bed and
implicit numerical dissipation. Mode-reconciliation impulse/work is measured
separately. Relative residuals use the old/new absolute inventory and the step
forces as their numerical scale; absolute residuals remain available. This is
not full mechanical-energy closure for pressure, horizontal transport, moving
geometry or the external subcycle. The existing external-mode subgrid momentum
approximation remains a scientific limitation when storage and conveying shape
differ. A small heat residual does not establish accurate circulation.

The reusable `fvqual-circulation-review` driver specifies tapered storage and
rectangular face-conveying width independently. Fourteen source-free heated
basin cases compare time steps, 10/20/40/80 layers, 10/20/40 columns, MUSCL-FCT
and upwind, separate mixing interventions, wind removal and conveying width.
Initial linear temperature uses an independently integrated storage moment.
Refinement is sampled at common physical x/depth positions within centre
support, not at equal layer indices. The W2 inferred two-lane experiment remains
a separate geometry sensitivity, not a formal horizontal-convergence sequence.

`circulation/qualified_results.json` and `circulation/qualified_review/` contain
the frozen-runtime diagnostics, figures, measured refinement ratios and native
recorded data. `qualified_seasonal/` contains the newly registered eight-case
March–December matrix. `original_seasonal/` retains the original solver baseline.
The earlier `seasonal/` and `corrected_seasonal/` folders preserve prototype
experiments; incomplete runs are marked superseded and excluded from acceptance.
The first prototype's near-zero incremental impulse scale is superseded by the
inventory-aware audit. It is not a temperature or velocity correction.

To reproduce in a fresh output folder on the current macOS development setup:

```sh
cmake --build build --parallel
python3 verification/degray/freeze_runtime.py build tests/output/my_phase5c/executed
tests/output/my_phase5c/executed/verification/fvqual-circulation-review tests/output/my_phase5c/basins.json
python3 verification/circulation_review.py tests/output/my_phase5c/basins.json tests/output/my_phase5c/basin_review --figures
python3 verification/degray/seasonal_validation.py --driver tests/output/my_phase5c/executed/apps/fvqual-degray --inputs <DeGray-InputFiles> --observed <prf_dam.npt> --output tests/output/my_phase5c/seasonal --cases baseline dt150 dt75 layers15 layers60 scalar_mixing_off viscosity_off no_wind_stress --jobs 4 --heat-audit
python3 verification/degray/review_artifacts.py tests/output/my_phase5c/seasonal --split-fields --caption-prefix "DeGray phase5c"
```

The converters need numpy/matplotlib only for figures. The local runtime freezer
copies FVQual drivers, core and sources and uses relative core lookup. Its SDK
runtime remains external, so this archive is not a release package. Per-variable
annual CSV artifacts keep field access from scanning unrelated diagnostics.
Every recorded run retains source/converter hashes and the original model clock.

### Evaluate in the rebuilt Composer

1. Open the delivered Bellinge phase5c review. In **Map**, use **Select** and
   click or box-select receiving cells. Right-click to plot the selected objects
   or a recorded variable. **Add from Map** retains existing curves. Inspect the
   nine outfall mappings and the recorded receiver budgets; the receiving mesh
   and fixed 18 °C inflow remain synthetic assumptions.
2. In **Results**, select a corrected DeGray run and its temperature field, then
   use **Show on Map** or **Show in 3D**. Hide Bellinge/basin layers when changing
   coordinates and calendars. Compare `scalar_mixing_off`, `viscosity_off` and
   `no_wind_stress` at matching times; use the profile view at the dam to compare
   physical depths. Different layer counts must not be compared by layer index.
3. Inspect `u_m_s`, `km_top_m2_s`, `kz_top_m2_s`, pressure and cumulative heat
   transport separately. Momentum work/residual scalars appear in the run tree.
   Compare the mode-correction work with vertical dissipation without treating
   the former as measured physical dissipation.
4. Open `circulation/qualified_review/basin_review.composition.json` for the
   controlled cases. Check wind/viscosity effects and profiles against
   `circulation_review.png` and `refinement.png`. Wheel zoom, middle-button pan,
   rectangle zoom, playback and CSV/PNG/PDF export should remain available.
5. Review annual profiles, seasonal/depth errors and common-support refinement
   before interpreting sensitivities. The retrospective autumn reservation is
   not independent validation. Scientific error tolerances and surveyed receiving
   geometry remain open; parameter bands must not conceal structural warm bias.

Actual-solver JVP/VJP is the subsequent implementation phase. Its verification
must cover directional differences, adjoint dot products and replay before
using sensitivities for uncertainty or calibration. Geometry, forcing,
parameters, numerical error and structural discrepancy remain distinct sources.

### October 13 qualified review results

All eight final seasonal experiments completed and passed the numerical gates.
On 991 common supported observations, the original solver has temperature RMSE
7.45 °C and the corrected solver 8.04 °C. Stable scalar mixing off gives 5.29 °C,
stable viscosity off 11.03 °C, and wind stress off 2.82 °C. These results identify
important wind/mixing sensitivities; they do not select calibrated settings or
resolve the structural warm bias. `review/diagnosis.json`,
`review/before_after_profiles.png` and `review/temperature_errors.png` retain
the direct original/corrected comparison. Reproduce the comparison, with both
registered baseline folders present, using
`python3 verification/degray/phase5c_review.py tests/output/phase5c`.

The complete FVQual run has 386 passing registered tests and one expected skip
for the disabled non-NetCDF branch. DeGray and the full Bellinge mesh are enabled;
eight observation/conversion checks pass. All 266,400 exported annual baseline
temperatures exactly equal the original solver CSV, including clock/layer axes.
Native graphical selection/plot/export checks pass 39 controlled-basin, 888
DeGray and 60 Bellinge selected-object values. The final nine-outfall Bellinge
live/replay check matches all 1,800 receiving-field comparison values exactly.

The delivered composition is
`tests/output/phase5c/qualified_swmm_coupling/bellinge/Bellinge_phase5c_review.composition.json`.
It opens three Bellinge curves and a 24-source recorded catalog. Its subsequent
Run writes separate review outputs. The local development app is
`HydroCoupleComposer/build/darwin-abi4/phase5c-review/gui/HydroCoupleComposer Phase 5c Review.app`.
This is a development review launch; full installation and assistive-technology
acceptance remain later gates.

## Phase 6a passive scalar sensitivity review

The first actual-solver derivative increment is qualified for passive upwind
advection and implicit vertical diffusion. JVP and VJP include initial cell
concentrations, gross inflow/withdrawal terms and a multiplier of the recorded
diffusion operator over multiple internal steps. The reverse solve uses the
true transpose on unequal-volume columns. Fixed observation/unit mappings
include repeated terms and their transpose. C/Python owned trace handles share
immutable records with solver checkpoints and can outlive their simulation.
See `python/README.md` for the callable surface and refusal matrix.

Hydraulics, geometry, source placement, Kz trajectory and time schedule are
frozen. This is not a thermal/hydrodynamic adjoint, SDK differentiable component
or SWMM parameter gradient. Temperature/salinity/age, MUSCL/FCT, horizontal
diffusion, reactions, settling, internal transfers and open boundaries are
refused. Logarithmic internal replay and full memory/performance qualification
remain open. Snapshot/restore is in memory, not a coordinated durable restart.

The October 13 Bellinge experiment keeps the registered nine-outfall water and
heat forcing, the synthetic three-column triangle/quad receiver and three
layers. It adds a **synthetic** conservative tracer: initial 0.0002 kg/m3 and
incoming 0.001 kg/m3. These are demonstration assumptions, not observed water
quality. Uncertainty is a Gaussian initial concentration (SD 0.00002 kg/m3) and
load multiplier (mean 1, SD 0.2), correlation -0.35, seed 60213. Negative
physical draws would fail and be counted; none of the 512 runs failed.

The local covariance and whole-run ensemble omit observation error,
structural discrepancy, geometry and numerical uncertainty. They are parameter
envelopes conditional on the recorded forcing. Do not interpret them as
validated predictive bounds or a calibrated posterior. The DeGray discrepancy
above remains unaccepted and unchanged.

Review in the launched **HydroCoupleComposer Sensitivity Review** app:

1. Results opens with baseline and lower/upper Gaussian parameter envelopes
   pinned for face 0, layer 1. Compare these with `sensitivity_review.png`, which
   also shows the actual whole-run ensemble percentiles. The GUI uses kg/m3;
   the scientific figure uses mg/L (multiply kg/m3 by 1000).
2. In Map, select a receiver face, Shift-select another, or draw a selection
   box. Right-click **Plot selected objects** to replace the pinned plot. Use
   **Add from Map** to add objects. Wheel zoom, middle-button pan and interactive
   Zoom tools retain their existing behavior.
3. Expand the recorded run and inspect `initial_sensitivity`, `load_sensitivity`
   and `parameter_sd`. Show them on Map/3D; change the vertical layer (0 is
   surface). Sensitivity fields show concentration response to one marginal
   parameter SD, before correlation. Select the variable before plotting it.
4. Use Profile/Section and animation to inspect the rising tracer after hour
   five. The flat receiver bed is synthetic. Keep water/heat/model uncertainty
   separate from the two tracer parameters shown here.
5. Export the data table and figure, save the composition, and reopen it. The
   native acceptance has already verified all 63 selected-object samples
   against the recorded field, their layer/time axes, plot zoom, playback and
   CSV/PNG/PDF exports.

Primary evidence under `tests/output/phase6a/`:

- `qualified_bellinge/qualification.json`, `ensemble_runs.json`, `run.json`,
  `sensitivity_review.png` and `.pdf`; all 512 runs are retained.
- `qualified_bellinge/Bellinge_sensitivity_review.composition.json` and sidecar;
  the launched copy is under Composer `build/darwin-abi4/phase6a-review/`.
- `full_tests.log`: 392 registered tests, 391 passing and one expected skip for
  the non-NetCDF failure branch in a NetCDF-enabled build. Actual DeGray and
  Bellinge mesh fixtures and the live SWMM Python bridge were enabled.
- `diffusion_qualified.json`: independently initialized cosine-column
  sensitivities converge toward the analytic no-flow diffusion reference with
  time/layer refinement. A 16-layer time sweep reduces RMS gradient error from
  0.0022043 (dt 10 s) to 0.0013029 (5 s) to 0.0008375 (2.5 s). Spatial and time
  errors coexist; the measurements do not establish reservoir field accuracy.
- At uniform diffusion uncertainty radii 0.05, 0.5 and 0.9, 512 actual runs each
  show maximum nonlinear departures 0.0000541, 0.0060907 and 0.0220446 kg/m3
  from the local tangent prediction. Local SD is a covariance approximation;
  these uniform distributions are not Gaussian 95% intervals.
- `temperature_feedback_diagnostic_all.json`: the 32-layer nominally uniform
  transported-temperature case changes Richardson/convective mixing branches
  under perturbation (Kz departure 0.090001 m2/s). Its whole-run derivative
  differs strongly from a frozen recorded-Kz derivative. This is a failed
  thermal-derivative regime diagnostic, not accepted thermal sensitivity.
- Bellinge whole-run directional differences use five epsilon values. Maximum
  error is 3.03e-13 kg/m3; the mapped adjoint dot residual is 4.07e-20. Forward
  replay and full-run ensemble/tangent agreement are within 6.51e-19 kg/m3.
  Tracer mass residuals are below 5.69e-14 kg. Water/temperature values are
  identical across tracer ensemble members. Ninety internal steps use an
  estimated 265740 bytes of trace payload. This is not peak process memory.
- Sample SD differs from local Gaussian SD by at most 13.56%, 8.85% and 3.81%
  for the common 32-, 64- and 512-member prefixes. Finite sample differences
  are a convergence diagnostic, not evidence of nonlinear tracer behavior.
- Composer `build/darwin-abi4/phase6a-native/` holds the native acceptance JSON,
  screenshots and exports. The local app still uses external SDK/GIS runtimes;
  it is not a phase 7 release installation.
- `qualified_bellinge/native_envelope_parity.json` verifies 63 reopened
  baseline/envelope values against the saved fields (CSV rounding below
  4.78e-16 kg/m3). Derived lower and upper views use the same physical `tracer`
  variable and units, preserving the plot compatibility guard. They are
  marginal parameter envelopes, not independently simulated trajectories.
  Versioned plot recipes and their generator hash are in `gui_metadata.json`;
  final native screenshots are under Composer
  `build/darwin-abi4/phase6a-bounds-acceptance/native/`.

`executed/` preserves the actual Bellinge executable/library/Python source
archive. The original diffusion fixture in it failed its thermal-feedback
finite-difference gate and was superseded by the explicit homogeneous-density
reference. `diffusion_executed/` preserves the corrected diffusion verifier and
the thermal diagnostic. `qualified_executed/` is an intermediate preserved
archive. No earlier source/evidence was overwritten. Archive SHA-256 hashes
identify the executed source and binaries; external SDK dependencies remain.
The GUI generator received a subsequent metadata correction; its executed
copy is preserved in `gui_generator/scalar_uncertainty_review.py`. Numerical
fields and the original run identity were retained when repairing the views.

From the FVQual root, reproduce into fresh folders:

```bash
PYTHONPATH=tests/output/phase6a/executed/sources/python python3 \
  tests/output/phase6a/gui_generator/scalar_uncertainty_review.py \
  tests/output/phase5c/qualified_swmm_coupling/bellinge \
  tests/output/phase6a/reproduced_bellinge \
  --library tests/output/phase6a/executed/src/libfvqual-core.dylib --samples 512

tests/output/phase6a/diffusion_executed/verification/fvqual-scalar-derivative-review \
  tests/output/phase6a/reproduced_diffusion.json
tests/output/phase6a/diffusion_executed/verification/fvqual-scalar-derivative-review \
  tests/output/phase6a/reproduced_temperature_diagnostic.json \
  --temperature-feedback-diagnostic
```

Use `--figures` with matplotlib/numpy to generate figures. Do not regenerate a
delivered manifest in place: run identities are part of saved plot provenance.

Next is phase 6b: qualify actual heat/momentum/mixing and reaction derivatives,
including initialization, thermal feedback and branch regimes. Then complete
geometry/withdrawal/bank and observation-weight derivatives, replay policy and
SDK integration gates before claiming phase 6 complete. Separate experiments
must retain DeGray structural discrepancy and compare local approximations
with full-run ensembles. SWMM parameter, feedback and advanced domain gradients
remain later gates; CPU packaging/platform installation is phase 7.

## Combined thermal geometry and replay review — October 15

The feasible portions of phases 6b, 6c and 6d are combined into one qualified
review package. A controlled native thermal column has analytic discrete
JVP/VJP for initial temperature/shear and eighteen meteorology, mixing,
optical, heat-capacity and geometry-scale parameters. Temperature, freshwater
density, Richardson mixing, wind shear and surface heat feedback are active.
Depth/area scaling and changing fixed-depth interpolation weights are included.
Blocked reverse replay uses one local step tape at a time. The same package
compares 256 full column runs and 256 full Bellinge receiving-model runs.

The column is one closed always-wet column, with upwind transport, Richardson
mixing and zero Manning/Coriolis. Its velocity outputs are departures from the
depth mean. It excludes sources, horizontal circulation, salinity, kinetics and
settling. Near-uniform-density/mixing/clamp boundaries refuse local gradients;
forward ensemble evaluation remains available. The full Bellinge thermal model
has explicitly named finite-difference sensitivity references and full-run
ensemble envelopes. It does **not** acquire an analytic horizontal-reservoir or
SWMM adjoint from the controlled-column implementation. The SDK differentiable
component interface remains unadvertised.

The rebuilt **HydroCoupleComposer Thermal Review** app opens the Bellinge
baseline and two ensemble envelopes, with an eight-source catalog including
the controlled column and the original/corrected DeGray reference runs. The
live review copy is under Composer `build/darwin-abi4/combined-thermal-review/`.
The saved primary composition is
`tests/output/phase6_combined/final/Thermal_review.composition.json`.

Evaluate it as follows:

1. In Results, compare the three pinned Bellinge curves and the data table.
   The bands are marginal sample 2.5/97.5 percentiles from stated parameter
   assumptions. They exclude observation/structural discrepancy and are not
   calibrated predictive bounds. The receiving mesh and inlet temperature are
   synthetic; the water forcing is the preserved Bellinge recording.
2. In Map, select faces or a selection box and right-click **Plot selected
   objects**. Use **Add from Map** for more faces, change the vertical layer
   (0 is surface), and check wheel zoom, middle-button pan and animation.
   Inspect `inlet_temperature_response`, `surface_heat_response` and
   `wind_stress_response` in the Bellinge run. These are scaled finite-difference
   responses, not analytic Bellinge gradients.
3. Select the controlled-column run and show its `temperature`, `wind_response`,
   `depth_response`, `parameter_sd` or shear fields. Hide the Bellinge map layer
   when switching to the column's local coordinates, then fit the view.
   Inspect its Profile or 3D layers. The footprint and prisms are synthetic;
   their constant displayed footprint does not depict the tapered storage
   cross-section. The standalone `Column_review.composition.json` also opens
   its selected-column series.
4. Compare the column's local covariance approximation with its native-run
   sample envelope in `thermal_review.png`/`.pdf`. The covariance-factor modes
   include initial sigma-profile offset SD .3 C, wind SD .3 m/s, air SD 1 C,
   wind/air correlation .35, depth-scale SD .025 and optical attenuation SD
   .03/m. Geometry sensitivities at fixed physical instrument depths include
   changing interpolation weights. The instrument observations here are
   illustrative model samples, not measured DeGray casts.
5. Export/save/reopen the review. Native acceptance verifies restored envelope
   values, recorded object curves and layers, zoom/playback, duplicate-safe
   Add from Map, and CSV/PNG/PDF exports for both the mixed Bellinge footprint
   and the single column. Do not overlay these synthetic experiments with the
   DeGray calendar/coordinates as if they were a common field realization.

The Bellinge ensemble uses initial 10 C, a synthetic inlet-temperature offset
SD 1 C, prescribed surface heat mean 100 W/m2 with SD 20, and signed x wind
stress mean .01 Pa with SD .002. These three illustrative parameters are
independent Gaussian inputs. The recorded water discharge/mapping, synthetic
receiver geometry and solver choices are held fixed. Every run is retained;
neither ensemble has failures. DeGray's previously measured warm bias remains
unaccepted; this package neither fits it nor changes its acceptance status.

Evidence under `tests/output/phase6_combined/`:

- `final/qualification.json`, `column_qualification.json` and
  `bellinge_qualification.json`; the individual files are in `final/`.
- `final/column_ensemble.json` and `bellinge_ensemble.json` retain all 512 runs.
- Column production/operator agreement is within 7.11e-15; a four-epsilon
  native directional sweep, including physical-depth sampling, has maximum
  error 5.03e-8 C. The mapped adjoint residual is about 1.09e-15.
- The 360-step example stores 361 states with every-step checkpoints and 20
  with blocked checkpoints. Parameter gradients agree within 4.55e-13. A local
  tape peaks at 981 nodes; the largest replay block holds 20 states. These are
  payload counts, not process peak memory, logarithmic checkpointing or a full
  memory/performance qualification.
- Bellinge's maximum heat-ledger residual is below .006 J. Its full-run
  finite-difference sweeps are preserved separately from the analytic column
  gates. `forward_parity.json` records exact before/after snapshots for a
  separate four-column native heat/wind case after production-kernel extraction.
- `final_full_tests.log`: 396 registered FVQual checks; 395 pass and one
  expected non-NetCDF failure-branch skip. Actual DeGray and regenerated
  Bellinge fixtures and the current SWMM live bridge were enabled.
- The regenerated Bellinge mesh now has 13861 vertices and 24295 triangular
  faces. Qualification independently reads its raw records, verifies all
  imported terrain coordinates and boundaries, and checks conservative wet
  execution with the explicit synthetic-bed override. Historical hard-coded
  counts were removed; prior review packages remain intact. This does not
  reinterpret terrain as surveyed receiving-reservoir bathymetry.
- All 70 Composer suites pass across the ordinary and required native
  graphics/local-network contexts. Logs are copied alongside this package.
  Native review artifacts are under Composer
  `build/darwin-abi4/combined-thermal-acceptance/native/` and
  `combined-column-native/`.
- `final/native_envelope_parity.json` verifies 63 Bellinge and 1083 column
  reopened baseline/envelope values, within 5e-11 C of saved data after CSV
  export rounding. Object selection checks independently verify 63 Bellinge
  and 361 selected-column samples.
- `qualified_runtime/` freezes the executed library, Python/driver sources
  and numerical source (197 hashed files). `executed/`, `development/` and
  `qualified/` preserve earlier development evidence. External SDK/GIS runtime
  dependencies still apply; this is not a clean release installation.

Reproduce from the FVQual root into a fresh directory:

```bash
PYTHONPATH=tests/output/phase6_combined/qualified_runtime/sources/python python3 \
  tests/output/phase6_combined/qualified_runtime/sources/verification/thermal_uncertainty_review.py \
  tests/output/phase5c/qualified_swmm_coupling/bellinge \
  tests/output/phase6_combined/reproduced \
  --library tests/output/phase6_combined/qualified_runtime/src/libfvqual-core.dylib \
  --samples 256
```

Add `--figures` with matplotlib/numpy. C/Python usage is documented in
`python/README.md`; `fvqual_column.h` states the exact restricted contract.

The remaining compatible derivative work is consolidated into one integration
package: general horizontal-reservoir and kinetics derivatives, remaining
storage/withdrawal/bank/event geometry, scalable general replay and SDK
differentiable-component integration, with conditional uncertainty and DeGray
discrepancy kept explicit. After that is the combined CPU packaging,
clean-install/platform, accessibility and release-evidence package (phase 7).
Feedback, coordinated durable restart, SWMM parameter gradients, wet/dry and
GPU work retain separate advanced gates. Plans stay local and uncommitted.

## Native reactive transport and SDK review — October 15

The next integration increment connects actual mixed-mesh upwind transport and
implicit vertical diffusion to nonlinear multispecies kinetics and the existing
HydroCouple differentiation contract. Initial pollutant/MSX concentrations and
reaction constants have analytic JVP/VJP. The cell-local reverse tape follows
the native adaptive Cash-Karp controller, including changes to accepted
internal step lengths. Expression domains and adaptive switching/scaling or
max-error ties refuse local gradients. Forward ensemble evaluation is separate
and remains available when a gradient refuses.

`FVQualComponent` implements `IDifferentiableModelComponent` under the explicit
`conditional_species_derivatives = 1` opt-in, which selects upwind transport.
Default components retain their transport scheme and capability declaration.
Supported species outputs also describe the state carried across updates;
`reaction_constants` exposes named constants in reaction-deck units. Inputs
not listed as differentiable are prescribed. Native SDK every-step and
logarithmic reverse workflows and cross-plugin RTTI/buffer tests exercise the
same contract. The SDK's existing checkpoint scheduler is reused.

This is **conditional species differentiation**. Hydraulics, temperature and
salinity, geometry, mixing, source discharge/composition and placement are
prescribed. Density-changing rate laws, self-shading, settling, horizontal
diffusion, internal transfers and open boundaries refuse the contract. The
controlled-column thermal derivatives from the preceding section are separate;
the combined work does not yet provide a horizontal thermal/momentum adjoint,
general geometry/withdrawal/bank derivatives or SWMM parameter gradients.

The Bellinge review uses the preserved phase5c water/heat recording, nine
outfall mappings and the synthetic three-face tri/quad receiving basin. It
does not use surveyed receiving-basin bathymetry or re-run the current full-size
Bellinge terrain mesh. Its synthetic A-to-B reaction conserves combined mass:
`dA/dt = -k*A*(1+B/cref)`, `dB/dt = -dA/dt`. Initial A/B are .0002/.00005
kg/m3, `k = .00002 s-1`, `cref = .001 kg/m3`; the synthetic inlet is .001
kg/m3 A and zero B. The 256 native ensemble runs all succeed. Gaussian SDs are
.00002 kg/m3 initial A, .00001 initial B and .000003 s-1 for k, with initial
A/k correlation -.35. Negative physical draws fail and are counted; none are
silently resampled. All sample parameters and outputs are retained.

The rebuilt **HydroCoupleComposer Reactive Review** app opens
`tests/output/phase6_reactive/final/Bellinge_reactive_review.composition.json`.
Its live review copy is in Composer `build/darwin-abi4/reactive-review/`.

Evaluate it as follows:

1. In Results, compare A's baseline and sample 2.5/97.5-percentile curves.
   The plotted values are kg/m3; `reactive_review.png`/`.pdf` uses mg/L. These
   are marginal parameter envelopes, with observation and structural error
   excluded. The three saved views retain A's physical identity and units.
2. In Map, select faces/box/Shift selection and right-click **Plot selected
   objects**. Use **Add from Map**, layer selection, wheel zoom, middle-button
   pan and animation. Layer 0 is surface. Native acceptance checks the selected
   series, zoom/playback, duplicate-safe additions, save/reopen and exports.
3. In the baseline run's catalog, inspect B, `kinetic_response_A`,
   `initial_response_A`, `parameter_sd_A` and the corresponding B fields.
   Responses describe changes for one stated parameter SD. They have
   concentration units; parameter SDs and their correlations are explicit in
   `qualification.json`. Show the fields in Map, Profile or 3D.
4. Compare the native ensemble with the local Gaussian approximation in
   `reactive_review.png`. `observation_design.json` supplies sensitivity SVD,
   confounding diagnostics and an illustrative Gaussian update from seven
   **synthetic** noisy observations. Its three prior-factor directions have
   numerical rank 3 and condition number about 4.43. This is a local survey
   diagnostic, not a field identifiability or calibration acceptance.
5. Inspect the retained withheld result: baseline RMSE .005371 mg/L becomes
   .006847 mg/L after the noisy update. The new posterior is evaluated through
   fresh native runs; the worse withheld error is not discarded. This example
   demonstrates why fitting observations and shrinking a local covariance do
   not establish predictive improvement. DeGray's warm bias remains unresolved.

The portable `examples/reactive_transport/composition.json` runs through the
ordinary CLI/GUI authoring path, with relative `@uri` reaction-deck resolution,
named constant overrides and process budgets. It has a separate synthetic
rectangular basin and inlet. Its A+B inventory closes independently after
save/reopen. Forward execution through a standard GUI workflow does not
automatically perform a reverse pass; the SDK review driver exercises that.

Evidence under `tests/output/phase6_reactive/`:

- `final/qualification.json`, `sdk_qualification.json`, `ensemble.json`,
  `observation_design.json`, `review_fields.json` and the GUI compositions.
- Four independent whole-native directional checks have maximum concentration
  error 2.10e-12 kg/m3; the mapped adjoint residual is below 9e-19. SDK and C ABI
  gradient routing agree exactly in this experiment. No conditional parameter
  perturbation changes water volume or temperature.
- SDK every-step/logarithmic storage is 20/6 composition checkpoints, with
  reverse peaks 21/10 and 19/40 replays. Gradients and native forward states
  agree. Latest-interval operator history is bounded separately; these counts
  are not a process peak-memory or performance qualification.
- `final_full_tests.log`: 405 registered FVQual checks, 404 pass and one
  expected non-NetCDF failure-branch skip. Actual DeGray, current regenerated
  Bellinge mesh and the live SWMM bridge are enabled. `optional_python_tests.log`
  separately qualifies the four new Python tests, including NumPy diagnostics.
- `gui_tests.log` records nine affected Composer suites passing. Native
  acceptance is in Composer `build/darwin-abi4/reactive-review/acceptance/native/`:
  three graphically selected faces and 63 samples, context menu, zoom/playback,
  reopening and CSV/PNG/PDF exports. `final/native_envelope_parity.json`
  independently verifies 63 reopened envelope values to within 5e-16 kg/m3
  after CSV rounding.
- `qualified_runtime/` freezes 212 hashed files: executed drivers/core,
  Python and numerical source. It requires the external SDK/GIS runtime and
  is an experiment archive, not a portable release installation.

Reproduce into a new directory from the FVQual root:

```sh
PYTHONPATH=tests/output/phase6_reactive/qualified_runtime/sources/python python3 \
  tests/output/phase6_reactive/qualified_runtime/sources/verification/reactive_uncertainty_review.py \
  tests/output/phase5c/qualified_swmm_coupling/bellinge \
  tests/output/phase6_reactive/reproduced \
  --library tests/output/phase6_reactive/qualified_runtime/src/libfvqual-core.dylib \
  --sdk-driver tests/output/phase6_reactive/qualified_runtime/verification/fvqual-reactive-sdk-review \
  --samples 256
```

Add `--figures` in an environment with NumPy/matplotlib to include the synthetic
observation update, SVD and figures. The library/driver/source and donor hashes
are retained. General reservoir thermal/momentum, settling, geometry and
coupled derivatives remain open, followed by clean CPU installation/platform
and accessibility/release qualification. Plans remain local and uncommitted.
