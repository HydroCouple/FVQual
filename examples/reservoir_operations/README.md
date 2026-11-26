# Scheduled reservoir operations

This synthetic 120-second example uses two triangles and one quad, three sigma
layers, and a −5 m bed. It exercises weather, a tributary, a pump/release group,
an exterior inflow edge, circulation, reactive decay, and particle settling
into an areal sediment store. It is a software qualification case.

From the FVQual repository:

```sh
build/apps/fvqual import-mesh examples/swmm_reservoir/reservoir.2dm \
  examples/swmm_reservoir/reservoir.mesh.json \
  --crs LOCAL --datum synthetic-bed --elevation-kind bed
build/apps/fvqual run examples/reservoir_operations/composition.json
python3 verification/reservoir_operations.py \
  --composer ../HydroCoupleComposer/build/darwin-abi4/HydroCoupleComposer.app/Contents/MacOS/HydroCoupleComposer
```

The qualification script independently integrates source, release, boundary
and solar inputs at the irregular 33.5 and 40.5 second control knots, checks
settling exchange and process residuals, and compares every recorded value
between the CLI and Composer. Outputs remain under ignored `tests/output`.

## Authoring in Composer

Open `composition.json`, select the reservoir, choose an argument in the
**Argument** dropdown, and choose **Basic** mode. **Edit…** opens contextual
controls for `meteorology`, `sources`, `structures`, `outlet_groups`,
`open_boundaries`, `internal_transfers`, and `reaction_deck`. Check a field to
include it; unchecked fields use documented solver defaults. Add/remove rows
and time records as needed. Apply validates a draft against the mesh, registered
species and other operations before changing the component. Cancel discards
unapplied edits. Advanced mode exposes the complete serialized argument.

Times are elapsed seconds from the reservoir epoch. Schedules start at zero,
increase strictly and hold values between knots; integration splits at those
knots. Weather blanks inherit previous values. Other scheduled controls need
a value at every record. An entirely blank control column uses its constant
value. Weather supports a uniform number or comma-separated per-column values.
Concentrations use species units declared in the registry/reaction deck;
columns and boundary endpoint nodes are zero-based mesh indices. Boundary
inflow discharge is per unit edge length in m²/s, sources use m³/s. Positive
source rates add water, negative rates withdraw ambient heat/species. Pump
capacity is negative for withdrawal; group demand is a positive magnitude.
Wind direction points towards an angle from +x, counterclockwise, in radians.
The latent heat term exchanges energy without removing evaporation water.
The operations weather deck cannot be combined with the legacy air-temperature
series or a wired air-temperature input. The inline reaction deck cannot be
combined with `kinetics_file`.

The component argument wrapper is `fvqual.operations/1.0` with `section` and
`value`. C API callers include `fvqual_configuration.h`, register their initial
species, and pass the corresponding top-level sections to
`fvq_sim_configure_json` before initialize. `reaction_deck` is optional text in
the same configuration and may declare further species. Query
`fvq_sim_budget_value` with the same `budget_*` names as component outputs.

## Budget and result review

Open the resulting `output/run.json` in **Runs**. Select a scalar item such as
`budget_water_residual`, `budget_heat_residual`, or `budget_tracer_residual`,
then **Plot series**; this plots and exports all recorded instants without a
map selection. Water is m³, heat is J with a 0 °C inventory reference, and
recognized mass concentration units convert to kg. Unrecognized units retain
their declared concentration × volume (or × area for benthic stores) basis.
Residuals compare the current inventory with the initial inventory plus signed
process amounts. Input/output columns are positive magnitudes; heat surface
terms and reaction/settling terms are signed. Internal transfers must have
matching input/output amounts. They do not increase the global inventory.

For `temperature`, choose **Show on Map** and select a column for a profile.
In layer properties, change **Rendering → Vertical layer (surface = 1)** to compare surface and bed
fields. The 3D view draws closed prisms for the selected layer at the recorded
interfaces. **Map → Section Line** cuts all layers using the actual recorded values;
drag across the mesh and release, then open the **Section** panel. Advance the recorded clock
to inspect changes. Results use the saved mesh/vertical geometry, independent
of a live model instance.

## Bellinge review

The current review composition is generated at
`tests/output/swmm_coupling_phase4/bellinge/review.composition.json`. Its
`inline_output/run.json` records the native six-hour Bellinge network event.
Select **openswmm → SWMM model and reservoir mapping → Basic → Edit…** to
review nine outfalls, 18 targets, conservative weights and the explicit fixed
18 °C inlet-temperature assumption. Receiver geometry is preserved; changing
the receiver mesh requires regenerating the geometry contract. Saving moves
the mapping into the composition with resolved model/output paths.

The receiving basin, location, mapping and temperature assumption are
synthetic. Bellinge is a coupling qualification, not a reservoir field
prediction. Use this standalone operations case to explore additional species
and settling. A live SWMM export must explicitly supply every registered
transported mass species; adding a species only to the receiver requires a
matching donor mapping or explicit zero-load declaration.

For the Bellinge run, **Runs → reservoir → water_volume → Plot series** should
show 50,000 to 62,162.10748 m³ across twenty records. The maximum water budget
residual is 8.19×10⁻¹⁰ m³; the heat residual is 0.490 J (1.11×10⁻¹³ of the
inventory). **Show on Map** for temperature reveals two triangles and one quad.
The final section spans a −5 m bed to an approximately 1.216 m stage. In
**Layer Properties → Symbology**, select graduated colours on **Value** to
distinguish the three temperature layers; a single symbol colours all layers
alike. The recorded 3D geometry uses the chosen vertical layer, while sections
show all three layers. The standalone operations example additionally verifies
all 3,384 CLI/Composer values exactly, including independently integrated
controls at 33.5 and 40.5 seconds.
