# DeGray Reservoir, March–December 1980

An end-to-end run of CE-QUAL-W2's DeGray example through FVQual: real
bathymetry, real meteorology, real inflows and outflows, real initial
conditions.

    fvqual-degray <w2-input-directory> <output-directory>

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
| Peak surface temperature | 31.4 °C, day 198 | ~30 °C in midsummer |
| Hypolimnion through summer | 6.5 → 7.4 °C | retains winter temperature |
| Maximum stratification | 24.3 °C | ~22 °C |
| Thermocline, midsummer | 5–7 m | 5–15 m |
| Thermocline, December | 35 m | turnover complete |
| Final profile | isothermal ~9.9 °C | isothermal ~9–10 °C |
| Volume change | +1.4% | small, inflow/outflow balance |

The reservoir starts nearly isothermal, stratifies through spring, holds a cold
hypolimnion all summer while the surface reaches 31 °C, then deepens its
thermocline through autumn and turns over by December. Surface temperature
tracks below the instantaneous equilibrium temperature throughout, which is what
thermal inertia does.

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
