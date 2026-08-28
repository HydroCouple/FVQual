# DeGray Reservoir, March–December 1980

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
