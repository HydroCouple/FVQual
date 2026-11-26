# Conditional reactive transport

This portable synthetic example converts A to B in a two-face, three-layer
basin with an inlet. Concentrations use kg/m3; `k` uses s-1 and `cref` kg/m3.
The reaction preserves A+B inventory. It is not a field chemistry model.

From the FVQual root:

```sh
build/apps/fvqual run examples/reactive_transport/composition.json
```

Open `composition.json` in Composer with the rebuilt FVQual component loaded.
Review A/B fields, layer selection, inventories and reaction/source budgets.
The example produces `output/run.json` and `output/fields.csv`.

`conditional_species_derivatives = 1` explicitly selects upwind transport and
enables the HydroCouple differentiable-component contract after a supported
prepare. The default is off. `reaction_constants` overrides existing deck
constants in their declared units; the engine does not infer units from an
arbitrary expression. Unknown names and nonfinite values are refused.

The component lists A/B outputs as differentiable states/outputs, and the
kinetic constants as differentiable arguments. Hydraulics, thermal forcing,
geometry, mixing and source composition are prescribed. SDK
`DifferentiableWorkflow` carries the state cotangents through time and supports
every-step and logarithmic checkpoints. Forward execution through the ordinary
GUI workflow does not itself request a reverse pass.

See `python/README.md` for owned reactive traces and `verification/degray/README.md`
for the Bellinge derivative, ensemble, SDK replay and observation-design review.
Settling, horizontal thermal/momentum derivatives, self-shading, general
geometry/withdrawal placement derivatives and SWMM parameter gradients remain
outside this contract.
