# Verification animations

Six benchmarks with closed-form or theoretical solutions, animated against the
model. Open [`index.html`](index.html) in a browser — it is self-contained, with
no server, network, or external assets.

| Scenario | Reference | What it exercises |
|---|---|---|
| Ritter dam break | Ritter's rarefaction, exact at every instant | External mode, dry front, shock capturing |
| Barotropic seiche | `eta = A cos(pi x/L) cos(2 pi t/T)`, `T = 2L/sqrt(gH)` | Surface wave speed and pressure gradient |
| Wind setup | `d(eta)/dx = tau/(rho g H)`, with an unforced control | Wind stress magnitude and steady balance |
| Transport schemes | Exact translated square pulse | MUSCL + FCT against first-order upwind |
| Internal seiche | Two-layer theory, `c = sqrt(g' h1 h2 / H)` | Mode splitting and baroclinic pressure |
| Lock exchange | Gravity-current front speed `0.5 sqrt(g' H)` | Baroclinic adjustment and density transport |

## Regenerating

```bash
cmake --build build --target fvqual-verify
./build/verification/fvqual-verify verification/data.json
python3 verification/build_page.py
```

`scenarios.cpp` runs each case through the real solver at its own adaptive time
step and decimates the result to animation frames; `build_page.py` inlines the
data into `page.html` to produce `index.html`. Both `data.json` and
`index.html` are committed so the page can be opened without a build.

Edit `page.html` for presentation changes and `scenarios.cpp` for the physics —
never `index.html`, which is generated.
