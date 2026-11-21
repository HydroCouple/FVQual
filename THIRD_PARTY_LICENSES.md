# Third-Party Licenses & Attribution

FVQual is MIT-licensed (see `LICENSE`). This file records third-party material
and the provenance of adapted formulations. Recording provenance is a standing
project policy, not a one-off: any formulation adapted from another model is
cited here with the source it came from and the nature of the adaptation.

## CE-QUAL-W2

CE-QUAL-W2 is released under the **MIT License**, Copyright (c) 2026
Environmental Modeling Team, Environmental Laboratory, Engineer Research and
Development Center (ERDC), U.S. Army Corps of Engineers.

- Upstream: <https://github.com/CE-QUAL-W2-ERDC/CE-QUAL-W2> (v4.5)
- Distribution: <https://www.cequalw2.org/>

**Policy.** FVQual reimplements W2 physics at the *formulation* level from the
published theory manual, with the manual's section/equation numbers cited in
code comments. Copyright does not extend to equations or algorithms, so
formulation-level reimplementation carries no license obligation; the citation
is an academic-attribution practice.

Where FVQual source is **substantially derived from W2 source code**, the
derived file carries the ERDC copyright notice and MIT permission text inline,
and is listed in the table below. Use the **v4.5 ERDC release** as the
provenance basis — earlier distributions (e.g. the v4.1.0 copy in the
HydroCouple tree) ship without license headers in their source files.

| FVQual file | W2 source consulted | Nature of adaptation |
|---|---|---|
| `src/forcing/sources.cpp` (selective-withdrawal branch) | `withdrawal.f90`, `DOWNSTREAM_WITHDRAWAL` (v5 distribution) | **Formulation level.** The governing equations are Bohan and Grace (1973) and Smith et al. (1987), published literature restated in the distributed W2 V5.0 technical note *"Updating Selective Withdrawal Algorithm and Computation"* (eqs. 1a–1d, 5, 6). The Fortran was read only to determine **which published variant the shipped model actually applies** — the classic zone march with a locally evaluated, implicit buoyancy frequency and `c_bi` of 1 (2 near a boundary), rather than the note's updated bisection form. No code was copied; the C++ shares no structure with the Fortran beyond the algorithm the papers describe. |

> **Note on the v5 tree.** The v5 sources in this working copy carry no
> per-file license headers, as the v4.1.0 copy in the HydroCouple tree does not.
> The equations implemented are from the peer-reviewed and technical literature
> cited above and are not themselves subject to copyright; the v4.5 ERDC MIT
> release remains the provenance basis should any file become substantially
> derived from W2 code.

## doxygen-awesome-css

Vendored under `docs/custom/` and used to theme the generated API
documentation. **MIT License**, Copyright (c) 2021–2023 jothepro.

- Upstream: <https://github.com/jothepro/doxygen-awesome-css>

The files are copied verbatim, with only the vendored `header.html` retargeted
from the theme's own branding to this project's. The full license text is
retained in the header of `docs/custom/css/doxygen-awesome.css`. The same theme
and Doxygen release are used by HydroCoupleSDK and the HydroCouple interface
definitions, so the three documentation sets render identically.

## openswmm.engine — shared lineage, not a dependency

`include/fvqual/kinetics/expression.h`, `src/kinetics/expression.cpp`,
`include/fvqual/kinetics/odesolver.h`, and `src/kinetics/odesolver.cpp` began as
`src/engine/math/MathExpr.*` and `src/engine/math/OdeSolver.*` in
[openswmm.engine](https://github.com/HydroCouple/openswmm.engine), which is
Apache-2.0.

They are reproduced here under MIT **by the copyright holder**, who holds the
rights to both projects. This is a relicensing of the author's own work, not a
redistribution under Apache-2.0, so FVQual remains a single-license project and
carries no NOTICE obligation for them.

The two copies are independent from here on but are expected to stay close
enough that a fix in one can be carried across by hand. One such fix already
exists only here: the expression parser accepted malformed input — `""`,
`"2 +"`, `"* 3"` all parsed as valid — and now rejects any postfix sequence that
does not evaluate to exactly one value. That hardening belongs upstream too.

## HydroCoupleSDK / HydroCouple

MIT License, Copyright (c) 2014-2026 Caleb Buahin. Consumed as a build
dependency; no source vendored.

## Build dependencies

Acquired through vcpkg at build time and not redistributed in this repository:

| Package | License |
|---|---|
| netcdf-c | BSD-3-Clause style (UCAR/Unidata) |
| HDF5 | BSD-3-Clause style (HDF Group) |
| GoogleTest | BSD-3-Clause |

## Example datasets

Validation datasets (DeGray Reservoir, Columbia Slough) are distributed with
CE-QUAL-W2. FVQual references them by download script rather than
redistributing copies; confirm redistribution terms before bundling any.
