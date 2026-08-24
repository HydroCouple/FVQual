# Contributing to FVQual

## Ground rules

**SI units internally, always.** Metres, seconds, kilograms, degrees
Celsius, W/m². Unit conversion happens only at IO boundaries — never inside a
kernel, a state store, or a physics routine. This rule exists because
US-customary internals were a recurring defect source in a sibling project.

**Every capability lands with its verification.** A pull request that adds
physics adds the test that pins it: an analytical benchmark, a conservation
identity, an invariant, or a comparison against a reference solution. "Builds
clean" is not verification.

**Record decisions where they bite.** Architectural choices belong in the
header that implements them, next to the code they constrain. If a change needs
to depart from one, update that comment in the same pull request rather than
leaving the explanation contradicting the code.

**Surgical changes.** Touch only what the change requires. Do not reformat or
"improve" adjacent code. Remove imports and helpers that your change orphaned;
leave pre-existing dead code alone and mention it instead.

**Reviewable file IO.** Tests and tools write to `tests/output/`, not to
temporary directories, so results can be inspected after a run.

## Code conventions

- C++20. No virtual dispatch in hot loops — select engines through function
  tables at setup time.
- Structure-of-arrays for all state and mesh data; keep hot arrays
  allocation-stable and device-mirrorable.
- Kernels are written device-portable-by-construction: fixed-size inline
  arrays, no allocation, no exceptions on the hot path. The core library never
  links Kokkos — GPU backends ship as a plugin behind a C ABI.
- Doxygen comments on public headers, in the HydroCoupleSDK style
  (`\brief`, `\param[in]`, `\threadsafety`).
- SPDX header on every source file: `// SPDX-License-Identifier: MIT`.

## Physics provenance

When implementing a formulation adapted from CE-QUAL-W2, cite the theory
manual section and equation number in a code comment. If you consult or adapt
W2 *source code* rather than the manual, record it in
`THIRD_PARTY_LICENSES.md` and carry the ERDC copyright notice in the derived
file. See `THIRD_PARTY_LICENSES.md`.

## Before opening a pull request

```bash
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Update `CHANGELOG.md` under `[Unreleased]` for anything user-visible.
