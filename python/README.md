# FVQual — Python bindings

**Status: not yet implemented.** This directory holds the design and
documentation for the forthcoming `fvqual` Python package. No extension module
is built here yet.

What exists today is the layer the bindings will sit on: a C89-clean C API with
opaque handles, in [`include/fvqual/capi/`](../include/fvqual/capi/), verified
by a pure-C test suite and a `-std=c89 -pedantic-errors` header compilation.

The structure follows
[openswmm.engine](https://github.com/HydroCouple/openswmm.engine): a C ABI
firewall, a Cython layer, and a thin Python package — one module per C API
header. See [`docs/design.rst`](docs/design.rst) for why, and for the runtime
control surface the bindings will expose.

## Building these docs

```bash
pip install -r docs/requirements.txt
make -C docs html
```

The rendered pages land in `docs/_build/html/`. Nothing here imports the
unbuilt package, so the docs build cleanly today.
