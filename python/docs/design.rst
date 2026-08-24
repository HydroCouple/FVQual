.. _bindings_design:

======================
Design of the bindings
======================

.. contents:: On this page
   :local:
   :depth: 2

The structure follows :openswmm:`openswmm.engine <>`: a C ABI in the middle, a
Cython layer on top of it, and a thin Python package over that. The reasons are
the same in both projects, and they are worth stating because the obvious
alternative — binding the C++ directly with pybind11 — looks simpler and is
not.


Why a C API sits in the middle
==============================

Binding C++ directly couples the wheel to the exact compiler, standard library,
and ABI that built it. Every refactor of a C++ header becomes a rebuild of the
extension, and every platform combination becomes a support matrix entry.

A C ABI with opaque handles is a firewall. The engine's C++ internals can be
restructured — as they were when the finite-volume sweeps were rewritten as
gathers for threading — without a compiled wheel noticing. The cost is that
every type crossing the boundary must be expressible in C, which is a
discipline rather than a limitation: it forces the runtime surface to be
explicit.

The C API lives in ``include/fvqual/capi/`` and is verified two ways: a
pure-C test suite compiled by a C compiler, and a header compilation under
``-std=c89 -pedantic-errors`` so the ABI promise is enforced rather than
stated.


Runtime forcing and operational control
=======================================

Scripted control **while the simulation runs** is a first-class requirement,
not an afterthought. The point of an embedded engine is to close a loop around
it: operate a gate from a rule, drive a release from a downstream target,
replay an observed schedule, or wrap the whole thing in an optimiser.

The engine side already exists. Every lever below is a plain setter designed to
be called between steps:

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - Engine call
     - What it controls
   * - ``SourceSet::setDischarge``
     - Inflow or withdrawal rate on a schedule
   * - ``StructureSet::setOpening``
     - Gate position in [0, 1] — the operator's actual lever
   * - ``StructureSet::setTailwater``
     - Downstream stage a structure discharges against
   * - ``StructureSet::setSchedule`` / ``clearSchedule``
     - Switch between a scheduled release, where the rating curve becomes a
       capacity limit, and the free rating curve
   * - ``OutletGroupSet::setDemand`` / ``setTarget``
     - Total release and target release temperature for a multi-port group
   * - ``TransferSet::setDischarge``
     - Rate of an internal transfer between columns

with the matching readers — ``discharge``, ``isRunning``, ``deliveredRelease``,
``achievedTemperature``, ``deliveredTemperature`` — so a control script can see
the consequence of what it just did before deciding the next step.

.. warning::

   These are **between-step** calls. ``advance()`` subcycles the barotropic
   mode internally and re-rates head-driven structures on every substep, so
   mutating a structure from inside a callback that fires mid-``advance()`` is
   undefined — for the same reason writing through a state view mid-step is.
   Drive the loop with ``steps()`` and set forcing at the top of each
   iteration.


Zero-copy state views
=====================

State is stored as flat, allocation-stable arrays for the lifetime of a
simulation, which is what makes a genuine zero-copy view possible: a NumPy
array can point straight at engine memory rather than copying it every step.

That is a sharp edge as well as a feature, so the contract is explicit:

- A view is valid until the simulation is destroyed or its mesh is replaced.
  Both bump a generation counter, and a stale view raises rather than reading
  freed memory.
- Writing through a view is legitimate for initial conditions and for forcing
  applied between steps.
- Writing through a view *during* a step is undefined. The engine is reading
  those same arrays.


Planned package layout
======================

One module per C API header, mirroring the openswmm rule:

.. code-block:: text

   python/
   └── fvqual/
       ├── __init__.py
       ├── _mesh.pyx          Mesh construction, UGRID input, geometry queries
       ├── _species.pyx       Species registry: temperature, age, constituents
       ├── _options.pyx       Solver options
       ├── _simulation.pyx    Lifecycle, steps()/run(), context manager
       ├── _state.pyx         Zero-copy views and bulk setters
       ├── _forcing.pyx       Sources, structures, outlet groups, transfers
       └── _output.pyx        UGRID writer registration

The simulation object is a context manager that guarantees teardown, and
``steps()`` yields control between internal steps — the iteration point a
control loop needs.
