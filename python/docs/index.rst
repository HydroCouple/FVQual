======
FVQual
======

**3D unstructured finite-volume hydrodynamic and water-quality engine —
Python documentation.**

FVQual solves three-dimensional hydrostatic free-surface flow and multispecies
water quality on an unstructured tri-quad horizontal mesh with
terrain-following sigma layers. It is the HydroCouple-native successor to the
CE-QUAL-W2 class of reservoir models, and it is built to be driven from a
script: opened, stepped, inspected, and steered while it runs.

.. admonition:: Status — initial C ABI package implemented
   :class: important

   ``python/fvqual`` provides owned meshes and simulations, interval transfer
   receipts, checkpoints, recorded/live surface-coupling examples, ensembles
   and reference-kernel, passive-scalar, controlled thermal-column and conditional
   reactive-transport JVP/VJP operations. It uses standard-library ctypes
   against the C89-clean ABI and requires a shared core library. See
   ``python/README.md`` and ``examples/swmm_reservoir/README.md`` for runnable
   examples and the exact supported scope.

   The design pages also describe future runtime controls and zero-copy views.
   Those capabilities, Cython extensions, binary wheels and full solver
   derivatives are not supplied by this initial package.

.. grid:: 2

    .. grid-item-card:: C++ / C API reference
        :link: https://hydrocouple.github.io/FVQual/

        The engine as it exists today. Doxygen reference for the solver
        classes, the mesh and state model, and the C ABI the Python package
        will bind.

    .. grid-item-card:: Design of the bindings
        :link: design
        :link-type: doc

        Why a C API sits between Python and C++, what runtime control the
        bindings expose, and the rules that keep a zero-copy view safe.


What the bindings are for
=========================

An engine that can only be run from a configuration file can be *executed*.
An engine that can be driven from a script can be *used*: closed inside a
control loop, wrapped in an optimiser, replayed against observed operations, or
stepped forward while a rule decides what the gates should do next.

That is the point of the binding work, and it is why the runtime control
surface — not the file formats — is what the design pages concentrate on.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   design
   api
   related


Related projects
================

FVQual belongs to the HydroCouple family and follows the same conventions for
build, testing, and documentation:

- :hydrocouple:`HydroCouple <>` — the component-based modeling interface
  definitions FVQual targets.
- :sdk:`HydroCoupleSDK <>` — the SDK providing the data plane, component
  lifecycle, and mesh IO.
- :openswmm:`openswmm.engine <>` — the SWMM engine whose C-API-plus-Cython
  structure these bindings follow, and whose multispecies reaction engine
  FVQual will reuse for water-quality kinetics.
