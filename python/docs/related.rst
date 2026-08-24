.. _related_projects:

================
Related projects
================

FVQual is one of a family of HydroCouple projects that share build, testing,
and documentation conventions. Each links to the others so the four sets of
documentation read as one.

.. list-table::
   :header-rows: 1
   :widths: 26 40 34

   * - Project
     - What it provides
     - Documentation
   * - `HydroCouple <https://github.com/HydroCouple/HydroCouple>`_
     - Component-based modeling interface definitions. FVQual is built to be a
       HydroCouple 2.0 model component.
     - `Interface reference <https://hydrocouple.github.io/HydroCouple/>`_
   * - `HydroCoupleSDK <https://github.com/HydroCouple/HydroCoupleSDK>`_
     - Concrete base implementations of the interfaces: the typed hyperslab
       data plane, component lifecycle, and mesh IO that FVQual builds on.
     - `SDK reference <https://hydrocouple.github.io/HydroCoupleSDK/>`_
   * - `openswmm.engine <https://github.com/HydroCouple/openswmm.engine>`_
     - The SWMM engine. Its C-API-plus-Cython binding structure is the model
       these bindings follow, and its multispecies reaction engine is what
       FVQual will reuse for water-quality kinetics.
     - `Python documentation <https://hydrocouple.github.io/openswmm.engine/>`_
   * - **FVQual**
     - This project: the 3D finite-volume hydrodynamic and water-quality
       engine.
     - `C++ API reference <https://hydrocouple.github.io/FVQual/>`_
