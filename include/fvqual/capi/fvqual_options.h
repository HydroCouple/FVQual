/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Caleb Buahin
 *
 * fvqual_options.h — physical and numerical settings of a simulation.
 *
 * A flat value struct mirroring the engine's option records. PY1 adds keyed
 * string access on top for scripting convenience; this is the typed surface
 * underneath it, and the one the C API guarantees.
 */

#ifndef FVQual_CAPI_OPTIONS_H
#define FVQual_CAPI_OPTIONS_H

#include "fvqual/capi/fvqual_engine.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief Advective discretization for transported scalars.
 */
typedef enum FVQ_AdvectionScheme
{
    FVQ_ADVECTION_UPWIND = 0, /*!< First-order upwind: bounded, diffusive. */
    FVQ_ADVECTION_MUSCL_FCT   /*!< MUSCL with Zalesak flux correction. */
} FVQ_AdvectionScheme;

/*!
 * \brief Simulation settings. SI units throughout.
 *
 * Fill with fvq_options_default() and override selectively.
 */
typedef struct FVQ_Options
{
    double gravity;          /*!< m/s². */
    double referenceDensity; /*!< Boussinesq reference density, kg/m³. */
    double coriolis;         /*!< Coriolis parameter, 1/s; zero for small basins. */

    /* Vertical mixing (Richardson-damped). */
    double neutralViscosity;      /*!< m²/s. */
    double backgroundViscosity;   /*!< m²/s. */
    double backgroundDiffusivity; /*!< m²/s; keep small or the thermocline erodes. */
    double convectiveMixing;      /*!< m²/s, applied where the column is unstable. */

    /* Horizontal mixing. */
    double horizontalViscosity;   /*!< m²/s. */
    double horizontalDiffusivity; /*!< m²/s. */

    FVQ_AdvectionScheme advectionScheme;

    /* Forcing. */
    double windStressX;     /*!< Pa. */
    double windStressY;     /*!< Pa. */
    double surfaceHeatFlux; /*!< W/m², positive warming. */
    double heatCapacity;    /*!< J/(kg·K). */

    double manningN; /*!< Bed roughness, s/m^(1/3). */

    /* Time stepping. */
    double cflNumber;   /*!< Courant number for the internal step, in (0, 1]. */
    double minTimeStep; /*!< s. */
    double maxTimeStep; /*!< s. */
} FVQ_Options;

/*!
 * \brief Fills \p options with the engine defaults.
 */
FVQual_CAPI void fvq_options_default(FVQ_Options *options);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* FVQual_CAPI_OPTIONS_H */
