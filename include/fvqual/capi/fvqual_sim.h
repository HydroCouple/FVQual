/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Caleb Buahin
 *
 * fvqual_sim.h — simulation lifecycle, stepping, and diagnostics.
 *
 * LIFECYCLE
 *
 *   create -> [add species] -> [set surface] -> initialize -> step/advance -> destroy
 *
 * The split between create and initialize exists because species must be
 * registered before state can be allocated: the state store is laid out
 * species-major over cells, and that layout has to be known before the first
 * byte is reserved. Making the caller state its species up front is what buys
 * the allocation-stable arrays that PY3's zero-copy views depend on.
 *
 * THREADING
 *
 * Nothing here touches global state, so independent simulations may run
 * concurrently on separate threads. fvq_sim_step and fvq_sim_advance perform
 * no callbacks and touch no interpreter, so a binding may release the GIL
 * around them.
 */

#ifndef FVQual_CAPI_SIM_H
#define FVQual_CAPI_SIM_H

#include "fvqual/capi/fvqual_engine.h"
#include "fvqual/capi/fvqual_mesh.h"
#include "fvqual/capi/fvqual_options.h"
#include "fvqual/capi/fvqual_sigma.h"

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief Opaque handle to a simulation.
 */
typedef struct FVQ_Simulation_ FVQ_Simulation;

/* ── Construction ──────────────────────────────────────────────────────── */

/*!
 * \brief Creates a simulation over a copy of \p mesh.
 *
 * The mesh is copied, so it may be destroyed immediately afterwards. The
 * simulation is not yet usable: register species, set the initial surface if
 * it is not uniform, then call fvq_sim_initialize.
 *
 * \param[in]  mesh    Horizontal mesh (copied).
 * \param[in]  sigma   Vertical layering.
 * \param[in]  options Settings, or NULL for defaults.
 * \param[out] out     Receives the simulation; destroy with fvq_sim_destroy.
 */
FVQual_CAPI FVQ_ErrorCode fvq_sim_create(const FVQ_Mesh *mesh, const FVQ_SigmaSpec *sigma,
                                         const FVQ_Options *options,
                                         FVQ_Simulation **out);

/*!
 * \brief Releases a simulation. Passing NULL is a no-op.
 */
FVQual_CAPI void fvq_sim_destroy(FVQ_Simulation *sim);

/*!
 * \brief Registers a transported species. Before initialize only.
 *
 * \param[in]  id           Unique identifier, non-empty.
 * \param[in]  units        Units string, may be NULL.
 * \param[in]  initialValue Uniform initial concentration.
 * \param[out] index        Receives the species index; may be NULL.
 */
FVQual_CAPI FVQ_ErrorCode fvq_sim_add_species(FVQ_Simulation *sim, const char *id,
                                              const char *units, double initialValue,
                                              int *index);

/*!
 * \brief Registers the reserved temperature species (deg C).
 *
 * Temperature is transported by the same machinery as every other scalar, not
 * by a special-purpose path — one advection scheme, one set of conservation
 * tests. Registering it is what makes the model baroclinic; without it the
 * water is homogeneous and only the barotropic mode does anything.
 */
FVQual_CAPI FVQ_ErrorCode fvq_sim_add_temperature(FVQ_Simulation *sim,
                                                  double initialCelsius, int *index);

/*!
 * \brief Registers the reserved water-age species (days).
 */
FVQual_CAPI FVQ_ErrorCode fvq_sim_add_age(FVQ_Simulation *sim, int *index);

/*!
 * \brief Sets the initial water-surface elevation per column (m).
 *
 * Before initialize only. Without this the surface is uniform at zero.
 *
 * \param[in] elevation One value per column.
 * \param[in] count     Must equal the mesh face count.
 */
FVQual_CAPI FVQ_ErrorCode fvq_sim_set_initial_surface(FVQ_Simulation *sim,
                                                      const double *elevation,
                                                      int64_t count);

/*!
 * \brief Allocates state and prepares the solver. Call once, after species.
 *
 * Fails with FVQ_ERR_DRY_COLUMN if any column has no water: three-dimensional
 * wetting and drying is not implemented, and the engine says so rather than
 * dividing by a vanishing thickness.
 */
FVQual_CAPI FVQ_ErrorCode fvq_sim_initialize(FVQ_Simulation *sim);

/* ── Stepping ──────────────────────────────────────────────────────────── */

/*!
 * \brief Largest stable internal time step for the current state (s).
 */
FVQual_CAPI double fvq_sim_stable_time_step(const FVQ_Simulation *sim);

/*!
 * \brief Advances exactly \p dt seconds, subcycling the external mode within.
 */
FVQual_CAPI FVQ_ErrorCode fvq_sim_step(FVQ_Simulation *sim, double dt);

/*!
 * \brief Advances \p duration seconds, substepping at the Courant limit and
 *        landing exactly on the target time.
 */
FVQual_CAPI FVQ_ErrorCode fvq_sim_advance(FVQ_Simulation *sim, double duration);

/*!
 * \brief Elapsed simulated time (s).
 */
FVQual_CAPI double fvq_sim_time(const FVQ_Simulation *sim);

/* ── Dimensions ────────────────────────────────────────────────────────── */

FVQual_CAPI int64_t fvq_sim_column_count(const FVQ_Simulation *sim);
FVQual_CAPI int     fvq_sim_layer_count(const FVQ_Simulation *sim);
FVQual_CAPI int64_t fvq_sim_cell_count(const FVQ_Simulation *sim);
FVQual_CAPI int     fvq_sim_species_count(const FVQ_Simulation *sim);

/* ── Diagnostics ───────────────────────────────────────────────────────── */

/*!
 * \brief Per-step diagnostics of the coupled solver.
 */
typedef struct FVQ_Diagnostics
{
    double timeStep;              /*!< Internal step last taken (s). */
    double maxHorizontalVelocity; /*!< Largest |u| over the cells (m/s). */
    double maxVerticalFlux;       /*!< Largest interface flux (m/s). */

    /*! Free-surface continuity residual (m/s). The vertical velocity is
     *  integrated up from a no-flux bed, so its value at the surface must
     *  vanish if the two modes agree. A non-negligible value means scalar
     *  transport can no longer be trusted, however plausible the fields look —
     *  worth asserting on in any long run. */
    double surfaceFluxResidual;

    double  totalVolume;   /*!< Domain water volume (m³). */
    double  maxRichardson; /*!< Largest gradient Richardson number. */
    int64_t externalSteps; /*!< External substeps in the last internal step. */
    int64_t steps;         /*!< Internal steps taken by the last advance. */
} FVQ_Diagnostics;

/*!
 * \brief Copies the current diagnostics.
 */
FVQual_CAPI FVQ_ErrorCode fvq_sim_diagnostics(const FVQ_Simulation *sim,
                                              FVQ_Diagnostics *out);

/*!
 * \brief Total water volume (m³).
 */
FVQual_CAPI double fvq_sim_total_volume(const FVQ_Simulation *sim);

/*!
 * \brief Volume-weighted total mass of a species (concentration x m³).
 */
FVQual_CAPI FVQ_ErrorCode fvq_sim_species_mass(const FVQ_Simulation *sim, int species,
                                               double *out);

/* Copies a transported scalar and its current geometry for observation/review.
   Cells are column-major, layers from surface to bed; interfaces have
   [columns, layers+1] elevations in metres. All counts must match. Errors leave
   outputs untouched. Copies remain valid after stepping/destroying the handle. */
FVQual_CAPI FVQ_ErrorCode fvq_sim_scalar_snapshot(const FVQ_Simulation *sim, int species,
    double *values, double *volumes, int64_t cells, double *interfaces, int64_t interfaceCount);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* FVQual_CAPI_SIM_H */
