/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Caleb Buahin
 * Reference kernels and restricted actual-solver passive scalar derivatives.
 */
#ifndef FVQUAL_CAPI_DIFFERENTIATION_H
#define FVQUAL_CAPI_DIFFERENTIATION_H
#include "fvqual/capi/fvqual_engine.h"
#include "fvqual/capi/fvqual_sim.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef enum FVQ_Kernel {
    /* x=[concentration, decay_rate_s-1, duration_s], y=C*exp(-k*dt). */
    FVQ_KERNEL_DECAY = 0,
    /* x=[temperature_C, heat_flux_W_m-2, duration_s, depth_m, density_kg_m-3, cp_J_kg_K].
       y=T+q*dt/(depth*density*cp). Fixed area, fully mixed column. */
    FVQ_KERNEL_COLUMN_HEAT = 1,
    /* x=[volume_m3, concentration, incoming_volume_m3, incoming_mass].
       y=(V*C+M)/(V+Vin). Compatible units for C and M. */
    FVQ_KERNEL_MIX = 2
} FVQ_Kernel;
/* Finite inputs; nonnegative decay rate/duration, positive depth/density/cp,
   nonnegative volumes with a positive combined volume. Invalid output is untouched. */
FVQual_CAPI FVQ_ErrorCode fvq_kernel_value(FVQ_Kernel kernel, const double *x, int count, double *value);
FVQual_CAPI FVQ_ErrorCode fvq_kernel_jvp(FVQ_Kernel kernel, const double *x,
    const double *direction, int count, double *product);
FVQual_CAPI FVQ_ErrorCode fvq_kernel_vjp(FVQ_Kernel kernel, const double *x,
    int count, double outputCotangent, double *inputCotangent);

/* Actual upwind advection + implicit vertical diffusion, passive pollutant
   only. Hydraulics, geometry, source placement, Kz and time schedule frozen.
   Refuses temperature/salinity/age, MUSCL, horizontal diffusion, reactions,
   settling, internal transfers and open boundaries. Source concentration
   derivatives include gross withdrawals. No SWMM parameter adjoint implied.
   Serial access to a simulation/trace. Recording is bounded to 1..10000 steps. */
typedef struct FVQ_ScalarTrace_ FVQ_ScalarTrace;
typedef struct FVQ_TraceSize {
    int64_t cells, steps, sources, storageBytes;
} FVQ_TraceSize;
FVQual_CAPI FVQ_ErrorCode fvq_sim_begin_scalar_trace(FVQ_Simulation *sim,
    int species, int64_t maxSteps);
FVQual_CAPI FVQ_ErrorCode fvq_sim_end_scalar_trace(FVQ_Simulation *sim);
/* Owned immutable snapshot, outlives simulation; destroy exactly once.
   Simulation checkpoint restore restores its recording history as well. */
FVQual_CAPI FVQ_ErrorCode fvq_sim_capture_scalar_trace(const FVQ_Simulation *sim,
    FVQ_ScalarTrace **out);
FVQual_CAPI void fvq_scalar_trace_destroy(FVQ_ScalarTrace *trace);
FVQual_CAPI FVQ_ErrorCode fvq_scalar_trace_size(const FVQ_ScalarTrace *trace, FVQ_TraceSize *out);
/* States: step-major [steps+1,cells], first block initial. Sources: concatenated
   per step (including inactive/withdrawal sources with zero input derivative).
   Times: [steps+1], seconds from simulation start. ID pointer owned by trace. */
FVQual_CAPI FVQ_ErrorCode fvq_scalar_trace_baseline(const FVQ_ScalarTrace *trace,
    double *states, int64_t stateCount, double *sources, int64_t sourceCount,
    double *times, int64_t timeCount);
FVQual_CAPI const char *fvq_scalar_trace_source_id(const FVQ_ScalarTrace *trace, int64_t source);
/* Input initial [cells], sources [sources], output states [steps+1,cells].
   Positive diffusion multiplier scales recorded vertical operator only.
   Input derivatives are evaluated at the recorded baseline (multiplier=1).
   Errors leave output arrays untouched; zero-length arrays may be NULL. */
FVQual_CAPI FVQ_ErrorCode fvq_scalar_trace_evaluate(const FVQ_ScalarTrace *trace,
    const double *initial, int64_t cells, const double *sources, int64_t sourceCount,
    double diffusionMultiplier, double *states, int64_t stateCount);
FVQual_CAPI FVQ_ErrorCode fvq_scalar_trace_jvp(const FVQ_ScalarTrace *trace,
    const double *initialDirection, int64_t cells, const double *sourceDirection, int64_t sourceCount,
    double diffusionDirection, double *stateDirection, int64_t stateCount);
FVQual_CAPI FVQ_ErrorCode fvq_scalar_trace_vjp(const FVQ_ScalarTrace *trace,
    const double *stateCotangents, int64_t stateCount, double *initialCotangents,
    int64_t cells, double *sourceCotangents, int64_t sourceCount, double *diffusionCotangent);
#ifdef __cplusplus
}
#endif
#endif
