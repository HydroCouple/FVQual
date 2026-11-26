/* SPDX-License-Identifier: MIT */
#ifndef FVQUAL_CAPI_REACTIVE_H
#define FVQUAL_CAPI_REACTIVE_H
#include "fvqual/capi/fvqual_sim.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Owned immutable native upwind/vertical-diffusion/multispecies-reaction trace.
 * Hydraulic and thermal trajectories, geometry, Kz, source discharge,
 * placement and composition are prescribed. Supports initial pollutant/MSX
 * concentrations and kinetic constants. No SWMM parameter or reservoir
 * thermal/momentum adjoint. Refuses density-changing rates, self-shading,
 * settling, horizontal diffusion, internal transfers and open boundaries.
 * Native adaptive kinetics is differentiated locally through its controller;
 * switching/domain boundaries and unconverged solves refuse derivatives.
 * Errors leave outputs untouched. Serial access per simulation/trace.
 */
typedef struct FVQ_ReactiveTrace_ FVQ_ReactiveTrace;
typedef struct FVQ_ReactiveSize { int64_t cells, species, steps, parameters; } FVQ_ReactiveSize;
FVQual_CAPI FVQ_ErrorCode fvq_sim_begin_reactive_trace(FVQ_Simulation *,int64_t maxSteps);
FVQual_CAPI FVQ_ErrorCode fvq_sim_end_reactive_trace(FVQ_Simulation *);
FVQual_CAPI FVQ_ErrorCode fvq_sim_capture_reactive_trace(const FVQ_Simulation *,FVQ_ReactiveTrace **out);
FVQual_CAPI void fvq_reactive_trace_destroy(FVQ_ReactiveTrace *);
FVQual_CAPI FVQ_ErrorCode fvq_reactive_trace_size(const FVQ_ReactiveTrace *,FVQ_ReactiveSize *);
/* JSON metadata pointer is owned by the trace. */
FVQual_CAPI const char *fvq_reactive_trace_metadata(const FVQ_ReactiveTrace *);
/* States are [steps+1,species,cells], species-major within each time. Times
 * are elapsed seconds. Constants follow metadata parameter_names order. */
FVQual_CAPI FVQ_ErrorCode fvq_reactive_trace_baseline(const FVQ_ReactiveTrace *,
    double *states,int64_t stateCount,double *parameters,int64_t parameterCount,double *times,int64_t timeCount);
FVQual_CAPI FVQ_ErrorCode fvq_reactive_trace_evaluate(const FVQ_ReactiveTrace *,
    const double *initial,int64_t initialCount,const double *parameters,int64_t parameterCount,double *states,int64_t stateCount);
FVQual_CAPI FVQ_ErrorCode fvq_reactive_trace_jvp(const FVQ_ReactiveTrace *,
    const double *initialDirection,int64_t initialCount,const double *parameterDirection,int64_t parameterCount,double *states,int64_t stateCount);
FVQual_CAPI FVQ_ErrorCode fvq_reactive_trace_vjp(const FVQ_ReactiveTrace *,
    const double *stateCotangents,int64_t stateCount,double *initialCotangents,int64_t initialCount,
    double *parameterCotangents,int64_t parameterCount,int64_t *peakCellTapeNodes);
/* Update an existing compiled constant; begin a new trace afterwards. */
FVQual_CAPI FVQ_ErrorCode fvq_sim_set_reaction_constant(FVQ_Simulation *,const char *name,double value);
/* Species ID pointer owned by the simulation; valid until reconfiguration or destruction. */
FVQual_CAPI const char *fvq_sim_species_id(const FVQ_Simulation *,int species);
#ifdef __cplusplus
}
#endif
#endif
