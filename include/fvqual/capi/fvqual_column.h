/* SPDX-License-Identifier: MIT
 * Controlled native thermal-column experiments and analytic discrete JVP/VJP.
 * One closed column, always wet, Richardson mixing, upwind, zero Manning and
 * Coriolis, no sources/kinetics/settling. Density/heat/wind/shear feedback and
 * depth/area scaling are differentiated. Outputs are temperature, shear_u,
 * shear_v, each surface-to-bed, at initial plus each prescribed step.
 * Not an arbitrary reservoir, SWMM or HydroCouple component adjoint.
 */
#ifndef FVQUAL_CAPI_COLUMN_H
#define FVQUAL_CAPI_COLUMN_H
#include "fvqual/capi/fvqual_engine.h"
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
typedef struct FVQ_ColumnTrace_ FVQ_ColumnTrace;
typedef struct FVQ_ColumnSize {
    int64_t stateSize, parameterCount, timeCount, checkpointCount;
    double forwardDifference, minimumBuoyancyMargin;
    int derivativeSupported;
} FVQ_ColumnSize;
typedef struct FVQ_ColumnAdjointStats {
    int64_t recomputedSteps, peakTapeNodes, peakBlockStates;
} FVQ_ColumnAdjointStats;
/* JSON accepts layers, depth_m, surface_area_m2, bed_area_fraction, dt_s,
   duration_s, checkpoint_stride, initial, parameters and weather_constants.
   Unknown settings are refused. Parameters/names/units and coverage are in
   metadata. Owned trace, destroy once; independent handles may run concurrently. */
FVQual_CAPI FVQ_ErrorCode fvq_column_create(const char *configuration, FVQ_ColumnTrace **out);
FVQual_CAPI void fvq_column_destroy(FVQ_ColumnTrace *trace);
FVQual_CAPI FVQ_ErrorCode fvq_column_size(const FVQ_ColumnTrace *trace, FVQ_ColumnSize *out);
/* Returned metadata UTF-8 pointer is owned until trace destruction. */
FVQual_CAPI const char *fvq_column_metadata(const FVQ_ColumnTrace *trace);
FVQual_CAPI FVQ_ErrorCode fvq_column_baseline(const FVQ_ColumnTrace *trace,
    double *initial, int64_t stateSize, double *parameters, int64_t parameterCount,
    double *times, int64_t timeCount);
/* Evaluate runs the actual solver, including perturbed mixing branches.
   JVP/VJP are local baseline derivatives and refuse branch/clamp boundaries.
   Arrays have exact lengths; outputs are untouched on error. */
FVQual_CAPI FVQ_ErrorCode fvq_column_evaluate(const FVQ_ColumnTrace *trace,
    const double *initial, int64_t stateSize, const double *parameters, int64_t parameterCount,
    double *states, int64_t outputCount);
FVQual_CAPI FVQ_ErrorCode fvq_column_jvp(const FVQ_ColumnTrace *trace,
    const double *initialDirection, int64_t stateSize, const double *parameterDirection, int64_t parameterCount,
    double *stateDirection, int64_t outputCount);
FVQual_CAPI FVQ_ErrorCode fvq_column_vjp(const FVQ_ColumnTrace *trace,
    const double *stateCotangents, int64_t outputCount, double *initialCotangents, int64_t stateSize,
    double *parameterCotangents, int64_t parameterCount, FVQ_ColumnAdjointStats *stats);
#ifdef __cplusplus
}
#endif
#endif
