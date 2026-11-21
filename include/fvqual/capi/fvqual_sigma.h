/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Caleb Buahin
 *
 * fvqual_sigma.h — vertical layering specification.
 *
 * A plain value struct rather than a handle: it carries no state and no
 * ownership, so a binding can build one on the stack and pass it by pointer
 * without any lifecycle to manage.
 */

#ifndef FVQual_CAPI_SIGMA_H
#define FVQual_CAPI_SIGMA_H

#include "fvqual/capi/fvqual_engine.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief How the normalized layer interfaces are distributed.
 */
typedef enum FVQ_SigmaDistribution
{
    FVQ_SIGMA_UNIFORM = 0,  /*!< Equal-thickness layers. */
    FVQ_SIGMA_STRETCHED,    /*!< Song-Haidvogel hyperbolic stretching. */
    FVQ_SIGMA_EXPLICIT      /*!< Caller-supplied normalized interfaces. */
} FVQ_SigmaDistribution;

/*!
 * \brief Vertical layering configuration.
 *
 * Fill with fvq_sigma_default() and then override; that way new fields added
 * in future versions acquire sensible values without every caller changing.
 */
typedef struct FVQ_SigmaSpec
{
    FVQ_SigmaDistribution distribution;

    /*! Number of layers (>= 1). Ignored when distribution is EXPLICIT. */
    int layerCount;

    /*! Stretching intensity (> 0) for STRETCHED; larger clusters more. */
    double theta;

    /*! Surface/bottom blend in [0, 1] for STRETCHED: 0 clusters toward the
     *  surface, which is where a stratified reservoir keeps its gradient. */
    double blend;

    /*! Minimum layer thickness (m, >= 0). The defense against vanishing
     *  layers in shallow arms. */
    double minLayerThickness;

    /*! Normalized interfaces for EXPLICIT: layerCount + 1 strictly increasing
     *  values, first 0 (surface), last 1 (bed). NULL otherwise. */
    const double *interfaces;

    /*! Number of entries in \c interfaces. */
    int interfaceCount;
} FVQ_SigmaSpec;

/*!
 * \brief Fills \p spec with defaults: 10 uniform layers, 0.05 m minimum
 *        thickness.
 */
FVQual_CAPI void fvq_sigma_default(FVQ_SigmaSpec *spec);

/*!
 * \brief Validates a specification without building anything.
 *
 * Useful for a binding that wants to reject bad input at the point the user
 * typed it rather than several calls later.
 */
FVQual_CAPI FVQ_ErrorCode fvq_sigma_validate(const FVQ_SigmaSpec *spec);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* FVQual_CAPI_SIGMA_H */
