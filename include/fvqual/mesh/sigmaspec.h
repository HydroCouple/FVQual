// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   sigmaspec.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  SigmaSpec — the vertical layering configuration of a domain
 *         (F1.2).
 *
 * Describes how a water column is divided vertically, independent of any
 * particular bathymetry: how many layers, how they are distributed between
 * surface and bed, and how thin a layer is allowed to become in shallow
 * water. Resolving a spec yields a normalized interface distribution on
 * [0, 1] that LayeredMesh maps onto real elevations per column.
 *
 * \threadsafety Value type; no shared state.
 */

#ifndef FVQual_MESH_SIGMASPEC_H
#define FVQual_MESH_SIGMASPEC_H

#include "fvqual/fvqual_export.h"

#include <string>
#include <vector>

namespace FVQual::Mesh {

/*!
 * \brief Vertical layering configuration.
 *
 * All lengths are meters (SI internally — see CONTRIBUTING).
 */
struct FVQual_EXPORT SigmaSpec
{
    /*!
     * \brief How the normalized interfaces are distributed.
     */
    enum class Distribution
    {
        Uniform,   //!< Equal-thickness layers.
        Stretched, //!< Song-Haidvogel-style hyperbolic stretching.
        Explicit   //!< Caller-supplied normalized interfaces.
    };

    Distribution distribution = Distribution::Uniform;

    /*!
     * \brief Number of vertical layers (>= 1). Ignored for Explicit, whose
     *        layer count is implied by the supplied interfaces.
     */
    int layerCount = 10;

    /*!
     * \brief Stretching intensity for Distribution::Stretched (> 0; larger
     *        clusters interfaces more strongly).
     */
    double theta = 0.0;

    /*!
     * \brief Surface/bottom blend in [0, 1] for Distribution::Stretched:
     *        0 clusters toward the surface, 1 toward the bed. Surface
     *        clustering is the usual choice for thermally stratified
     *        reservoirs, where the gradient lives in the epilimnion.
     */
    double b = 0.0;

    /*!
     * \brief Minimum layer thickness (m, >= 0). Thin layers are lifted to
     *        this thickness and the surplus taken proportionally from
     *        thicker ones — the defense against vanishing sigma layers in
     *        shallow arms.
     */
    double minLayerThickness = 0.05;

    /*!
     * \brief Normalized interfaces for Distribution::Explicit: layerCount + 1
     *        strictly increasing values, first 0 (surface), last 1 (bed).
     */
    std::vector<double> explicitInterfaces;

    /*!
     * \brief Resolves the spec into a normalized interface distribution.
     * \param[out] interfaces Normalized interfaces on [0, 1], surface first.
     * \param[out] message    Failure description when returning false.
     */
    [[nodiscard]] bool resolve(std::vector<double> &interfaces,
                               std::string &message) const;

    /*!
     * \brief Layer count the spec will actually produce.
     */
    [[nodiscard]] int resolvedLayerCount() const noexcept;
};

} // namespace FVQual::Mesh

#endif // FVQual_MESH_SIGMASPEC_H
