// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   ilayeredmeshcomponentdataitem.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  ILayeredMeshComponentDataItem — proposed HydroCouple interface for a
 *         layered 2-D mesh (F8.2).
 *
 * **The gap.** HydroCouple describes an unstructured mesh through
 * `IPolyhedralSurfaceComponentDataItem`, whose values are indexed by one entity
 * dimension — patch, edge, or vertex. A layered model has values per patch *per
 * layer*, and there is no interface for that. `IRegularGrid3DComponentDataItem`
 * carries three dimensions but only for a structured grid, so it cannot describe
 * an unstructured surface extruded downward.
 *
 * A stratified reservoir engine whose coupling surface cannot express
 * stratification has given up most of what it is for, so this is the interface
 * that has to exist before `FVQualComponent` can exchange anything worth
 * exchanging.
 *
 * **Why it lives here and not in HydroCouple.** It is a proposal, and a proposal
 * for an interface is worth more with a working implementation and a real
 * consumer behind it than as a header on its own. It is written in FVQual's
 * namespace so that adopting it upstream is a deliberate act rather than a
 * silent dependency, and so that nothing in HydroCouple has to change before
 * this can be built and measured.
 *
 * **The part that is easy to get wrong.** A consumer handed temperature on
 * `{patch, layer}` cannot interpret it. Layer 7 of column 3 is at *some*
 * elevation, and under a terrain-following sigma grid that elevation moves every
 * step as the free surface rises and falls. Values without their vertical
 * coordinate are not data, they are an array.
 *
 * The layered UGRID writer already solves this for files, by declaring the CF
 * `ocean_sigma_coordinate` and publishing its formula terms so any CF-aware
 * reader reconstructs true elevations without the writer's help. This interface
 * takes the same position: it publishes the sigma coordinate and the two terms
 * that close it, so a coupling partner can compute
 *
 *     z(patch, k) = eta(patch) + sigma(k) * (depth(patch) + eta(patch))
 *
 * for itself, at the time the values are for. An interface that made the
 * consumer ask the producer for elevations out of band would be one where the
 * two can silently disagree about which step they are describing.
 *
 * **Sign conventions**, stated because there are several in use and picking
 * silently is how a coupled model ends up upside down:
 *
 * - Elevations are positive up, in the mesh's own vertical datum.
 * - `sigmaInterfaces()` runs 0 at the free surface to −1 at the bed, which is
 *   the CF convention rather than FVQual's internal one. The implementation
 *   converts.
 * - `depth()` is positive down, also per CF.
 */

#ifndef FVQual_COMPONENT_ILAYEREDMESHCOMPONENTDATAITEM_H
#define FVQual_COMPONENT_ILAYEREDMESHCOMPONENTDATAITEM_H

#include "fvqual/fvqual_export.h"

#include "hydrocouplespatial.h"

#include <cstdint>

namespace FVQual::Component {

/*!
 * \brief A polyhedral surface whose values are resolved into vertical layers.
 *
 * Proposed for `HydroCouple::Spatial`. Extends the polyhedral-surface item
 * rather than replacing it: a layered mesh *is* a surface, so a partner that
 * only understands surfaces can still discover the patch topology through the
 * base interface and take the surface layer.
 *
 * **Shape.** Values are indexed `{patch, layer}`, layer 0 at the free surface
 * and increasing downward — the same ordering `LayeredMesh` uses internally, so
 * no reindexing happens at the boundary and none can be got wrong there.
 */
class FVQual_EXPORT ILayeredMeshComponentDataItem
    : public virtual HydroCouple::Spatial::IPolyhedralSurfaceComponentDataItem
{
public:
    virtual ~ILayeredMeshComponentDataItem() = default;

    /*!
     * \brief The IDimension of the layers (the last dimension of shape()).
     */
    [[nodiscard]] virtual HydroCouple::IDimension *layerDimension() const = 0;

    /*!
     * \brief Number of layers, which is `layerDimension()->length()`.
     */
    [[nodiscard]] virtual int layerCount() const = 0;

    /*!
     * \brief Normalized sigma at each layer *interface*, `layerCount() + 1`
     *        values, 0 at the free surface to −1 at the bed (CF convention).
     *
     * Constant in time: the grid stretches with the surface, it does not
     * re-partition. Together with eta and depth this is the whole vertical
     * coordinate.
     */
    [[nodiscard]] virtual const double *sigmaInterfaces() const = 0;

    /*!
     * \brief Free-surface elevation of each patch (m, positive up), at the time
     *        the current values describe.
     *
     * A formula term, not a diagnostic. It moves every step, and it must come
     * from the same step as the values or the elevations a consumer computes
     * will describe a mesh that no longer exists.
     */
    [[nodiscard]] virtual const double *surfaceElevation() const = 0;

    /*!
     * \brief Still-water depth of each patch (m, positive down), per CF.
     */
    [[nodiscard]] virtual const double *depth() const = 0;

    /*!
     * \brief Elevation of interface \p interfaceIndex in \p patchIndex (m).
     *
     * A convenience over the formula terms, provided because the arithmetic is
     * easy to get subtly wrong — the sign of sigma and the sign of depth pull in
     * opposite directions — and because a wrong elevation is not something a
     * coupled run reports. Implementations must return exactly what the formula
     * terms imply.
     *
     * \param[in] patchIndex     Patch, in the surface's own ordering.
     * \param[in] interfaceIndex 0 at the free surface through `layerCount()` at
     *                           the bed.
     */
    [[nodiscard]] virtual double interfaceElevation(int64_t patchIndex,
                                                    int interfaceIndex) const = 0;
};

} // namespace FVQual::Component

#endif // FVQual_COMPONENT_ILAYEREDMESHCOMPONENTDATAITEM_H
