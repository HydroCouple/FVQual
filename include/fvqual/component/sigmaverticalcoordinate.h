// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   sigmaverticalcoordinate.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  FVQual's sigma grid, presented as HydroCouple::Spatial::IVerticalCoordinate
 *         (F8.2e).
 *
 * **What this replaces.** An earlier header in this directory proposed
 * `ILayeredMeshComponentDataItem` in FVQual's own namespace, carrying the sigma
 * coordinate as three formula terms — `sigmaInterfaces()`, `surfaceElevation()`
 * and `depth()` — for a consumer to combine itself. That interface has since
 * been adopted upstream into `hydrocouplespatial.h`, and the vertical coordinate
 * was separated out of it, so this file is now an *implementation* rather than a
 * proposal.
 *
 * **Why the formula terms moved out of the interface.** Publishing terms is right
 * for a file, where a reader implements a handful of CF formulae once and reuses
 * them. It is wrong for a coupling boundary: it obliges every partner to
 * reimplement the arithmetic, the sign of sigma and the sign of depth pull in
 * opposite directions, and a partner that gets it wrong produces a plausible
 * profile and no error at all. FVQual knows where its own layers are, so it
 * answers with elevations. The terms are still published to the NetCDF output,
 * where the CF convention applies and a reader expects them.
 *
 * **Conventions.** Elevations positive up in the mesh's datum; interface 0 at the
 * free surface descending to `layerCount()` at the bed. FVQual's internal sigma
 * runs 0 at the surface to 1 at the bed, which is the opposite sign from CF's
 * `ocean_sigma_coordinate`; that conversion stays inside the UGRID writer and
 * does not appear here, because this interface trades in elevations and never
 * exposes a sigma whose sign could be misread.
 *
 * \threadsafety Const after construction apart from \ref refresh(), which must
 *               not run concurrently with readers.
 */

#ifndef FVQual_COMPONENT_SIGMAVERTICALCOORDINATE_H
#define FVQual_COMPONENT_SIGMAVERTICALCOORDINATE_H

#include "fvqual/fvqual_export.h"

#include "hydrocouplespatial.h"

#include <cstdint>

namespace FVQual::Mesh {
class LayeredMesh;
}

namespace FVQual::Component {

/*!
 * \brief Terrain-following sigma, as a HydroCouple vertical coordinate.
 *
 * Holds a reference to the live `LayeredMesh` rather than a copy of its
 * elevations, so a consumer reading during an exchange sees the surface the
 * model currently has. \ref refresh() exists for the case where the mesh object
 * itself is replaced.
 */
class FVQual_EXPORT SigmaVerticalCoordinate final
    : public virtual HydroCouple::Spatial::IVerticalCoordinate
{
public:
    /*!
     * \brief Binds to a mesh.
     * \param[in] mesh Layered mesh; must outlive this object.
     */
    explicit SigmaVerticalCoordinate(const FVQual::Mesh::LayeredMesh &mesh) noexcept;

    ~SigmaVerticalCoordinate() override = default;

    //! Always VerticalCoordinateKind::Sigma.
    [[nodiscard]] HydroCouple::Spatial::VerticalCoordinateKind
    kind() const override;

    [[nodiscard]] int64_t layerCount() const override;
    [[nodiscard]] int64_t columnCount() const override;
    [[nodiscard]] uint64_t geometryEpoch() const override;

    /*!
     * \brief Always true.
     *
     * Sigma layers stretch with the free surface, so every interface elevation
     * changes whenever the stage does. A consumer must not cache them.
     */
    [[nodiscard]] bool isTimeVarying() const override;

    /*!
     * \brief Elevation of one interface (m, positive up).
     *
     * Out-of-range indices return a quiet NaN rather than reading past the end,
     * because a coupling partner's indexing error should surface as an
     * unmistakable value and not as a plausible elevation from a neighbouring
     * column.
     */
    [[nodiscard]] double interfaceElevation(int64_t cellIndex,
                                            int64_t interfaceIndex) const override;

    [[nodiscard]] std::span<const double> interfaceElevations() const override;

    //! Rebinds to a different mesh.
    void refresh(const FVQual::Mesh::LayeredMesh &mesh) noexcept;

private:
    const FVQual::Mesh::LayeredMesh *m_mesh = nullptr;
    uint64_t m_epochOffset = 0;
};

} // namespace FVQual::Component

#endif // FVQual_COMPONENT_SIGMAVERTICALCOORDINATE_H
