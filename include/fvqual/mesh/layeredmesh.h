// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   layeredmesh.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  LayeredMesh — the 3-D computational mesh: a horizontal tri-quad
 *         mesh times terrain-following sigma layers.
 *
 * The 3-D mesh is never stored as general polyhedral connectivity. It is the
 * product of a 2-D mesh and a layer count:
 *
 *     cell(column, k) = column * layerCount + k
 *
 * with k = 0 at the surface descending to k = layerCount - 1 at the bed. The
 * consequences are the point of the design:
 *
 * - Vertical neighbours are +/-1, so a column's cells are contiguous and the
 *   implicit vertical solves (diffusion, settling, light attenuation) are
 *   cache-friendly tridiagonal sweeps — one thread per column, the pattern
 *   that maps well to a GPU.
 * - Horizontal neighbours come from the 2-D topology, so a face-flux loop is
 *   the fixed 2-D edge list iterated per layer.
 * - It is exactly the layered-2D UGRID convention used for output.
 *
 * Layers are terrain-following: interface elevations are regenerated from the
 * free surface every time it moves, so no layer is ever added or removed (the
 * bookkeeping CE-QUAL-W2 needs with z-layers disappears).
 *
 * \note Distinct from HydroCouple::SDK::Tools::LayeredMesh, which is a
 *       node-extruded prism/hexahedron mesh for visualization and export.
 *       This class is the finite-volume computational mesh.
 *
 * \threadsafety Concurrent reads are safe; updateSurface() mutates.
 */

#ifndef FVQual_MESH_LAYEREDMESH_H
#define FVQual_MESH_LAYEREDMESH_H

#include "fvqual/fvqual_export.h"
#include "fvqual/mesh/horizontalmesh.h"
#include "fvqual/mesh/sigmaspec.h"

#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace FVQual::Mesh {

/*!
 * \brief Outcome of the hydrostatic-consistency (Haney) mesh screening.
 *
 * Terrain-following coordinates commit a well-known sin: the horizontal
 * pressure gradient becomes the small difference of two large terms wherever
 * layers slope steeply, and truncation error shows up as spurious currents in
 * a fluid that should be at rest. The bathymetric steepness ratio
 *
 *     r = |H_L - H_R| / (H_L + H_R)
 *
 * across each interior edge predicts where that happens. r <= 0.2 is the
 * conventional target for sigma models; beyond it, smooth the bathymetry or
 * use a hybrid z-sigma vertical coordinate below the transition depth.
 *
 * This is a screening diagnostic, not a correctness gate: the model runs with
 * a failing mesh, but the Phase 3 quiescent-stratified-basin test is expected
 * to expose the consequences.
 */
struct FVQual_EXPORT HaneyReport
{
    double  maxRatio     = 0.0;  //!< Largest steepness ratio found.
    int64_t worstEdge    = -1;   //!< Edge where it occurred (-1 when none).
    int64_t exceedances  = 0;    //!< Interior edges above the limit.
    double  limit        = 0.2;  //!< Limit the report was computed against.

    /*!
     * \brief True when no interior edge exceeded the limit.
     */
    [[nodiscard]] bool passed() const noexcept { return exceedances == 0; }

    /*!
     * \brief One-line human-readable summary.
     */
    [[nodiscard]] std::string summary() const;
};

/*!
 * \brief Terrain-following layered finite-volume mesh.
 */
class FVQual_EXPORT LayeredMesh
{
public:
    LayeredMesh() = default;

    /*!
     * \brief Builds the layered mesh over a horizontal mesh.
     *
     * Bed elevations come from the horizontal mesh's per-face bed. Surface
     * elevations are the initial free surface, one per column.
     *
     * \param[in]  horizontal       Horizontal mesh (must outlive this object).
     * \param[in]  spec             Vertical layering configuration.
     * \param[in]  surfaceElevation Initial free surface per column (m).
     * \param[out] message          Failure description when returning false.
     */
    [[nodiscard]] bool build(const HorizontalMesh &horizontal, const SigmaSpec &spec,
                             std::span<const double> surfaceElevation,
                             std::string &message);

    /*!
     * \brief Convenience overload: a spatially uniform initial water surface.
     */
    [[nodiscard]] bool build(const HorizontalMesh &horizontal, const SigmaSpec &spec,
                             double uniformSurfaceElevation, std::string &message);

    /*!
     * \brief Regenerates interface elevations, thicknesses, and volumes for a
     *        moved free surface. Layer count and topology never change.
     *
     * \param[in]  surfaceElevation New free surface per column (m).
     * \param[out] message          Failure description when returning false.
     */
    [[nodiscard]] bool updateSurface(std::span<const double> surfaceElevation,
                                     std::string &message);

    // ── Dimensions ────────────────────────────────────────────────────────

    [[nodiscard]] int layerCount() const noexcept { return m_layerCount; }

    /*!
     * \brief Number of water columns — one per horizontal face.
     */
    [[nodiscard]] int64_t columnCount() const noexcept { return m_columnCount; }

    [[nodiscard]] int64_t cellCount() const noexcept
    {
        return m_columnCount * static_cast<int64_t>(m_layerCount);
    }

    /*!
     * \brief Number of vertical interfaces per column (layerCount + 1).
     */
    [[nodiscard]] int interfaceCount() const noexcept { return m_layerCount + 1; }

    /*!
     * \brief Total interface slots across all columns.
     */
    [[nodiscard]] int64_t interfaceSlotCount() const noexcept
    {
        return m_columnCount * static_cast<int64_t>(interfaceCount());
    }

    // ── Indexing ──────────────────────────────────────────────────────────

    /*!
     * \brief Flat index of cell k in \p column (k = 0 at the surface).
     */
    [[nodiscard]] int64_t cell(int64_t column, int k) const noexcept
    {
        return column * static_cast<int64_t>(m_layerCount) + k;
    }

    /*!
     * \brief Flat index of interface k in \p column.
     */
    [[nodiscard]] int64_t interfaceSlot(int64_t column, int k) const noexcept
    {
        return column * static_cast<int64_t>(interfaceCount()) + k;
    }

    // ── Geometry ──────────────────────────────────────────────────────────

    /*!
     * \brief Interface elevations (m, positive up), column-major by
     *        interface: strictly decreasing within a column.
     */
    [[nodiscard]] const std::vector<double> &interfaceZ() const noexcept
    {
        return m_interfaceZ;
    }

    /*!
     * \brief Elevation of interface k in \p column (m).
     */
    [[nodiscard]] double z(int64_t column, int k) const noexcept
    {
        return m_interfaceZ[static_cast<size_t>(interfaceSlot(column, k))];
    }

    /*!
     * \brief Layer thicknesses (m), indexed by cell.
     */
    [[nodiscard]] const std::vector<double> &thickness() const noexcept
    {
        return m_thickness;
    }

    /*!
     * \brief Cell volumes (m³), indexed by cell. Prism: plan area x thickness.
     */
    [[nodiscard]] const std::vector<double> &cellVolume() const noexcept
    {
        return m_cellVolume;
    }

    /*!
     * \brief Total water depth per column (m).
     */
    [[nodiscard]] const std::vector<double> &columnDepth() const noexcept
    {
        return m_columnDepth;
    }

    /*!
     * \brief Free-surface elevation per column (m).
     */
    [[nodiscard]] const std::vector<double> &surfaceElevation() const noexcept
    {
        return m_surfaceElevation;
    }

    /*!
     * \brief Total water volume of the domain (m³).
     */
    [[nodiscard]] double totalVolume() const noexcept;

    /*!
     * \brief The horizontal mesh this layering sits on.
     */
    [[nodiscard]] const HorizontalMesh &horizontal() const noexcept
    {
        return *m_horizontal;
    }

    /*!
     * \brief Normalized sigma interfaces this mesh was built with.
     */
    [[nodiscard]] const std::vector<double> &sigmaInterfaces() const noexcept
    {
        return m_sigma;
    }

    /*!
     * \brief Screens the mesh for hydrostatic consistency (see HaneyReport).
     * \param[in] limit Steepness-ratio limit; 0.2 by convention.
     */
    [[nodiscard]] HaneyReport haneyReport(double limit = 0.2) const;

private:
    [[nodiscard]] bool regenerate(std::span<const double> surfaceElevation,
                                  std::string &message);

    const HorizontalMesh *m_horizontal = nullptr;

    int     m_layerCount  = 0;
    int64_t m_columnCount = 0;
    double  m_minLayerThickness = 0.0;

    std::vector<double> m_sigma;
    std::vector<double> m_interfaceZ;
    std::vector<double> m_thickness;
    std::vector<double> m_cellVolume;
    std::vector<double> m_columnDepth;
    std::vector<double> m_surfaceElevation;
};

} // namespace FVQual::Mesh

#endif // FVQual_MESH_LAYEREDMESH_H
