// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   statestore.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  StateStore — the model's prognostic state in structure-of-arrays
 *         form (F1.3, §5.2).
 *
 * Every hot array is a flat `std::vector<double>` sized once at allocate()
 * and never resized during a run. That is what makes the state
 * device-mirrorable: a Kokkos view can wrap each block and stay valid for the
 * whole simulation. Public members rather than accessors, deliberately
 * — kernels index these directly and the indirection would buy nothing.
 *
 * Indexing conventions:
 * - Column arrays are indexed by column (the external/barotropic mode).
 * - Cell arrays are indexed by `mesh.cell(column, k)`.
 * - Interface arrays are indexed by `mesh.interfaceSlot(column, k)` and hold
 *   layerCount + 1 values per column (surface and bed included).
 * - Species are species-major: species s of cell i lives at
 *   `species[s * cellCount + i]`, so one species over the whole domain is a
 *   contiguous block — the layout the transport kernels and the reaction
 *   integrator both want.
 *
 * Vertical velocity `w` is carried as a full prognostic-capable array even
 * though the hydrostatic core only diagnoses it from continuity. That is the
 * non-hydrostatic seam: adding a pressure-correction module
 * later needs no change to the state layout.
 *
 * \threadsafety Plain data; callers coordinate access.
 */

#ifndef FVQual_STATE_STATESTORE_H
#define FVQual_STATE_STATESTORE_H

#include "fvqual/fvqual_export.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/state/speciesregistry.h"

#include "hydrocouplesdk/io/snapshot.h"

#include <cstdint>
#include <span>
#include <string>
#include <vector>

namespace FVQual::State {

/*!
 * \brief Prognostic and derived state of one model instance.
 */
class FVQual_EXPORT StateStore
{
public:
    StateStore() = default;

    /*!
     * \brief Sizes every block for a mesh and species set and applies each
     *        species' uniform initial value.
     *
     * \param[in] mesh     Layered mesh (must outlive this store).
     * \param[in] species  Registered species.
     */
    void allocate(const Mesh::LayeredMesh &mesh, const SpeciesRegistry &species);

    // ── Dimensions ────────────────────────────────────────────────────────

    [[nodiscard]] int64_t columnCount() const noexcept { return m_columnCount; }
    [[nodiscard]] int64_t cellCount() const noexcept { return m_cellCount; }
    [[nodiscard]] int64_t interfaceSlotCount() const noexcept { return m_interfaceSlots; }
    [[nodiscard]] int     speciesCount() const noexcept { return m_speciesCount; }
    [[nodiscard]] int     layerCount() const noexcept { return m_layerCount; }

    // ── External (barotropic) mode — per column ───────────────────────────

    std::vector<double> eta;  //!< Free-surface elevation (m).
    std::vector<double> ubar; //!< Depth-averaged x velocity (m/s).
    std::vector<double> vbar; //!< Depth-averaged y velocity (m/s).

    // ── Internal (baroclinic) mode — per cell ─────────────────────────────

    std::vector<double> u;       //!< x velocity (m/s).
    std::vector<double> v;       //!< y velocity (m/s).
    std::vector<double> density; //!< Water density (kg/m³).

    // ── Per interface ─────────────────────────────────────────────────────

    std::vector<double> w;  //!< Vertical velocity (m/s) — see the seam above.
    std::vector<double> kz; //!< Vertical eddy diffusivity for scalars (m²/s).
    std::vector<double> km; //!< Vertical eddy viscosity for momentum (m²/s).

    // ── Transported scalars — species-major over cells ────────────────────

    /*!
     * \brief All species concentrations: `species[s * cellCount + cell]`.
     */
    std::vector<double> species;

    /*!
     * \brief Mutable contiguous block of species \p s over every cell.
     */
    [[nodiscard]] std::span<double> speciesBlock(int s) noexcept
    {
        return {species.data() + static_cast<size_t>(s) * static_cast<size_t>(m_cellCount),
                static_cast<size_t>(m_cellCount)};
    }

    /*!
     * \brief Read-only contiguous block of species \p s over every cell.
     */
    [[nodiscard]] std::span<const double> speciesBlock(int s) const noexcept
    {
        return {species.data() + static_cast<size_t>(s) * static_cast<size_t>(m_cellCount),
                static_cast<size_t>(m_cellCount)};
    }

    /*!
     * \brief Reference to species \p s in \p cell.
     */
    [[nodiscard]] double &speciesAt(int s, int64_t cell) noexcept
    {
        return species[static_cast<size_t>(s) * static_cast<size_t>(m_cellCount) +
                       static_cast<size_t>(cell)];
    }

    [[nodiscard]] double speciesAt(int s, int64_t cell) const noexcept
    {
        return species[static_cast<size_t>(s) * static_cast<size_t>(m_cellCount) +
                       static_cast<size_t>(cell)];
    }

    // ── Diagnostics ───────────────────────────────────────────────────────

    /*!
     * \brief Volume-weighted total mass of species \p s (concentration x m³).
     *
     * The quantity the Eulerian engine must conserve to machine precision;
     * the transport cross-check gate is written against it.
     *
     * \param[in] s    Species index.
     * \param[in] mesh Mesh supplying the cell volumes.
     */
    [[nodiscard]] double totalMass(int s, const Mesh::LayeredMesh &mesh) const noexcept;

    /*!
     * \brief Volume-weighted mean of species \p s.
     */
    [[nodiscard]] double volumeWeightedMean(int s,
                                            const Mesh::LayeredMesh &mesh) const noexcept;

    // ── Snapshot capture ──────────────────────────────────────────────────

    /*!
     * \brief Materializes a per-cell field as a layered snapshot slice shaped
     *        { columnCount, layerCount } — the shape LayeredUgridWriter
     *        expects for a face-located layered variable.
     *
     * \param[in] itemId Field identity (also the output variable name).
     * \param[in] values Cell-indexed values; must have cellCount() entries.
     */
    [[nodiscard]] HydroCouple::SDK::IO::FieldSlice
    layeredSlice(const std::string &itemId, std::span<const double> values) const;

    /*!
     * \brief Materializes a per-column field as a snapshot slice shaped
     *        { columnCount }.
     */
    [[nodiscard]] HydroCouple::SDK::IO::FieldSlice
    columnSlice(const std::string &itemId, std::span<const double> values) const;

private:
    int64_t m_columnCount     = 0;
    int64_t m_cellCount       = 0;
    int64_t m_interfaceSlots  = 0;
    int     m_speciesCount    = 0;
    int     m_layerCount      = 0;
};

} // namespace FVQual::State

#endif // FVQual_STATE_STATESTORE_H
