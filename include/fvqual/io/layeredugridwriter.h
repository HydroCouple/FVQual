// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   layeredugridwriter.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  LayeredUgridWriter — CF-1.11 / UGRID-1.0 output with a vertical
 *         layer dimension.
 *
 * The SDK's UGRID writers serialize one value per node, edge, or face. A 3-D
 * layered model needs one value per face *per layer*, which UGRID expresses
 * through the layered-2D convention: 2-D mesh topology, plus a vertical
 * dimension, plus CF vertical-coordinate metadata tying the two together.
 * This is what Delft3D-FM writes and what xugrid and ugrid-checks understand,
 * and it is why FVQual does not attempt UGRID's full 3-D volume topology,
 * whose tool support is thin.
 *
 * The vertical coordinate is declared as CF `ocean_sigma_coordinate`:
 *
 *     z(t, k, f) = eta(t, f) + sigma(k) * (depth(f) + eta(t, f))
 *
 * so any CF-aware reader can reconstruct true cell elevations at every output
 * time without FVQual's help. Note the sign convention change: CF sigma runs
 * 0 at the surface to -1 at the bed and `depth` is positive down, whereas
 * FVQual's internal sigma runs 0 to 1 downward with elevations positive up.
 * The writer converts.
 *
 * Implements the SDK `IOutputWriter` contract, so it is driven by the SDK
 * `IOThread` — the simulation never blocks on IO and the writer never sees
 * concurrent calls.
 *
 * \threadsafety Driven exclusively by the IOThread.
 */

#ifndef FVQual_IO_LAYEREDUGRIDWRITER_H
#define FVQual_IO_LAYEREDUGRIDWRITER_H

#include "fvqual/fvqual_export.h"

#include "hydrocouplesdk/io/meshdefinition.h"
#include "hydrocouplesdk/io/outputwriter.h"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace FVQual::IO {

/*!
 * \brief Vertical-coordinate metadata needed to write CF formula terms.
 */
struct FVQual_EXPORT LayeredUgridMetadata
{
    /*!
     * \brief Normalized sigma interfaces in FVQual convention: layerCount + 1
     *        values, 0 at the surface increasing to 1 at the bed.
     */
    std::vector<double> sigmaInterfaces;

    /*!
     * \brief Bed elevation per face (m, positive up). Written as CF `depth`
     *        (positive down) after negation.
     */
    std::vector<double> bedElevation;

    /*!
     * \brief Layer count implied by the interfaces.
     */
    [[nodiscard]] int layerCount() const noexcept
    {
        return sigmaInterfaces.size() >= 2
                   ? static_cast<int>(sigmaInterfaces.size()) - 1
                   : 0;
    }
};

/*!
 * \brief Layered CF/UGRID NetCDF-4 writer.
 */
class FVQual_EXPORT LayeredUgridWriter final : public HydroCouple::SDK::IO::IOutputWriter
{
public:
    /*!
     * \brief Constructs the writer.
     * \param[in] path     Destination .nc file; parent directories created on
     *                     initialize().
     * \param[in] mesh     Horizontal UGRID mesh to serialize.
     * \param[in] vertical Sigma interfaces and bed elevations.
     */
    LayeredUgridWriter(std::filesystem::path path,
                       HydroCouple::SDK::IO::MeshDefinition mesh,
                       LayeredUgridMetadata vertical);

    ~LayeredUgridWriter() override;

    LayeredUgridWriter(const LayeredUgridWriter &)            = delete;
    LayeredUgridWriter &operator=(const LayeredUgridWriter &) = delete;

    /*!
     * \brief Registers a layered face field, shaped { face, layer } in the
     *        snapshot and written as { time, face, layer }.
     *
     * Must be called before initialize().
     *
     * \param[in] itemId       FieldSlice::itemId to match (also the variable
     *                         name).
     * \param[in] units        CF `units` ("" omits the attribute).
     * \param[in] standardName CF `standard_name` ("" omits the attribute).
     */
    void addLayeredField(const std::string &itemId, const std::string &units = "",
                         const std::string &standardName = "");

    /*!
     * \brief Registers a per-column (depth-independent) face field, shaped
     *        { face } in the snapshot and written as { time, face }.
     */
    void addColumnField(const std::string &itemId, const std::string &units = "",
                        const std::string &standardName = "");

    /*!
     * \brief Sets which registered column field carries the free surface used
     *        in the CF `formula_terms`. Defaults to "eta"; the field must be
     *        registered with addColumnField() or initialize() fails.
     */
    void setSurfaceFieldId(const std::string &itemId);

    // ── IOutputWriter ─────────────────────────────────────────────────────

    [[nodiscard]] std::string name() const override;
    [[nodiscard]] int         initialize() override;
    [[nodiscard]] int         write(const HydroCouple::SDK::IO::Snapshot &snapshot) override;
    [[nodiscard]] int         finalize() override;
    [[nodiscard]] std::string errorMessage() const override;

    [[nodiscard]] std::vector<HydroCouple::SDK::IO::ResultEntry> catalog() const override;

    /*!
     * \brief The destination path.
     */
    [[nodiscard]] const std::filesystem::path &path() const noexcept;

    /*!
     * \brief Whether this build can write NetCDF at all.
     */
    [[nodiscard]] static bool available() noexcept;

private:
    struct Impl; //!< Hides the NetCDF C API from the public surface.
    std::unique_ptr<Impl> m_impl;
};

} // namespace FVQual::IO

#endif // FVQual_IO_LAYEREDUGRIDWRITER_H
