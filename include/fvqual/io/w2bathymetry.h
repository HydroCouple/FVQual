// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   w2bathymetry.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  Reading CE-QUAL-W2 bathymetry into a FVQual mesh (F8.1).
 *
 * W2 is the reason FVQual exists, and its bathymetry files are where thirty
 * years of reservoir schematisations live. Being able to read one is what turns
 * "the engine passes its analytical benchmarks" into "the engine reproduces a
 * real reservoir", which is a different and much harder claim.
 *
 * **The format.** A CSV with a `$` marker, a segment-index row, then rows of
 * `DLX` (segment lengths), `ELWS` (initial water surface), `PHI0` (orientation),
 * and `FRIC` (friction), followed by a `LAYERH` header and one row per layer:
 * the layer thickness, then that layer's width in every segment, then the layer
 * number. Layer 1 is the top of the grid.
 *
 * **The datum is not in the file.** The layer rows carry a `K` column and the
 * header names an `ELEV` column, but no example in the W2 distribution
 * populates it — elevations come from the control file's location card, as
 * `EBOT`, the elevation of the bottom of the grid. So it is a required
 * parameter here rather than something guessed from the data. For DeGray,
 * `EBOT = 66.35 m` with 36 two-meter layers puts the top of the grid at
 * 138.35 m and the initial surface at 123.80 m, giving 7.29e8 m³ over
 * 5.09e7 m² — which matches the reservoir's published figures and is how the
 * interpretation was checked.
 *
 * **How a W2 grid becomes a FVQual mesh.** W2 is laterally averaged: a segment
 * is one cell wide and its width varies by layer. FVQual is three-dimensional,
 * and since F1.5 its columns carry a plan area that varies with elevation — so
 * a segment maps to one column whose hypsographic curve is `width(z) x DLX`.
 * That preserves the storage the schematisation encodes without pretending the
 * reservoir is a box.
 *
 * **Where the curve is sampled matters.** W2's width is a step function of
 * elevation — constant within a layer — while a hypsographic curve is
 * piecewise linear. Sampling at layer *midpoints* is what makes the two agree:
 * the integral of a linear function over a layer equals its midpoint value
 * times the thickness, so each layer's volume comes out at exactly
 * `B x DLX x H`. Sampling at layer boundaries instead would average adjacent
 * widths and lose volume wherever the reservoir narrows.
 */

#ifndef FVQual_IO_W2BATHYMETRY_H
#define FVQual_IO_W2BATHYMETRY_H

#include "fvqual/fvqual_export.h"
#include "fvqual/mesh/layeredmesh.h"

#include "hydrocouplesdk/io/meshdefinition.h"

#include <iosfwd>
#include <string>
#include <vector>

namespace FVQual::IO {

/*!
 * \brief A CE-QUAL-W2 bathymetry file, and the mesh it describes.
 */
class FVQual_EXPORT W2Bathymetry
{
public:
    /*!
     * \brief Reads \p path.
     *
     * \param[in]  path            Bathymetry CSV.
     * \param[in]  bottomElevation `EBOT` from the control file's location card
     *                             (m). Required: it is not in the bathymetry.
     * \param[out] message         Failure description, with a line number.
     */
    [[nodiscard]] bool read(const std::string &path, double bottomElevation,
                            std::string &message);

    //! As above, from an open stream.
    [[nodiscard]] bool parse(std::istream &input, const std::string &origin,
                             double bottomElevation, std::string &message);

    // ── What the file said ────────────────────────────────────────────────

    [[nodiscard]] int segmentCount() const noexcept
    {
        return static_cast<int>(m_length.size());
    }

    [[nodiscard]] int layerCount() const noexcept
    {
        return static_cast<int>(m_layerHeight.size());
    }

    //! Segment lengths, `DLX` (m).
    [[nodiscard]] const std::vector<double> &segmentLength() const noexcept
    {
        return m_length;
    }

    //! Segment orientations, `PHI0` (radians clockwise from north).
    [[nodiscard]] const std::vector<double> &orientation() const noexcept
    {
        return m_orientation;
    }

    //! Bed friction, `FRIC` — Manning's n or a Chezy coefficient, per the model.
    [[nodiscard]] const std::vector<double> &friction() const noexcept
    {
        return m_friction;
    }

    //! Initial water-surface elevation per segment, `ELWS` (m).
    [[nodiscard]] const std::vector<double> &initialSurface() const noexcept
    {
        return m_initialSurface;
    }

    //! Layer thicknesses, top of the grid first (m).
    [[nodiscard]] const std::vector<double> &layerHeight() const noexcept
    {
        return m_layerHeight;
    }

    //! Width of \p layer in \p segment (m). Layer 0 is the top of the grid.
    [[nodiscard]] double width(int layer, int segment) const;

    //! Elevation of the bottom of \p layer (m).
    [[nodiscard]] double layerBottom(int layer) const;

    //! Elevation of the bottom of the grid, as supplied (m).
    [[nodiscard]] double bottomElevation() const noexcept { return m_bottom; }

    // ── Derived quantities, for checking against the source model ─────────

    /*!
     * \brief Water volume at a given surface elevation (m³).
     *
     * Computed the way W2 computes it — the sum of `B x DLX x H` over wetted
     * cells, with the layer containing the surface counted only up to it — so
     * it is an independent reference to gate the mesh against rather than a
     * restatement of it.
     */
    [[nodiscard]] double volumeAt(double surfaceElevation) const;

    //! Water-surface area at a given elevation (m²).
    [[nodiscard]] double surfaceAreaAt(double surfaceElevation) const;

    /*!
     * \brief Segments that carry water anywhere in the grid.
     *
     * W2 pads a branch with zero-width segments at each end. They are not part
     * of the reservoir and become no column here.
     */
    [[nodiscard]] const std::vector<int> &activeSegments() const noexcept
    {
        return m_active;
    }

    // ── Building a FVQual mesh ────────────────────────────────────────────

    /*!
     * \brief Lays the active segments out as a chain of quadrilateral columns.
     *
     * Segment centers follow the polyline implied by `DLX` and `PHI0`, and each
     * column is as wide as its segment is at the initial water surface. The
     * width at other elevations is carried by the hypsography rather than the
     * plan geometry.
     *
     * \param[out] definition Mesh for HorizontalMesh::build().
     * \param[out] bed        Bed elevation per column (m).
     * \param[out] message    Failure description when returning false.
     */
    [[nodiscard]] bool toMeshDefinition(HydroCouple::SDK::IO::MeshDefinition &definition,
                                        std::vector<double> &bed,
                                        std::string &message) const;

    /*!
     * \brief Gives every column its storage curve, sampled at layer midpoints.
     *
     * Call after LayeredMesh::build(), then updateSurface() to apply — see the
     * sequence note on LayeredMesh::setHypsography().
     */
    [[nodiscard]] bool applyHypsography(Mesh::LayeredMesh &mesh,
                                        std::string &message) const;

    /*!
     * \brief Lays the same reservoir out \p lateralCells wide instead of one.
     *
     * **Why this exists.** CE-QUAL-W2 is laterally averaged, so a model built
     * from it is one segment wide and every return path is forced through the
     * vertical. DeGray's residual interior warm bias of +3 to +4 °C is suspected
     * to be that schematisation rather than the physics: a synthetic basin under
     * the same wind transports 39% less heat below 20 m at eight cells wide.
     * That test used a made-up basin. This one uses DeGray.
     *
     * **No new bathymetry is invented.** `width(layer, segment)` already carries
     * the lateral information — it is the section's width at each elevation.
     * Cell `l` of a segment is given the lateral span
     *
     *     [ -H0 + 2 H0 l/L , -H0 + 2 H0 (l+1)/L ]
     *
     * where `H0` is the segment's half-width at its initial surface, and its
     * plan area at elevation `z` is `DLX` times that span's overlap with
     * `[-W(z)/2, +W(z)/2]`. The overlaps sum to `W(z)` exactly, so the expanded
     * mesh holds **the same volume at every stage** as the one-column chain. The
     * cross-section shape falls out of the curve: cells near the centerline stay
     * wet to the bed, outer cells dry out partway up, which is what a real
     * section does and what a single averaged column cannot represent.
     *
     * A column's bed is the lowest elevation at which its span still overlaps
     * the section at all.
     *
     * Columns are ordered segment-major: `segment * lateralCells + lateral`, so
     * `applyLateralHypsography()` must be given the same \p lateralCells.
     *
     * \param[in]  lateralCells Cells across the channel (>= 1). 1 reproduces
     *                          toMeshDefinition()'s topology.
     * \param[out] definition   Mesh for HorizontalMesh::build().
     * \param[out] bed          Bed elevation per column (m).
     * \param[out] message      Failure description when returning false.
     */
    [[nodiscard]] bool toLateralMeshDefinition(
        int lateralCells, HydroCouple::SDK::IO::MeshDefinition &definition,
        std::vector<double> &bed, std::string &message) const;

    /*!
     * \brief The storage curves for a mesh built by toLateralMeshDefinition().
     *
     * \param[in,out] mesh         Layered mesh built on that definition.
     * \param[in]     lateralCells The same count passed to the builder.
     * \param[out]    message      Failure description when returning false.
     */
    [[nodiscard]] bool applyLateralHypsography(Mesh::LayeredMesh &mesh,
                                               int lateralCells,
                                               std::string &message) const;

    /*!
     * \brief How many lanes each segment gets, and where its columns start.
     *
     * **Why a plan object rather than two parameters.** The mesh builder and the
     * hypsography have to agree exactly about lane counts and column ordering.
     * Passing the same two numbers to both invites them to disagree — and this
     * work has already produced one bug of precisely that shape. A plan computed
     * once and handed to both cannot.
     */
    struct LateralPlan
    {
        std::vector<int>     lanes;  //!< Lanes per active segment.
        std::vector<int64_t> offset; //!< First column of each segment.
        int64_t              columnCount = 0;

        [[nodiscard]] int laneCount(size_t segment) const
        {
            return lanes[segment];
        }
    };

    /*!
     * \brief Chooses lanes per segment so no cell is narrower than \p minCellWidth.
     *
     * **The problem this solves.** Refining every segment equally refines the
     * wrong place. DeGray's narrowest active segment is 79 m across while its
     * widest is 5530 m, so four lanes everywhere makes 19.8 m cells in the
     * tributary arms — against kilometer-long segments — and the barotropic
     * Courant limit collapses on cells that contribute nothing to the question
     * being asked. The overturning cell under investigation lives in the wide
     * water near the dam.
     *
     * A minimum cell width states the constraint directly: refine a segment only
     * as far as it can afford. Narrow arms keep one lane; the wide lower
     * reservoir gets the resolution.
     *
     * The resulting mesh is non-uniform but still **conforming**: where the lane
     * count changes, the shared boundary carries the union of both sides' lane
     * positions and the coarser cell is a polygon that walks all of them, so
     * every edge is shared by exactly two faces and no node hangs.
     *
     * \param[in] maxLanes     Ceiling on lanes per segment (>= 1).
     * \param[in] minCellWidth Narrowest lateral cell to allow (m).
     */
    [[nodiscard]] LateralPlan planLateralCells(int maxLanes, double minCellWidth) const;

    //! Builds the mesh a \ref LateralPlan describes.
    [[nodiscard]] bool toPlannedMeshDefinition(
        const LateralPlan &plan, HydroCouple::SDK::IO::MeshDefinition &definition,
        std::vector<double> &bed, std::string &message) const;

    //! The storage curves for a mesh built by toPlannedMeshDefinition().
    [[nodiscard]] bool applyPlannedHypsography(Mesh::LayeredMesh &mesh,
                                               const LateralPlan &plan,
                                               std::string &message) const;

private:
    /*!
     * \brief Width of \p segment at its initial water surface (m).
     *
     * The plan geometry the model runs at. Falls back to the widest layer when
     * the surface sits above every wetted one.
     */
    [[nodiscard]] double surfaceWidth(int segment) const;

    std::vector<double> m_length;
    std::vector<double> m_orientation;
    std::vector<double> m_friction;
    std::vector<double> m_initialSurface;
    std::vector<double> m_layerHeight;

    //! Width by [layer][segment], layer 0 at the top of the grid.
    std::vector<std::vector<double>> m_width;

    std::vector<int> m_active;

    double m_bottom = 0.0;
};

} // namespace FVQual::IO

#endif // FVQual_IO_W2BATHYMETRY_H
