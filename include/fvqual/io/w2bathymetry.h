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
 * `EBOT = 66.35 m` with 36 two-metre layers puts the top of the grid at
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
     * Segment centres follow the polyline implied by `DLX` and `PHI0`, and each
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

private:
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
