// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   sources.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  Inflows and withdrawals, and the question of what depth they act at
 *         (F3.4 and Phase 4).
 *
 * Until now the model has been a closed box. A reservoir is not: it has
 * tributaries at the head, outlet works at the dam, and the whole point of
 * operating one is choosing *which water* leaves. That choice is vertical, and
 * it is why this file is mostly about placement rather than about bookkeeping.
 *
 * **Inflows plunge.** A cold tributary entering a warm reservoir does not mix
 * across the surface; it sinks until it finds water of its own density and
 * travels along that horizon as an interflow. Placing inflow at the surface
 * because that is where the river is would put a winter flood's cold, turbid
 * water directly into the epilimnion, which is both wrong and consequential —
 * it is the difference between a turbidity plume that reaches an intake and
 * one that slides underneath it.
 *
 * **Withdrawals draw from a layer, not a point.** Stratification confines the
 * water approaching an outlet to a horizontal band whose thickness depends on
 * the discharge and the strength of the density gradient. Selective withdrawal
 * — choosing an outlet elevation to control release temperature — is the
 * single most-used operational feature of CE-QUAL-W2, and it exists precisely
 * because that band is narrow enough to aim.
 *
 * Both are volumetric source terms on the continuity equation, so the
 * mass-consistency contract of the mode splitting extends to cover them: the
 * external mode's column source and the internal mode's cell sources must sum
 * to the same thing, or the free-surface continuity residual stops vanishing
 * and scalar transport silently stops conserving mass.
 *
 * \threadsafety Immutable after configuration except through distribute().
 */

#ifndef FVQual_FORCING_SOURCES_H
#define FVQual_FORCING_SOURCES_H

#include "fvqual/fvqual_export.h"
#include "fvqual/mesh/layeredmesh.h"

#include <span>
#include <string>
#include <vector>

namespace FVQual::Forcing {

/*!
 * \brief How a source's discharge is spread over the water column.
 */
enum class Placement
{
    /*!
     * \brief Spread over the full depth in proportion to layer thickness.
     */
    Uniform,

    /*!
     * \brief Concentrated in the surface layer — evaporation, precipitation,
     *        or a spillway.
     */
    Surface,

    /*!
     * \brief Concentrated in the bottom layer — a bottom outlet or seepage.
     */
    Bed,

    /*!
     * \brief Concentrated at a fixed elevation, regardless of what the density
     *        field is doing.
     */
    Elevation,

    /*!
     * \brief Inflow finds its own level: it descends until the ambient density
     *        matches its own, then spreads there. The physically right default
     *        for a tributary.
     */
    DensityMatched,

    /*!
     * \brief Withdrawal draws through an envelope centred on the outlet, whose
     *        half-height follows from the discharge and the local
     *        stratification. Outflows only.
     */
    SelectiveWithdrawal
};

/*!
 * \brief One inflow or withdrawal.
 */
struct FVQual_EXPORT Source
{
    /*!
     * \brief Column the source acts on.
     */
    int64_t column = 0;

    /*!
     * \brief Volumetric discharge (m³/s). Positive adds water, negative
     *        removes it.
     */
    double discharge = 0.0;

    Placement placement = Placement::Uniform;

    /*!
     * \brief Outlet or inlet elevation (m, positive up) for Elevation and
     *        SelectiveWithdrawal.
     */
    double elevation = 0.0;

    /*!
     * \brief Outlet width (m) for SelectiveWithdrawal. Zero or negative treats
     *        the outlet as a point sink; positive treats it as a line sink,
     *        which is what a long spillway crest or a slot gate behaves like.
     */
    double outletWidth = 0.0;

    /*!
     * \brief Effective horizontal angle of withdrawal (radians).
     *
     * The near-field geometry the outlet presents to the pool. A port on a wide
     * dam with a vertical face sees a half-space, θ = π, which is the default
     * and the case the W2 technical note documents; a port recessed in a
     * channel sees less. It enters the point-sink zone limit directly
     * (Smith et al. 1987, eq. 1a), so it is the lever that gives selective
     * withdrawal a weak two-dimensionality.
     */
    double withdrawalAngle = 3.14159265358979323846;

    /*!
     * \brief Inflow temperature (deg C), used both to set the incoming water's
     *        properties and to find its neutral-buoyancy depth. Ignored for
     *        withdrawals, which carry the ambient values by construction.
     */
    double temperature = 10.0;

    /*!
     * \brief Inflow concentration of each registered species, indexed by
     *        species. Entries beyond its length are treated as zero. For
     *        withdrawals this is ignored.
     */
    std::vector<double> concentration;

    /*!
     * \brief Identifier for reporting.
     */
    std::string id;
};

/*!
 * \brief The set of sources acting on a model, and their vertical placement.
 */
class FVQual_EXPORT SourceSet
{
public:
    /*!
     * \brief Registers a source.
     * \param[in]  source  The inflow or withdrawal.
     * \param[out] message Failure description when returning -1.
     * \return Its index, or -1 when the specification is invalid.
     */
    int add(const Source &source, std::string &message);

    [[nodiscard]] int count() const noexcept { return static_cast<int>(m_sources.size()); }

    [[nodiscard]] const Source &source(int index) const noexcept
    {
        return m_sources[static_cast<size_t>(index)];
    }

    /*!
     * \brief Updates the discharge of an existing source, which is how an
     *        operating schedule is applied mid-run.
     */
    [[nodiscard]] bool setDischarge(int index, double discharge, std::string &message);

    /*!
     * \brief Clears every source.
     */
    void clear() noexcept;

    /*!
     * \brief Total discharge added to the domain (m³/s), inflows less
     *        withdrawals — the quantity a volume budget must account for.
     */
    [[nodiscard]] double netDischarge() const noexcept;

    /*!
     * \brief Distributes every source over the cells of its column.
     *
     * Must be called before each step, because placement depends on the
     * present density field: a plunging inflow follows the thermocline as it
     * moves, and a withdrawal envelope widens as stratification weakens.
     *
     * \param[in]  mesh             Layered mesh with current geometry.
     * \param[in]  density          Cell densities (kg/m³).
     * \param[in]  gravity          Gravitational acceleration (m/s²).
     * \param[in]  referenceDensity Boussinesq reference density (kg/m³).
     * \param[out] cellDischarge    Per-cell volumetric source (m³/s), resized
     *                              and overwritten.
     * \param[out] columnDischarge  Per-column total (m³/s), resized and
     *                              overwritten. Equals the per-cell sum by
     *                              construction — the property the
     *                              mass-consistency contract rests on.
     */
    void distribute(const Mesh::LayeredMesh &mesh, std::span<const double> density,
                    double gravity, double referenceDensity,
                    std::vector<double> &cellDischarge,
                    std::vector<double> &columnDischarge);

    /*!
     * \brief Refreshes only the per-column totals from the current discharges.
     *
     * The external mode consumes column totals alone — the vertical placement
     * matters only to transport, which runs once per internal step. Structures
     * re-rated every barotropic substep therefore need this and not the full
     * distribute(), whose envelope and density work is both unnecessary at that
     * cadence and far too expensive to repeat a hundred times a step.
     *
     * \param[out] columnDischarge Per-column net discharge (m³/s); overwritten.
     */
    void columnTotals(std::span<double> columnDischarge) const;

    /*!
     * \brief Flux-weighted value a withdrawal is drawing, from cached weights.
     *
     * Uses the weights the last distribute() actually applied, rather than
     * recomputing them. That distinction is what makes an internal transfer
     * exactly conservative: the withdrawal removes mass through those weights,
     * so the composition handed to the receiving inflow has to be the mean
     * through the *same* weights. Recomputing against newer geometry would
     * leave a small residual that no tolerance could justify.
     *
     * \param[in] mesh   The layered mesh.
     * \param[in] index  Source index.
     * \param[in] values Per-cell field to average.
     * \return The weighted mean, or zero if the source moved no water.
     */
    [[nodiscard]] double weightedMean(const Mesh::LayeredMesh &mesh, int index,
                                      std::span<const double> values) const;

    /*!
     * \brief Sets the composition an inflow delivers for one species.
     *
     * \param[in]  index            Source index.
     * \param[in]  species          Species index.
     * \param[in]  temperatureIndex Index of temperature, or -1.
     * \param[in]  value            Value the inflow carries.
     * \param[out] message          Failure description when returning false.
     */
    [[nodiscard]] bool setInflowValue(int index, int species, int temperatureIndex,
                                      double value, std::string &message);

    /*!
     * \brief Flux-weighted value a withdrawal would draw at a stated discharge.
     *
     * The question multi-port operation turns on: if this port passed
     * \p probeDischarge, what temperature (or concentration) would come out of
     * it? It has to be asked of ports that are currently *shut*, which is why
     * the discharge is supplied rather than read from the source — a closed
     * outlet has no withdrawal envelope, so its own discharge cannot answer it.
     *
     * \param[in] mesh             The layered mesh.
     * \param[in] index            Source index.
     * \param[in] probeDischarge   Discharge to size the envelope with (m³/s).
     * \param[in] density          Cell densities (kg/m³).
     * \param[in] gravity          Gravitational acceleration (m/s²).
     * \param[in] referenceDensity Boussinesq reference density (kg/m³).
     * \param[in] values           Per-cell field to average.
     * \return The envelope-weighted mean of \p values.
     */
    [[nodiscard]] double withdrawalMean(const Mesh::LayeredMesh &mesh, int index,
                                        double probeDischarge,
                                        std::span<const double> density, double gravity,
                                        double referenceDensity,
                                        std::span<const double> values) const;

    /*!
     * \brief Discharge-weighted inflow concentration of one species in each
     *        cell, for the cells inflow reaches; zero elsewhere.
     *
     * Valid after distribute(). This is what the transport needs to mix an
     * arriving load in correctly: the source term in the scalar equation is
     * `S * c_source`, and for an outflow `c_source` is the ambient value,
     * which is why only inflows appear here.
     *
     * \param[in]  mesh             Layered mesh.
     * \param[in]  species          Species index.
     * \param[in]  temperatureIndex Index of the reserved temperature species,
     *                              or -1; temperature is carried on the
     *                              source's own field rather than in its
     *                              concentration vector.
     * \param[out] out              Per-cell concentration, resized.
     */
    void inflowConcentration(const Mesh::LayeredMesh &mesh, int species,
                             int temperatureIndex, std::vector<double> &out) const;

    /*!
     * \brief Half-height of the withdrawal envelope computed for each source
     *        on the last distribute() (m); zero for sources that are not
     *        selective withdrawals. Exposed because it is the number an
     *        operator reasons about when choosing a port.
     */
    [[nodiscard]] const std::vector<double> &envelopeHalfHeight() const noexcept
    {
        return m_envelope;
    }

    /*!
     * \brief Elevation each source was actually placed at on the last
     *        distribute() (m) — for a plunging inflow this is the depth it
     *        found, which is a diagnostic worth watching.
     */
    [[nodiscard]] const std::vector<double> &placedElevation() const noexcept
    {
        return m_placed;
    }

private:
    /*!
     * \brief Fills \p weights (one per layer, summing to one) for one source.
     */
    /*!
     * \brief Vertical weights for a source, free of any recorded state.
     */
    void weightsCore(const Source &source, const Mesh::LayeredMesh &mesh,
                     std::span<const double> density, double gravity,
                     double referenceDensity, std::vector<double> &weights,
                     double &envelope, double &placed) const;

    void weightsFor(const Source &source, const Mesh::LayeredMesh &mesh,
                    std::span<const double> density, double gravity,
                    double referenceDensity, int index, std::vector<double> &weights);

    std::vector<Source> m_sources;
    std::vector<double> m_envelope;
    std::vector<double> m_placed;
    std::vector<double> m_weights;
    std::vector<double> m_sourceWeights;
};

} // namespace FVQual::Forcing

#endif // FVQual_FORCING_SOURCES_H
