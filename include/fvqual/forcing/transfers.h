// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   transfers.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  Internal transfers between columns (F4.3).
 *
 * Water moved from one place in the model to another: a pumped-storage return,
 * a bypass around a dam segment, a destratification pump lifting hypolimnetic
 * water to the surface, a withdrawal tower discharging into another branch.
 *
 * **What makes a transfer more than two independent sources.** An inflow
 * registered by hand carries a composition the modeller states. A transfer
 * cannot: the water arriving at the receiving column *is* the water that left
 * the donating one, so its temperature and every constituent it carries are
 * determined by where in the donor's column it was drawn from. Pump 6 °C water
 * out of a hypolimnion and 6 °C water is what arrives, whatever the receiving
 * column happens to be doing. That coupling is the entire content of this file
 * — the discharges are bookkeeping.
 *
 * **Conservation is exact, not approximate.** The withdrawal removes mass
 * through the vertical weights the source distribution applied,
 * `dt · Q · Σ_k w_k c_k`. The inflow adds `dt · Q · c_in`. Setting `c_in` to the
 * mean through those *same cached* weights makes the two cancel identically, so
 * an internal transfer moves constituent mass without creating or destroying
 * any. Recomputing the weights against newer geometry instead would leave a
 * residual that no tolerance could justify — which is why SourceSet exposes
 * weightedMean() (cached) separately from withdrawalMean() (probed).
 *
 * \threadsafety One set per solver; the update methods mutate the paired
 *               sources.
 */

#ifndef FVQual_FORCING_TRANSFERS_H
#define FVQual_FORCING_TRANSFERS_H

#include "fvqual/forcing/sources.h"
#include "fvqual/fvqual_export.h"
#include "fvqual/mesh/layeredmesh.h"

#include <span>
#include <string>
#include <vector>

namespace FVQual::Forcing {

/*!
 * \brief One internal transfer between two columns.
 */
struct FVQual_EXPORT Transfer
{
    std::string id;

    /*!
     * \brief Column the water is taken from.
     */
    int64_t fromColumn = 0;

    /*!
     * \brief Column the water is delivered to.
     */
    int64_t toColumn = 0;

    /*!
     * \brief Transfer rate (m³/s, magnitude).
     */
    double discharge = 0.0;

    /*!
     * \brief How the withdrawal draws from the donating column. Defaults to
     *        selective withdrawal, since the intake is an outlet like any other
     *        and which layers it pulls from is what sets the delivered water.
     */
    Placement fromPlacement = Placement::SelectiveWithdrawal;

    /*!
     * \brief Intake elevation in the donating column (m).
     */
    double fromElevation = 0.0;

    /*!
     * \brief Intake width for a line-sink envelope (m); zero is a point sink.
     */
    double outletWidth = 0.0;

    /*!
     * \brief How the inflow enters the receiving column. Defaults to
     *        density-matched, because transferred water finds its own level:
     *        cold hypolimnetic water pumped into a warm column plunges rather
     *        than floating on top of it.
     */
    Placement toPlacement = Placement::DensityMatched;

    /*!
     * \brief Discharge elevation in the receiving column (m), for placements
     *        that use one.
     */
    double toElevation = 0.0;
};

/*!
 * \brief The internal transfers acting on a model.
 */
class FVQual_EXPORT TransferSet
{
public:
    /*!
     * \brief Registers a transfer and the pair of sources that carry it.
     *
     * \param[in]     transfer The transfer.
     * \param[in,out] sources  Source set to register the pair in.
     * \param[out]    message  Failure description when returning -1.
     * \return Its index, or -1 when the specification is invalid.
     */
    int add(const Transfer &transfer, SourceSet &sources, std::string &message);

    [[nodiscard]] int count() const noexcept
    {
        return static_cast<int>(m_transfers.size());
    }

    [[nodiscard]] const Transfer &transfer(int index) const noexcept
    {
        return m_transfers[static_cast<size_t>(index)];
    }

    /*!
     * \brief Sets a transfer rate (m³/s, magnitude).
     */
    [[nodiscard]] bool setDischarge(int index, double discharge, std::string &message);

    /*!
     * \brief Writes the transfer rates into the paired sources.
     *
     * Called once per internal step, before the sources are distributed, so the
     * withdrawal and the inflow are placed against the same state.
     */
    [[nodiscard]] bool updateDischarges(SourceSet &sources, std::string &message);

    /*!
     * \brief Gives each receiving inflow the water its donor is drawing.
     *
     * Called once per species at transport time, after distribute() has cached
     * the vertical weights and before the species is advected — the values must
     * be the ones the withdrawal will remove mass at, which are the current,
     * pre-advection ones.
     *
     * Cannot fail: the source indices it writes through are this class's own,
     * fixed when the transfer was registered. It is deliberately not given a
     * failure path, because the only honest thing a caller could do with one is
     * abandon a half-transported species.
     *
     * \param[in]     mesh             The layered mesh.
     * \param[in]     species          Species index being transported.
     * \param[in]     temperatureIndex Index of temperature, or -1.
     * \param[in]     values           Current per-cell values of that species.
     * \param[in,out] sources          Sources to write the composition into.
     */
    void carry(const Mesh::LayeredMesh &mesh, int species, int temperatureIndex,
               std::span<const double> values, SourceSet &sources);

    /*!
     * \brief Temperature the transfer last delivered (deg C).
     *
     * The diagnostic an operator of a destratification pump actually wants: not
     * how much water was moved, but how cold it was when it arrived.
     */
    [[nodiscard]] double deliveredTemperature(int index) const noexcept
    {
        return m_deliveredTemperature[static_cast<size_t>(index)];
    }

    /*!
     * \brief Index of the withdrawal source of a transfer.
     */
    [[nodiscard]] int intakeSource(int index) const noexcept
    {
        return m_intake[static_cast<size_t>(index)];
    }

    /*!
     * \brief Index of the inflow source of a transfer.
     */
    [[nodiscard]] int outfallSource(int index) const noexcept
    {
        return m_outfall[static_cast<size_t>(index)];
    }

private:
    std::vector<Transfer> m_transfers;
    std::vector<int>      m_intake;
    std::vector<int>      m_outfall;
    std::vector<double>   m_deliveredTemperature;
};

} // namespace FVQual::Forcing

#endif // FVQual_FORCING_TRANSFERS_H
