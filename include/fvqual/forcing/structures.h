// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   structures.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  Outlet works: weirs, gates, pipes, and pumps (F4.1).
 *
 * A `Source` says how much water leaves and from which depth. A structure says
 * *why* that number is what it is. The distinction matters because a real
 * outlet is not a prescribed discharge — it is a hole in a dam, and how much
 * goes through it depends on how deep the water is above it. Give a reservoir
 * a fixed release and it will happily drain below its own spillway crest;
 * give it a weir and it draws down and stops, because the head that drives the
 * flow disappears as the flow removes it.
 *
 * That feedback is the whole point. Structures compute a discharge from the
 * current water level each step and write it into the source set, which then
 * decides what depth the water actually comes from. The two compose: a gate
 * knows its rating curve, and selective withdrawal knows which layers feed it.
 *
 * **Direction.** Weirs, orifices, and pipes are outlet works and release only.
 * If the tailwater rises above the driving head they close rather than
 * reversing — a real submerged gate would admit flow backwards, but modelling
 * that means deciding what the returning water is made of, and a structure
 * that silently invents composition is worse than one that admits it stops.
 * Pumps may run in either direction, because a pumped inflow has a named
 * source whose properties the caller states.
 *
 * \threadsafety One set per solver; update() mutates the recorded discharges.
 */

#ifndef FVQual_FORCING_STRUCTURES_H
#define FVQual_FORCING_STRUCTURES_H

#include "fvqual/forcing/sources.h"
#include "fvqual/fvqual_export.h"
#include "fvqual/mesh/layeredmesh.h"

#include <span>
#include <string>
#include <vector>

namespace FVQual::Forcing {

/*!
 * \brief The rating law a structure follows.
 */
enum class StructureType
{
    /*!
     * \brief Sharp-crested weir: `Q = C L H^(3/2)`, reduced by Villemonte's
     *        submergence factor when the tailwater rises over the crest.
     */
    Weir,

    /*!
     * \brief Orifice or gate: `Q = C A g_open sqrt(2 g dH)`, where the driving
     *        head is measured to the orifice centre when the outlet discharges
     *        freely and to the tailwater when it is submerged.
     */
    Orifice,

    /*!
     * \brief Pipe or culvert: the same square-root law, always driven by the
     *        head difference, with the coefficient absorbing entrance and
     *        friction losses.
     */
    Pipe,

    /*!
     * \brief Pump: a prescribed capacity, switched on and off by water level.
     *        The only structure whose discharge does not depend on head.
     */
    Pump
};

/*!
 * \brief One hydraulic structure.
 */
struct FVQual_EXPORT Structure
{
    std::string   id;
    StructureType type   = StructureType::Weir;
    int64_t       column = 0;

    /*!
     * \brief Weir crest, orifice centre, or pipe invert elevation (m).
     */
    double crestElevation = 0.0;

    /*!
     * \brief Discharge coefficient. Zero selects a conventional default for
     *        the type: 1.84 for a sharp-crested weir, 0.62 for an orifice,
     *        0.80 for a pipe.
     */
    double coefficient = 0.0;

    /*!
     * \brief Weir crest length (m). Weirs only.
     */
    double width = 0.0;

    /*!
     * \brief Flow area (m²). Orifices and pipes only.
     */
    double area = 0.0;

    /*!
     * \brief Gate opening as a fraction in [0, 1]. Orifices and pipes only;
     *        this is the operational control an operator actually moves.
     */
    double opening = 1.0;

    /*!
     * \brief Downstream water-surface elevation (m). Below the crest the
     *        structure discharges freely.
     */
    double tailwaterElevation = -1.0e30;

    /*!
     * \brief Pump capacity (m³/s), signed: negative removes water from the
     *        reservoir, positive adds it. Pumps only.
     */
    double capacity = 0.0;

    /*!
     * \brief Release \c scheduledDischarge instead of the free rating.
     *
     * A structure left to its rating curve passes whatever the head drives
     * through it — correct for an ungated spillway, wrong for almost anything
     * an operator touches. Under a schedule the rating curve changes role
     * entirely: it stops being the discharge and becomes the *capacity*, the
     * most the head can deliver through that opening. The release is then the
     * scheduled demand, or the capacity, whichever is smaller. This is what
     * lets a multi-port outlet be operated to a target rather than merely
     * observed (see OutletGroupSet).
     */
    bool scheduled = false;

    /*!
     * \brief Demanded release magnitude when \c scheduled (m³/s).
     */
    double scheduledDischarge = 0.0;

    /*!
     * \brief Switch the structure on and off by water-surface elevation.
     */
    bool controlled = false;

    /*!
     * \brief Level at or above which a controlled structure starts (m).
     */
    double startElevation = 0.0;

    /*!
     * \brief Level at or below which it stops (m).
     *
     * Kept distinct from \c startElevation on purpose: with a single
     * threshold, a pump that lowers the level far enough to switch itself off
     * in one step will switch back on in the next and chatter, cycling every
     * step for the rest of the run. The gap between the two is the deadband
     * that prevents it, exactly as it does on real plant.
     */
    double stopElevation = 0.0;

    /*!
     * \brief How the discharge is spread over the column. Defaults to
     *        selective withdrawal, which is what an outlet does.
     */
    Placement placement = Placement::SelectiveWithdrawal;

    /*!
     * \brief Outlet width for a line-sink withdrawal envelope (m); zero
     *        treats it as a point sink.
     */
    double outletWidth = 0.0;

    /*!
     * \brief Temperature of pumped-in water (deg C). Pumps with positive
     *        capacity only.
     */
    double temperature = 10.0;

    /*!
     * \brief Concentration of pumped-in water, by species index.
     */
    std::vector<double> concentration;
};

/*!
 * \brief The structures acting on a model, and the sources they drive.
 */
class FVQual_EXPORT StructureSet
{
public:
    /*!
     * \brief Registers a structure and the source it will drive.
     *
     * The structure computes discharge; the source it creates in \p sources
     * carries that discharge into the solver at the right depth. They are
     * registered together so the pairing cannot come apart.
     *
     * \param[in]     structure The outlet works.
     * \param[in,out] sources   Source set to register the matching source in.
     * \param[out]    message   Failure description when returning -1.
     * \return Its index, or -1 when the specification is invalid.
     */
    int add(const Structure &structure, SourceSet &sources, std::string &message);

    [[nodiscard]] int count() const noexcept
    {
        return static_cast<int>(m_structures.size());
    }

    [[nodiscard]] const Structure &structure(int index) const noexcept
    {
        return m_structures[static_cast<size_t>(index)];
    }

    /*!
     * \brief Index, in the paired SourceSet, of the source this structure drives.
     */
    [[nodiscard]] int sourceIndex(int index) const noexcept
    {
        return m_sourceIndex[static_cast<size_t>(index)];
    }

    /*!
     * \brief Sets a gate opening in [0, 1] — the operational control.
     */
    [[nodiscard]] bool setOpening(int index, double opening, std::string &message);

    /*!
     * \brief Sets the downstream water level a structure discharges against.
     */
    [[nodiscard]] bool setTailwater(int index, double elevation, std::string &message);

    /*!
     * \brief Puts a structure on a scheduled release of \p discharge (m³/s).
     *
     * The rating curve becomes the capacity limit; see Structure::scheduled.
     * A negative discharge is refused rather than silently made positive,
     * because a caller asking an outlet to release −50 m³/s has confused a
     * magnitude with a direction and should hear about it.
     */
    [[nodiscard]] bool setSchedule(int index, double discharge, std::string &message);

    /*!
     * \brief Returns a structure to its free rating curve.
     */
    [[nodiscard]] bool clearSchedule(int index, std::string &message);

    /*!
     * \brief Recomputes every structure's discharge from the current water
     *        levels and writes the results into \p sources.
     *
     * Called by the solver on every *barotropic* substep, not once per internal
     * step. That placement is the whole of what makes the feedback work. A
     * free-surface-driven release is a barotropic process: the head that drives
     * it changes on the gravity-wave timescale, which mode splitting resolves
     * in the external mode at a few seconds per step while the internal step
     * may be several minutes long. Rating a spillway once per internal step
     * freezes its discharge at the head the step began with, and a reservoir
     * released for four minutes at its initial head draws down *below its own
     * crest* — a state the rating curve says is impossible. Re-rating each
     * substep lets the head fall as the water leaves, so the release decays
     * with it and the level asymptotes to the crest on its own, with no
     * limiter needed to hold it there.
     *
     * \param[in]     surfaceElevation Water-surface elevation per column (m).
     * \param[in,out] sources          Source set to update.
     * \param[out]    message          Failure description when returning false.
     */
    [[nodiscard]] bool update(std::span<const double> surfaceElevation,
                              SourceSet &sources, std::string &message);

    /*!
     * \brief Discharge computed for a structure on the last update (m³/s),
     *        signed as the source is: negative releases.
     */
    [[nodiscard]] double discharge(int index) const noexcept
    {
        return m_discharge[static_cast<size_t>(index)];
    }

    /*!
     * \brief Whether a controlled structure is currently running.
     */
    [[nodiscard]] bool isRunning(int index) const noexcept
    {
        return m_running[static_cast<size_t>(index)];
    }

    /*!
     * \brief Total discharge through all structures (m³/s).
     */
    [[nodiscard]] double totalDischarge() const noexcept;

    /*!
     * \brief Evaluates a rating law in isolation, with no state involved.
     *
     * Exposed because a rating curve is a thing engineers check against a
     * table, and being able to ask the model for `Q(H)` directly is worth more
     * than making them run a simulation to find out.
     *
     * \param[in] structure The structure.
     * \param[in] upstream  Upstream water-surface elevation (m).
     * \return Discharge magnitude (m³/s), always non-negative.
     */
    [[nodiscard]] static double rate(const Structure &structure, double upstream);

private:
    std::vector<Structure> m_structures;
    std::vector<int>       m_sourceIndex;
    std::vector<double>    m_discharge;
    std::vector<char>      m_running;
};

} // namespace FVQual::Forcing

#endif // FVQual_FORCING_STRUCTURES_H
