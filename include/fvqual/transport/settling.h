// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   settling.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  Downward settling of particulates, with a handoff to the bed (F6.2c).
 *
 * Algae sink, detritus sinks, suspended solids sink. In a reservoir this is not
 * a detail: it is how the epilimnion clears, how the light field recovers after
 * a bloom, and how organic matter arrives at the sediment to exert the oxygen
 * demand that F6.2b made expressible.
 *
 * **Why this is not part of the vertical advection field.** The mode-splitting
 * machinery diagnoses one vertical velocity per interface from continuity, and
 * the whole consistency argument — the residual that must vanish at the free
 * surface — depends on that field being shared by every species. A settling
 * velocity is *per species*, so it cannot live there. It is a separate operator,
 * which is also how CE-QUAL-W2 treats it.
 *
 * **Why it is solved implicitly.** Settling is downward-only, so upwinding takes
 * every cell's inflow from the cell above it, and the resulting matrix is lower
 * triangular. That means a single forward substitution from the surface down
 * solves it exactly — the same cost as an explicit update, with three
 * differences that all matter:
 *
 * - **No CFL constraint.** An explicit scheme needs `w dt < h`, and the binding
 *   case is a thin sigma layer at the surface rather than a fast settling
 *   velocity: a 0.1 m top layer and 10 m/d caps the step near 900 s, which is a
 *   constraint the hydrodynamics does not otherwise have.
 * - **Positivity, unconditionally.** Every coefficient in the substitution is
 *   positive, so no dt can drive a concentration negative. An explicit scheme
 *   that overshoots produces a negative particulate, which the kinetics will
 *   then turn into nonsense.
 * - **Exact conservation.** What leaves the water column telescopes to
 *   the material handed to the bed. For a prism this is `w dt C_bed` per
 *   unit area. Depth-varying storage also intercepts particles on projected
 *   banks; bank and bottom deposition share the column's aggregate bed store,
 *   expressed per unit plan-polygon area. Bank-specific sediment compartments
 *   and sloped-bed reaction exposure are not represented by that aggregate.
 *
 * The cost is numerical diffusion: first-order upwind in an implicit form is
 * more diffusive than the same scheme explicit. For a process whose velocity is
 * of order a metre a day acting on a smooth profile that is a small price, and
 * it buys away a stability limit and a whole class of negative concentrations.
 *
 * **Units at the bed.** Settled material accumulates in a benthic species
 * measured per unit **area** (g/m²), matching the areal convention that
 * RateScope::Bed uses. Sediment is a surface, not a volume, and a benthic store
 * expressed as a concentration in the bottom cell would change value whenever
 * the layering changed.
 *
 * \threadsafety One set per solver; apply() is internally parallel over
 *               columns, and each column's substitution is independent.
 */

#ifndef FVQual_TRANSPORT_SETTLING_H
#define FVQual_TRANSPORT_SETTLING_H

#include "fvqual/fvqual_export.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/state/speciesregistry.h"
#include "fvqual/state/statestore.h"

#include <string>
#include <vector>

namespace FVQual::Transport {

/*!
 * \brief Per-species settling velocities and their benthic destinations.
 */
class FVQual_EXPORT SettlingSet
{
public:
    /*!
     * \brief Sizes the set for \p species.
     */
    void initialize(const State::SpeciesRegistry &species);

    /*!
     * \brief Sets a species' settling velocity (m/s, downward).
     *
     * Zero — the default — means the species does not settle, so it costs
     * nothing to leave dissolved constituents alone.
     *
     * \param[in]  species  Species index; must be a transported species.
     * \param[in]  velocity Settling velocity (m/s), non-negative.
     * \param[out] message  Failure description when returning false.
     */
    [[nodiscard]] bool setVelocity(int species, double velocity,
                                   const State::SpeciesRegistry &registry,
                                   std::string &message);

    /*!
     * \brief Directs a settling species' bed flux into a benthic species.
     *
     * \param[in]  species Settling species index.
     * \param[in]  benthic Receiving species; must be SpeciesKind::MsxBenthic,
     *                     because a transported receiver would be advected
     *                     away again and the sediment would flow downstream.
     * \param[out] message Failure description when returning false.
     */
    [[nodiscard]] bool setBenthicReceiver(int species, int benthic,
                                          const State::SpeciesRegistry &registry,
                                          std::string &message);

    //! Settling velocity of \p species (m/s), or 0.
    [[nodiscard]] double velocity(int species) const;

    //! Receiving benthic species for \p species, or -1.
    [[nodiscard]] int benthicReceiver(int species) const;

    //! Whether any species settles.
    [[nodiscard]] bool any() const noexcept;

    /*!
     * \brief Settles every species by \p dt.
     *
     * \param[in]     mesh  Layered mesh, for thicknesses and plan areas.
     * \param[in,out] state State store; settling species and their benthic
     *                      receivers are updated.
     * \param[in]     dt    Interval (s).
     */
    void apply(const Mesh::LayeredMesh &mesh, State::StateStore &state, double dt);

    /*!
     * \brief Mass delivered to the bed for \p species since the last reset (g).
     *
     * Accumulated whether or not there is a benthic receiver. Without one the
     * mass leaves the model, and an unreported loss is indistinguishable from a
     * conservation bug — so it is reported.
     */
    [[nodiscard]] double deliveredToBed(int species) const;

    //! Mass delivered to the bed with nowhere to go (g).
    [[nodiscard]] double lostToBed(int species) const;

    void resetTotals();

private:
    std::vector<double> m_velocity;
    std::vector<int>    m_receiver;
    std::vector<double> m_delivered;
    std::vector<double> m_lost;
};

} // namespace FVQual::Transport

#endif // FVQual_TRANSPORT_SETTLING_H
