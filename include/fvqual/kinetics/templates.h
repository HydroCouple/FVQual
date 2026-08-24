// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   templates.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  Ready-made constituent models built on the reaction network (F6.4).
 *
 * Everything in Phase 6 up to this point is machinery: expressions, scopes,
 * light, settling. This is what the machinery is for. A modeller who wants a
 * dissolved-oxygen model should get one, not a blank network and a list of
 * things that would need writing.
 *
 * **These are builders, not a second engine.** Every template does exactly what
 * a user could do by hand with `addConstant`, `addRate` and `setVelocity` — the
 * rate laws it installs are readable strings, inspectable afterwards, and
 * removable. Nothing here is privileged, and no template can express anything a
 * hand-written network could not. That is the point: a template is a starting
 * position and an worked example, not a black box.
 *
 * **Two phases, because the lifecycle has two phases.** Species must exist
 * before the state store is allocated, and rate laws can only be compiled once
 * the solver is initialized. So each template registers its species first, then
 * configures the network after the solver exists. Trying to do both at once
 * would mean either allocating state twice or resolving names against a
 * registry that is still growing.
 *
 * **Units.** Time is seconds everywhere in the engine, but rate coefficients are
 * quoted per day in every reference a modeller will consult. Each template's
 * parameters are therefore **per day**, converted internally. The alternative —
 * SI in the struct — invites a silent factor of 86400, which is exactly the kind
 * of error that produces a plausible-looking wrong answer.
 */

#ifndef FVQual_KINETICS_TEMPLATES_H
#define FVQual_KINETICS_TEMPLATES_H

#include "fvqual/fvqual_export.h"
#include "fvqual/kinetics/reactionnetwork.h"
#include "fvqual/state/speciesregistry.h"

#include <string>

namespace FVQual::Hydro {
class InternalModeSolver;
}

namespace FVQual::Kinetics {

//! Seconds in a day, the conversion every template applies to its parameters.
inline constexpr double kSecondsPerDay = 86400.0;

/*!
 * \brief Water age: time since the parcel was last at the surface (days).
 *
 * Age is the cheapest diagnostic a reservoir model can carry and often the most
 * revealing — it shows residence time, short-circuiting, and the reach of an
 * inflow, without any water-quality parameters to argue about.
 *
 * **The surface condition is a relaxation, not a reset**, and that is a
 * deliberate choice rather than a shortcut. The definition says age is zero at
 * the surface; imposing that as a hard assignment inside a split operator would
 * make the surface layer's contribution depend on its thickness, and would hand
 * the integrator a discontinuity to resolve every step. A first-order
 * relaxation over `surfaceContactTime` is smooth, stable, and leaves the surface
 * age at approximately that timescale — so choose it short compared with the
 * ages of interest, which for a reservoir means hours against weeks.
 */
struct FVQual_EXPORT WaterAgeTemplate
{
    /*!
     * \brief Timescale over which surface water is reset toward zero age (days).
     *
     * The default of an hour leaves surface age at about 0.04 d, which is
     * negligible against residence times but far longer than a time step.
     */
    double surfaceContactTime = 1.0 / 24.0;

    //! Registers `__AGE__`. Call before the state store is allocated.
    [[nodiscard]] bool registerSpecies(State::SpeciesRegistry &registry,
                                       std::string &message) const;

    //! Installs the rate laws. Call after the solver is initialized.
    [[nodiscard]] bool configure(Hydro::InternalModeSolver &solver,
                                 std::string &message) const;
};

/*!
 * \brief A single constituent decaying at first order, with Arrhenius scaling.
 *
 * The workhorse: coliform die-off, a decaying tracer, a conservative substance
 * when the rate is zero. Temperature dependence is written as a rate law over
 * the temperature species rather than built in, which is what makes it visible
 * and editable afterwards.
 */
struct FVQual_EXPORT FirstOrderDecayTemplate
{
    //! Constituent id to register and decay.
    std::string id = "TRACER";

    //! Units recorded in the registry.
    std::string units = "g/m3";

    //! Initial concentration everywhere.
    double initialValue = 0.0;

    //! Decay coefficient at 20 °C (1/day). Zero makes the constituent
    //! conservative — expressed by installing no rate law at all, so the
    //! transport engine carries it untouched.
    double decayRatePerDay = 0.0;

    //! Arrhenius temperature multiplier per degree. 1.0 disables the
    //! dependence, and then the rate law does not reference temperature.
    double temperatureTheta = 1.047;

    //! Settling velocity (m/day); zero for a dissolved constituent.
    double settlingVelocityPerDay = 0.0;

    /*!
     * \brief Id of a benthic species to receive settled material.
     *
     * Empty means the settled mass leaves the model, which SettlingSet reports
     * rather than hides. Only meaningful when \ref settlingVelocityPerDay is
     * positive.
     */
    std::string benthicId;

    [[nodiscard]] bool registerSpecies(State::SpeciesRegistry &registry,
                                       std::string &message) const;

    [[nodiscard]] bool configure(Hydro::InternalModeSolver &solver,
                                 std::string &message) const;
};

} // namespace FVQual::Kinetics

#endif // FVQual_KINETICS_TEMPLATES_H
