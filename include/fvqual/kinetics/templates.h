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

/*!
 * \brief Dissolved oxygen with carbonaceous demand, nutrients and algae.
 *
 * The standard eutrophication model, in the form QUAL2K and CE-QUAL-W2 share:
 * oxygen driven by reaeration at the surface, carbonaceous decay and
 * nitrification in the water column, sediment demand at the bed, and algae that
 * produce oxygen while growing and consume it while respiring.
 *
 * **Stoichiometry is the whole game.** The reason this template exists rather
 * than a page of documentation telling people how to write it is that the
 * coefficients have to agree with each other or the budgets do not close.
 * Nitrogen is the clearest case: growth removes it from ammonium and nitrate,
 * respiration and mortality put it back, and nitrification moves it between the
 * two pools — and every one of those terms has to carry the same `g N per g
 * algae` or the model quietly manufactures nitrogen. The gates check the
 * budgets rather than the coefficients, because a budget that closes is
 * evidence and a coefficient that looks right is not.
 *
 * **Nutrient uptake is split without a regulariser.** Algae draw nitrogen from
 * ammonium and nitrate in proportion to their availability. Written naively
 * that needs `NH4 / (NH4 + NO3)`, which is 0/0 in water containing neither.
 * Instead the split reuses the half-saturation denominator that limitation
 * already carries — `NH4 / (Ksn + NH4 + NO3)` — so the two uptake terms sum to
 * exactly the nitrogen the growth term requires, the denominator can never
 * vanish, and nitrogen conservation is exact rather than exact-to-within-an-
 * epsilon.
 *
 * **What is deliberately not here.** Oxygen saturation is a constant, not a
 * function of temperature and elevation; organic nitrogen and phosphorus have
 * no separate pools, so mortality returns nutrients directly to the inorganic
 * ones; and light limitation is Monod rather than Steele, so there is no
 * photoinhibition. Each is expressible by hand on top of what the template
 * installs, and each would be a poor thing to hide inside a default.
 */
struct FVQual_EXPORT DissolvedOxygenTemplate
{
    // ── Identifiers ───────────────────────────────────────────────────────

    std::string oxygenId     = "DO";
    std::string cbodId       = "CBOD";
    std::string ammoniumId   = "NH4";
    std::string nitrateId    = "NO3";
    std::string phosphateId  = "PO4";
    std::string algaeId      = "ALGAE";

    /*!
     * \brief Benthic species to receive settled algae; empty to let it leave.
     *
     * SettlingSet reports what it loses either way, so an empty id is a
     * modelling choice rather than a silent sink.
     */
    std::string algaeBenthicId;

    // ── What to include ───────────────────────────────────────────────────

    //! Algae, and with them light limitation and oxygen production.
    bool includeAlgae = true;

    /*!
     * \brief Nitrogen and phosphorus pools.
     *
     * Without them the model is Streeter-Phelps with a sediment demand, which
     * is a legitimate and often sufficient model. Algal growth is then limited
     * by light alone.
     */
    bool includeNutrients = true;

    // ── Oxygen ────────────────────────────────────────────────────────────

    //! Saturation concentration (g O₂/m³). Constant; see the class note.
    double saturation = 9.0;

    /*!
     * \brief Surface transfer velocity for oxygen (m/day).
     *
     * A velocity rather than a first-order rate, because gas exchange is a
     * property of the surface: quoting it as `1/day` builds in an assumed
     * depth, and the model then reaerates differently for no physical reason
     * when the reservoir is drawn down.
     */
    double reaerationVelocityPerDay = 2.0;

    //! Sediment oxygen demand (g O₂/m²/day), applied at the bed.
    double sedimentOxygenDemandPerDay = 1.0;

    // ── Carbonaceous demand ───────────────────────────────────────────────

    //! CBOD decay coefficient at 20 °C (1/day).
    double cbodDecayPerDay = 0.23;

    // ── Nitrogen ──────────────────────────────────────────────────────────

    //! Nitrification, ammonium to nitrate, at 20 °C (1/day).
    double nitrificationPerDay = 0.1;

    //! Oxygen consumed per unit nitrogen nitrified (g O₂/g N).
    double oxygenPerNitrogen = 4.57;

    // ── Algae ─────────────────────────────────────────────────────────────

    double algalGrowthPerDay      = 2.0; //!< Maximum growth at 20 °C (1/day).
    double algalRespirationPerDay = 0.1; //!< Respiration at 20 °C (1/day).
    double algalMortalityPerDay   = 0.1; //!< Mortality to CBOD at 20 °C (1/day).
    double algalSettlingPerDay    = 0.1; //!< Settling velocity (m/day).

    //! Half-saturation irradiance for growth (W/m²).
    double lightHalfSaturation = 100.0;

    //! Half-saturation nitrogen, as ammonium plus nitrate (g N/m³).
    double nitrogenHalfSaturation = 0.014;

    //! Half-saturation phosphorus (g P/m³).
    double phosphorusHalfSaturation = 0.002;

    //! Specific attenuation of light by algae (m²/g), for self-shading.
    double algalLightAttenuation = 0.02;

    // ── Stoichiometry ─────────────────────────────────────────────────────

    double oxygenPerAlgae     = 1.4;   //!< g O₂ per g algal biomass.
    double nitrogenPerAlgae   = 0.08;  //!< g N per g algal biomass.
    double phosphorusPerAlgae = 0.011; //!< g P per g algal biomass.

    // ── Temperature ───────────────────────────────────────────────────────

    /*!
     * \brief Arrhenius multiplier per degree, applied to every rate.
     *
     * 1.0 removes the dependence entirely, and then no rate law references
     * temperature — which is what makes the template usable in a model that
     * does not register it.
     */
    double temperatureTheta = 1.047;

    // ── Initial conditions ────────────────────────────────────────────────

    double initialOxygen     = 9.0;
    double initialCbod       = 0.0;
    double initialAmmonium   = 0.0;
    double initialNitrate    = 0.0;
    double initialPhosphate  = 0.0;
    double initialAlgae      = 0.0;

    /*!
     * \brief Registers every constituent the configuration calls for.
     *
     * Call before the state store is allocated. Temperature is *not* registered
     * here: it belongs to the hydrodynamics, and a template that registered it
     * would fight with a model that already had.
     */
    [[nodiscard]] bool registerSpecies(State::SpeciesRegistry &registry,
                                       std::string &message) const;

    /*!
     * \brief Installs the rate laws, settling and light attenuation.
     */
    [[nodiscard]] bool configure(Hydro::InternalModeSolver &solver,
                                 std::string &message) const;
};

} // namespace FVQual::Kinetics

#endif // FVQual_KINETICS_TEMPLATES_H
