// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/kinetics/templates.h"

#include "fvqual/hydro/internalmode.h"

#include <sstream>
#include <vector>

namespace FVQual::Kinetics {

namespace {

/*!
 * \brief Formats \p value with enough digits to survive the round trip.
 *
 * Rate coefficients become text and are parsed back, so the default six
 * significant figures would quietly perturb them.
 */
std::string number(double value)
{
    std::ostringstream out;
    out.precision(17);
    out << value;
    return out.str();
}

} // namespace

// ── Water age ─────────────────────────────────────────────────────────────

bool WaterAgeTemplate::registerSpecies(State::SpeciesRegistry &registry,
                                       std::string &message) const
{
    return registry.addAge(message) >= 0;
}

bool WaterAgeTemplate::configure(Hydro::InternalModeSolver &solver,
                                 std::string &message) const
{
    if (!(surfaceContactTime > 0.0))
    {
        message = "the surface contact timescale must be positive";
        return false;
    }

    ReactionNetwork &network = solver.reactions();

    const std::string age = State::SpeciesRegistry::kAge;

    // One day of age per day elapsed. In SI that is 1/86400 per second, and
    // writing it as the reciprocal rather than a decimal keeps the intent
    // legible in the rate law itself.
    if (!network.addRate(age, "1 / " + number(kSecondsPerDay), message))
        return false;

    // Relaxation toward zero at the surface. Areal, like every surface law, so
    // the engine divides by the surface cell's thickness — which means the
    // written rate has to be multiplied by that thickness to leave a
    // thickness-independent relaxation. The alternative is a reset whose
    // strength depends on how the column happens to be sliced.
    const double rate = 1.0 / (surfaceContactTime * kSecondsPerDay);

    const std::string relaxation =
        "0 - " + number(rate) + " * " + age + " * " + ReactionNetwork::kThickness;

    return network.addRate(age, relaxation, RateScope::Surface, message);
}

// ── First-order decay ─────────────────────────────────────────────────────

bool FirstOrderDecayTemplate::registerSpecies(State::SpeciesRegistry &registry,
                                              std::string &message) const
{
    if (id.empty())
    {
        message = "a constituent needs an id";
        return false;
    }

    State::SpeciesDefinition definition;
    definition.id           = id;
    definition.units        = units;
    definition.initialValue = initialValue;

    if (registry.add(definition, message) < 0)
        return false;

    if (benthicId.empty())
        return true;

    State::SpeciesDefinition bed;
    bed.id    = benthicId;
    bed.kind  = State::SpeciesKind::MsxBenthic;
    bed.units = "g/m2";

    return registry.add(bed, message) >= 0;
}

bool FirstOrderDecayTemplate::configure(Hydro::InternalModeSolver &solver,
                                        std::string &message) const
{
    if (decayRatePerDay < 0.0 || settlingVelocityPerDay < 0.0)
    {
        message = "decay and settling rates must be non-negative";
        return false;
    }

    if (decayRatePerDay > 0.0)
    {
        const double rate = decayRatePerDay / kSecondsPerDay;

        // Temperature dependence only appears when it is switched on. A
        // theta of exactly one would multiply by 1 forever, and worse, would
        // make every such model depend on temperature being registered.
        std::string law = "0 - " + number(rate) + " * " + id;

        if (temperatureTheta != 1.0)
        {
            law = "0 - " + number(rate) + " * " + number(temperatureTheta) + " ^ (" +
                  State::SpeciesRegistry::kTemperature + " - 20) * " + id;
        }

        if (!solver.reactions().addRate(id, law, message))
            return false;
    }

    if (settlingVelocityPerDay <= 0.0)
        return true;

    // The registry is the solver's, so the indices resolved here are the ones
    // the state store was allocated with.
    const int species = solver.speciesIndex(id);

    if (species < 0)
    {
        message = "constituent '" + id + "' is not registered";
        return false;
    }

    if (!solver.settling().setVelocity(species,
                                       settlingVelocityPerDay / kSecondsPerDay,
                                       solver.species(), message))
    {
        return false;
    }

    if (benthicId.empty())
        return true;

    const int bed = solver.speciesIndex(benthicId);

    if (bed < 0)
    {
        message = "benthic species '" + benthicId + "' is not registered";
        return false;
    }

    return solver.settling().setBenthicReceiver(species, bed, solver.species(), message);
}

// ── Dissolved oxygen suite ────────────────────────────────────────────────

namespace {

/*!
 * \brief Registers one constituent, or reports why it could not be.
 */
bool addConstituent(State::SpeciesRegistry &registry, const std::string &id,
                    const std::string &units, double initialValue,
                    std::string &message)
{
    if (id.empty())
    {
        message = "a constituent needs an id";
        return false;
    }

    State::SpeciesDefinition definition;
    definition.id           = id;
    definition.units        = units;
    definition.initialValue = initialValue;

    return registry.add(definition, message) >= 0;
}

//! Wraps \p text in parentheses, so terms compose without precedence surprises.
std::string group(const std::string &text)
{
    return "(" + text + ")";
}

//! Sums a list of signed terms into one rate expression, or "0" if empty.
std::string sum(const std::vector<std::string> &terms)
{
    if (terms.empty())
        return "0";

    std::string out = terms.front();

    for (size_t i = 1; i < terms.size(); ++i)
        out += " + " + terms[i];

    return out;
}

} // namespace

bool DissolvedOxygenTemplate::registerSpecies(State::SpeciesRegistry &registry,
                                              std::string &message) const
{
    if (!addConstituent(registry, oxygenId, "g/m3", initialOxygen, message))
        return false;

    if (!addConstituent(registry, cbodId, "g/m3", initialCbod, message))
        return false;

    if (includeNutrients)
    {
        if (!addConstituent(registry, ammoniumId, "g/m3", initialAmmonium, message))
            return false;

        if (!addConstituent(registry, nitrateId, "g/m3", initialNitrate, message))
            return false;

        if (!addConstituent(registry, phosphateId, "g/m3", initialPhosphate, message))
            return false;
    }

    if (includeAlgae)
    {
        if (!addConstituent(registry, algaeId, "g/m3", initialAlgae, message))
            return false;

        if (!algaeBenthicId.empty())
        {
            State::SpeciesDefinition bed;
            bed.id    = algaeBenthicId;
            bed.kind  = State::SpeciesKind::MsxBenthic;
            bed.units = "g/m2";

            if (registry.add(bed, message) < 0)
                return false;
        }
    }

    message.clear();
    return true;
}

bool DissolvedOxygenTemplate::configure(Hydro::InternalModeSolver &solver,
                                        std::string &message) const
{
    if (saturation < 0.0 || reaerationVelocityPerDay < 0.0 ||
        sedimentOxygenDemandPerDay < 0.0 || cbodDecayPerDay < 0.0 ||
        nitrificationPerDay < 0.0 || algalGrowthPerDay < 0.0 ||
        algalRespirationPerDay < 0.0 || algalMortalityPerDay < 0.0 ||
        algalSettlingPerDay < 0.0)
    {
        message = "rate coefficients must be non-negative";
        return false;
    }

    if (includeAlgae && (lightHalfSaturation <= 0.0 ||
                         (includeNutrients && (nitrogenHalfSaturation <= 0.0 ||
                                               phosphorusHalfSaturation <= 0.0))))
    {
        // Positive half-saturation constants are what keep every limitation
        // denominator away from zero, so this is a well-posedness requirement
        // rather than a taste in parameters.
        message = "half-saturation constants must be positive";
        return false;
    }

    ReactionNetwork &network = solver.reactions();

    // The Arrhenius factor, or nothing at all when the dependence is switched
    // off — so a model without a temperature species still compiles.
    const std::string arrhenius =
        temperatureTheta == 1.0
            ? std::string()
            : " * " + number(temperatureTheta) + " ^ " +
                  group(std::string(State::SpeciesRegistry::kTemperature) + " - 20");

    const double decay   = cbodDecayPerDay / kSecondsPerDay;
    const double nitrify = nitrificationPerDay / kSecondsPerDay;

    // ── Shared sub-expressions ────────────────────────────────────────────

    // Oxygen demand exerted by carbonaceous material, in g O₂/m³/s.
    const std::string cbodOxidation = number(decay) + " * " + cbodId + arrhenius;

    // Nitrogen oxidised, in g N/m³/s.
    const std::string nitrification =
        includeNutrients ? number(nitrify) + " * " + ammoniumId + arrhenius
                         : std::string();

    // Algal terms. `base` is growth without its nitrogen limitation, because
    // the nitrogen uptake split below needs that factor separated out — see the
    // class note.
    std::string base;
    std::string growth;
    std::string respiration;
    std::string mortality;
    std::string nitrogenDenominator;

    if (includeAlgae)
    {
        const double mu   = algalGrowthPerDay / kSecondsPerDay;
        const double resp = algalRespirationPerDay / kSecondsPerDay;
        const double mort = algalMortalityPerDay / kSecondsPerDay;

        const std::string light =
            group(std::string(ReactionNetwork::kLight) + " / " +
                  group(number(lightHalfSaturation) + " + " + ReactionNetwork::kLight));

        base = number(mu) + arrhenius + " * " + light;

        if (includeNutrients)
        {
            nitrogenDenominator = group(number(nitrogenHalfSaturation) + " + " +
                                        ammoniumId + " + " + nitrateId);

            const std::string phosphorus =
                group(phosphateId + " / " +
                      group(number(phosphorusHalfSaturation) + " + " + phosphateId));

            base += " * " + phosphorus;
        }

        base += " * " + algaeId;

        growth = includeNutrients ? group(base) + " * " +
                                        group(group(ammoniumId + " + " + nitrateId) +
                                              " / " + nitrogenDenominator)
                                  : base;

        respiration = number(resp) + " * " + algaeId + arrhenius;
        mortality   = number(mort) + " * " + algaeId + arrhenius;
    }

    // ── Dissolved oxygen ──────────────────────────────────────────────────

    {
        std::vector<std::string> terms;

        terms.push_back("0 - " + group(cbodOxidation));

        if (includeNutrients)
        {
            terms.push_back("0 - " + number(oxygenPerNitrogen) + " * " +
                            group(nitrification));
        }

        if (includeAlgae)
        {
            terms.push_back(number(oxygenPerAlgae) + " * " + group(growth));
            terms.push_back("0 - " + number(oxygenPerAlgae) + " * " +
                            group(respiration));
        }

        if (!network.addRate(oxygenId, sum(terms), message))
            return false;

        if (reaerationVelocityPerDay > 0.0)
        {
            // Areal, so the velocity is used directly: m/s x g/m³ = g/m²/s, and
            // the engine divides by the surface cell's thickness.
            const std::string exchange =
                number(reaerationVelocityPerDay / kSecondsPerDay) + " * " +
                group(number(saturation) + " - " + oxygenId);

            if (!network.addRate(oxygenId, exchange, RateScope::Surface, message))
                return false;
        }

        if (sedimentOxygenDemandPerDay > 0.0)
        {
            const std::string demand =
                "0 - " + number(sedimentOxygenDemandPerDay / kSecondsPerDay) + arrhenius;

            if (!network.addRate(oxygenId, demand, RateScope::Bed, message))
                return false;
        }
    }

    // ── Carbonaceous demand ───────────────────────────────────────────────

    {
        std::vector<std::string> terms;

        terms.push_back("0 - " + group(cbodOxidation));

        if (includeAlgae)
        {
            // Dead algae become oxygen demand, in the oxygen equivalents CBOD
            // is measured in.
            terms.push_back(number(oxygenPerAlgae) + " * " + group(mortality));
        }

        if (!network.addRate(cbodId, sum(terms), message))
            return false;
    }

    // ── Nutrients ─────────────────────────────────────────────────────────

    if (includeNutrients)
    {
        std::vector<std::string> ammonium;
        std::vector<std::string> nitrate;
        std::vector<std::string> phosphate;

        ammonium.push_back("0 - " + group(nitrification));
        nitrate.push_back(group(nitrification));

        if (includeAlgae)
        {
            // The uptake split. Both terms share the limitation denominator, so
            // they sum to exactly `nitrogenPerAlgae * growth` — which is what
            // makes the nitrogen budget close to round-off rather than to the
            // size of a regulariser.
            const std::string fromAmmonium = number(nitrogenPerAlgae) + " * " +
                                             group(base) + " * " +
                                             group(ammoniumId + " / " +
                                                   nitrogenDenominator);

            const std::string fromNitrate = number(nitrogenPerAlgae) + " * " +
                                            group(base) + " * " +
                                            group(nitrateId + " / " +
                                                  nitrogenDenominator);

            ammonium.push_back("0 - " + group(fromAmmonium));
            nitrate.push_back("0 - " + group(fromNitrate));

            // Respiration and mortality return nitrogen to the ammonium pool.
            // There are no organic pools here, so it returns directly; see the
            // class note.
            ammonium.push_back(number(nitrogenPerAlgae) + " * " +
                               group(respiration + " + " + mortality));

            phosphate.push_back("0 - " + number(phosphorusPerAlgae) + " * " +
                                group(growth));
            phosphate.push_back(number(phosphorusPerAlgae) + " * " +
                                group(respiration + " + " + mortality));
        }

        if (!network.addRate(ammoniumId, sum(ammonium), message))
            return false;

        if (!network.addRate(nitrateId, sum(nitrate), message))
            return false;

        if (!phosphate.empty() && !network.addRate(phosphateId, sum(phosphate), message))
            return false;
    }

    // ── Algae ─────────────────────────────────────────────────────────────

    if (includeAlgae)
    {
        std::vector<std::string> terms;

        terms.push_back(group(growth));
        terms.push_back("0 - " + group(respiration));
        terms.push_back("0 - " + group(mortality));

        if (!network.addRate(algaeId, sum(terms), message))
            return false;

        const int algae = solver.speciesIndex(algaeId);

        if (algae < 0)
        {
            message = "algal species '" + algaeId + "' is not registered";
            return false;
        }

        if (!solver.light().setSpecificAttenuation(algae, algalLightAttenuation,
                                                   message))
        {
            return false;
        }

        if (algalSettlingPerDay > 0.0)
        {
            if (!solver.settling().setVelocity(algae,
                                               algalSettlingPerDay / kSecondsPerDay,
                                               solver.species(), message))
            {
                return false;
            }

            if (!algaeBenthicId.empty())
            {
                const int bed = solver.speciesIndex(algaeBenthicId);

                if (bed < 0)
                {
                    message = "benthic species '" + algaeBenthicId +
                              "' is not registered";
                    return false;
                }

                if (!solver.settling().setBenthicReceiver(algae, bed, solver.species(),
                                                          message))
                {
                    return false;
                }
            }
        }
    }

    message.clear();
    return true;
}

} // namespace FVQual::Kinetics
