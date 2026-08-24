// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/kinetics/templates.h"

#include "fvqual/hydro/internalmode.h"

#include <sstream>

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

} // namespace FVQual::Kinetics
