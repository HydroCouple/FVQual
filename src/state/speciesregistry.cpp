// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/state/speciesregistry.h"

namespace FVQual::State {

int SpeciesRegistry::add(const SpeciesDefinition &definition, std::string &message)
{
    if (definition.id.empty())
    {
        message = "species id must not be empty";
        return -1;
    }

    if (indexOf(definition.id) >= 0)
    {
        message = "species '" + definition.id + "' is already registered";
        return -1;
    }

    if (definition.kind == SpeciesKind::Temperature && m_temperatureIndex >= 0)
    {
        message = "temperature is already registered";
        return -1;
    }

    if (definition.kind == SpeciesKind::Age && m_ageIndex >= 0)
    {
        message = "water age is already registered";
        return -1;
    }

    const int index = static_cast<int>(m_definitions.size());
    m_definitions.push_back(definition);

    if (definition.kind == SpeciesKind::Temperature)
        m_temperatureIndex = index;
    else if (definition.kind == SpeciesKind::Age)
        m_ageIndex = index;

    message.clear();
    return index;
}

int SpeciesRegistry::addTemperature(double initialCelsius, std::string &message)
{
    SpeciesDefinition definition;
    definition.id           = kTemperature;
    definition.kind         = SpeciesKind::Temperature;
    definition.units        = "degC";
    definition.standardName = "sea_water_temperature";
    definition.initialValue = initialCelsius;

    return add(definition, message);
}

int SpeciesRegistry::addAge(std::string &message)
{
    SpeciesDefinition definition;
    definition.id           = kAge;
    definition.kind         = SpeciesKind::Age;
    definition.units        = "d";
    definition.standardName = "sea_water_age_since_surface_contact";
    definition.initialValue = 0.0;

    return add(definition, message);
}

int SpeciesRegistry::indexOf(const std::string &id) const noexcept
{
    for (size_t i = 0; i < m_definitions.size(); ++i)
    {
        if (m_definitions[i].id == id)
            return static_cast<int>(i);
    }

    return -1;
}

} // namespace FVQual::State
