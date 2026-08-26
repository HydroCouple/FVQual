// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/forcing/forcingfield.h"

#include <utility>

namespace FVQual::Forcing {

bool ForcingField::setPerColumn(std::vector<double> values, int64_t columns,
                                std::string &message)
{
    if (columns <= 0)
    {
        message = "forcing field: the mesh reports no columns";
        return false;
    }

    if (static_cast<int64_t>(values.size()) != columns)
    {
        message = "forcing field: given " + std::to_string(values.size()) +
                  " values for " + std::to_string(columns) + " columns";
        return false;
    }

    // Commit only on success, so a refused call leaves the previous forcing
    // in place rather than a half-applied mixture of old and new.
    m_values = std::move(values);

    return true;
}

void MeteorologyField::setUniform(const Hydro::Meteorology &weather)
{
    base = weather;

    airTemperature.setUniform(weather.airTemperature);
    dewPoint.setUniform(weather.dewPoint);
    windSpeed.setUniform(weather.windSpeed);
    windDirection.setUniform(weather.windDirection);
    cloudCover.setUniform(weather.cloudCover);
    shortwave.setUniform(weather.shortwave);
    pressure.setUniform(weather.pressure);
    albedo.setUniform(weather.albedo);
}

Hydro::Meteorology MeteorologyField::at(int64_t column) const
{
    // Start from the coefficients, then overwrite everything spatial. Listing
    // the spatial members here rather than copying selectively is deliberate:
    // adding a variable to Meteorology and forgetting it here leaves it at the
    // base value, which is a defensible default rather than a wrong one.
    Hydro::Meteorology weather = base;

    weather.airTemperature = airTemperature.at(column);
    weather.dewPoint       = dewPoint.at(column);
    weather.windSpeed      = windSpeed.at(column);
    weather.windDirection  = windDirection.at(column);
    weather.cloudCover     = cloudCover.at(column);
    weather.shortwave      = shortwave.at(column);
    weather.pressure       = pressure.at(column);
    weather.albedo         = albedo.at(column);

    // surfaceAbsorption is intentionally left at its base value — see the note
    // in the header.
    return weather;
}

bool MeteorologyField::isUniform() const noexcept
{
    return perColumnCount() == 0;
}

int MeteorologyField::perColumnCount() const noexcept
{
    const ForcingField *const fields[] = {&airTemperature, &dewPoint,
                                          &windSpeed,      &windDirection,
                                          &cloudCover,     &shortwave,
                                          &pressure,       &albedo};

    int count = 0;

    for (const ForcingField *field : fields)
        count += field->isUniform() ? 0 : 1;

    return count;
}

} // namespace FVQual::Forcing
