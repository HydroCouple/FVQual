// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/hydro/heatbudget.h"
#include "fvqual/hydro/columnphysics.h"

#include <algorithm>
#include <cmath>

namespace FVQual::Hydro {

double saturationVapourPressure(double celsius) noexcept
{
    return ColumnPhysics::vapour(celsius);
}

HeatFluxTerms surfaceHeatFlux(const Meteorology &weather,
                              double surfaceTemperature) noexcept
{
    const auto terms=ColumnPhysics::surface(weather,surfaceTemperature,weather.windSpeed,
        weather.airTemperature,weather.dewPoint,weather.shortwave,weather.cloudCover,
        weather.evaporationA,weather.evaporationB);
    return {terms.netShortwave,terms.surfaceSolar,terms.penetratingSolar,terms.longwaveDown,
            terms.longwaveUp,terms.latent,terms.sensible};
}

double WindStress::magnitude() const noexcept
{
    return std::sqrt(x * x + y * y);
}

WindStress windStress(const Meteorology &weather) noexcept
{
    const auto stress=ColumnPhysics::wind(weather,weather.windSpeed,weather.windDirection,weather.dragCoefficient);
    return {stress.first,stress.second};
}

namespace {

constexpr double kPi = 3.14159265358979323846;

[[nodiscard]] double radians(double degrees) noexcept
{
    return degrees * kPi / 180.0;
}

} // namespace

SolarPosition solarPosition(double julianDay, double latitude) noexcept
{
    SolarPosition position;

    const double day = std::floor(julianDay);

    // Cooper's declination. Good to about half a degree, which is far inside
    // the uncertainty of the cloud and transmittance factors it feeds.
    position.declination =
        radians(23.45) * std::sin(2.0 * kPi * (284.0 + day) / 365.0);

    // Hour angle, zero at solar noon and 15 degrees per hour. The equation of
    // time is neglected: it shifts solar noon by at most a quarter of an hour,
    // which moves daily total radiation by well under a percent.
    const double hours = (julianDay - day) * 24.0;

    position.hourAngle = radians(15.0 * (hours - 12.0));

    const double phi = radians(latitude);

    const double cosZenith = std::sin(phi) * std::sin(position.declination) +
                             std::cos(phi) * std::cos(position.declination) *
                                 std::cos(position.hourAngle);

    position.elevation = std::asin(std::clamp(cosZenith, -1.0, 1.0));

    return position;
}

double dayLength(double julianDay, double latitude) noexcept
{
    const SolarPosition noon = solarPosition(std::floor(julianDay) + 0.5, latitude);

    const double phi = radians(latitude);

    // cos(omega_sunset) = -tan(phi) tan(delta), clamped for the polar day and
    // polar night where the sun never sets or never rises.
    const double argument = -std::tan(phi) * std::tan(noon.declination);

    if (argument <= -1.0)
        return 24.0;

    if (argument >= 1.0)
        return 0.0;

    return 2.0 * std::acos(argument) * 24.0 / (2.0 * kPi);
}

double clearSkySolarRadiation(double julianDay, double latitude,
                              double cloudCover) noexcept
{
    const SolarPosition position = solarPosition(julianDay, latitude);

    const double sinElevation = std::sin(position.elevation);

    if (sinElevation <= 0.0)
        return 0.0; // Below the horizon.

    const double day = std::floor(julianDay);

    // Solar constant, modulated by the Earth's orbital eccentricity.
    const double extraterrestrial =
        1367.0 * (1.0 + 0.033 * std::cos(2.0 * kPi * day / 365.0));

    // Air mass, and the usual empirical clear-sky transmittance. The air mass
    // is capped because 1/sin(elevation) runs away at the horizon, where the
    // flux is negligible anyway and the formula is not valid.
    const double airMass = std::min(1.0 / sinElevation, 38.0);

    const double transmittance = std::pow(0.7, std::pow(airMass, 0.678));

    const double clear = extraterrestrial * sinElevation * transmittance;

    const double cloud = std::clamp(cloudCover, 0.0, 1.0);

    return clear * (1.0 - 0.65 * cloud * cloud);
}

double equilibriumTemperature(const Meteorology &weather) noexcept
{
    // The balance is monotone decreasing in surface temperature: raising it
    // increases back radiation, evaporation and conduction, and changes nothing
    // else. So a sign change brackets the root and bisection converges on it.
    double low  = -50.0;
    double high = 100.0;

    // Guard against a bracket that does not contain a root — which would mean
    // the physics had changed, not that the solver had failed.
    if (surfaceHeatFlux(weather, low).total() < 0.0)
        return low;

    if (surfaceHeatFlux(weather, high).total() > 0.0)
        return high;

    for (int i = 0; i < 200; ++i)
    {
        const double middle = 0.5 * (low + high);

        if (surfaceHeatFlux(weather, middle).total() > 0.0)
            low = middle;
        else
            high = middle;
    }

    return 0.5 * (low + high);
}

} // namespace FVQual::Hydro
