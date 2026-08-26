// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/hydro/heatbudget.h"

#include <algorithm>
#include <cmath>

namespace FVQual::Hydro {

namespace {

//! Stefan-Boltzmann constant (W/(m²·K⁴)).
constexpr double kStefanBoltzmann = 5.670374419e-8;

//! Absolute zero, for the radiation terms (deg C).
constexpr double kAbsoluteZero = -273.15;

[[nodiscard]] double kelvin(double celsius) noexcept
{
    return celsius - kAbsoluteZero;
}

} // namespace

double saturationVapourPressure(double celsius) noexcept
{
    return 610.78 * std::exp(17.27 * celsius / (celsius + 237.3));
}

HeatFluxTerms surfaceHeatFlux(const Meteorology &weather,
                              double surfaceTemperature) noexcept
{
    HeatFluxTerms terms;

    // ── Shortwave ─────────────────────────────────────────────────────────

    terms.netShortwave = (1.0 - weather.albedo) * weather.shortwave;

    terms.surfaceSolar = weather.surfaceAbsorption * terms.netShortwave;
    terms.penetratingSolar = terms.netShortwave - terms.surfaceSolar;

    // ── Longwave in ───────────────────────────────────────────────────────

    // Swinbank's emissivity from air temperature alone, with the usual cloud
    // enhancement. It is an empirical fit, not a derivation, and it is the term
    // most worth replacing with a measured incoming-longwave record when one
    // exists.
    const double airK      = kelvin(weather.airTemperature);
    const double emissivity = 0.937e-5 * airK * airK;

    // Clamped, because clearSkySolarRadiation clamps and these two must agree
    // about what a cloud fraction is. They did not, and it cost a whole DeGray
    // run: CE-QUAL-W2 records cloud in tenths, so a value of 10 arrived here
    // and `1 + 0.17 C^2` turned the atmospheric longwave into eighteen times
    // itself — an equilibrium temperature of 64 C under 11 C air.
    const double cover = std::clamp(weather.cloudCover, 0.0, 1.0);

    const double cloud = 1.0 + 0.17 * cover * cover;

    terms.longwaveDown = (1.0 - weather.longwaveReflection) * emissivity *
                         kStefanBoltzmann * airK * airK * airK * airK * cloud;

    // ── Longwave out ──────────────────────────────────────────────────────

    const double surfaceK = kelvin(surfaceTemperature);

    terms.longwaveUp = -weather.waterEmissivity * kStefanBoltzmann * surfaceK *
                       surfaceK * surfaceK * surfaceK;

    // ── Latent and sensible ───────────────────────────────────────────────

    // The wind function, shared by both: they differ only in what drives them.
    const double speed = std::max(0.0, weather.windSpeed);

    const double wind =
        weather.evaporationA +
        weather.evaporationB *
            (weather.windExponent == 1.0 ? speed
                                         : std::pow(speed, weather.windExponent));

    const double saturated = saturationVapourPressure(surfaceTemperature);

    // Vapour pressure of the air is the saturation pressure at its dew point —
    // that is what a dew point *is*, and it avoids needing relative humidity
    // and air temperature to agree.
    const double ambient = saturationVapourPressure(weather.dewPoint);

    // Negative when the surface is wetter than the air, which is the usual
    // direction: evaporation cools.
    terms.latent = -wind * (saturated - ambient);

    // Conduction, through the Bowen ratio. Scaled by pressure because the
    // psychrometric constant is: at altitude the same temperature difference
    // moves less heat.
    const double gamma = weather.psychrometric * weather.pressure / 101325.0;

    terms.sensible = -wind * gamma * (surfaceTemperature - weather.airTemperature);

    return terms;
}

double WindStress::magnitude() const noexcept
{
    return std::sqrt(x * x + y * y);
}

WindStress windStress(const Meteorology &weather) noexcept
{
    WindStress stress;

    const double speed = std::max(0.0, weather.windSpeed);

    const double magnitude =
        weather.airDensity * weather.dragCoefficient * speed * speed;

    stress.x = magnitude * std::cos(weather.windDirection);
    stress.y = magnitude * std::sin(weather.windDirection);

    return stress;
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
