// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/hydro/heatbudget.h"

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

    const double cloud = 1.0 + 0.17 * weather.cloudCover * weather.cloudCover;

    terms.longwaveDown = (1.0 - weather.longwaveReflection) * emissivity *
                         kStefanBoltzmann * airK * airK * airK * airK * cloud;

    // ── Longwave out ──────────────────────────────────────────────────────

    const double surfaceK = kelvin(surfaceTemperature);

    terms.longwaveUp = -weather.waterEmissivity * kStefanBoltzmann * surfaceK *
                       surfaceK * surfaceK * surfaceK;

    // ── Latent and sensible ───────────────────────────────────────────────

    // The wind function, shared by both: they differ only in what drives them.
    const double wind = weather.evaporationA +
                        weather.evaporationB * std::max(0.0, weather.windSpeed);

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
