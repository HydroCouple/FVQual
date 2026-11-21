// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   heatbudget.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  Surface heat exchange from meteorology (F3.5).
 *
 * Until now the surface heat flux has been a number the modeler supplies. That
 * is fine for a controlled experiment and useless for a reservoir, because the
 * flux is not an input — it depends on the water temperature it is trying to
 * set. Warm water radiates and evaporates harder, which is the negative
 * feedback that gives a lake an equilibrium temperature at all. A prescribed
 * flux has no such feedback and will happily drive a column past boiling.
 *
 * **The terms**, all in W/m², positive into the water:
 *
 * - **Net shortwave** `(1 - albedo) * incident`. Most of it does not stop at
 *   the surface.
 * - **Atmospheric longwave**, from air temperature with a cloud correction, less
 *   the small fraction reflected.
 * - **Back radiation**, `epsilon sigma T_s^4`, the largest single loss and the
 *   one that makes the budget self-limiting.
 * - **Latent heat**, an evaporative loss driven by the vapor-pressure deficit
 *   between the water surface and the air, through a wind function.
 * - **Sensible heat**, conduction, tied to the latent term by the psychrometric
 *   constant.
 *
 * **Shortwave is not a surface flux, and that is the whole point.** Longwave,
 * latent and sensible exchange happen in the top millimeters and belong to the
 * surface layer. Sunlight does not: it penetrates, and the depth over which it
 * is absorbed is what builds a thermocline. Handing all of it to the surface
 * layer produces a model that stratifies for the wrong reason and mixes down
 * far too readily. So the penetrating part is distributed through the water
 * column using the same Beer-Lambert march the kinetics already use for light —
 * the absorbed energy in a cell is the difference in irradiance across it,
 * which F6.2a already computes and already gates as an energy budget.
 *
 * **Equilibrium temperature** falls out of the formulation rather than being
 * imposed: the temperature at which the terms cancel. It is the single most
 * useful number for sanity-checking a heat budget, so it is exposed and gated.
 *
 * \threadsafety The evaluation is a pure function of meteorology and surface
 *               temperature.
 */

#ifndef FVQual_HYDRO_HEATBUDGET_H
#define FVQual_HYDRO_HEATBUDGET_H

#include "fvqual/fvqual_export.h"

namespace FVQual::Hydro {

/*!
 * \brief Meteorological forcing, uniform over the water body.
 *
 * Spatially uniform because a reservoir is small against a weather system. Set
 * these from a met record each step; they are ordinary settable state, not
 * configuration.
 */
struct FVQual_EXPORT Meteorology
{
    double airTemperature = 20.0;  //!< Dry-bulb air temperature (deg C).
    double dewPoint       = 10.0;  //!< Dew-point temperature (deg C).
    double windSpeed      = 2.0;   //!< Wind speed at the reference height (m/s).

    /*!
     * \brief Direction the wind blows *towards* (radians, zero along +x,
     *        increasing counter-clockwise).
     *
     * FVQual's own convention, stated because there are at least three in use
     * and picking silently is how a reservoir ends up mixing along the wrong
     * axis. Meteorological records overwhelmingly give the direction wind comes
     * *from*, and CE-QUAL-W2's `PHI` is one of them — converting is the
     * caller's job, at the boundary where the two conventions meet. That is the
     * same rule cloud cover now follows, and for the same reason: the last time
     * a unit convention was resolved anywhere but the boundary it cost a whole
     * DeGray run.
     */
    double windDirection = 0.0;
    double cloudCover     = 0.0;   //!< Fraction of sky covered, 0 to 1.
    double shortwave      = 0.0;   //!< Incident shortwave irradiance (W/m²).
    double pressure       = 101325.0; //!< Atmospheric pressure (Pa).

    //! Shortwave reflected at the surface (-).
    double albedo = 0.06;

    /*!
     * \brief Fraction of net shortwave absorbed at the surface (-).
     *
     * CE-QUAL-W2's BETA. The remainder penetrates and is absorbed through the
     * water column — see the file note on why that distinction matters.
     */
    double surfaceAbsorption = 0.45;

    //! Emissivity of the water surface (-), for back radiation.
    double waterEmissivity = 0.97;

    //! Fraction of incoming longwave reflected at the surface (-).
    double longwaveReflection = 0.03;

    /*!
     * \brief Wind function coefficients for the latent heat flux.
     *
     * `f(U) = evaporationA + evaporationB * U^windExponent`, in W/(m²·Pa), so
     * that the latent flux is `f(U) * (e_s(T_s) - e_a)`.
     *
     * The exponent exists because calibrated wind functions are routinely
     * quadratic — CE-QUAL-W2's DeGray application uses `9.2 + 0.46 W²` in
     * W/(m²·mmHg) — and a linear form cannot represent one at both low and high
     * wind.
     *
     * The defaults are exactly that W2 form converted to SI by dividing by
     * 133.322 Pa/mmHg. They replace an earlier guess that was roughly eight
     * times too small: it gave a latent flux of 10 W/m² at a 1000 Pa vapor
     * deficit where the physical value is nearer 80, and since evaporation is
     * the dominant cooling term at high temperature, that pushed the
     * equilibrium temperature tens of degrees too high. Nothing in the unit
     * tests caught it, because they all computed their expectations from the
     * same function; a whole-reservoir run did, immediately.
     *
     * This is the least universal part of any heat budget and the first thing
     * recalibrated per site, which is why all three are exposed.
     */
    double evaporationA = 6.901e-2;
    double evaporationB = 3.450e-3;
    double windExponent = 2.0;

    /*!
     * \brief Height the stored \ref windSpeed refers to (m), and the height the
     *        evaporative wind function is calibrated for (m).
     *
     * **Why two heights for one wind.** Momentum stress and the evaporative wind
     * function are calibrated against different reference heights, and CE-QUAL-W2
     * carries both: it scales the measured wind to 10 m for the surface stress
     * and to 2 m for the wind function, through the logarithmic profile
     *
     *     U(z) = U(h) * ln(z / z0) / ln(h / z0).
     *
     * FVQual drove both from one wind. With DeGray's anemometer at 10 m and a
     * roughness length of 0.01 m the 2 m wind is 0.767 of it, so applying
     * coefficients calibrated at 2 m to the 10 m wind inflates the quadratic
     * term by 1/0.767^2, about 1.7x. Excess evaporation cools the surface, and
     * the measured surface bias is 1.9 °C too cold.
     *
     * Setting \ref evaporationHeight equal to \ref windHeight disables the
     * conversion, which is what a model with no stated heights should do.
     */
    double windHeight        = 10.0;
    double evaporationHeight = 10.0;

    //! Aerodynamic roughness length for that profile (m).
    double windRoughness = 0.01;

    /*!
     * \brief Psychrometric constant (Pa/K) at standard pressure.
     *
     * Ties the sensible flux to the latent one: the Bowen ratio is
     * `gamma (T_s - T_a) / (e_s - e_a)`. Scaled by the actual pressure inside
     * the evaluation, so an alpine reservoir gets the right answer.
     */
    double psychrometric = 66.1;

    /*!
     * \brief Air density (kg/m³) and surface drag coefficient (-), for the
     *        momentum the wind carries.
     *
     * Separate from the wind function above, which governs the *heat* the wind
     * carries. Both are driven by `windSpeed`, and keeping them apart is what
     * lets one be wired and the other forgotten — which is exactly what
     * happened: the DeGray run set wind for the heat budget and left the
     * momentum stress at zero, so wind evaporated water and moved none. See
     * \ref windStress, which exists so that supplying wind supplies both.
     */
    double airDensity      = 1.225;
    double dragCoefficient = 1.3e-3;
};

/*!
 * \brief Surface wind stress components (Pa).
 */
struct FVQual_EXPORT WindStress
{
    double x = 0.0;
    double y = 0.0;

    //! Magnitude (Pa).
    [[nodiscard]] double magnitude() const noexcept;
};

/*!
 * \brief Wind stress from wind speed and direction, `tau = rho_a C_D W²`.
 *
 * Derived here rather than supplied by the caller alongside the wind speed,
 * deliberately. A model that takes speed for its heat budget and stress for its
 * momentum has two inputs describing one wind, and nothing checks that they
 * agree — or that the second was set at all. Deriving it means specifying wind
 * once reaches both consumers.
 *
 * The quadratic drag law with a constant coefficient is the standard first
 * approximation. `dragCoefficient` is exposed because it is genuinely
 * site-dependent — sheltered reservoirs are routinely calibrated well below the
 * open-water value.
 */
[[nodiscard]] FVQual_EXPORT WindStress windStress(const Meteorology &weather) noexcept;

/*!
 * \brief The individual terms of the surface exchange (W/m²), positive into
 *        the water.
 */
struct FVQual_EXPORT HeatFluxTerms
{
    double netShortwave  = 0.0; //!< After albedo; sum of the two below.
    double surfaceSolar  = 0.0; //!< The part absorbed in the surface layer.
    double penetratingSolar = 0.0; //!< The part that enters the water column.
    double longwaveDown   = 0.0; //!< Atmospheric, less reflection.
    double longwaveUp     = 0.0; //!< Back radiation. Negative.
    double latent         = 0.0; //!< Evaporative. Negative when evaporating.
    double sensible       = 0.0; //!< Conductive.

    /*!
     * \brief Everything that acts *at* the surface (W/m²).
     *
     * Excludes the penetrating shortwave, which is a distributed source and is
     * applied through the water column instead. Adding it here would be the
     * classic error the file note describes.
     */
    [[nodiscard]] double surfaceNet() const noexcept
    {
        return surfaceSolar + longwaveDown + longwaveUp + latent + sensible;
    }

    //! Everything, including the penetrating shortwave.
    [[nodiscard]] double total() const noexcept
    {
        return surfaceNet() + penetratingSolar;
    }
};

/*!
 * \brief Saturation vapor pressure over water (Pa) at \p celsius.
 *
 * Magnus-Tetens. Accurate to a few tenths of a percent between -30 and 50 °C,
 * which is the range a reservoir model sees.
 */
[[nodiscard]] FVQual_EXPORT double saturationVapourPressure(double celsius) noexcept;

/*!
 * \brief Evaluates the surface exchange for a given water-surface temperature.
 *
 * \param[in] weather            Current meteorology.
 * \param[in] surfaceTemperature Water-surface temperature (deg C).
 */
[[nodiscard]] FVQual_EXPORT HeatFluxTerms surfaceHeatFlux(const Meteorology &weather,
                                                          double surfaceTemperature)
    noexcept;

/*!
 * \brief Water temperature at which the surface terms cancel (deg C).
 *
 * The temperature a shallow, well-mixed body would drift to and stay at under
 * unchanging weather. Solved by bisection, because the balance is monotone in
 * surface temperature — back radiation and evaporation both rise with it — so a
 * root exists and is unique, and bisection cannot run away where a Newton step
 * on a quartic could.
 *
 * Uses total(), penetrating shortwave included: a shallow body absorbs all of
 * it, so the equilibrium is the one the whole budget implies.
 *
 * \param[in] weather Current meteorology.
 */
[[nodiscard]] FVQual_EXPORT double equilibriumTemperature(
    const Meteorology &weather) noexcept;

/*!
 * \brief Where the sun is, for a place and an instant.
 */
struct FVQual_EXPORT SolarPosition
{
    double declination = 0.0; //!< Solar declination (radians).
    double hourAngle   = 0.0; //!< Hour angle, zero at solar noon (radians).
    double elevation   = 0.0; //!< Elevation above the horizon (radians), <0 at night.
};

/*!
 * \brief Solar position for a Julian day of year and a location.
 *
 * \param[in] julianDay Day of year with a fractional part, 1.0 at midnight on
 *                      1 January. Local time, not UTC — a reservoir model runs
 *                      on the clock its meteorology was recorded in.
 * \param[in] latitude  Degrees north.
 */
[[nodiscard]] FVQual_EXPORT SolarPosition solarPosition(double julianDay,
                                                        double latitude) noexcept;

/*!
 * \brief Shortwave irradiance on a horizontal surface (W/m²).
 *
 * Needed because CE-QUAL-W2 computes rather than reads solar radiation when its
 * `SROC` switch is off, which is how the DeGray application is configured — its
 * meteorology file has no shortwave column at all. A heat budget that can only
 * be driven by measured radiation cannot run that model.
 *
 * Extraterrestrial irradiance is reduced by an air-mass transmittance and then
 * by cloud. Clear-sky transmittance follows the usual `0.7^(m^0.678)` form; the
 * cloud reduction is `1 - 0.65 C²`, which is what W2 uses. Both are empirical
 * fits, and neither is a substitute for a measured record where one exists.
 *
 * \param[in] julianDay  Day of year with a fractional part.
 * \param[in] latitude   Degrees north.
 * \param[in] cloudCover Fraction of sky covered, 0 to 1.
 */
[[nodiscard]] FVQual_EXPORT double clearSkySolarRadiation(double julianDay,
                                                          double latitude,
                                                          double cloudCover) noexcept;

/*!
 * \brief Hours between sunrise and sunset at a latitude and day of year.
 *
 * Exposed because it is the cleanest check on the solar geometry: twelve hours
 * at the equator all year, and the polar day and night at high latitude.
 */
[[nodiscard]] FVQual_EXPORT double dayLength(double julianDay,
                                             double latitude) noexcept;

} // namespace FVQual::Hydro

#endif // FVQual_HYDRO_HEATBUDGET_H
