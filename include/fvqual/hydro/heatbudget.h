// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   heatbudget.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  Surface heat exchange from meteorology (F3.5).
 *
 * Until now the surface heat flux has been a number the modeller supplies. That
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
 * - **Latent heat**, an evaporative loss driven by the vapour-pressure deficit
 *   between the water surface and the air, through a wind function.
 * - **Sensible heat**, conduction, tied to the latent term by the psychrometric
 *   constant.
 *
 * **Shortwave is not a surface flux, and that is the whole point.** Longwave,
 * latent and sensible exchange happen in the top millimetres and belong to the
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
     * `f(U) = evaporationA + evaporationB * U`, in W/(m²·Pa), so that the
     * latent flux is `f(U) * (e_s(T_s) - e_a)`. The defaults are a
     * Ryan-Harleman-style fit in SI; they are exposed because the wind function
     * is the least universal part of any heat budget and is routinely
     * recalibrated per site.
     */
    double evaporationA = 3.0e-3;
    double evaporationB = 1.8e-3;

    /*!
     * \brief Psychrometric constant (Pa/K) at standard pressure.
     *
     * Ties the sensible flux to the latent one: the Bowen ratio is
     * `gamma (T_s - T_a) / (e_s - e_a)`. Scaled by the actual pressure inside
     * the evaluation, so an alpine reservoir gets the right answer.
     */
    double psychrometric = 66.1;
};

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
 * \brief Saturation vapour pressure over water (Pa) at \p celsius.
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

} // namespace FVQual::Hydro

#endif // FVQual_HYDRO_HEATBUDGET_H
