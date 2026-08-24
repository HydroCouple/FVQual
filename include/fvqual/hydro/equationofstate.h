// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   equationofstate.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  Water density as a function of temperature, dissolved solids, and
 *         suspended solids (F3.2).
 *
 * Density is what makes a reservoir model three-dimensional in any
 * interesting sense: it sets the stratification that suppresses vertical
 * mixing, steers an inflow to the depth where it matches, and drives the
 * baroclinic circulation the internal mode solves for. A percent error here
 * is not a percent error in the answer — it is the difference between a lake
 * that stratifies and one that does not.
 *
 * \threadsafety Pure functions.
 */

#ifndef FVQual_HYDRO_EQUATIONOFSTATE_H
#define FVQual_HYDRO_EQUATIONOFSTATE_H

#include "fvqual/fvqual_export.h"

namespace FVQual::Hydro {

/*!
 * \brief Density of water from temperature, salinity, and suspended solids.
 *
 * SI units throughout: degrees Celsius, kg/m³, and g/m³ for solids.
 */
class FVQual_EXPORT EquationOfState
{
public:
    /*!
     * \brief Density of pure water at one atmosphere (kg/m³).
     *
     * The standard fifth-order polynomial in temperature, which reproduces
     * the density maximum near 4 °C — the feature responsible for winter
     * inverse stratification in deep temperate reservoirs, and the reason a
     * linear approximation is not good enough.
     *
     * \param[in] temperature Water temperature (°C).
     */
    [[nodiscard]] static double freshwaterDensity(double temperature) noexcept;

    /*!
     * \brief Density including dissolved and suspended solids (kg/m³).
     *
     * The salinity terms follow the UNESCO one-atmosphere international
     * equation of state. Suspended solids add their excess mass over the
     * water they displace, using a particle density of 2650 kg/m³ (quartz).
     *
     * \param[in] temperature      Water temperature (°C).
     * \param[in] salinity         Salinity or total dissolved solids (kg/m³,
     *                             numerically equal to practical salinity for
     *                             the dilute range reservoirs live in).
     * \param[in] suspendedSolids  Inorganic suspended solids (kg/m³).
     */
    [[nodiscard]] static double density(double temperature, double salinity = 0.0,
                                        double suspendedSolids = 0.0) noexcept;

    /*!
     * \brief Particle density used for the suspended-solids term (kg/m³).
     */
    static constexpr double kParticleDensity = 2650.0;

    /*!
     * \brief Reference density for the Boussinesq approximation (kg/m³).
     */
    static constexpr double kReferenceDensity = 1000.0;
};

} // namespace FVQual::Hydro

#endif // FVQual_HYDRO_EQUATIONOFSTATE_H
