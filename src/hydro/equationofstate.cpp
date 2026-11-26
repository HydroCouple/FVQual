// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/hydro/equationofstate.h"
#include "fvqual/hydro/columnphysics.h"

#include <algorithm>
#include <cmath>

namespace FVQual::Hydro {

double EquationOfState::freshwaterDensity(double temperature) noexcept
{
    return ColumnPhysics::freshwaterDensity(temperature);
}

double EquationOfState::density(double temperature, double salinity,
                                double suspendedSolids) noexcept
{
    const double t = std::clamp(temperature, -2.0, 100.0);
    const double s = std::max(0.0, salinity);

    double result = freshwaterDensity(t);

    if (s > 0.0)
    {
        const double a = 0.824493 - 4.0899e-3 * t + 7.6438e-5 * t * t -
                         8.2467e-7 * t * t * t + 5.3875e-9 * t * t * t * t;

        const double b = -5.72466e-3 + 1.0227e-4 * t - 1.6546e-6 * t * t;

        result += a * s + b * s * std::sqrt(s) + 4.8314e-4 * s * s;
    }

    if (suspendedSolids > 0.0)
    {
        // Suspended particles displace water, so only their excess density
        // counts toward the mixture.
        result += suspendedSolids * (1.0 - result / kParticleDensity);
    }

    return result;
}

} // namespace FVQual::Hydro
