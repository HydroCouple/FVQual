// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/hydro/equationofstate.h"

#include <gtest/gtest.h>

using FVQual::Hydro::EquationOfState;

TEST(EquationOfState, MatchesReferenceFreshwaterDensities)
{
    // Standard tabulated values for pure water at one atmosphere.
    EXPECT_NEAR(EquationOfState::freshwaterDensity(0.0), 999.84, 0.01);
    EXPECT_NEAR(EquationOfState::freshwaterDensity(4.0), 999.97, 0.01);
    EXPECT_NEAR(EquationOfState::freshwaterDensity(10.0), 999.70, 0.01);
    EXPECT_NEAR(EquationOfState::freshwaterDensity(20.0), 998.21, 0.01);
    EXPECT_NEAR(EquationOfState::freshwaterDensity(30.0), 995.65, 0.01);
}

TEST(EquationOfState, ReproducesTheDensityMaximumNearFourDegrees)
{
    // The density maximum near 4 °C is what drives winter inverse
    // stratification in deep temperate reservoirs; a linear equation of state
    // would miss it entirely.
    double peakTemperature = 0.0;
    double peakDensity     = 0.0;

    for (int i = 0; i <= 200; ++i)
    {
        const double t = 0.05 * i;
        const double d = EquationOfState::freshwaterDensity(t);

        if (d > peakDensity)
        {
            peakDensity     = d;
            peakTemperature = t;
        }
    }

    EXPECT_NEAR(peakTemperature, 4.0, 0.3);
    EXPECT_GT(EquationOfState::freshwaterDensity(4.0),
              EquationOfState::freshwaterDensity(0.0));
    EXPECT_GT(EquationOfState::freshwaterDensity(4.0),
              EquationOfState::freshwaterDensity(8.0));
}

TEST(EquationOfState, DensityFallsMonotonicallyAboveTheMaximum)
{
    for (double t = 5.0; t < 35.0; t += 1.0)
        EXPECT_GT(EquationOfState::density(t), EquationOfState::density(t + 1.0))
            << "at " << t << " °C";
}

TEST(EquationOfState, SalinityAndSolidsIncreaseDensity)
{
    const double fresh = EquationOfState::density(20.0);

    EXPECT_GT(EquationOfState::density(20.0, 1.0), fresh);
    EXPECT_GT(EquationOfState::density(20.0, 35.0), EquationOfState::density(20.0, 1.0));

    // Seawater at 20 °C, salinity 35 is close to 1024.8 kg/m³.
    EXPECT_NEAR(EquationOfState::density(20.0, 35.0), 1024.8, 0.5);

    // Suspended solids contribute only their excess over displaced water.
    const double withSolids = EquationOfState::density(20.0, 0.0, 1.0);
    EXPECT_GT(withSolids, fresh);
    EXPECT_LT(withSolids - fresh, 1.0);
}

TEST(EquationOfState, StratificationContrastIsPhysicallySized)
{
    // A 25 °C epilimnion over a 5 °C hypolimnion is a typical summer
    // reservoir. The contrast is small in absolute terms — under 0.3% — which
    // is precisely why the baroclinic pressure gradient must be computed
    // carefully rather than as the difference of two large numbers.
    const double warm = EquationOfState::density(25.0);
    const double cold = EquationOfState::density(5.0);

    EXPECT_GT(cold - warm, 2.0);
    EXPECT_LT(cold - warm, 4.0);
}
