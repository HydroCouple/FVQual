// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_forcingfield.cpp
 * \brief Forcing that is uniform or per column, and wind that reaches both of
 *        its consumers (F3.7).
 *
 * Two things need gating here, and only one of them is the new feature.
 *
 * **The feature:** a per-column field must actually be read. The failure mode
 * is not a crash but silence — a field set and never consulted looks exactly
 * like a field that happens to be uniform, and every conservation gate still
 * passes. So the tests set columns to *different* values and require the
 * outputs to differ, which is the ablation stated as an assertion.
 *
 * **The regression:** wind must reach the momentum equation as well as the heat
 * budget. The DeGray run set `windSpeed` and left `windStressX/Y` at zero, so
 * the wind evaporated water and moved none — and because that zeroed the shear,
 * it also collapsed Richardson-damped mixing to molecular diffusivity for an
 * entire simulated year. 254 unit tests passed over it, because each consumer
 * was correct in isolation and nothing asserted they were fed from one source.
 * `windStress` exists so that specifying wind specifies both, and the gates
 * below pin the arithmetic against numbers computed outside the code.
 */

#include "fvqual/forcing/forcingfield.h"
#include "fvqual/hydro/heatbudget.h"

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

using FVQual::Forcing::ForcingField;
using FVQual::Forcing::MeteorologyField;
using FVQual::Hydro::Meteorology;
using FVQual::Hydro::surfaceHeatFlux;
using FVQual::Hydro::windStress;

namespace {

//! Mild overcast conditions, so the radiative terms stand alone.
Meteorology mild()
{
    Meteorology weather;
    weather.airTemperature = 15.0;
    weather.dewPoint       = 10.0;
    weather.windSpeed      = 3.0;
    weather.cloudCover     = 0.2;
    weather.shortwave      = 400.0;
    return weather;
}

} // namespace

// ── ForcingField ──────────────────────────────────────────────────────────

TEST(ForcingField, AUniformFieldReadsTheSameEverywhere)
{
    const ForcingField field(7.25);

    EXPECT_TRUE(field.isUniform());

    for (int64_t c = 0; c < 50; ++c)
        EXPECT_EQ(field.at(c), 7.25);
}

TEST(ForcingField, PerColumnValuesAreReadPerColumn)
{
    // The gate the feature exists for: a field that is set but never consulted
    // is indistinguishable from a uniform one, so require the values to differ.
    ForcingField field;
    std::string  message;

    ASSERT_TRUE(field.setPerColumn({1.0, 2.0, 3.0}, 3, message)) << message;

    EXPECT_FALSE(field.isUniform());

    EXPECT_EQ(field.at(0), 1.0);
    EXPECT_EQ(field.at(1), 2.0);
    EXPECT_EQ(field.at(2), 3.0);

    EXPECT_NE(field.at(0), field.at(2));
}

TEST(ForcingField, AConstantPerColumnFieldEqualsTheUniformPathExactly)
{
    // Not "close": the same value through a different code path must produce
    // bit-identical arithmetic, which is what lets the existing heat-budget
    // gates keep their meaning once the supplier changes.
    const ForcingField uniform(3.75);

    ForcingField spread;
    std::string  message;

    ASSERT_TRUE(spread.setPerColumn({3.75, 3.75, 3.75, 3.75}, 4, message)) << message;

    for (int64_t c = 0; c < 4; ++c)
        EXPECT_EQ(spread.at(c), uniform.at(c));
}

TEST(ForcingField, AWrongLengthArrayIsRefusedAndChangesNothing)
{
    ForcingField field(9.0);
    std::string  message;

    EXPECT_FALSE(field.setPerColumn({1.0, 2.0}, 5, message));

    // The message must name both counts; "invalid input" sends the reader back
    // to the source to find out what was actually wrong.
    EXPECT_NE(message.find('2'), std::string::npos) << message;
    EXPECT_NE(message.find('5'), std::string::npos) << message;

    // Commit-on-success: the previous forcing must survive a refused call.
    EXPECT_TRUE(field.isUniform());
    EXPECT_EQ(field.at(0), 9.0);
}

TEST(ForcingField, SettingUniformDiscardsPerColumnValues)
{
    ForcingField field;
    std::string  message;

    ASSERT_TRUE(field.setPerColumn({1.0, 2.0}, 2, message)) << message;
    ASSERT_FALSE(field.isUniform());

    field.setUniform(4.5);

    EXPECT_TRUE(field.isUniform());
    EXPECT_EQ(field.at(0), 4.5);
    EXPECT_EQ(field.at(1), 4.5);
}

// ── MeteorologyField ──────────────────────────────────────────────────────

TEST(MeteorologyField, UniformSetupReproducesItsSourceStructExactly)
{
    const Meteorology weather = mild();

    MeteorologyField field;
    field.setUniform(weather);

    EXPECT_TRUE(field.isUniform());
    EXPECT_EQ(field.perColumnCount(), 0);

    for (int64_t c = 0; c < 8; ++c)
    {
        const Meteorology at = field.at(c);

        EXPECT_EQ(at.airTemperature, weather.airTemperature);
        EXPECT_EQ(at.dewPoint, weather.dewPoint);
        EXPECT_EQ(at.windSpeed, weather.windSpeed);
        EXPECT_EQ(at.cloudCover, weather.cloudCover);
        EXPECT_EQ(at.shortwave, weather.shortwave);
        EXPECT_EQ(at.albedo, weather.albedo);

        // And the flux, which is what actually matters downstream.
        EXPECT_EQ(surfaceHeatFlux(at, 18.0).total(),
                  surfaceHeatFlux(weather, 18.0).total());
    }
}

TEST(MeteorologyField, CoefficientsSurviveWhereObservationsVary)
{
    // The wind function is a property of the formulation, not of a place, so
    // it must not be disturbed by giving wind speed a per-column value.
    MeteorologyField field;
    field.setUniform(mild());

    field.base.evaporationB = 5.0e-3;
    field.base.windExponent = 2.0;

    std::string message;
    ASSERT_TRUE(field.windSpeed.setPerColumn({1.0, 9.0}, 2, message)) << message;

    EXPECT_EQ(field.at(0).evaporationB, 5.0e-3);
    EXPECT_EQ(field.at(1).evaporationB, 5.0e-3);
    EXPECT_EQ(field.at(1).windExponent, 2.0);
}

TEST(MeteorologyField, ShadingHalfABasinCoolsOnlyThatHalf)
{
    // The ablation, stated as an assertion. If the shortwave field were
    // ignored, both columns would return the same flux and this fails.
    MeteorologyField field;
    field.setUniform(mild());

    std::string message;
    ASSERT_TRUE(field.shortwave.setPerColumn({800.0, 0.0}, 2, message)) << message;

    const double sunlit = surfaceHeatFlux(field.at(0), 18.0).total();
    const double shaded = surfaceHeatFlux(field.at(1), 18.0).total();

    EXPECT_GT(sunlit, shaded) << "the shaded column received as much heat as the sunlit one";

    // The difference is the whole net shortwave, which is arithmetic we can
    // state independently: 800 W/m² less the albedo.
    EXPECT_NEAR(sunlit - shaded, 800.0 * (1.0 - field.at(0).albedo), 1.0e-9);
}

TEST(MeteorologyField, PerColumnCountReportsWhatVaries)
{
    MeteorologyField field;
    field.setUniform(mild());

    std::string message;

    ASSERT_TRUE(field.windSpeed.setPerColumn({1.0, 2.0}, 2, message)) << message;
    EXPECT_EQ(field.perColumnCount(), 1);

    ASSERT_TRUE(field.cloudCover.setPerColumn({0.1, 0.9}, 2, message)) << message;
    EXPECT_EQ(field.perColumnCount(), 2);
    EXPECT_FALSE(field.isUniform());
}

// ── Wind stress ───────────────────────────────────────────────────────────

TEST(WindStress, MagnitudeMatchesTheDragLawComputedByHand)
{
    // Against outside arithmetic rather than against the function itself: at
    // 10 m/s with rho 1.225 and C_D 1.3e-3, tau = 1.225 * 0.0013 * 100
    // = 0.15925 Pa. That is also a sanity check on the order of magnitude —
    // open-water stress at 10 m/s is a fraction of a pascal, not pascals.
    Meteorology weather;
    weather.windSpeed     = 10.0;
    weather.windDirection = 0.0;

    const auto stress = windStress(weather);

    EXPECT_NEAR(stress.magnitude(), 0.15925, 1.0e-6);
    EXPECT_NEAR(stress.x, 0.15925, 1.0e-6);
    EXPECT_NEAR(stress.y, 0.0, 1.0e-12);
}

TEST(WindStress, StressIsQuadraticInWindSpeed)
{
    Meteorology weather;
    weather.windDirection = 0.0;

    weather.windSpeed = 4.0;
    const double slow = windStress(weather).magnitude();

    weather.windSpeed = 8.0;
    const double fast = windStress(weather).magnitude();

    EXPECT_NEAR(fast / slow, 4.0, 1.0e-9) << "drag is not quadratic in speed";
}

TEST(WindStress, DirectionRotatesTheStressAndPreservesItsMagnitude)
{
    Meteorology weather;
    weather.windSpeed = 6.0;

    const double expected = 1.225 * 1.3e-3 * 36.0;

    weather.windDirection = 0.0;
    const auto east = windStress(weather);

    weather.windDirection = 2.0 * std::atan(1.0); // pi/2
    const auto north = windStress(weather);

    EXPECT_NEAR(east.x, expected, 1.0e-9);
    EXPECT_NEAR(east.y, 0.0, 1.0e-12);

    EXPECT_NEAR(north.x, 0.0, 1.0e-12);
    EXPECT_NEAR(north.y, expected, 1.0e-9);

    EXPECT_NEAR(east.magnitude(), north.magnitude(), 1.0e-12);
}

TEST(WindStress, StillAirExertsNoStress)
{
    Meteorology weather;
    weather.windSpeed = 0.0;

    EXPECT_EQ(windStress(weather).magnitude(), 0.0);
}

TEST(WindStress, TheSameWindDrivesBothHeatAndMomentum)
{
    // The regression gate. A wind speed that cools the water must also push it;
    // the defect was two consumers of one variable with only one of them wired,
    // and no test asserted the pairing.
    Meteorology calm = mild();
    calm.windSpeed   = 0.0;

    Meteorology blowing = mild();
    blowing.windSpeed   = 10.0;

    const double calmLatent    = surfaceHeatFlux(calm, 25.0).latent;
    const double blowingLatent = surfaceHeatFlux(blowing, 25.0).latent;

    EXPECT_LT(blowingLatent, calmLatent) << "wind did not increase evaporative cooling";

    EXPECT_EQ(windStress(calm).magnitude(), 0.0);
    EXPECT_GT(windStress(blowing).magnitude(), 0.0)
        << "wind cooled the water but exerted no stress — the DeGray defect";
}
