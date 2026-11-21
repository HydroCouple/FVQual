// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_heatbudget.cpp
 * \brief Verification of the meteorological surface heat budget (F3.5).
 *
 * A prescribed heat flux cannot be wrong in an interesting way — it is whatever
 * the modeler typed. A computed one can, and the ways it goes wrong are
 * characteristic: a sign slip that heats at night, a wind function off by a
 * factor that makes a lake evaporate away, or sunlight dumped into the surface
 * layer so the model stratifies for the wrong reason.
 *
 * So the gates here are about behavior rather than arithmetic. The budget must
 * be **self-limiting** — warm water must lose heat faster, which is what gives a
 * lake an equilibrium temperature — and a column left under steady weather must
 * actually find that temperature. Energy must be conserved: what the flux says
 * entered has to appear in the water. And the sunlight has to go *into the
 * column*, because that is the difference between a thermocline and a hot skin.
 */

#include "fvqual/hydro/heatbudget.h"
#include "fvqual/hydro/internalmode.h"
#include "fvqual/mesh/horizontalmesh.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/mesh/sigmaspec.h"
#include "fvqual/state/speciesregistry.h"
#include "fvqual/state/statestore.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <string>

using FVQual::Hydro::equilibriumTemperature;
using FVQual::Hydro::HeatFluxTerms;
using FVQual::Hydro::InternalModeOptions;
using FVQual::Hydro::InternalModeSolver;
using FVQual::Hydro::Meteorology;
using FVQual::Hydro::saturationVapourPressure;
using FVQual::Hydro::WindStress;
using FVQual::Hydro::windStress;
using FVQual::Hydro::surfaceHeatFlux;
using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;
using FVQual::State::SpeciesRegistry;
using FVQual::State::StateStore;
namespace Testing = FVQual::Testing;

namespace {

//! Mild, still, overcast night: no sun, so the radiative terms stand alone.
Meteorology night()
{
    Meteorology weather;
    weather.airTemperature = 15.0;
    weather.dewPoint       = 10.0;
    weather.windSpeed      = 2.0;
    weather.cloudCover     = 0.0;
    weather.shortwave      = 0.0;
    return weather;
}

//! A still basin with temperature, for the solver-level gates.
struct Basin
{
    HorizontalMesh     horizontal;
    LayeredMesh        layered;
    SpeciesRegistry    registry;
    StateStore         state;
    InternalModeSolver solver;

    int temperature = -1;

    double depth     = 0.0;
    double thickness = 0.0;
    int    layers    = 0;

    void build(int layerCount, double totalDepth, double initial)
    {
        std::string message;

        ASSERT_TRUE(horizontal.build(Testing::channelStrip(3, 300.0, 40.0), message))
            << message;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = -totalDepth;

        SigmaSpec spec;
        spec.layerCount        = layerCount;
        spec.minLayerThickness = 0.0;

        ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

        ASSERT_GE(registry.addTemperature(initial, message), 0) << message;
        temperature = registry.temperatureIndex();

        state.allocate(layered, registry);

        for (int64_t cell = 0; cell < layered.cellCount(); ++cell)
            state.speciesAt(temperature, cell) = initial;

        InternalModeOptions options;
        options.gravity   = 9.80665;
        options.manningN  = 0.0;
        options.cflNumber = 0.5;

        // No vertical mixing, so heat stays where it was put and the vertical
        // structure the gates measure is the budget's doing and nothing else.
        options.neutralViscosity      = 0.0;
        options.backgroundViscosity   = 0.0;
        options.backgroundDiffusivity = 0.0;
        options.convectiveMixing      = 0.0;

        ASSERT_TRUE(solver.initialize(layered, state, registry, options, message))
            << message;

        depth     = totalDepth;
        layers    = layerCount;
        thickness = totalDepth / static_cast<double>(layerCount);
    }

    [[nodiscard]] double at(int k) const
    {
        return state.speciesAt(temperature, layered.cell(0, k));
    }

    //! Column heat content relative to 0 °C, per unit area (J/m²).
    [[nodiscard]] double heatPerArea() const
    {
        const double capacity = 1000.0 * 4186.0;

        double total = 0.0;

        for (int k = 0; k < layers; ++k)
            total += at(k) * thickness * capacity;

        return total;
    }

    void run(double dt, int steps)
    {
        std::string message;

        for (int i = 0; i < steps; ++i)
            ASSERT_TRUE(solver.step(dt, message)) << message;
    }
};

} // namespace

// ── Wind reaches the momentum equation, not only the heat budget ──────────

TEST(HeatBudget, WindDrivenBySunAndAirAlsoDrivesTheWater)
{
    // The regression gate for the defect a whole-reservoir comparison found and
    // 254 unit tests did not: the DeGray driver set `windSpeed` for the heat
    // budget, `windStressX/Y` stayed at their default zero, and wind evaporated
    // water while moving none. Nothing failed, because each consumer was correct
    // in isolation and no test asserted they were fed from one source.
    //
    // The second-order consequence was worse than the missing current. Vertical
    // mixing is Richardson-damped, `Ri = N²/max(shear, 1e-12)`, so with nothing
    // driving the surface there is no shear, Ri is enormous, and the eddy
    // diffusivity collapses to its 1e-7 background — the mixing scheme runs every
    // step and does nothing. A parameterization switched off by an unset input in
    // another module is invisible to every budget check.
    Basin still;
    still.build(/*layers=*/6, /*totalDepth=*/12.0, 15.0);

    Meteorology calm = night();
    calm.windSpeed   = 0.0;

    still.solver.meteorology() = calm;
    still.solver.enableMeteorology(/*attenuation=*/0.45);

    still.run(60.0, 40);

    Basin blown;
    blown.build(/*layers=*/6, /*totalDepth=*/12.0, 15.0);

    Meteorology blowing  = night();
    blowing.windSpeed    = 12.0;
    blowing.windDirection = 0.0;

    blown.solver.meteorology() = blowing;
    blown.solver.enableMeteorology(/*attenuation=*/0.45);

    blown.run(60.0, 40);

    const int64_t surface = blown.layered.cell(0, 0);

    const double stillSpeed = std::abs(still.state.u[static_cast<size_t>(surface)]);
    const double blownSpeed = std::abs(blown.state.u[static_cast<size_t>(surface)]);

    EXPECT_NEAR(stillSpeed, 0.0, 1.0e-12) << "still air moved the water";

    EXPECT_GT(blownSpeed, 1.0e-4)
        << "wind cooled the water but did not move it — the DeGray defect";

    // Direction, not just magnitude: this wind blows along +x.
    EXPECT_GT(blown.state.u[static_cast<size_t>(surface)], 0.0)
        << "the surface current opposes the wind driving it";
}

// ── The terms ─────────────────────────────────────────────────────────────

TEST(HeatBudget, SaturationVapourPressureMatchesKnownValues)
{
    // Anchored against the textbook values, because everything evaporative
    // scales directly off this and an error here would be invisible elsewhere.
    EXPECT_NEAR(saturationVapourPressure(0.0), 611.0, 3.0);
    EXPECT_NEAR(saturationVapourPressure(10.0), 1228.0, 10.0);
    EXPECT_NEAR(saturationVapourPressure(20.0), 2338.0, 20.0);
    EXPECT_NEAR(saturationVapourPressure(30.0), 4243.0, 40.0);

    // Monotone, which the bisection in equilibriumTemperature relies on.
    for (double t = -20.0; t < 45.0; t += 1.0)
        EXPECT_GT(saturationVapourPressure(t + 1.0), saturationVapourPressure(t));
}

TEST(HeatBudget, TheBudgetIsSelfLimiting)
{
    // The property that makes a computed flux better than a prescribed one:
    // warming the water increases its losses, so the flux cannot run away. Back
    // radiation, evaporation and conduction all rise with surface temperature
    // and nothing else depends on it.
    const Meteorology weather = night();

    double previous = surfaceHeatFlux(weather, -5.0).total();

    for (double t = -4.0; t <= 40.0; t += 1.0)
    {
        const double flux = surfaceHeatFlux(weather, t).total();

        EXPECT_LT(flux, previous) << "flux did not fall as the surface warmed, at " << t;

        previous = flux;
    }
}

TEST(HeatBudget, WarmWaterUnderAColdNightSkyLosesHeat)
{
    // The sign convention, pinned where it is unambiguous. Warm water, cool air,
    // no sun: every term but the atmospheric longwave is a loss, and the total
    // must be negative.
    Meteorology weather = night();
    weather.airTemperature = 5.0;
    weather.dewPoint       = 0.0;

    const HeatFluxTerms terms = surfaceHeatFlux(weather, 20.0);

    EXPECT_LT(terms.total(), 0.0) << "warm water gained heat on a cold night";

    EXPECT_LT(terms.longwaveUp, 0.0);
    EXPECT_LT(terms.latent, 0.0) << "evaporation warmed the water";
    EXPECT_LT(terms.sensible, 0.0) << "conduction ran from cold air into warm water";
    EXPECT_GT(terms.longwaveDown, 0.0);
    EXPECT_EQ(terms.netShortwave, 0.0);
}

TEST(HeatBudget, SunlightIsSplitBetweenTheSurfaceAndTheColumn)
{
    Meteorology weather = night();
    weather.shortwave         = 800.0;
    weather.albedo            = 0.06;
    weather.surfaceAbsorption = 0.45;

    const HeatFluxTerms terms = surfaceHeatFlux(weather, 15.0);

    const double net = 800.0 * (1.0 - 0.06);

    EXPECT_NEAR(terms.netShortwave, net, 1.0e-9);
    EXPECT_NEAR(terms.surfaceSolar, 0.45 * net, 1.0e-9);
    EXPECT_NEAR(terms.penetratingSolar, 0.55 * net, 1.0e-9);

    // The penetrating part is *not* in the surface flux — that is the whole
    // distinction, and the thing a careless implementation gets wrong.
    EXPECT_NEAR(terms.total() - terms.surfaceNet(), terms.penetratingSolar, 1.0e-9);
}

TEST(HeatBudget, EquilibriumTemperatureIsWhereTheFluxVanishes)
{
    for (const Meteorology &weather :
         {night(), [] {
              Meteorology w = night();
              w.shortwave   = 600.0;
              return w;
          }(),
          [] {
              Meteorology w  = night();
              w.windSpeed    = 8.0;
              w.airTemperature = 30.0;
              w.dewPoint     = 25.0;
              return w;
          }()})
    {
        const double equilibrium = equilibriumTemperature(weather);

        EXPECT_NEAR(surfaceHeatFlux(weather, equilibrium).total(), 0.0, 1.0e-6)
            << "equilibrium " << equilibrium;

        // Either side of it, the flux pushes back toward it.
        EXPECT_GT(surfaceHeatFlux(weather, equilibrium - 1.0).total(), 0.0);
        EXPECT_LT(surfaceHeatFlux(weather, equilibrium + 1.0).total(), 0.0);
    }

    // Sunlight raises it, which is the sanity check anyone would apply.
    Meteorology sunny = night();
    sunny.shortwave   = 600.0;

    EXPECT_GT(equilibriumTemperature(sunny), equilibriumTemperature(night()) + 3.0);
}

// ── Through the solver ────────────────────────────────────────────────────

TEST(HeatBudgetSolver, AShallowColumnDriftsToTheEquilibriumTemperature)
{
    // The behavioral gate. A shallow, unstratified column under unchanging
    // weather must approach the equilibrium temperature the budget implies and
    // then stay there — that is what "equilibrium" has to mean for the number
    // to be worth reporting.
    Meteorology weather = night();
    weather.shortwave   = 300.0;
    weather.windSpeed   = 4.0;

    const double target = equilibriumTemperature(weather);

    // One layer, and that is not a shortcut: equilibrium temperature is defined
    // for a well-mixed body, and this fixture has vertical mixing switched off
    // so the vertical gates below can attribute structure to the budget alone.
    // With two layers and no mixing the layers are thermally independent — the
    // lower one receives sunlight and has no loss term at all, so it heats
    // without bound and no whole-column equilibrium exists for it to find. A
    // single layer *is* the well-mixed body the number describes.
    //
    // Shallow, so the thermal time constant (about 27 hours here) is short
    // enough to reach equilibrium within a test.
    Basin basin;
    basin.build(/*layers=*/1, /*totalDepth=*/0.25, target - 4.0);

    basin.solver.meteorology() = weather;

    // Nearly transparent water, so the whole column absorbs the sunlight and
    // behaves as the well-mixed body the equilibrium assumes.
    basin.solver.enableMeteorology(/*attenuation=*/0.05);

    ASSERT_LT(basin.at(0), target - 3.0);

    basin.run(600.0, 500); // ~83 hours, three thermal time constants

    for (int k = 0; k < basin.layers; ++k)
    {
        EXPECT_NEAR(basin.at(k), target, 0.5)
            << "layer " << k << " settled at " << basin.at(k) << " against " << target;
    }
}

TEST(HeatBudgetSolver, TheColumnDoesNotOvershootFromAboveEither)
{
    // Approached from the warm side, which is the direction a runaway would
    // show up in: if the losses were understated the column would sit above
    // equilibrium indefinitely.
    Meteorology weather = night();
    weather.shortwave   = 300.0;
    weather.windSpeed   = 4.0;

    const double target = equilibriumTemperature(weather);

    Basin basin;
    basin.build(/*layers=*/1, /*totalDepth=*/0.25, target + 4.0);

    basin.solver.meteorology() = weather;
    basin.solver.enableMeteorology(0.05);

    basin.run(600.0, 500);

    EXPECT_NEAR(basin.at(0), target, 0.5);
}

TEST(HeatBudgetSolver, EnergyEnteringTheWaterMatchesTheFlux)
{
    // Conservation. Over a short interval the surface temperature barely moves,
    // so the flux is nearly constant and the heat gained must equal flux times
    // time — with the penetrating sunlight included, since it is absorbed
    // within the column rather than lost.
    Meteorology weather = night();
    weather.shortwave   = 500.0;

    Basin basin;
    basin.build(/*layers=*/6, /*totalDepth=*/6.0, 15.0);

    basin.solver.meteorology() = weather;

    // Strong attenuation, so all the sunlight is absorbed well above the bed
    // and the bed term is not carrying the budget.
    basin.solver.enableMeteorology(/*attenuation=*/2.0);

    const double before = basin.heatPerArea();
    const double flux   = surfaceHeatFlux(weather, 15.0).total();

    constexpr double dt    = 30.0;
    constexpr int    steps = 20;

    basin.run(dt, steps);

    const double gained   = basin.heatPerArea() - before;
    const double expected = flux * dt * steps;

    ASSERT_GT(std::abs(expected), 0.0);

    // 1% — the surface temperature does drift a little over ten minutes, which
    // changes the flux, and that is physics rather than error.
    EXPECT_NEAR(gained, expected, 0.01 * std::abs(expected))
        << "gained " << gained << " J/m², flux implies " << expected;
}

TEST(HeatBudgetSolver, SunlightHeatsTheColumnAndNotJustTheSkin)
{
    // The reason the penetrating shortwave is distributed. With absorption
    // spread over a few meters the interior warms too; if all the sunlight were
    // handed to the surface layer, everything below it would be untouched, and
    // the model would stratify because of a numerical choice rather than an
    // optical one.
    Meteorology weather = night();
    weather.shortwave   = 800.0;
    weather.windSpeed   = 0.0;

    Basin basin;
    basin.build(/*layers=*/8, /*totalDepth=*/8.0, 15.0);

    basin.solver.meteorology() = weather;
    basin.solver.enableMeteorology(/*attenuation=*/0.4);

    basin.run(60.0, 60); // one hour

    // Every layer above the bed gained heat, not only the top one.
    for (int k = 1; k < basin.layers - 1; ++k)
    {
        EXPECT_GT(basin.at(k), 15.0)
            << "layer " << k << " was not reached by the sun";
    }

    // And the profile is stably ordered, warmest at the top: absorption falls
    // off with depth, so the heating must too.
    for (int k = 1; k < basin.layers - 1; ++k)
    {
        EXPECT_LT(basin.at(k), basin.at(k - 1))
            << "layer " << k << " warmed more than the one above it";
    }

    // The surface gained most, and by a clear margin — the gate would pass on a
    // nearly uniform warming otherwise.
    EXPECT_GT(basin.at(0) - 15.0, 2.0 * (basin.at(basin.layers - 2) - 15.0));
}

TEST(HeatBudgetSolver, TheBudgetIsInertUntilItIsEnabled)
{
    // An existing model must not change behavior because this code landed.
    // Without enableMeteorology(), the prescribed flux is still what acts.
    Meteorology weather = night();
    weather.shortwave   = 800.0;

    Basin basin;
    basin.build(/*layers=*/4, /*totalDepth=*/4.0, 15.0);

    basin.solver.meteorology() = weather;

    EXPECT_FALSE(basin.solver.meteorologyEnabled());

    basin.run(60.0, 20);

    for (int k = 0; k < basin.layers; ++k)
        EXPECT_NEAR(basin.at(k), 15.0, 1.0e-9) << "layer " << k;
}

// ── Solar geometry and radiation ──────────────────────────────────────────

TEST(SolarRadiation, DeclinationMatchesTheSolsticesAndEquinoxes)
{
    // Solar geometry is the one part of a heat budget with exact reference
    // values, so it is checked against them rather than against itself. The
    // Earth's tilt is 23.44 degrees, and the declination reaches it at the
    // solstices and passes through zero at the equinoxes.
    using FVQual::Hydro::solarPosition;

    constexpr double kDegrees = 57.29577951308232;

    // Day numbers: equinoxes near 80 and 266, solstices near 172 and 355.
    EXPECT_NEAR(solarPosition(80.5, 0.0).declination * kDegrees, 0.0, 1.0);
    EXPECT_NEAR(solarPosition(266.5, 0.0).declination * kDegrees, 0.0, 1.5);

    EXPECT_NEAR(solarPosition(172.5, 0.0).declination * kDegrees, 23.45, 0.5);
    EXPECT_NEAR(solarPosition(355.5, 0.0).declination * kDegrees, -23.45, 0.5);

    // And it never leaves the tropics.
    for (double day = 1.0; day < 366.0; day += 1.0)
    {
        const double declination = solarPosition(day + 0.5, 0.0).declination * kDegrees;

        EXPECT_LE(std::abs(declination), 23.46) << "day " << day;
    }
}

TEST(SolarRadiation, DayLengthIsTwelveHoursAtTheEquatorAndPolarAtThePoles)
{
    using FVQual::Hydro::dayLength;

    // The equator gets twelve hours all year — the cleanest check there is on
    // the geometry, and one that fails for almost any sign or units error.
    for (double day : {1.0, 80.0, 172.0, 266.0, 355.0})
        EXPECT_NEAR(dayLength(day + 0.5, 0.0), 12.0, 0.05) << "day " << day;

    // Midsummer above the Arctic circle: the sun does not set.
    EXPECT_NEAR(dayLength(172.5, 78.0), 24.0, 1.0e-9);

    // Midwinter there: it does not rise.
    EXPECT_NEAR(dayLength(355.5, 78.0), 0.0, 1.0e-9);

    // Mid-latitude summer is long and winter is short, and they are symmetric
    // about twelve.
    const double summer = dayLength(172.5, 45.0);
    const double winter = dayLength(355.5, 45.0);

    EXPECT_GT(summer, 15.0);
    EXPECT_LT(winter, 9.0);
    EXPECT_NEAR(summer + winter, 24.0, 0.3);
}

TEST(SolarRadiation, RadiationIsZeroAtNightAndPeaksAtSolarNoon)
{
    using FVQual::Hydro::clearSkySolarRadiation;

    constexpr double latitude = 34.2; // DeGray.

    // Midnight and the small hours: the sun is below the horizon.
    for (double hour : {0.0, 2.0, 22.0})
    {
        const double day = 172.0 + hour / 24.0;

        EXPECT_EQ(clearSkySolarRadiation(day, latitude, 0.0), 0.0)
            << "hour " << hour;
    }

    // Noon is the maximum, and the curve is symmetric about it.
    const double noon = clearSkySolarRadiation(172.5, latitude, 0.0);

    EXPECT_GT(noon, 700.0) << "midsummer noon at " << latitude << " deg";
    EXPECT_LT(noon, 1100.0);

    for (double offset : {1.0, 2.0, 3.0})
    {
        const double before = clearSkySolarRadiation(172.5 - offset / 24.0, latitude, 0.0);
        const double after  = clearSkySolarRadiation(172.5 + offset / 24.0, latitude, 0.0);

        EXPECT_LT(before, noon) << "offset " << offset;
        EXPECT_NEAR(before, after, 1.0e-9 * std::max(before, 1.0))
            << "the day is not symmetric about noon at offset " << offset;
    }
}

TEST(SolarRadiation, SummerDeliversMoreThanWinterAndCloudReducesIt)
{
    using FVQual::Hydro::clearSkySolarRadiation;

    constexpr double latitude = 34.2;

    // Daily totals, by trapezoid over the day.
    const auto dailyTotal = [&](double day, double cloud) {
        double total = 0.0;

        for (int i = 0; i < 24 * 6; ++i)
        {
            const double t = day + static_cast<double>(i) / (24.0 * 6.0);

            total += clearSkySolarRadiation(t, latitude, cloud) * 600.0;
        }

        return total; // J/m²
    };

    const double summer = dailyTotal(172.0, 0.0);
    const double winter = dailyTotal(355.0, 0.0);

    EXPECT_GT(summer, 2.0 * winter)
        << "summer " << summer << " J/m², winter " << winter;

    // Magnitudes a solar engineer would recognize: a clear midsummer day at
    // this latitude delivers of order 25-35 MJ/m².
    EXPECT_GT(summer, 20.0e6);
    EXPECT_LT(summer, 40.0e6);

    // Overcast roughly halves it, per the 1 - 0.65 C^2 reduction.
    const double overcast = dailyTotal(172.0, 1.0);

    EXPECT_NEAR(overcast / summer, 0.35, 0.02);
}

TEST(SolarRadiation, TheWindFunctionCanBeQuadratic)
{
    // DeGray's calibration is 9.2 + 0.46 W^2, and a linear form cannot follow
    // one at both low and high wind — at 8 m/s the quadratic term is six times
    // the constant.
    // 6 m/s rather than 8: with the W2 coefficients the quadratic-to-linear
    // ratio at 8 m/s comes out at exactly 3.0000, which sat on the wrong side
    // of a `> 3.0` bound by four hundredths of a watt. A coincidence of the
    // numbers, and a reminder to leave a gate room to breathe.
    Meteorology linear = night();
    linear.windSpeed    = 6.0;
    linear.windExponent = 1.0;

    Meteorology quadratic = linear;
    quadratic.windExponent = 2.0;

    const double linearLatent    = std::abs(surfaceHeatFlux(linear, 25.0).latent);
    const double quadraticLatent = std::abs(surfaceHeatFlux(quadratic, 25.0).latent);

    EXPECT_GT(quadraticLatent, 2.0 * linearLatent)
        << "linear " << linearLatent << " W/m², quadratic " << quadraticLatent;

    // At 1 m/s the two forms agree, since U and U² are both 1 there — which is
    // what makes the exponent a genuine shape change rather than a scale factor.
    Meteorology calm = linear;
    calm.windSpeed = 1.0;

    Meteorology calmQuadratic = calm;
    calmQuadratic.windExponent = 2.0;

    EXPECT_NEAR(surfaceHeatFlux(calm, 25.0).latent,
                surfaceHeatFlux(calmQuadratic, 25.0).latent, 1.0e-12);
}

TEST(SolarRadiation, ComputedRadiationDrivesTheBudgetItIsMeantTo)
{
    // The reason this exists: CE-QUAL-W2 computes shortwave when SROC is off,
    // which is how DeGray is configured — its meteorology file has no shortwave
    // column. So the budget has to be drivable from date and latitude alone.
    using FVQual::Hydro::clearSkySolarRadiation;

    Meteorology weather = night();
    weather.airTemperature = 28.0;
    weather.dewPoint       = 20.0;

    // Midsummer noon at DeGray's latitude.
    weather.shortwave = clearSkySolarRadiation(172.5, 34.2, 0.0);

    ASSERT_GT(weather.shortwave, 0.0);

    const double warm = equilibriumTemperature(weather);

    // The same weather at midnight has no sun at all.
    weather.shortwave = clearSkySolarRadiation(172.0, 34.2, 0.0);

    ASSERT_EQ(weather.shortwave, 0.0);

    const double cool = equilibriumTemperature(weather);

    EXPECT_GT(warm, cool + 5.0)
        << "sunlit equilibrium " << warm << " against night " << cool;
}

// ── Magnitudes, against the world rather than against ourselves ───────────

TEST(HeatBudget, TheTermsAreTheRightSizeAndNotJustSelfConsistent)
{
    // The gate that was missing, and whose absence let a wind function eight
    // times too weak pass 247 tests. Everything else here computes its
    // expectation from the same function it is testing, so a coefficient wrong
    // by an order of magnitude stays invisible: the budget is still
    // self-limiting, still conserves energy, still reaches "its" equilibrium.
    // Only a comparison against numbers from outside the code finds it.
    Meteorology weather;
    weather.airTemperature = 20.0;
    weather.dewPoint       = 12.0;
    weather.windSpeed      = 2.0;
    weather.cloudCover     = 0.0;
    weather.shortwave      = 0.0;

    const HeatFluxTerms terms = surfaceHeatFlux(weather, 22.0);

    // Back radiation from water near 20 C is around 400 W/m2 — Stefan-Boltzmann
    // at 295 K times an emissivity of 0.97.
    EXPECT_GT(-terms.longwaveUp, 350.0);
    EXPECT_LT(-terms.longwaveUp, 450.0);

    // Atmospheric longwave is a little less, since the sky is colder and its
    // emissivity below one.
    EXPECT_GT(terms.longwaveDown, 250.0);
    EXPECT_LT(terms.longwaveDown, 400.0);

    // Evaporation at a ~1.2 kPa deficit and a light breeze: tens of W/m2, and
    // emphatically not single digits.
    const double deficit = saturationVapourPressure(22.0) -
                           saturationVapourPressure(weather.dewPoint);

    EXPECT_GT(deficit, 1000.0);
    EXPECT_GT(-terms.latent, 40.0)
        << "evaporation of only " << -terms.latent
        << " W/m2 at a " << deficit << " Pa deficit is far too weak";
    EXPECT_LT(-terms.latent, 250.0);

    // The wind function itself, stated as a number a reader can check: the
    // latent flux per Pascal of deficit.
    const double perPascal = -terms.latent / deficit;

    EXPECT_GT(perPascal, 0.04) << "wind function " << perPascal << " W/(m2.Pa)";
    EXPECT_LT(perPascal, 0.20);
}

TEST(HeatBudget, EquilibriumTemperatureStaysInPhysicalBounds)
{
    // A lake is not far from the air above it *in the dark*: with no sun, the
    // radiative and evaporative terms tie the equilibrium to air temperature
    // within a few degrees. Sunlight raises it, and by a lot — 700 W/m2 held
    // steady is peak noon sustained, and it genuinely drives equilibrium 20 C
    // above a cold air mass. So the claim has to be split, which my first
    // version of this gate did not do: it bounded everything at +/-15 C and
    // duly failed on sunlit cases that were perfectly physical.
    for (double airTemperature : {5.0, 15.0, 25.0, 35.0})
    {
        Meteorology dark;
        dark.airTemperature = airTemperature;
        dark.dewPoint       = airTemperature - 5.0;
        dark.windSpeed      = 3.0;
        dark.shortwave      = 0.0;

        const double unlit = equilibriumTemperature(dark);

        EXPECT_GT(unlit, airTemperature - 12.0)
            << "air " << airTemperature << " C gave a dark equilibrium of " << unlit;
        EXPECT_LT(unlit, airTemperature + 5.0)
            << "air " << airTemperature << " C gave a dark equilibrium of " << unlit;

        // Sun raises it, monotonically, and never to anything absurd.
        double previous = unlit;

        for (double sun : {200.0, 400.0, 700.0})
        {
            Meteorology weather = dark;
            weather.shortwave   = sun;

            const double equilibrium = equilibriumTemperature(weather);

            EXPECT_GT(equilibrium, previous)
                << "more sun did not warm the equilibrium, at " << sun << " W/m2";

            EXPECT_LT(equilibrium, airTemperature + 40.0)
                << "air " << airTemperature << " C and " << sun
                << " W/m2 gave " << equilibrium;

            previous = equilibrium;
        }
    }
}

TEST(HeatBudget, AnOutOfRangeCloudFractionCannotRunAway)
{
    // Cloud cover is a fraction here, but the world supplies it in eighths,
    // tenths and percent. The cloud factor is quadratic, so a value of 10 where
    // 1.0 was meant multiplies the atmospheric longwave by eighteen — which is
    // exactly what happened on the first DeGray run, and produced an
    // equilibrium temperature of 64 C beneath 11 C air.
    //
    // The conversion belongs at the boundary where the conventions meet, and it
    // is done there. This is the second line: whatever arrives, the term stays
    // physical.
    Meteorology weather;
    weather.airTemperature = 11.0;
    weather.dewPoint       = 7.0;
    weather.windSpeed      = 5.0;
    weather.shortwave      = 0.0;

    // The reference is fully overcast — cloudCover = 1.0. My first version left
    // it at the default of zero and so compared a clear sky against a clouded
    // one, which failed for the entirely correct reason that cloud raises
    // atmospheric longwave.
    weather.cloudCover = 1.0;

    const double overcast = surfaceHeatFlux(weather, 10.0).longwaveDown;

    for (double bad : {1.0, 8.0, 10.0, 100.0})
    {
        weather.cloudCover = bad;

        const HeatFluxTerms terms = surfaceHeatFlux(weather, 10.0);

        EXPECT_LE(terms.longwaveDown, overcast * 1.0001)
            << "cloud cover " << bad << " gave " << terms.longwaveDown << " W/m2";

        // And the equilibrium stays somewhere a lake could actually sit.
        EXPECT_LT(equilibriumTemperature(weather), weather.airTemperature + 20.0)
            << "cloud cover " << bad;
    }
}

TEST(HeatBudget, EqualReferenceHeightsLeaveTheWindFunctionUntouched)
{
    // The default, and what every existing model gets. If this moved, the fix
    // for one dataset would have silently rewritten every other result.
    Meteorology weather;
    weather.airTemperature = 18.0;
    weather.dewPoint       = 8.0;
    weather.windSpeed      = 5.0;

    const HeatFluxTerms plain = surfaceHeatFlux(weather, 20.0);

    weather.windHeight        = 10.0;
    weather.evaporationHeight = 10.0;
    weather.windRoughness     = 0.01;

    const HeatFluxTerms stated = surfaceHeatFlux(weather, 20.0);

    EXPECT_DOUBLE_EQ(plain.latent, stated.latent);
    EXPECT_DOUBLE_EQ(plain.sensible, stated.sensible);
}

TEST(HeatBudget, TheEvaporativeWindIsBroughtToItsOwnReferenceHeight)
{
    // CE-QUAL-W2 scales the measured wind to 10 m for the surface stress and to
    // 2 m for the evaporative wind function. DeGray's anemometer is at 10 m, so
    // driving the wind function at 10 m applies coefficients calibrated at 2 m
    // to a wind that is 1/0.767 too fast — and since the function is quadratic,
    // the wind-driven part of it comes out about 1.7x too large.
    Meteorology weather;
    weather.airTemperature = 18.0;
    weather.dewPoint       = 8.0;
    weather.windSpeed      = 5.0;
    weather.windRoughness  = 0.01;

    weather.windHeight        = 10.0;
    weather.evaporationHeight = 10.0;
    const HeatFluxTerms atTen = surfaceHeatFlux(weather, 20.0);

    weather.evaporationHeight = 2.0;
    const HeatFluxTerms atTwo = surfaceHeatFlux(weather, 20.0);

    // Evaporating, so latent is negative; the corrected one must cool less.
    ASSERT_LT(atTen.latent, 0.0);
    EXPECT_GT(atTwo.latent, atTen.latent) << "the 2 m wind must evaporate less";

    // Against the profile computed here rather than read back from the code.
    const double ratio = std::log(2.0 / 0.01) / std::log(10.0 / 0.01);
    EXPECT_NEAR(ratio, 0.767, 0.002);

    const double expected =
        weather.evaporationA + weather.evaporationB * std::pow(5.0 * ratio, 2.0);
    const double drove =
        weather.evaporationA + weather.evaporationB * std::pow(5.0, 2.0);

    // The fluxes are proportional to the wind function, so their ratio is it.
    EXPECT_NEAR(atTwo.latent / atTen.latent, expected / drove, 1.0e-9);

    // And it must be a real reduction, not a rounding difference.
    EXPECT_LT(expected / drove, 0.85);
}

TEST(HeatBudget, MomentumStressIgnoresTheEvaporationHeight)
{
    // The stress is calibrated at the height the wind is stored at. If the
    // evaporation correction leaked into it, this fix would have quietly
    // weakened the wind stress too — and wind stress is the one lever that
    // moves DeGray's interior bias.
    Meteorology weather;
    weather.windSpeed     = 6.0;
    weather.windRoughness = 0.01;
    weather.windHeight    = 10.0;

    weather.evaporationHeight = 10.0;
    const WindStress full = windStress(weather);

    weather.evaporationHeight = 2.0;
    const WindStress reduced = windStress(weather);

    EXPECT_DOUBLE_EQ(full.magnitude(), reduced.magnitude());
}
