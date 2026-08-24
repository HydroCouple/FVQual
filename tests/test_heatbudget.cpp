// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_heatbudget.cpp
 * \brief Verification of the meteorological surface heat budget (F3.5).
 *
 * A prescribed heat flux cannot be wrong in an interesting way — it is whatever
 * the modeller typed. A computed one can, and the ways it goes wrong are
 * characteristic: a sign slip that heats at night, a wind function off by a
 * factor that makes a lake evaporate away, or sunlight dumped into the surface
 * layer so the model stratifies for the wrong reason.
 *
 * So the gates here are about behaviour rather than arithmetic. The budget must
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
    // The behavioural gate. A shallow, unstratified column under unchanging
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
    // spread over a few metres the interior warms too; if all the sunlight were
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
    // An existing model must not change behaviour because this code landed.
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
