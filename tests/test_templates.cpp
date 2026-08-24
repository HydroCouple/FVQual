// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_templates.cpp
 * \brief Verification of the constituent template library (F6.4).
 *
 * A template's job is to produce a model that is *right*, so these gates check
 * the physics it configures rather than that the calls succeeded. Water age has
 * to advance at a day per day; a decay template has to reproduce its closed
 * form; a temperature-corrected rate has to be corrected by the factor the
 * literature quotes.
 *
 * The most useful thing a template can get wrong is a unit conversion, because
 * a rate that is 86400 times too small produces output that looks entirely
 * plausible. Every gate below is therefore written in the units a modeller
 * would use — days — and compares against a value computed independently in
 * those units.
 */

#include "fvqual/hydro/internalmode.h"
#include "fvqual/kinetics/templates.h"
#include "fvqual/mesh/horizontalmesh.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/mesh/sigmaspec.h"
#include "fvqual/state/speciesregistry.h"
#include "fvqual/state/statestore.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <string>

using FVQual::Hydro::InternalModeOptions;
using FVQual::Hydro::InternalModeSolver;
using FVQual::Kinetics::FirstOrderDecayTemplate;
using FVQual::Kinetics::kSecondsPerDay;
using FVQual::Kinetics::WaterAgeTemplate;
using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;
using FVQual::State::SpeciesRegistry;
using FVQual::State::StateStore;
namespace Testing = FVQual::Testing;

namespace {

/*!
 * \brief A still, unstratified basin — the setting in which a template's rate
 *        laws are the only thing changing the answer.
 */
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

    //! Registers species, then builds. Templates must register before this.
    void build(int layers, double totalDepth, double waterTemperature = 20.0)
    {
        std::string message;

        ASSERT_TRUE(horizontal.build(Testing::channelStrip(3, 300.0, 40.0), message))
            << message;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = -totalDepth;

        SigmaSpec spec;
        spec.layerCount        = layers;
        spec.minLayerThickness = 0.0;

        ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

        state.allocate(layered, registry);

        temperature = registry.temperatureIndex();

        if (temperature >= 0)
        {
            for (int64_t cell = 0; cell < layered.cellCount(); ++cell)
                state.speciesAt(temperature, cell) = waterTemperature;
        }

        InternalModeOptions options;
        options.gravity               = 9.80665;
        options.manningN              = 0.0;
        options.cflNumber             = 0.5;
        options.surfaceHeatFlux       = 0.0;
        options.neutralViscosity      = 0.0;
        options.backgroundViscosity   = 0.0;
        options.backgroundDiffusivity = 0.0;
        options.convectiveMixing      = 0.0;

        ASSERT_TRUE(solver.initialize(layered, state, registry, options, message))
            << message;

        depth     = totalDepth;
        thickness = totalDepth / static_cast<double>(layers);
    }

    void setUniform(int species, double value)
    {
        for (int64_t cell = 0; cell < layered.cellCount(); ++cell)
            state.speciesAt(species, cell) = value;
    }

    [[nodiscard]] double at(int species, int k) const
    {
        return state.speciesAt(species, layered.cell(0, k));
    }

    void run(double dt, int steps)
    {
        std::string message;

        for (int i = 0; i < steps; ++i)
            ASSERT_TRUE(solver.step(dt, message)) << message;
    }
};

} // namespace

// ── Water age ─────────────────────────────────────────────────────────────

TEST(Templates, WaterAgeAdvancesOneDayPerDay)
{
    // The unit gate. Age is quoted in days and the engine runs in seconds, so
    // a template that forgot the conversion would report ages 86400 times too
    // small — and "0.0002 days" is not obviously wrong at a glance.
    WaterAgeTemplate age;

    Basin basin;

    std::string message;
    ASSERT_TRUE(age.registerSpecies(basin.registry, message)) << message;

    basin.build(/*layers=*/5, /*totalDepth=*/20.0);

    const int index = basin.solver.speciesIndex(SpeciesRegistry::kAge);
    ASSERT_GE(index, 0);

    ASSERT_TRUE(age.configure(basin.solver, message)) << message;
    ASSERT_TRUE(basin.solver.compileReactions(message)) << message;

    // A tenth of a day, exactly: 10 x 864 s. Short enough that the suite stays
    // usable, long enough that a missing 86400 would be off by four orders of
    // magnitude rather than hidden in a tolerance.
    constexpr double dt    = 864.0;
    constexpr int    steps = 10;
    constexpr double days  = dt * steps / kSecondsPerDay;

    basin.run(dt, steps);

    ASSERT_NEAR(days, 0.1, 1.0e-12);

    // An interior cell, away from the surface relaxation, must have aged
    // exactly the elapsed time.
    for (int k = 1; k < 5; ++k)
    {
        EXPECT_NEAR(basin.at(index, k), days, 1.0e-9)
            << "layer " << k << " aged " << basin.at(index, k) << " days";
    }
}

TEST(Templates, WaterAgeIsHeldNearZeroAtTheSurface)
{
    // The surface condition, and an honest statement of what it is: a
    // relaxation, so the surface age settles near the contact timescale rather
    // than at exactly zero. The gate is written against that number, not
    // against zero, because claiming zero would be claiming something the
    // implementation does not do.
    WaterAgeTemplate age;
    age.surfaceContactTime = 0.2 / 24.0; // 12 minutes, so the run reaches steady state

    Basin basin;

    std::string message;
    ASSERT_TRUE(age.registerSpecies(basin.registry, message)) << message;

    basin.build(/*layers=*/5, /*totalDepth=*/20.0);

    const int index = basin.solver.speciesIndex(SpeciesRegistry::kAge);

    ASSERT_TRUE(age.configure(basin.solver, message)) << message;
    ASSERT_TRUE(basin.solver.compileReactions(message)) << message;

    basin.run(864.0, 10);

    // Steady state of dA/dt = 1 - A/tau is A = tau, in days.
    const double expected = age.surfaceContactTime;

    EXPECT_NEAR(basin.at(index, 0), expected, 0.05 * expected)
        << "surface age settled at " << basin.at(index, 0) << " days";

    // And it is far below the interior, which is the whole point.
    EXPECT_LT(basin.at(index, 0), 0.1 * basin.at(index, 1));
}

TEST(Templates, TheSurfaceResetDoesNotDependOnLayerThickness)
{
    // The reason the relaxation carries an explicit __THICKNESS__ factor. A
    // Surface-scoped law is areal and the engine divides by the surface cell's
    // thickness, so a rate written without that factor would relax harder in a
    // thin layer — the reset would silently depend on the sigma grid.
    double previous = 0.0;

    for (int layers : {4, 8, 16})
    {
        WaterAgeTemplate age;

        Basin basin;

        std::string message;
        ASSERT_TRUE(age.registerSpecies(basin.registry, message)) << message;

        basin.build(layers, /*totalDepth=*/20.0);

        const int index = basin.solver.speciesIndex(SpeciesRegistry::kAge);

        ASSERT_TRUE(age.configure(basin.solver, message)) << message;
        ASSERT_TRUE(basin.solver.compileReactions(message)) << message;

        basin.run(864.0, 10);

        const double surface = basin.at(index, 0);

        if (previous > 0.0)
        {
            EXPECT_NEAR(surface, previous, 1.0e-6 * previous)
                << layers << " layers gave a surface age of " << surface
                << " against " << previous;
        }

        previous = surface;
    }
}

// ── First-order decay ─────────────────────────────────────────────────────

TEST(Templates, DecayMatchesItsClosedFormInTheUnitsItWasGivenIn)
{
    constexpr double perDay = 5.0;
    constexpr double c0     = 10.0;

    FirstOrderDecayTemplate decay;
    decay.id               = "COLI";
    decay.decayRatePerDay  = perDay;
    decay.temperatureTheta = 1.0; // isolate the rate from the correction
    decay.initialValue     = c0;

    Basin basin;

    std::string message;
    ASSERT_TRUE(decay.registerSpecies(basin.registry, message)) << message;

    basin.build(/*layers=*/4, /*totalDepth=*/12.0);

    const int index = basin.solver.speciesIndex("COLI");
    ASSERT_GE(index, 0);

    basin.setUniform(index, c0);

    ASSERT_TRUE(decay.configure(basin.solver, message)) << message;
    ASSERT_TRUE(basin.solver.compileReactions(message)) << message;

    constexpr double dt    = 864.0;
    constexpr int    steps = 10;
    const double     days  = dt * steps / kSecondsPerDay;

    basin.run(dt, steps);

    ASSERT_NEAR(days, 0.1, 1.0e-12);

    const double exact = c0 * std::exp(-perDay * days);

    for (int k = 0; k < 4; ++k)
        EXPECT_NEAR(basin.at(index, k), exact, 1.0e-6 * c0) << "layer " << k;

    // Worth stating in absolute terms: a tenth of a day at 5/d leaves 61%. A
    // rate divided by 86400 twice, or not at all, would be nowhere near this.
    EXPECT_NEAR(exact / c0, 0.6065, 0.0005);
}

TEST(Templates, TheTemperatureCorrectionIsTheOneTheLiteratureQuotes)
{
    // theta^(T-20) with theta = 1.047 is the standard BOD correction, and it
    // roughly doubles the rate over 15 degrees. Checked against the factor
    // computed here, and against a reference model run at 20 °C so the
    // correction is isolated from everything else.
    constexpr double perDay = 0.3;
    constexpr double theta  = 1.047;
    constexpr double warm   = 30.0;
    constexpr double c0     = 8.0;
    constexpr double dt     = 864.0;
    constexpr int    steps  = 10; // a tenth of a day

    const auto finalValue = [&](double waterTemperature) {
        FirstOrderDecayTemplate decay;
        decay.id               = "CBOD";
        decay.decayRatePerDay  = perDay;
        decay.temperatureTheta = theta;

        Basin basin;

        std::string message;
        EXPECT_TRUE(decay.registerSpecies(basin.registry, message)) << message;
        EXPECT_GE(basin.registry.addTemperature(waterTemperature, message), 0) << message;

        basin.build(/*layers=*/4, /*totalDepth=*/12.0, waterTemperature);

        const int index = basin.solver.speciesIndex("CBOD");
        basin.setUniform(index, c0);

        EXPECT_TRUE(decay.configure(basin.solver, message)) << message;
        EXPECT_TRUE(basin.solver.compileReactions(message)) << message;

        basin.run(dt, steps);

        return basin.at(index, 1);
    };

    const double atWarm = finalValue(warm);
    const double days   = dt * steps / kSecondsPerDay;

    const double correctedRate = perDay * std::pow(theta, warm - 20.0);
    const double exact         = c0 * std::exp(-correctedRate * days);

    EXPECT_NEAR(atWarm, exact, 1.0e-5 * c0);

    // The correction must actually bite: at 30 °C the rate is ~1.58x the 20 °C
    // value, so the warm case must decay visibly further.
    const double atReference = c0 * std::exp(-perDay * days);

    EXPECT_LT(atWarm, 0.985 * atReference)
        << "warm " << atWarm << " against reference " << atReference;
}

TEST(Templates, AZeroRateInstallsNoLawAtAll)
{
    // How a conservative tracer is expressed. Not a rate law that evaluates to
    // zero — no law, so the transport engine carries the constituent and the
    // integrator is never called for it.
    FirstOrderDecayTemplate tracer;
    tracer.id              = "CL";
    tracer.decayRatePerDay = 0.0;

    Basin basin;

    std::string message;
    ASSERT_TRUE(tracer.registerSpecies(basin.registry, message)) << message;

    basin.build(/*layers=*/4, /*totalDepth=*/12.0);

    const int index = basin.solver.speciesIndex("CL");
    basin.setUniform(index, 5.0);

    ASSERT_TRUE(tracer.configure(basin.solver, message)) << message;

    EXPECT_TRUE(basin.solver.reactions().empty())
        << "a conservative constituent installed a rate law";

    basin.run(864.0, 10);

    for (int k = 0; k < 4; ++k)
        EXPECT_NEAR(basin.at(index, k), 5.0, 1.0e-9) << "layer " << k;
}

TEST(Templates, SettlingAndItsBenthicReceiverAreWiredThrough)
{
    constexpr double perDay = 20.0;
    constexpr double c0     = 6.0;

    FirstOrderDecayTemplate solids;
    solids.id                     = "TSS";
    solids.settlingVelocityPerDay = perDay;
    solids.benthicId              = "TSS_BED";

    Basin basin;

    std::string message;
    ASSERT_TRUE(solids.registerSpecies(basin.registry, message)) << message;

    basin.build(/*layers=*/6, /*totalDepth=*/12.0);

    const int index = basin.solver.speciesIndex("TSS");
    const int bed   = basin.solver.speciesIndex("TSS_BED");

    ASSERT_GE(index, 0);
    ASSERT_GE(bed, 0);

    basin.setUniform(index, c0);

    ASSERT_TRUE(solids.configure(basin.solver, message)) << message;

    EXPECT_NEAR(basin.solver.settling().velocity(index), perDay / kSecondsPerDay,
                1.0e-18);
    EXPECT_EQ(basin.solver.settling().benthicReceiver(index), bed);

    const auto totalMass = [&]() {
        double water = 0.0;

        for (int64_t cell = 0; cell < basin.layered.cellCount(); ++cell)
        {
            water += basin.state.speciesAt(index, cell) *
                     basin.layered.cellVolume()[static_cast<size_t>(cell)];
        }

        double onBed = 0.0;

        for (int64_t c = 0; c < basin.layered.columnCount(); ++c)
        {
            onBed += basin.state.speciesAt(bed, basin.layered.cell(c, 5)) *
                     basin.horizontal.faceArea()[static_cast<size_t>(c)];
        }

        return water + onBed;
    };

    const double initial = totalMass();

    basin.run(864.0, 10); // a tenth of a day

    // 20 m/d for a tenth of a day is 2 m of settling through a 12 m column, so
    // a sixth of the mass should have reached the bed.
    EXPECT_NEAR(totalMass(), initial, 1.0e-9 * initial);
    EXPECT_GT(basin.solver.settling().deliveredToBed(index), 0.05 * initial);
}

TEST(Templates, MisconfigurationIsRefused)
{
    Basin basin;
    std::string message;

    FirstOrderDecayTemplate nameless;
    nameless.id.clear();
    EXPECT_FALSE(nameless.registerSpecies(basin.registry, message));

    FirstOrderDecayTemplate negative;
    negative.id              = "X";
    negative.decayRatePerDay = -1.0;

    ASSERT_TRUE(negative.registerSpecies(basin.registry, message)) << message;
    basin.build(/*layers=*/3, /*totalDepth=*/9.0);
    EXPECT_FALSE(negative.configure(basin.solver, message));

    WaterAgeTemplate age;
    age.surfaceContactTime = 0.0;
    EXPECT_FALSE(age.configure(basin.solver, message));
}
