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

// ── Dissolved oxygen suite ────────────────────────────────────────────────

namespace {

using FVQual::Kinetics::DissolvedOxygenTemplate;

/*!
 * \brief Builds a still basin carrying a configured DO suite.
 *
 * Returns by pointer because the solver is not movable, and the tests need the
 * whole assembly to outlive the helper.
 */
struct Suite
{
    Basin                   basin;
    DissolvedOxygenTemplate model;

    int oxygen = -1, cbod = -1, ammonium = -1, nitrate = -1, phosphate = -1, algae = -1;

    void build(int layers, double depth, double waterTemperature = 20.0)
    {
        std::string message;

        ASSERT_TRUE(model.registerSpecies(basin.registry, message)) << message;

        if (model.temperatureTheta != 1.0)
            ASSERT_GE(basin.registry.addTemperature(waterTemperature, message), 0)
                << message;

        basin.build(layers, depth, waterTemperature);

        oxygen = basin.solver.speciesIndex(model.oxygenId);
        cbod   = basin.solver.speciesIndex(model.cbodId);

        if (model.includeNutrients)
        {
            ammonium  = basin.solver.speciesIndex(model.ammoniumId);
            nitrate   = basin.solver.speciesIndex(model.nitrateId);
            phosphate = basin.solver.speciesIndex(model.phosphateId);
        }

        if (model.includeAlgae)
            algae = basin.solver.speciesIndex(model.algaeId);

        ASSERT_TRUE(model.configure(basin.solver, message)) << message;
        ASSERT_TRUE(basin.solver.compileReactions(message)) << message;
    }

    //! Total nitrogen in a cell: both inorganic pools plus what the algae hold.
    [[nodiscard]] double nitrogenAt(int k) const
    {
        return basin.at(ammonium, k) + basin.at(nitrate, k) +
               model.nitrogenPerAlgae * basin.at(algae, k);
    }

    [[nodiscard]] double phosphorusAt(int k) const
    {
        return basin.at(phosphate, k) + model.phosphorusPerAlgae * basin.at(algae, k);
    }
};

} // namespace

TEST(DissolvedOxygen, StreeterPhelpsMatchesTheAnalyticalDeficit)
{
    // The classic benchmark, and the reason the suite can be trusted before any
    // of the biology is switched on. In a well-mixed column with only
    // carbonaceous decay and reaeration, the oxygen deficit has a closed form:
    //
    //   D(t) = kd L0 / (ka - kd) (e^-kd t - e^-ka t) + D0 e^-ka t
    //
    // A single layer, because reaeration is a surface process: in a stratified
    // column only the top layer exchanges, and Streeter-Phelps assumes the whole
    // depth does.
    constexpr double depth      = 3.0;
    constexpr double kd         = 0.3;  // 1/day
    constexpr double velocity   = 1.5;  // m/day
    constexpr double saturation = 9.0;
    constexpr double l0         = 12.0;

    // A small initial deficit, so the oxygen sag the gate is checking is most
    // of what the answer contains. At d0 = 3 the curve only rose by 13% over
    // this interval and the "the deficit actually moved" check below — which is
    // there to stop the analytic comparison passing on a flat line — failed on
    // the fixture rather than on the model.
    constexpr double d0         = 1.0;

    Suite suite;
    suite.model.includeAlgae               = false;
    suite.model.includeNutrients           = false;
    suite.model.temperatureTheta           = 1.0;
    suite.model.cbodDecayPerDay            = kd;
    suite.model.reaerationVelocityPerDay   = velocity;
    suite.model.sedimentOxygenDemandPerDay = 0.0;
    suite.model.saturation                 = saturation;

    suite.build(/*layers=*/1, depth);

    suite.basin.setUniform(suite.oxygen, saturation - d0);
    suite.basin.setUniform(suite.cbod, l0);

    // The first-order reaeration coefficient this configuration implies.
    const double ka = velocity / depth;

    constexpr double dt    = 864.0;
    constexpr int    steps = 20;
    const double     days  = dt * steps / kSecondsPerDay;

    suite.basin.run(dt, steps);

    const double exact = kd * l0 / (ka - kd) * (std::exp(-kd * days) - std::exp(-ka * days)) +
                         d0 * std::exp(-ka * days);

    const double deficit = saturation - suite.basin.at(suite.oxygen, 0);

    EXPECT_NEAR(deficit, exact, 1.0e-6 * saturation)
        << "deficit " << deficit << " against " << exact;

    // The CBOD itself is a plain exponential, which pins the decay coefficient
    // separately from the reaeration.
    EXPECT_NEAR(suite.basin.at(suite.cbod, 0), l0 * std::exp(-kd * days),
                1.0e-7 * l0);

    // And the deficit is well away from both its initial value and zero, so the
    // gate is exercising the curve rather than an endpoint.
    EXPECT_GT(deficit, 1.2 * d0);
}

TEST(DissolvedOxygen, OxygenMinusCarbonaceousDemandIsConservedWithoutExchange)
{
    // With reaeration and sediment demand off, oxidising a unit of CBOD costs
    // exactly a unit of oxygen — so DO - CBOD is an identity of the
    // stoichiometry. It is also the statement that the ultimate deficit equals
    // the initial demand, which is the whole content of Streeter-Phelps.
    Suite suite;
    suite.model.includeAlgae               = false;
    suite.model.includeNutrients           = false;
    suite.model.temperatureTheta           = 1.0;
    suite.model.cbodDecayPerDay            = 2.0;
    suite.model.reaerationVelocityPerDay   = 0.0;
    suite.model.sedimentOxygenDemandPerDay = 0.0;

    suite.build(/*layers=*/4, 12.0);

    suite.basin.setUniform(suite.oxygen, 9.0);
    suite.basin.setUniform(suite.cbod, 5.0);

    const double invariant = 9.0 - 5.0;

    for (int i = 0; i < 10; ++i)
    {
        suite.basin.run(864.0, 1);

        for (int k = 0; k < 4; ++k)
        {
            EXPECT_NEAR(suite.basin.at(suite.oxygen, k) - suite.basin.at(suite.cbod, k),
                        invariant, 1.0e-9)
                << "step " << i << " layer " << k;
        }
    }

    // And it moved: most of the demand should be gone after a day at 2/d.
    EXPECT_LT(suite.basin.at(suite.cbod, 0), 0.85 * 5.0);
}

TEST(DissolvedOxygen, NitrificationCostsTheStoichiometricOxygen)
{
    // 4.57 g O₂ per g N is the number every reference quotes, and it has to
    // appear as the ratio between what leaves ammonium and what leaves oxygen.
    constexpr double ratio = 4.57;

    Suite suite;
    suite.model.includeAlgae               = false;
    suite.model.temperatureTheta           = 1.0;
    suite.model.cbodDecayPerDay            = 0.0;
    suite.model.reaerationVelocityPerDay   = 0.0;
    suite.model.sedimentOxygenDemandPerDay = 0.0;
    suite.model.nitrificationPerDay        = 1.0;
    suite.model.oxygenPerNitrogen          = ratio;

    suite.build(/*layers=*/3, 9.0);

    constexpr double nh0 = 0.5;
    constexpr double do0 = 9.0;

    suite.basin.setUniform(suite.oxygen, do0);
    suite.basin.setUniform(suite.ammonium, nh0);
    suite.basin.setUniform(suite.nitrate, 0.0);
    suite.basin.setUniform(suite.phosphate, 0.0);
    suite.basin.setUniform(suite.cbod, 0.0);

    suite.basin.run(864.0, 10);

    const double nitrified = nh0 - suite.basin.at(suite.ammonium, 1);
    const double consumed  = do0 - suite.basin.at(suite.oxygen, 1);

    ASSERT_GT(nitrified, 0.02) << "nothing was nitrified";

    EXPECT_NEAR(consumed / nitrified, ratio, 1.0e-6 * ratio);

    // Nitrogen moves between pools rather than disappearing.
    EXPECT_NEAR(suite.basin.at(suite.ammonium, 1) + suite.basin.at(suite.nitrate, 1),
                nh0, 1.0e-9);
}

TEST(DissolvedOxygen, NitrogenAndPhosphorusBudgetsClose)
{
    // The gate the whole template exists for. Growth, respiration, mortality and
    // nitrification each move nutrients between pools, and every one of those
    // terms has to carry the same stoichiometric coefficient or the model
    // quietly manufactures nitrogen. Checked every step, and with settling off
    // so nothing legitimately leaves the water column.
    Suite suite;
    suite.model.temperatureTheta       = 1.0;
    suite.model.algalSettlingPerDay    = 0.0;
    suite.model.algalGrowthPerDay      = 2.0;
    suite.model.algalRespirationPerDay = 0.3;
    suite.model.algalMortalityPerDay   = 0.2;
    suite.model.nitrificationPerDay    = 0.5;

    suite.build(/*layers=*/4, 8.0);

    // Enough light for growth to be limited by nutrients rather than switched
    // off entirely.
    suite.basin.solver.lightOptions().incidentShortwave     = 400.0;
    suite.basin.solver.lightOptions().surfaceAbsorption     = 0.0;
    suite.basin.solver.lightOptions().backgroundAttenuation = 0.2;

    suite.basin.setUniform(suite.oxygen, 9.0);
    suite.basin.setUniform(suite.cbod, 1.0);
    suite.basin.setUniform(suite.ammonium, 0.4);
    suite.basin.setUniform(suite.nitrate, 0.2);
    suite.basin.setUniform(suite.phosphate, 0.05);
    suite.basin.setUniform(suite.algae, 2.0);

    const double nitrogen   = suite.nitrogenAt(1);
    const double phosphorus = suite.phosphorusAt(1);

    ASSERT_GT(nitrogen, 0.0);
    ASSERT_GT(phosphorus, 0.0);

    for (int i = 0; i < 10; ++i)
    {
        suite.basin.run(864.0, 1);

        for (int k = 0; k < 4; ++k)
        {
            EXPECT_NEAR(suite.nitrogenAt(k), nitrogen, 1.0e-8 * nitrogen)
                << "nitrogen, step " << i << " layer " << k;
            EXPECT_NEAR(suite.phosphorusAt(k), phosphorus, 1.0e-8 * phosphorus)
                << "phosphorus, step " << i << " layer " << k;
        }
    }

    // The budgets must be closing over real activity, not over a model that sat
    // still: the algae should have grown measurably.
    EXPECT_GT(suite.basin.at(suite.algae, 0), 1.05 * 2.0)
        << "algae did not grow; the budget gate is vacuous";
}

TEST(DissolvedOxygen, GrowthIsTheProductOfItsLimitationFactors)
{
    // Each limitation factor computed independently here and multiplied, against
    // what the configured network actually does over one short step. This is
    // where a mis-parenthesised expression shows up — the budgets above would
    // still close if growth were limited by the wrong product.
    constexpr double muMax  = 2.0;
    constexpr double ksLight = 100.0;
    constexpr double ksN     = 0.02;
    constexpr double ksP     = 0.004;
    constexpr double nh4     = 0.3;
    constexpr double no3     = 0.1;
    constexpr double po4     = 0.02;
    constexpr double algae0  = 1.0;
    constexpr double incident = 500.0;
    constexpr double kd       = 0.25;
    constexpr double depth    = 8.0;
    constexpr int    layers   = 4;

    Suite suite;
    suite.model.temperatureTheta          = 1.0;
    suite.model.algalGrowthPerDay         = muMax;
    suite.model.algalRespirationPerDay    = 0.0;
    suite.model.algalMortalityPerDay      = 0.0;
    suite.model.algalSettlingPerDay       = 0.0;
    suite.model.algalLightAttenuation     = 0.0; // no self-shading, so light is known
    suite.model.lightHalfSaturation       = ksLight;
    suite.model.nitrogenHalfSaturation    = ksN;
    suite.model.phosphorusHalfSaturation  = ksP;
    suite.model.nitrificationPerDay       = 0.0;
    suite.model.cbodDecayPerDay           = 0.0;
    suite.model.reaerationVelocityPerDay  = 0.0;
    suite.model.sedimentOxygenDemandPerDay = 0.0;

    suite.build(layers, depth);

    suite.basin.solver.lightOptions().incidentShortwave     = incident;
    suite.basin.solver.lightOptions().surfaceAbsorption     = 0.0;
    suite.basin.solver.lightOptions().backgroundAttenuation = kd;

    suite.basin.setUniform(suite.oxygen, 9.0);
    suite.basin.setUniform(suite.cbod, 0.0);
    suite.basin.setUniform(suite.ammonium, nh4);
    suite.basin.setUniform(suite.nitrate, no3);
    suite.basin.setUniform(suite.phosphate, po4);
    suite.basin.setUniform(suite.algae, algae0);

    // One short step, so the concentrations the factors depend on have barely
    // moved and the linearisation is accurate to well inside the tolerance.
    constexpr double dt = 0.5;

    suite.basin.run(dt, 1);

    const double h   = depth / static_cast<double>(layers);
    const double tau = kd * h;

    for (int k = 0; k < layers; ++k)
    {
        // The layer-mean irradiance, from Beer-Lambert, computed here.
        const double top   = incident * std::exp(-kd * static_cast<double>(k) * h);
        const double light = top * (1.0 - std::exp(-tau)) / tau;

        const double fL = light / (ksLight + light);
        const double fN = (nh4 + no3) / (ksN + nh4 + no3);
        const double fP = po4 / (ksP + po4);

        const double rate = muMax / kSecondsPerDay * fL * fN * fP;

        // The closed form, not `algae0 + rate * dt`. Growth is exponential, and
        // the linear estimate carries a truncation error of (r dt)^2/2 that a
        // tolerance scaled to the increment sits below — the same mistake this
        // suite has now made twice, so it is written out rather than fixed
        // quietly.
        const double expected = algae0 * std::exp(rate * dt);

        EXPECT_NEAR(suite.basin.at(suite.algae, k), expected, 1.0e-9)
            << "layer " << k << ": fL " << fL << " fN " << fN << " fP " << fP;
    }

    // Oxygen is produced in proportion, which pins the photosynthetic
    // stoichiometry against the same growth.
    const double grown = suite.basin.at(suite.algae, 0) - algae0;

    EXPECT_NEAR(suite.basin.at(suite.oxygen, 0) - 9.0,
                suite.model.oxygenPerAlgae * grown, 1.0e-6 * suite.model.oxygenPerAlgae * grown);
}

TEST(DissolvedOxygen, WithoutNutrientsGrowthIsLimitedByLightAlone)
{
    // The reduced configuration has to be a real model rather than a broken
    // one: no nutrient species exist, so no rate law may reference them.
    Suite suite;
    suite.model.includeNutrients = false;
    suite.model.temperatureTheta = 1.0;

    suite.build(/*layers=*/3, 9.0);

    EXPECT_LT(suite.ammonium, 0);
    EXPECT_GE(suite.algae, 0);

    suite.basin.solver.lightOptions().incidentShortwave = 300.0;

    suite.basin.setUniform(suite.algae, 1.0);
    suite.basin.setUniform(suite.oxygen, 9.0);

    suite.basin.run(864.0, 5);

    EXPECT_GT(suite.basin.at(suite.algae, 0), 1.0) << "algae did not grow on light alone";
}

TEST(DissolvedOxygen, InTheDarkAlgaeOnlyLose)
{
    // Zero irradiance means zero growth, so respiration and mortality are all
    // that is left. Worth pinning because the light factor sits in a
    // denominator, and a formulation that read `Ks/(Ks+I)` by mistake would
    // grow fastest at night.
    Suite suite;
    suite.model.temperatureTheta       = 1.0;
    suite.model.algalSettlingPerDay    = 0.0;
    suite.model.algalRespirationPerDay = 0.5;
    suite.model.algalMortalityPerDay   = 0.5;

    suite.build(/*layers=*/3, 9.0);

    suite.basin.solver.lightOptions().incidentShortwave = 0.0;

    suite.basin.setUniform(suite.oxygen, 9.0);
    suite.basin.setUniform(suite.ammonium, 0.3);
    suite.basin.setUniform(suite.nitrate, 0.1);
    suite.basin.setUniform(suite.phosphate, 0.02);
    suite.basin.setUniform(suite.algae, 3.0);

    suite.basin.run(864.0, 10);

    EXPECT_LT(suite.basin.at(suite.algae, 1), 3.0) << "algae grew in the dark";

    // Respiration consumes oxygen; with no photosynthesis the water column can
    // only lose it.
    EXPECT_LT(suite.basin.at(suite.oxygen, 1), 9.0);
}

TEST(DissolvedOxygen, MisconfigurationIsRefused)
{
    Suite suite;
    suite.model.temperatureTheta = 1.0;
    suite.build(/*layers=*/3, 9.0);

    std::string message;

    DissolvedOxygenTemplate bad;
    bad.temperatureTheta = 1.0;
    bad.cbodDecayPerDay  = -1.0;
    EXPECT_FALSE(bad.configure(suite.basin.solver, message));

    DissolvedOxygenTemplate zeroHalfSaturation;
    zeroHalfSaturation.temperatureTheta    = 1.0;
    zeroHalfSaturation.lightHalfSaturation = 0.0;
    EXPECT_FALSE(zeroHalfSaturation.configure(suite.basin.solver, message));
}
