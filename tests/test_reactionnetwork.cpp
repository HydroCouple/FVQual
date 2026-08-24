// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_reactionnetwork.cpp
 * \brief Verification of operator-split reaction kinetics (F6.1b).
 *
 * Two kinds of claim, and they are not the same kind.
 *
 * The engine's own claims are testable as identities: a rate law integrates to
 * its closed form, an unresolvable name is refused at configuration rather than
 * silently read as zero, and a species with no rate law is left alone.
 *
 * The claim a *network* makes is different. Whether constituent mass is
 * conserved depends entirely on whether the rate laws sum to zero across
 * species, which is the modeller's arithmetic and not the engine's. So the
 * conservation gate here does not test the engine — it tests that a network
 * written to conserve actually does, and that one written not to is *allowed*
 * to lose mass rather than being silently corrected. Both are legitimate models,
 * and an engine that quietly fixed the second would be lying about the first.
 *
 * **Ablations these gates were checked against**, because a gate that cannot
 * fail is decoration:
 *
 * - Reaction operator removed from `step()` → three solver gates fail, and the
 *   `profile().reactions` check catches the case where the fixture's rate
 *   constant is what is really zero.
 * - Pre-step compile guard removed → the refusal gate fails, because the step
 *   transports before discovering it cannot react.
 * - `bound != used` check removed from `compile()` → the unresolvable-name gate
 *   fails, and the typo becomes a rate law that evaluates as zero.
 * - Per-thread workspace slicing replaced by a shared slot, at 4 threads →
 *   four gates fail. The threading is verified, not assumed.
 */

#include "fvqual/hydro/internalmode.h"
#include "fvqual/kinetics/reactionnetwork.h"
#include "fvqual/mesh/horizontalmesh.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/mesh/sigmaspec.h"
#include "fvqual/state/speciesregistry.h"
#include "fvqual/state/statestore.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

using FVQual::Hydro::InternalModeOptions;
using FVQual::Hydro::InternalModeSolver;
using FVQual::Kinetics::ReactionNetwork;
using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;
using FVQual::State::SpeciesDefinition;
using FVQual::State::SpeciesRegistry;
using FVQual::State::StateStore;
namespace Testing = FVQual::Testing;

namespace {

/*!
 * \brief A small mesh with a species set, enough to integrate a network on.
 */
struct Box
{
    HorizontalMesh  horizontal;
    LayeredMesh     layered;
    SpeciesRegistry registry;
    StateStore      state;

    void build(int layers = 4)
    {
        std::string message;

        ASSERT_TRUE(horizontal.build(Testing::channelStrip(3, 300.0, 40.0), message))
            << message;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = -10.0;

        SigmaSpec spec;
        spec.layerCount        = layers;
        spec.minLayerThickness = 0.0;

        ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;
    }

    int addSpecies(const std::string &id, double initial)
    {
        SpeciesDefinition definition;
        definition.id           = id;
        definition.units        = "g/m3";
        definition.initialValue = initial;

        std::string message;
        const int   index = registry.add(definition, message);

        EXPECT_GE(index, 0) << message;
        return index;
    }

    void allocate() { state.allocate(layered, registry); }

    void setUniform(int species, double value)
    {
        for (int64_t i = 0; i < layered.cellCount(); ++i)
            state.speciesAt(species, i) = value;
    }

    [[nodiscard]] double at(int species, int64_t cell) const
    {
        return state.speciesAt(species, cell);
    }

    [[nodiscard]] int64_t cells() const { return layered.cellCount(); }
};

} // namespace

// ── Configuration ─────────────────────────────────────────────────────────

TEST(ReactionNetwork, MalformedConfigurationIsRefused)
{
    Box box;
    box.build();
    box.addSpecies("A", 1.0);

    ReactionNetwork network;
    std::string     message;

    // A rate law that is not a valid expression.
    EXPECT_FALSE(network.addRate("A", "0 - k *", message));
    EXPECT_FALSE(message.empty());

    // A species with no name.
    EXPECT_FALSE(network.addRate("", "0 - A", message));

    // Two laws for one species: ambiguous, not additive.
    ASSERT_TRUE(network.addRate("A", "0 - 0.1 * A", message)) << message;
    EXPECT_FALSE(network.addRate("A", "0 - 0.2 * A", message));

    // A duplicate constant.
    ASSERT_TRUE(network.addConstant("k", 0.1, message)) << message;
    EXPECT_FALSE(network.addConstant("k", 0.2, message));

    // Integrating before compiling is a programming error, not a silent no-op.
    box.allocate();
    EXPECT_FALSE(network.integrate(box.registry, box.state, box.cells(), 1.0, message));
}

TEST(ReactionNetwork, AnUnresolvableNameIsCaughtAtCompileTime)
{
    // The failure this check exists for: `bind_variables` leaves an unknown
    // name unbound, and an unbound variable evaluates as zero. A typo would
    // otherwise turn "0 - kdecay * A" into "0", producing a conservative
    // tracer that the modeller believes is decaying.
    Box box;
    box.build();
    box.addSpecies("A", 1.0);

    ReactionNetwork network;
    std::string     message;

    ASSERT_TRUE(network.addRate("A", "0 - kdecay * A", message)) << message;

    EXPECT_FALSE(network.compile(box.registry, message))
        << "an undefined constant was accepted";
    EXPECT_NE(message.find("kdecay"), std::string::npos) << message;

    // Declaring it resolves the network.
    ASSERT_TRUE(network.addConstant("kdecay", 0.15, message)) << message;
    EXPECT_TRUE(network.compile(box.registry, message)) << message;
}

TEST(ReactionNetwork, ARateLawForAnUnregisteredSpeciesIsRefused)
{
    Box box;
    box.build();
    box.addSpecies("A", 1.0);

    ReactionNetwork network;
    std::string     message;

    ASSERT_TRUE(network.addRate("B", "0 - 0.1 * B", message)) << message;

    EXPECT_FALSE(network.compile(box.registry, message));
    EXPECT_NE(message.find('B'), std::string::npos) << message;
}

// ── Integration against closed forms ──────────────────────────────────────

TEST(ReactionNetwork, FirstOrderDecayMatchesItsClosedForm)
{
    constexpr double k  = 0.2;
    constexpr double c0 = 8.0;
    constexpr double dt = 300.0;

    Box box;
    box.build();
    const int a = box.addSpecies("A", c0);
    box.allocate();
    box.setUniform(a, c0);

    ReactionNetwork network;
    std::string     message;

    ASSERT_TRUE(network.addConstant("k", k, message)) << message;
    ASSERT_TRUE(network.addRate("A", "0 - k * A", message)) << message;
    ASSERT_TRUE(network.compile(box.registry, message)) << message;

    // Stepped repeatedly, as the solver will step it, so a splitting error
    // would accumulate rather than cancel.
    constexpr int steps = 8;

    for (int i = 0; i < steps; ++i)
    {
        ASSERT_TRUE(network.integrate(box.registry, box.state, box.cells(), dt, message))
            << message;
    }

    const double exact = c0 * std::exp(-k * dt * steps);

    for (int64_t cell = 0; cell < box.cells(); ++cell)
    {
        EXPECT_NEAR(box.at(a, cell), exact, 1.0e-8 * c0)
            << "cell " << cell << " reached " << box.at(a, cell) << " against " << exact;
    }

    EXPECT_EQ(network.diagnostics().negativeCells, 0);
    EXPECT_EQ(network.diagnostics().unconvergedCells, 0);
}

TEST(ReactionNetwork, ASpeciesWithNoRateLawIsLeftAlone)
{
    // How a conservative tracer is expressed: by saying nothing about it. If
    // the network touched it, every inert constituent would need a "0" rate law
    // and forgetting one would be a silent decay.
    Box box;
    box.build();
    const int a     = box.addSpecies("A", 5.0);
    const int inert = box.addSpecies("Inert", 3.0);
    box.allocate();
    box.setUniform(a, 5.0);
    box.setUniform(inert, 3.0);

    ReactionNetwork network;
    std::string     message;

    ASSERT_TRUE(network.addRate("A", "0 - 0.5 * A", message)) << message;
    ASSERT_TRUE(network.compile(box.registry, message)) << message;
    ASSERT_TRUE(network.integrate(box.registry, box.state, box.cells(), 100.0, message))
        << message;

    EXPECT_LT(box.at(a, 0), 5.0) << "the reacting species should have decayed";

    for (int64_t cell = 0; cell < box.cells(); ++cell)
        EXPECT_EQ(box.at(inert, cell), 3.0) << "the inert species moved";
}

TEST(ReactionNetwork, ARateLawMayDependOnAnotherSpecies)
{
    // Temperature is a transported species here, so Arrhenius dependence is
    // something the modeller writes rather than something the engine special
    // cases. This is the mechanism that makes that true.
    constexpr double k20   = 0.3;
    constexpr double theta = 1.047;
    constexpr double c0    = 6.0;
    constexpr double dt    = 200.0;

    Box box;
    box.build();
    const int a = box.addSpecies("A", c0);
    const int t = box.addSpecies("T", 25.0);
    box.allocate();
    box.setUniform(a, c0);
    box.setUniform(t, 25.0);

    ReactionNetwork network;
    std::string     message;

    ASSERT_TRUE(network.addConstant("k20", k20, message)) << message;
    ASSERT_TRUE(network.addConstant("theta", theta, message)) << message;
    ASSERT_TRUE(network.addRate("A", "0 - k20 * theta ^ (T - 20) * A", message)) << message;
    ASSERT_TRUE(network.compile(box.registry, message)) << message;
    ASSERT_TRUE(network.integrate(box.registry, box.state, box.cells(), dt, message))
        << message;

    const double rate  = k20 * std::pow(theta, 25.0 - 20.0);
    const double exact = c0 * std::exp(-rate * dt);

    EXPECT_NEAR(box.at(a, 0), exact, 1.0e-8 * c0);

    // Temperature has no rate law, so it must be untouched — the dependence is
    // one-way.
    EXPECT_EQ(box.at(t, 0), 25.0);
}

// ── Conservation: a property of the network, not of the engine ─────────────

TEST(ReactionNetwork, AClosedNetworkConservesTheSumItWasWrittenToConserve)
{
    // A -> B -> C with the rate laws summing to zero. The invariant belongs to
    // the network's arithmetic, and the engine's job is only not to spoil it.
    constexpr double k1 = 0.02;
    constexpr double k2 = 0.008;
    constexpr double a0 = 10.0;

    Box box;
    box.build();
    const int a = box.addSpecies("A", a0);
    const int b = box.addSpecies("B", 0.0);
    const int c = box.addSpecies("C", 0.0);
    box.allocate();
    box.setUniform(a, a0);
    box.setUniform(b, 0.0);
    box.setUniform(c, 0.0);

    ReactionNetwork network;
    std::string     message;

    ASSERT_TRUE(network.addConstant("k1", k1, message)) << message;
    ASSERT_TRUE(network.addConstant("k2", k2, message)) << message;
    ASSERT_TRUE(network.addRate("A", "0 - k1 * A", message)) << message;
    ASSERT_TRUE(network.addRate("B", "k1 * A - k2 * B", message)) << message;
    ASSERT_TRUE(network.addRate("C", "k2 * B", message)) << message;
    ASSERT_TRUE(network.compile(box.registry, message)) << message;

    constexpr double dt    = 60.0;
    constexpr int    steps = 20;

    for (int i = 0; i < steps; ++i)
    {
        ASSERT_TRUE(network.integrate(box.registry, box.state, box.cells(), dt, message))
            << message;

        // Checked every step, not just at the end: a drift that cancels by
        // coincidence at the final time would otherwise pass.
        const double total = box.at(a, 0) + box.at(b, 0) + box.at(c, 0);

        EXPECT_NEAR(total, a0, 1.0e-9 * a0) << "after step " << i;
    }

    // And the trajectory is right, not merely conservative — a network that
    // moved nothing at all would also conserve.
    const double t      = dt * steps;
    const double exactA = a0 * std::exp(-k1 * t);
    const double exactB = a0 * k1 * (std::exp(-k1 * t) - std::exp(-k2 * t)) / (k2 - k1);

    EXPECT_NEAR(box.at(a, 0), exactA, 1.0e-7 * a0);
    EXPECT_NEAR(box.at(b, 0), exactB, 1.0e-7 * a0);
}

TEST(ReactionNetwork, AnOpenNetworkIsAllowedToLoseMass)
{
    // The counterpart, and the reason the engine does not police conservation.
    // Simple first-order decay to nothing is a perfectly ordinary model: BOD
    // exerted, a tracer lost to the sediment. If the engine "corrected" this it
    // would be refusing to run the model it was given.
    Box box;
    box.build();
    const int a = box.addSpecies("A", 4.0);
    box.allocate();
    box.setUniform(a, 4.0);

    ReactionNetwork network;
    std::string     message;

    ASSERT_TRUE(network.addRate("A", "0 - 0.01 * A", message)) << message;
    ASSERT_TRUE(network.compile(box.registry, message)) << message;
    ASSERT_TRUE(network.integrate(box.registry, box.state, box.cells(), 600.0, message))
        << message;

    EXPECT_LT(box.at(a, 0), 4.0) << "the decay was suppressed";
    EXPECT_GT(box.at(a, 0), 0.0);
}

TEST(ReactionNetwork, AZeroOrderSinkGoesNegativeAndSaysSo)
{
    // A rate law that does not stop at zero will cross it. The engine does not
    // clamp: a clamp restores plausibility while destroying the mass balance,
    // and a modeller who cannot see it happening cannot fix the law that caused
    // it. What the engine owes is a report.
    Box box;
    box.build();
    const int a = box.addSpecies("A", 1.0);
    box.allocate();
    box.setUniform(a, 1.0);

    ReactionNetwork network;
    std::string     message;

    // Constant removal, independent of concentration — it runs straight past
    // zero.
    ASSERT_TRUE(network.addRate("A", "0 - 0.01", message)) << message;
    ASSERT_TRUE(network.compile(box.registry, message)) << message;
    ASSERT_TRUE(network.integrate(box.registry, box.state, box.cells(), 500.0, message))
        << message;

    EXPECT_LT(box.at(a, 0), 0.0) << "expected this fixture to cross zero";

    EXPECT_LT(network.diagnostics().minimumValue, 0.0)
        << "the excursion was not reported";
    EXPECT_EQ(network.diagnostics().negativeCells, box.cells())
        << "every cell should have been counted";
}

// ── Through the full solver ───────────────────────────────────────────────

namespace {

/*!
 * \brief A closed stratified basin with a decaying tracer.
 *
 * Deliberately a *whole* model rather than the network in isolation: the
 * question the gates below ask is whether the split operator is placed
 * correctly in the step, and that cannot be answered without the step.
 */
struct Basin
{
    HorizontalMesh     horizontal;
    LayeredMesh        layered;
    SpeciesRegistry    registry;
    StateStore         state;
    InternalModeSolver solver;

    int temperature = -1;
    int dye         = -1;

    /*!
     * \brief Settable before build(), so a test can ask for a still column.
     *
     * The per-layer light gates need vertical mixing switched off: a neutrally
     * stratified column gets the *maximum* Pacanowski-Philander diffusivity,
     * which would blend the differential growth the gate is measuring.
     */
    InternalModeOptions options = [] {
        InternalModeOptions defaults;
        defaults.gravity         = 9.80665;
        defaults.manningN        = 0.0;
        defaults.cflNumber       = 0.5;
        defaults.surfaceHeatFlux = 0.0;
        return defaults;
    }();

    void suppressVerticalMixing()
    {
        options.neutralViscosity      = 0.0;
        options.backgroundViscosity   = 0.0;
        options.backgroundDiffusivity = 0.0;
        options.convectiveMixing      = 0.0;
    }

    //! \param front Initial temperature as a function of x — a step drives motion.
    template <typename TemperatureFunction, typename DyeFunction>
    void build(int columns, int layers, double length, double width, double depth,
               TemperatureFunction front, DyeFunction dyeField)
    {
        std::string message;

        ASSERT_TRUE(horizontal.build(Testing::channelStrip(columns, length, width),
                                     message))
            << message;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = -depth;

        SigmaSpec spec;
        spec.layerCount        = layers;
        spec.minLayerThickness = 0.0;

        ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

        ASSERT_GE(registry.addTemperature(10.0, message), 0) << message;
        temperature = registry.temperatureIndex();

        SpeciesDefinition definition;
        definition.id           = "DYE";
        definition.units        = "g/m3";
        definition.initialValue = 0.0;

        ASSERT_GE(registry.add(definition, message), 0) << message;
        dye = registry.indexOf("DYE");

        state.allocate(layered, registry);

        for (int64_t c = 0; c < layered.columnCount(); ++c)
        {
            const double x = horizontal.faceCentroidX()[static_cast<size_t>(c)];

            for (int k = 0; k < layered.layerCount(); ++k)
            {
                const int64_t cell = layered.cell(c, k);

                state.speciesAt(temperature, cell) = front(x);
                state.speciesAt(dye, cell)         = dyeField(x);
            }
        }

        ASSERT_TRUE(solver.initialize(layered, state, registry, options, message))
            << message;
    }

    [[nodiscard]] double peakSpeed() const
    {
        double peak = 0.0;

        for (int64_t i = 0; i < state.cellCount(); ++i)
        {
            peak = std::max(peak, std::hypot(state.u[static_cast<size_t>(i)],
                                             state.v[static_cast<size_t>(i)]));
        }

        return peak;
    }
};

//! A horizontally uniform field, whatever x is.
auto uniform(double value)
{
    return [value](double) { return value; };
}

} // namespace

TEST(ReactionNetworkSolver, AnUncompiledNetworkRefusesTheStepBeforeChangingAnything)
{
    // Two claims, and the second is the one that needs a fixture with motion in
    // it. Refusing is easy — the network itself refuses to integrate when it has
    // not been compiled. Refusing *before* the step does any work is the point:
    // otherwise the caller is handed an error alongside a state that has been
    // transported but not reacted, which is worse than either outcome alone.
    //
    // So the basin carries a temperature front. A step here demonstrably moves
    // water and dye, which means an unchanged state after the refusal is
    // evidence rather than a tautology.
    Basin basin;
    basin.build(20, 6, 200.0, 40.0, 10.0,
                [](double x) { return x < 100.0 ? 8.0 : 18.0; },
                [](double x) { return x < 100.0 ? 2.0 : 0.0; });

    std::string message;
    ASSERT_TRUE(basin.solver.reactions().addRate("DYE", "0 - 0.001 * DYE", message))
        << message;

    const std::vector<double> etaBefore = basin.state.eta;

    std::vector<double> dyeBefore;
    for (int64_t cell = 0; cell < basin.state.cellCount(); ++cell)
        dyeBefore.push_back(basin.state.speciesAt(basin.dye, cell));

    EXPECT_FALSE(basin.solver.step(10.0, message));
    EXPECT_FALSE(message.empty());

    EXPECT_EQ(basin.state.eta, etaBefore) << "the refused step moved the free surface";

    for (int64_t cell = 0; cell < basin.state.cellCount(); ++cell)
    {
        ASSERT_EQ(basin.state.speciesAt(basin.dye, cell),
                  dyeBefore[static_cast<size_t>(cell)])
            << "the refused step transported cell " << cell;
    }

    // Compiled, the same step runs — and does change things, which is what
    // makes the comparison above meaningful.
    ASSERT_TRUE(basin.solver.compileReactions(message)) << message;
    ASSERT_TRUE(basin.solver.step(10.0, message)) << message;

    EXPECT_NE(basin.state.eta, etaBefore) << "this fixture does not move; the "
                                             "unchanged-state check above proves nothing";
}

TEST(ReactionNetworkSolver, DecayThroughTheFullSolverMatchesItsClosedForm)
{
    // Unstratified and at rest, so transport moves nothing and the tracer's
    // trajectory is the reaction's alone. Any deviation from exp(-kt) is the
    // splitting or the placement, not the flow.
    constexpr double k     = 0.004;
    constexpr double c0    = 5.0;
    constexpr double dt    = 20.0;
    constexpr int    steps = 30;

    Basin basin;
    basin.build(6, 5, 300.0, 40.0, 12.0, uniform(12.0), uniform(c0));

    std::string message;

    ASSERT_TRUE(basin.solver.reactions().addConstant("k", k, message)) << message;
    ASSERT_TRUE(basin.solver.reactions().addRate("DYE", "0 - k * DYE", message))
        << message;
    ASSERT_TRUE(basin.solver.compileReactions(message)) << message;

    for (int i = 0; i < steps; ++i)
        ASSERT_TRUE(basin.solver.step(dt, message)) << message;

    const double exact = c0 * std::exp(-k * dt * steps);

    for (int64_t cell = 0; cell < basin.state.cellCount(); ++cell)
    {
        EXPECT_NEAR(basin.state.speciesAt(basin.dye, cell), exact, 1.0e-6 * c0)
            << "cell " << cell;
    }

    // The reactions were actually reached — without this the test above would
    // also pass on a step() that skipped them and a fixture with k = 0.
    EXPECT_GT(basin.solver.profile().reactions, 0)
        << "no time was spent in the reaction operator";

    // Temperature carries no rate law, so it must be exactly where it started.
    EXPECT_NEAR(basin.state.speciesAt(basin.temperature, 0), 12.0, 1.0e-9);
}

TEST(ReactionNetworkSolver, LinearDecayInAStirredBasinFollowsTheBulkLaw)
{
    // The gate that transport and reaction actually compose. For a linear rate
    // law, d/dt of the *volume integral* is -k times that integral no matter
    // what the flow does, because advection only rearranges. So the total mass
    // has a closed form even though no individual cell does — and a split that
    // double-counted or dropped a step would break it.
    constexpr double k = 0.003;

    Basin basin;

    // A temperature front drives a gravity current; the dye starts on one side
    // so it is genuinely advected rather than sitting still.
    basin.build(40, 10, 400.0, 40.0, 10.0,
                [](double x) { return x < 200.0 ? 8.0 : 18.0; },
                [](double x) { return x < 200.0 ? 4.0 : 0.0; });

    std::string message;

    ASSERT_TRUE(basin.solver.reactions().addConstant("k", k, message)) << message;
    ASSERT_TRUE(basin.solver.reactions().addRate("DYE", "0 - k * DYE", message))
        << message;
    ASSERT_TRUE(basin.solver.compileReactions(message)) << message;

    const double initialMass = basin.solver.totalSpeciesMass(basin.dye);
    ASSERT_GT(initialMass, 0.0);

    constexpr double duration = 300.0;

    ASSERT_TRUE(basin.solver.advance(duration, message)) << message;

    // The fixture has to be stirring, or this gate degenerates into the
    // previous one.
    EXPECT_GT(basin.peakSpeed(), 0.01)
        << "the basin never moved; this is not a transport-reaction test";

    const double expected = initialMass * std::exp(-k * duration);
    const double mass     = basin.solver.totalSpeciesMass(basin.dye);

    EXPECT_NEAR(mass, expected, 1.0e-5 * initialMass)
        << "mass " << mass << " against " << expected;

    // And the dye did move. Measured as the mass that crossed into the
    // initially clean half, which starts at exactly zero — so any positive
    // value is transport and nothing else.
    //
    // The first version of this check asked for a positive concentration in
    // *every* cell, which failed: 300 s of a gravity current does not carry dye
    // 195 m to the far wall, so the far cells were legitimately still clean.
    // That was a wrong expectation about the fixture rather than a defect.
    double crossed = 0.0;

    for (int64_t c = 0; c < basin.layered.columnCount(); ++c)
    {
        if (basin.horizontal.faceCentroidX()[static_cast<size_t>(c)] < 200.0)
            continue;

        for (int k = 0; k < basin.layered.layerCount(); ++k)
        {
            const int64_t cell = basin.layered.cell(c, k);

            crossed += basin.state.speciesAt(basin.dye, cell) *
                       basin.layered.cellVolume()[static_cast<size_t>(cell)];
        }
    }

    EXPECT_GT(crossed, 1.0e-3 * mass)
        << "no dye crossed the front; the bulk law was satisfied by a still basin";
}

// ── Derived environmental fields ──────────────────────────────────────────

TEST(ReactionNetwork, AReferencedFieldThatWasNotSuppliedIsRefused)
{
    // The same discipline as an unbound name, one level up. An empty light
    // field would read as zero irradiance, which is not an error anyone
    // notices — it is a model in which photosynthesis quietly never happens.
    Box box;
    box.build();
    box.addSpecies("A", 1.0);
    box.allocate();

    ReactionNetwork network;
    std::string     message;

    ASSERT_TRUE(network.addRate("A", "0 - 0.001 * __LIGHT__ * A", message)) << message;
    ASSERT_TRUE(network.compile(box.registry, message)) << message;

    EXPECT_TRUE(network.usesLight());

    EXPECT_FALSE(network.integrate(box.registry, box.state, box.cells(), 1.0, message));
    EXPECT_NE(message.find("__LIGHT__"), std::string::npos) << message;

    // Supplied, it runs.
    const std::vector<double> light(static_cast<size_t>(box.cells()), 200.0);

    FVQual::Kinetics::DerivedFields fields;
    fields.light = light;

    EXPECT_TRUE(network.integrate(box.registry, box.state, box.cells(), 1.0, message,
                                  fields))
        << message;

    // A field shorter than the mesh is refused too — it would be a stale array
    // from a different mesh, read past its end.
    const std::vector<double> tooShort(static_cast<size_t>(box.cells()) - 1, 200.0);

    fields.light = tooShort;

    EXPECT_FALSE(network.integrate(box.registry, box.state, box.cells(), 1.0, message,
                                  fields));
}

TEST(ReactionNetwork, AnUnusedFieldNeedNotBeSupplied)
{
    // The counterpart: a network that never mentions light must not demand it.
    Box box;
    box.build();
    box.addSpecies("A", 1.0);
    box.allocate();
    box.setUniform(0, 1.0);

    ReactionNetwork network;
    std::string     message;

    ASSERT_TRUE(network.addRate("A", "0 - 0.001 * A", message)) << message;
    ASSERT_TRUE(network.compile(box.registry, message)) << message;

    EXPECT_FALSE(network.usesLight());
    EXPECT_TRUE(network.integrate(box.registry, box.state, box.cells(), 1.0, message))
        << message;
}

TEST(ReactionNetwork, SpeciesDerivedFieldsAndConstantsCoexistInOneExpression)
{
    // The structural gate. The variable table is species, then the three
    // derived fields, then constants, and the constants moved when the derived
    // block was inserted between them. An off-by-three in that offset would
    // feed a rate law an irradiance where it expected a coefficient — a wrong
    // answer, not a crash, so it needs a test that names all three kinds at
    // once and checks the arithmetic.
    Box box;
    box.build();
    const int a = box.addSpecies("A", 2.0);
    const int b = box.addSpecies("B", 3.0);
    box.allocate();
    box.setUniform(a, 2.0);
    box.setUniform(b, 3.0);

    ReactionNetwork network;
    std::string     message;

    ASSERT_TRUE(network.addConstant("k", 0.5, message)) << message;
    ASSERT_TRUE(network.addConstant("m", 7.0, message)) << message;

    // Every category, and the constants last so a shifted offset shows up.
    ASSERT_TRUE(network.addRate("A", "0 - k * A * B * __LIGHT__ / m / __THICKNESS__",
                                message))
        << message;
    ASSERT_TRUE(network.compile(box.registry, message)) << message;

    constexpr double irradiance = 400.0;
    constexpr double thickness  = 2.5;

    const std::vector<double> light(static_cast<size_t>(box.cells()), irradiance);
    const std::vector<double> thick(static_cast<size_t>(box.cells()), thickness);

    FVQual::Kinetics::DerivedFields fields;
    fields.light     = light;
    fields.thickness = thick;

    constexpr double dt = 0.01;

    ASSERT_TRUE(network.integrate(box.registry, box.state, box.cells(), dt, message,
                                 fields))
        << message;

    // B carries no rate law, so the coefficient multiplying A is constant and
    // the law is first-order with a closed form. The first version of this gate
    // compared against `A0 + rate * dt` instead, with a tolerance a thousand
    // times below that linearisation's own truncation error — it failed by
    // exactly `A0 c² dt²/2`, which was my arithmetic and not a defect.
    const double coefficient = 0.5 * 3.0 * irradiance / 7.0 / thickness;
    const double expected    = 2.0 * std::exp(-coefficient * dt);

    // The tolerance is the integrator's own: it is asked for 1e-8 relative
    // accuracy, so demanding better than that of it would be testing the
    // controller's luck rather than this table's layout. Measured error 1e-9.
    EXPECT_NEAR(box.at(a, 0), expected, 1.0e-8 * expected)
        << "coefficient " << coefficient;

    // The change is large, so a shifted table offset could not slip through as
    // a small discrepancy.
    EXPECT_LT(box.at(a, 0), 0.8 * 2.0);

    EXPECT_EQ(box.at(b, 0), 3.0) << "B has no rate law and must not have moved";
}

TEST(ReactionNetwork, ASpeciesNamedLikeADerivedFieldIsRefused)
{
    Box box;
    box.build();
    box.addSpecies("__LIGHT__", 1.0);
    box.allocate();

    ReactionNetwork network;
    std::string     message;

    ASSERT_TRUE(network.addRate("__LIGHT__", "0 - 0.1 * __LIGHT__", message)) << message;

    EXPECT_FALSE(network.compile(box.registry, message));
    EXPECT_NE(message.find("__LIGHT__"), std::string::npos) << message;
}

TEST(ReactionNetworkSolver, LightLimitedGrowthUsesEachLayersOwnIrradiance)
{
    // What the whole of F6.2a exists for. Growth is light-limited, light falls
    // off with depth, so each layer grows at its own rate — and every rate has a
    // closed form, because with no self-shading the irradiance is constant in
    // time.
    //
    // This also gates that the light reaching the kinetics is per-cell. A single
    // column-averaged scalar, or the surface value used everywhere, would pass a
    // "growth happened" check and fail this one.
    constexpr double incident   = 400.0;
    constexpr double kd         = 0.5;
    constexpr double halfSat    = 150.0; // W/m², the Monod half-saturation
    constexpr double muMax      = 0.002; // 1/s
    constexpr double c0         = 1.0;
    constexpr double totalDepth = 12.0;
    constexpr int    layers     = 6;
    constexpr double dt         = 10.0;
    constexpr int    steps      = 20;

    Basin basin;
    basin.suppressVerticalMixing();
    basin.build(4, layers, 200.0, 40.0, totalDepth, uniform(12.0), uniform(c0));

    basin.solver.lightOptions().incidentShortwave     = incident;
    basin.solver.lightOptions().surfaceAbsorption     = 0.0;
    basin.solver.lightOptions().backgroundAttenuation = kd;

    std::string message;

    ASSERT_TRUE(basin.solver.reactions().addConstant("mu", muMax, message)) << message;
    ASSERT_TRUE(basin.solver.reactions().addConstant("ks", halfSat, message)) << message;
    ASSERT_TRUE(basin.solver.reactions().addRate(
        "DYE", "mu * __LIGHT__ / (ks + __LIGHT__) * DYE", message))
        << message;
    ASSERT_TRUE(basin.solver.compileReactions(message)) << message;

    for (int i = 0; i < steps; ++i)
        ASSERT_TRUE(basin.solver.step(dt, message)) << message;

    const double h       = totalDepth / static_cast<double>(layers);
    const double tau     = kd * h;
    const double elapsed = dt * steps;

    double previous = 0.0;

    for (int k = 0; k < layers; ++k)
    {
        // The layer-mean irradiance, recomputed here from Beer-Lambert rather
        // than read back from the light field, so the two are independent.
        const double top        = incident * std::exp(-kd * static_cast<double>(k) * h);
        const double irradiance = top * (1.0 - std::exp(-tau)) / tau;

        const double rate  = muMax * irradiance / (halfSat + irradiance);
        const double exact = c0 * std::exp(rate * elapsed);

        const double actual = basin.state.speciesAt(basin.dye, basin.layered.cell(0, k));

        EXPECT_NEAR(actual, exact, 1.0e-6 * exact) << "layer " << k;

        // Strictly decreasing with depth. Without this the gate would pass on a
        // uniform field that happened to match the mean.
        if (k > 0)
            EXPECT_LT(actual, previous) << "layer " << k << " grew as fast as the one above";

        previous = actual;
    }

    // And the spread is large, so "each layer its own rate" is a real claim
    // here rather than a distinction within the tolerance.
    const double top    = basin.state.speciesAt(basin.dye, basin.layered.cell(0, 0));
    const double bottom = basin.state.speciesAt(basin.dye, basin.layered.cell(0, layers - 1));

    EXPECT_GT(top / bottom, 1.2) << "top " << top << " bottom " << bottom;
}

TEST(ReactionNetworkSolver, DepthAndThicknessComeFromTheMesh)
{
    // Both are read straight from the geometry, so the gate is that the values a
    // rate law sees are the ones the mesh holds — depth positive downward from
    // the moving surface, thickness the cell's own.
    constexpr double totalDepth = 20.0;
    constexpr int    layers     = 5;
    constexpr double dt         = 1.0e-3;

    Basin basin;
    basin.suppressVerticalMixing();
    basin.build(4, layers, 200.0, 40.0, totalDepth, uniform(12.0), uniform(0.0));

    std::string message;

    // The ratio, not the product: multiplication is commutative, so a rate law
    // written as `__DEPTH__ * __THICKNESS__` would pass this gate with the two
    // fields wired to each other's slots. Dividing makes the swap visible.
    //
    // depth/thickness at the centre of layer k is exactly k + 1/2, which is also
    // a different value in every layer and zero in none.
    ASSERT_TRUE(basin.solver.reactions().addRate("DYE", "__DEPTH__ / __THICKNESS__",
                                                 message))
        << message;
    ASSERT_TRUE(basin.solver.compileReactions(message)) << message;
    ASSERT_TRUE(basin.solver.step(dt, message)) << message;

    for (int k = 0; k < layers; ++k)
    {
        const double expected = (static_cast<double>(k) + 0.5) * dt;

        EXPECT_NEAR(basin.state.speciesAt(basin.dye, basin.layered.cell(0, k)), expected,
                    1.0e-9 * expected)
            << "layer " << k;
    }
}

TEST(ReactionNetworkSolver, TemperatureIsAddressableByItsRegistryId)
{
    // Backs the header's claim that temperature is not a special case. It is
    // registered under a reserved id, and that id has to be usable in a rate
    // law or Arrhenius dependence would need engine support after all.
    constexpr double k20   = 0.01;
    constexpr double theta = 1.047;
    constexpr double water = 24.0;
    constexpr double c0    = 2.0;
    constexpr double dt    = 25.0;
    constexpr int    steps = 12;

    Basin basin;
    basin.build(4, 4, 200.0, 40.0, 10.0, uniform(water), uniform(c0));

    std::string message;

    ASSERT_TRUE(basin.solver.reactions().addConstant("k20", k20, message)) << message;
    ASSERT_TRUE(basin.solver.reactions().addConstant("theta", theta, message)) << message;
    ASSERT_TRUE(basin.solver.reactions().addRate(
        "DYE", "0 - k20 * theta ^ (__TEMPERATURE__ - 20) * DYE", message))
        << message;
    ASSERT_TRUE(basin.solver.compileReactions(message)) << message;

    for (int i = 0; i < steps; ++i)
        ASSERT_TRUE(basin.solver.step(dt, message)) << message;

    const double rate  = k20 * std::pow(theta, water - 20.0);
    const double exact = c0 * std::exp(-rate * dt * steps);

    EXPECT_NEAR(basin.state.speciesAt(basin.dye, 0), exact, 1.0e-6 * c0);

    // The rate must be the warm one, not the 20 °C one — otherwise the
    // temperature term could be silently evaluating as zero.
    const double atTwenty = c0 * std::exp(-k20 * dt * steps);

    EXPECT_LT(basin.state.speciesAt(basin.dye, 0), atTwenty)
        << "the temperature dependence had no effect";
}
