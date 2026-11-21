// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_structures.cpp
 * \brief Verification of outlet works (F4.1).
 *
 * Two layers of claim. The rating laws are arithmetic and are checked against
 * their closed forms directly — an engineer comparing a gate against a
 * manufacturer's table should get the same number the model does.
 *
 * The claim that matters more is the *feedback*: a structure's discharge falls
 * as the reservoir it drains draws down, and a model without that loop will
 * cheerfully empty a lake through its own spillway. The drawdown test is the
 * one that would catch that, and it is the closest thing here to a real
 * operating scenario.
 */

#include "fvqual/forcing/structures.h"
#include "fvqual/hydro/internalmode.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

using FVQual::Forcing::Placement;
using FVQual::Forcing::SourceSet;
using FVQual::Forcing::Structure;
using FVQual::Forcing::StructureSet;
using FVQual::Forcing::StructureType;
using FVQual::Hydro::InternalModeOptions;
using FVQual::Hydro::InternalModeSolver;
using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;
using FVQual::State::SpeciesRegistry;
using FVQual::State::StateStore;
namespace Testing = FVQual::Testing;

namespace {

constexpr double kGravity = 9.80665;

Structure weir(double crest, double width, double coefficient = 1.84)
{
    Structure structure;
    structure.id             = "spillway";
    structure.type           = StructureType::Weir;
    structure.crestElevation = crest;
    structure.width          = width;
    structure.coefficient    = coefficient;
    return structure;
}

Structure orifice(double centre, double area, double coefficient = 0.62)
{
    Structure structure;
    structure.id             = "gate";
    structure.type           = StructureType::Orifice;
    structure.crestElevation = centre;
    structure.area           = area;
    structure.coefficient    = coefficient;
    return structure;
}

/*!
 * \brief A flat-bedded reservoir with structures attached.
 */
struct Reservoir
{
    HorizontalMesh     horizontal;
    LayeredMesh        layered;
    SpeciesRegistry    registry;
    StateStore         state;
    InternalModeSolver solver;

    void build(int columns, int layers, double length, double width, double bed,
               double surface)
    {
        std::string message;

        ASSERT_TRUE(horizontal.build(Testing::channelStrip(columns, length, width),
                                     message))
            << message;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = bed;

        SigmaSpec spec;
        spec.layerCount        = layers;
        spec.minLayerThickness = 0.0;

        ASSERT_TRUE(layered.build(horizontal, spec, surface, message)) << message;
        ASSERT_GE(registry.addTemperature(10.0, message), 0) << message;

        state.allocate(layered, registry);

        InternalModeOptions options;
        options.gravity   = kGravity;
        options.manningN  = 0.0;
        options.cflNumber = 0.5;

        ASSERT_TRUE(solver.initialize(layered, state, registry, options, message))
            << message;
    }

    [[nodiscard]] double meanSurface() const
    {
        double total = 0.0;

        for (int64_t c = 0; c < layered.columnCount(); ++c)
            total += solver.external().surfaceElevation()[static_cast<size_t>(c)];

        return total / static_cast<double>(layered.columnCount());
    }
};

} // namespace

// ── Rating laws ───────────────────────────────────────────────────────────

TEST(Structures, WeirFollowsItsRatingCurve)
{
    const Structure spillway = weir(/*crest=*/-2.0, /*width=*/15.0, /*C=*/1.84);

    // Q = C L H^(3/2)
    for (double head : {0.25, 0.5, 1.0, 2.0})
    {
        const double expected = 1.84 * 15.0 * std::pow(head, 1.5);

        EXPECT_NEAR(StructureSet::rate(spillway, -2.0 + head), expected, 1.0e-9)
            << "at head " << head;
    }

    // Below the crest there is nothing to drive it.
    EXPECT_EQ(StructureSet::rate(spillway, -2.0), 0.0);
    EXPECT_EQ(StructureSet::rate(spillway, -3.0), 0.0);
}

TEST(Structures, SubmergenceThrottlesAWeir)
{
    Structure spillway = weir(-2.0, 15.0);

    const double upstream = 0.0;                       // 2 m of head
    const double free     = StructureSet::rate(spillway, upstream);

    ASSERT_GT(free, 0.0);

    // Tailwater below the crest changes nothing.
    spillway.tailwaterElevation = -3.0;
    EXPECT_NEAR(StructureSet::rate(spillway, upstream), free, 1.0e-12);

    // Rising tailwater drowns it progressively.
    spillway.tailwaterElevation = -1.0;
    const double partly = StructureSet::rate(spillway, upstream);

    spillway.tailwaterElevation = -0.5;
    const double mostly = StructureSet::rate(spillway, upstream);

    EXPECT_LT(partly, free);
    EXPECT_LT(mostly, partly);

    // Level tailwater means no head difference and no flow.
    spillway.tailwaterElevation = upstream;
    EXPECT_EQ(StructureSet::rate(spillway, upstream), 0.0);
}

TEST(Structures, OrificeFollowsItsRatingCurveAndItsGate)
{
    Structure gate = orifice(/*centre=*/-8.0, /*area=*/2.5, /*C=*/0.62);

    const double head     = 8.0; // surface at zero, center at -8
    const double expected = 0.62 * 2.5 * std::sqrt(2.0 * kGravity * head);

    EXPECT_NEAR(StructureSet::rate(gate, 0.0), expected, 1.0e-9);

    // The gate opening is a linear multiplier — the operator's actual lever.
    gate.opening = 0.5;
    EXPECT_NEAR(StructureSet::rate(gate, 0.0), 0.5 * expected, 1.0e-9);

    gate.opening = 0.0;
    EXPECT_EQ(StructureSet::rate(gate, 0.0), 0.0);

    // Submerging the orifice switches the driving head to the difference
    // across it, which is less.
    gate.opening            = 1.0;
    gate.tailwaterElevation = -4.0;

    const double submerged = 0.62 * 2.5 * std::sqrt(2.0 * kGravity * 4.0);
    EXPECT_NEAR(StructureSet::rate(gate, 0.0), submerged, 1.0e-9);
    EXPECT_LT(submerged, expected);
}

TEST(Structures, OutletWorksDoNotRunBackwards)
{
    // A tailwater above the reservoir would physically drive flow inwards.
    // The model closes instead, because admitting water means deciding what it
    // is made of, and inventing that silently is worse than stopping.
    Structure gate = orifice(-8.0, 2.5);
    gate.tailwaterElevation = 5.0;

    EXPECT_EQ(StructureSet::rate(gate, 0.0), 0.0);

    Structure spillway = weir(-2.0, 15.0);
    spillway.tailwaterElevation = 5.0;

    EXPECT_EQ(StructureSet::rate(spillway, 0.0), 0.0);
}

TEST(Structures, MalformedStructuresAreRefused)
{
    SourceSet   sources;
    StructureSet structures;
    std::string message;

    Structure noCrest;
    noCrest.type = StructureType::Weir;
    EXPECT_LT(structures.add(noCrest, sources, message), 0);
    EXPECT_FALSE(message.empty());

    Structure noArea;
    noArea.type = StructureType::Orifice;
    EXPECT_LT(structures.add(noArea, sources, message), 0);

    Structure idlePump;
    idlePump.type = StructureType::Pump;
    EXPECT_LT(structures.add(idlePump, sources, message), 0);

    // An inverted deadband would make a controlled structure chatter, so it is
    // rejected at configuration rather than discovered at runtime.
    Structure inverted     = weir(-2.0, 10.0);
    inverted.controlled    = true;
    inverted.startElevation = 0.0;
    inverted.stopElevation  = 1.0;
    EXPECT_LT(structures.add(inverted, sources, message), 0);
}

// ── Control ───────────────────────────────────────────────────────────────

TEST(Structures, PumpHysteresisDoesNotChatter)
{
    SourceSet    sources;
    StructureSet structures;
    std::string  message;

    Structure pump;
    pump.id             = "release pump";
    pump.type           = StructureType::Pump;
    pump.capacity       = -12.0; // removes water
    pump.controlled     = true;
    pump.startElevation = 1.0;
    pump.stopElevation  = 0.0;
    pump.placement      = Placement::Uniform;

    const int index = structures.add(pump, sources, message);
    ASSERT_GE(index, 0) << message;

    const auto levelIs = [&](double level) {
        const std::vector<double> surface(1, level);
        EXPECT_TRUE(structures.update(surface, sources, message)) << message;
        return structures.isRunning(index);
    };

    EXPECT_FALSE(levelIs(0.5)) << "should start idle inside the deadband";
    EXPECT_TRUE(levelIs(1.2)) << "should start above the start level";

    // The defining property: inside the deadband it holds whatever state it
    // was in. A single-threshold controller would switch here and cycle every
    // step for the rest of the run.
    EXPECT_TRUE(levelIs(0.5)) << "should keep running inside the deadband";
    EXPECT_TRUE(levelIs(0.1));

    EXPECT_FALSE(levelIs(-0.1)) << "should stop below the stop level";
    EXPECT_FALSE(levelIs(0.5)) << "should stay off inside the deadband";

    // And the discharge follows the state.
    (void)levelIs(2.0);
    EXPECT_NEAR(structures.discharge(index), -12.0, 1.0e-12);

    (void)levelIs(-1.0);
    EXPECT_EQ(structures.discharge(index), 0.0);
}

TEST(Structures, StructureDrivesItsSource)
{
    SourceSet    sources;
    StructureSet structures;
    std::string  message;

    const int index = structures.add(weir(-2.0, 10.0), sources, message);
    ASSERT_GE(index, 0) << message;
    ASSERT_EQ(sources.count(), 1);

    const std::vector<double> surface(1, 0.0);
    ASSERT_TRUE(structures.update(surface, sources, message)) << message;

    const double expected = -1.84 * 10.0 * std::pow(2.0, 1.5);

    EXPECT_NEAR(structures.discharge(index), expected, 1.0e-9);
    EXPECT_NEAR(sources.source(0).discharge, expected, 1.0e-9)
        << "the source the structure drives must carry its discharge";

    // Releases are negative, as a withdrawal is.
    EXPECT_LT(sources.netDischarge(), 0.0);
}

// ── The feedback that matters ─────────────────────────────────────────────

TEST(Structures, ReservoirDrawsDownToItsWeirCrestAndStops)
{
    // The integration test, and the one a prescribed discharge would fail. As
    // the reservoir falls, the head on the crest falls with it, so the release
    // decays and the level asymptotes to the crest instead of continuing to
    // the bed.
    constexpr double kCrest = -3.0;

    Reservoir reservoir;
    reservoir.build(12, 8, 600.0, 60.0, /*bed=*/-15.0, /*surface=*/0.0);

    std::string message;

    Structure spillway  = weir(kCrest, /*width=*/40.0);
    spillway.column     = 11;
    spillway.placement  = Placement::SelectiveWithdrawal;

    const int index = reservoir.solver.structures().add(
        spillway, reservoir.solver.sources(), message);
    ASSERT_GE(index, 0) << message;

    const double start = reservoir.meanSurface();
    ASSERT_NEAR(start, 0.0, 1.0e-9);

    // Sampled early: with 3 m of head over a 40 m crest this pool empties to
    // the crest in roughly five minutes, so a coarse first sample would miss
    // the decay entirely and see only the finished state.
    ASSERT_TRUE(reservoir.solver.advance(30.0, message)) << message;

    const double afterFirst   = reservoir.meanSurface();
    const double firstRelease = std::abs(reservoir.solver.structures().discharge(index));

    EXPECT_LT(afterFirst, start) << "the reservoir should be drawing down";
    EXPECT_GT(firstRelease, 0.0);

    ASSERT_TRUE(reservoir.solver.advance(90.0, message)) << message;

    const double afterSecond   = reservoir.meanSurface();
    const double secondRelease = std::abs(reservoir.solver.structures().discharge(index));

    // The feedback: less head, less flow.
    EXPECT_LT(afterSecond, afterFirst);
    EXPECT_LT(secondRelease, firstRelease)
        << "release " << secondRelease << " should have fallen from " << firstRelease;

    // The level must not pass the crest, however long it runs — the physical
    // statement a prescribed release cannot make.
    ASSERT_TRUE(reservoir.solver.advance(20000.0, message)) << message;

    EXPECT_GT(reservoir.meanSurface(), kCrest - 1.0e-6)
        << "drew down past the crest to " << reservoir.meanSurface();

    EXPECT_LT(reservoir.meanSurface(), afterSecond + 1.0e-9);
    EXPECT_LT(std::abs(reservoir.solver.structures().discharge(index)), firstRelease);

    // Mode-splitting consistency has to survive a structure-driven release.
    EXPECT_LT(reservoir.solver.diagnostics().surfaceFluxResidual, 1.0e-10);
}

TEST(Structures, DrawdownFollowsTheClosedFormWeirLaw)
{
    // The quantitative version of the test above. A weir draining a reservoir
    // of constant plan area has an exact solution: A dh/dt = -C L h^(3/2)
    // integrates to
    //
    //     h(t) = h0 / (1 + C L sqrt(h0) t / (2 A))^2
    //
    // so the level approaches the crest as 1/t^2 and never reaches it. Testing
    // against this rather than against a tolerance is what distinguishes a
    // model that drains at the right *rate* from one that merely ends up in the
    // right place, and it is the check that fails loudly if the release is ever
    // rated against a stale head again.
    constexpr double kCrest  = -3.0;
    constexpr double kLength = 600.0;
    constexpr double kWidth  = 60.0;
    constexpr double kCrestLength = 40.0;
    constexpr double kC      = 1.84;
    constexpr double kH0     = 3.0;

    Reservoir reservoir;
    reservoir.build(12, 8, kLength, kWidth, /*bed=*/-15.0, /*surface=*/0.0);

    std::string message;

    Structure spillway = weir(kCrest, kCrestLength, kC);
    spillway.column    = 11;
    spillway.placement = Placement::SelectiveWithdrawal;

    ASSERT_GE(reservoir.solver.structures().add(spillway, reservoir.solver.sources(),
                                                message),
              0)
        << message;

    const double area = kLength * kWidth;
    const double k    = kC * kCrestLength * std::sqrt(kH0) / (2.0 * area);

    // The closed form assumes the pool is level at every instant. A real basin
    // is not: a gravity wave needs L/sqrt(g H) — about 50 s here — to carry the
    // drawdown from the spillway to the far end, and until it has, the column
    // at the weir sits below the level-pool height and releases less. The
    // comparison therefore starts after two basin transits, and the early
    // deficit is asserted separately below as the physical effect it is.
    double time = 0.0;

    for (double target : {100.0, 300.0, 1000.0, 3000.0, 10000.0})
    {
        ASSERT_TRUE(reservoir.solver.advance(target - time, message)) << message;
        time = target;

        const double modelled =
            reservoir.solver.external().surfaceElevation()[11] - kCrest;
        const double exact = kH0 / std::pow(1.0 + k * time, 2.0);

        EXPECT_NEAR(modelled, exact, 0.10 * exact)
            << "at t = " << time << " s the head was " << modelled
            << " m against the closed-form " << exact << " m";

        // The rating curve must still hold at every instant, on the head the
        // structure actually sees — the discharge and the level cannot drift
        // apart the way they do when a stale head is used.
        EXPECT_NEAR(std::abs(reservoir.solver.structures().discharge(0)),
                    kC * kCrestLength * std::pow(modelled, 1.5),
                    1.0e-6 + 0.02 * kC * kCrestLength * std::pow(modelled, 1.5))
            << "discharge and head disagree at t = " << time << " s";
    }
}

TEST(Structures, TheBasinTakesTimeToLevelBeforeTheLevelPoolLawApplies)
{
    // The counterpart to the exclusion above: before a gravity wave has crossed
    // the basin, the model should *not* match the level-pool law, and should
    // release less than it, because the water at the far end has not yet
    // arrived. A model that matched the closed form from t = 0 would be
    // reporting a basin that responds instantaneously over 600 m.
    constexpr double kCrest = -3.0;

    Reservoir reservoir;
    reservoir.build(12, 8, 600.0, 60.0, -15.0, 0.0);

    std::string message;

    Structure spillway = weir(kCrest, 40.0);
    spillway.column    = 11;
    spillway.placement = Placement::SelectiveWithdrawal;

    ASSERT_GE(reservoir.solver.structures().add(spillway, reservoir.solver.sources(),
                                                message),
              0)
        << message;

    // One quarter of a basin transit.
    ASSERT_TRUE(reservoir.solver.advance(10.0, message)) << message;

    const double k = 1.84 * 40.0 * std::sqrt(3.0) / (2.0 * 600.0 * 60.0);
    const double levelPool = 3.0 / std::pow(1.0 + k * 10.0, 2.0);
    const double atWeir = reservoir.solver.external().surfaceElevation()[11] - kCrest;

    EXPECT_LT(atWeir, levelPool)
        << "the column at the spillway should have drawn down further than a "
           "level pool, not less";

    // The far end should barely have moved yet.
    EXPECT_GT(reservoir.solver.external().surfaceElevation()[0], kCrest + 2.9)
        << "the far end of the basin responded faster than a gravity wave";
}

TEST(Structures, GateClosureStopsTheRelease)
{
    // The operational lever, end to end: closing the gate stops the water.
    Reservoir reservoir;
    reservoir.build(10, 6, 400.0, 50.0, -20.0, 0.0);

    std::string message;

    Structure gate = orifice(/*centre=*/-15.0, /*area=*/3.0);
    gate.column    = 9;

    const int index =
        reservoir.solver.structures().add(gate, reservoir.solver.sources(), message);
    ASSERT_GE(index, 0) << message;

    ASSERT_TRUE(reservoir.solver.advance(300.0, message)) << message;

    const double openLevel = reservoir.meanSurface();
    EXPECT_LT(openLevel, 0.0) << "an open gate should have released water";
    EXPECT_LT(reservoir.solver.structures().discharge(index), 0.0);

    ASSERT_TRUE(reservoir.solver.structures().setOpening(index, 0.0, message)) << message;
    ASSERT_TRUE(reservoir.solver.advance(300.0, message)) << message;

    EXPECT_EQ(reservoir.solver.structures().discharge(index), 0.0);
    EXPECT_NEAR(reservoir.meanSurface(), openLevel, 1.0e-6)
        << "the level should hold once the gate is shut";
}
