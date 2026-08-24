// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_outletgroups.cpp
 * \brief Verification of multi-port outlet operation (F4.2).
 *
 * The claim under test is an operational one: given ports at different
 * elevations in a stratified reservoir, the model can hit a target release
 * temperature that no single port reaches, and the split it chooses is the one
 * heat conservation dictates rather than a plausible-looking guess.
 *
 * The gates that matter most are the ones about *failure* to hit the target.
 * A blending model that always reports success is useless to an operator: the
 * questions worth answering are what happens when the target is outside
 * anything the reservoir holds, and what happens when the ports cannot pass the
 * flow demanded. Both must be visible rather than papered over.
 */

#include "fvqual/forcing/outletgroups.h"
#include "fvqual/hydro/internalmode.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

using FVQual::Forcing::OutletGroup;
using FVQual::Forcing::OutletGroupSet;
using FVQual::Forcing::Placement;
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

constexpr double kGravity   = 9.80665;
constexpr int    kColumns   = 10;
constexpr int    kLayers    = 12;
constexpr double kLength    = 500.0;
constexpr double kWidth     = 80.0;
constexpr double kBed       = -24.0;
constexpr int    kDamColumn = 9;

/*!
 * \brief A thermally stratified reservoir with outlet works at the dam.
 *
 * Temperature falls linearly with depth from \p surfaceTemperature at the free
 * surface to \p bedTemperature at the bed, which is the summer profile the
 * whole exercise is about. A uniform profile is used by one test to show the
 * operator's lever disappearing.
 */
struct StratifiedReservoir
{
    HorizontalMesh     horizontal;
    LayeredMesh        layered;
    SpeciesRegistry    registry;
    StateStore         state;
    InternalModeSolver solver;
    int                temperature = -1;

    void build(double surfaceTemperature, double bedTemperature)
    {
        std::string message;

        ASSERT_TRUE(horizontal.build(
            Testing::channelStrip(kColumns, kLength, kWidth), message))
            << message;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = kBed;

        SigmaSpec spec;
        spec.layerCount        = kLayers;
        spec.minLayerThickness = 0.0;

        ASSERT_TRUE(layered.build(horizontal, spec, /*surface=*/0.0, message)) << message;

        temperature = registry.addTemperature(surfaceTemperature, message);
        ASSERT_GE(temperature, 0) << message;

        state.allocate(layered, registry);

        for (int64_t c = 0; c < layered.columnCount(); ++c)
        {
            for (int k = 0; k < kLayers; ++k)
            {
                const int64_t cell = layered.cell(c, k);

                const double centre =
                    0.5 * (layered.z(c, k) + layered.z(c, k + 1));
                const double fraction = centre / kBed; // 0 at surface, 1 at bed

                state.speciesAt(temperature, cell) =
                    surfaceTemperature +
                    fraction * (bedTemperature - surfaceTemperature);
            }
        }

        InternalModeOptions options;
        options.gravity   = kGravity;
        options.manningN  = 0.0;
        options.cflNumber = 0.5;

        ASSERT_TRUE(solver.initialize(layered, state, registry, options, message))
            << message;
    }

    /*!
     * \brief Adds a scheduled port at \p elevation, wide enough not to be the
     *        binding constraint unless a test means it to be.
     */
    int addPort(const std::string &id, double elevation, double area = 30.0)
    {
        Structure port;
        port.id             = id;
        port.type           = StructureType::Orifice;
        port.column         = kDamColumn;
        port.crestElevation = elevation;
        port.area           = area;
        port.placement      = Placement::SelectiveWithdrawal;
        port.scheduled      = true;

        std::string message;
        const int   index =
            solver.structures().add(port, solver.sources(), message);

        EXPECT_GE(index, 0) << message;
        return index;
    }

    /*!
     * \brief Registers a blending group over \p ports.
     */
    int addGroup(const std::vector<int> &ports, double demand, double target)
    {
        OutletGroup group;
        group.id                = "outlet works";
        group.ports             = ports;
        group.demand            = demand;
        group.targetTemperature = target;

        std::string message;
        const int   index =
            solver.outletGroups().add(group, solver.structures(), message);

        EXPECT_GE(index, 0) << message;
        return index;
    }

    /*!
     * \brief Runs the group decision without advancing time.
     */
    void operate()
    {
        std::string message;

        ASSERT_TRUE(solver.outletGroups().update(
            layered, solver.sources(), state.density,
            state.speciesBlock(temperature), kGravity, 998.2, solver.structures(),
            message))
            << message;
    }
};

} // namespace

// ── Configuration ─────────────────────────────────────────────────────────

TEST(OutletGroups, MalformedGroupsAreRefused)
{
    StructureSet   structures;
    OutletGroupSet groups;
    FVQual::Forcing::SourceSet sources;
    std::string    message;

    Structure port;
    port.type           = StructureType::Orifice;
    port.crestElevation = -10.0;
    port.area           = 4.0;

    const int a = structures.add(port, sources, message);
    ASSERT_GE(a, 0) << message;

    OutletGroup empty;
    empty.ports = {};
    EXPECT_LT(groups.add(empty, structures, message), 0);
    EXPECT_FALSE(message.empty());

    OutletGroup missing;
    missing.ports = {a, 7};
    EXPECT_LT(groups.add(missing, structures, message), 0);

    // The one that would silently release twice what was asked.
    OutletGroup duplicated;
    duplicated.ports = {a, a};
    EXPECT_LT(groups.add(duplicated, structures, message), 0);

    OutletGroup negative;
    negative.ports  = {a};
    negative.demand = -5.0;
    EXPECT_LT(groups.add(negative, structures, message), 0);

    OutletGroup fine;
    fine.ports = {a};
    EXPECT_GE(groups.add(fine, structures, message), 0) << message;
}

TEST(OutletGroups, AScheduledStructureReleasesItsScheduleNotItsRating)
{
    // The change of role the whole feature rests on: under a schedule the
    // rating curve stops being the discharge and becomes the capacity.
    StructureSet structures;
    FVQual::Forcing::SourceSet sources;
    std::string  message;

    Structure gate;
    gate.type           = StructureType::Orifice;
    gate.crestElevation = -20.0;
    gate.area           = 8.0;

    const int index = structures.add(gate, sources, message);
    ASSERT_GE(index, 0) << message;

    const std::vector<double> surface(1, 0.0);
    const double              capacity = StructureSet::rate(gate, 0.0);

    ASSERT_GT(capacity, 50.0) << "fixture needs headroom above the scheduled 40";

    // Free rating: it passes what the head drives.
    ASSERT_TRUE(structures.update(surface, sources, message)) << message;
    EXPECT_NEAR(std::abs(structures.discharge(index)), capacity, 1.0e-9);

    // Scheduled below capacity: it passes the schedule.
    ASSERT_TRUE(structures.setSchedule(index, 40.0, message)) << message;
    ASSERT_TRUE(structures.update(surface, sources, message)) << message;
    EXPECT_NEAR(std::abs(structures.discharge(index)), 40.0, 1.0e-9);

    // Scheduled above capacity: it passes the capacity, not the wish.
    ASSERT_TRUE(structures.setSchedule(index, 10.0 * capacity, message)) << message;
    ASSERT_TRUE(structures.update(surface, sources, message)) << message;
    EXPECT_NEAR(std::abs(structures.discharge(index)), capacity, 1.0e-9);

    // And it goes back.
    ASSERT_TRUE(structures.clearSchedule(index, message)) << message;
    ASSERT_TRUE(structures.update(surface, sources, message)) << message;
    EXPECT_NEAR(std::abs(structures.discharge(index)), capacity, 1.0e-9);

    // A negative schedule is a confused magnitude, and is refused.
    EXPECT_FALSE(structures.setSchedule(index, -5.0, message));
    EXPECT_FALSE(message.empty());
}

// ── The blend ─────────────────────────────────────────────────────────────

TEST(OutletGroups, PortsAtDifferentDepthsDeliverDifferentTemperatures)
{
    // The precondition for everything else: if the ports all delivered the same
    // water there would be nothing to blend. This is also the gate that fails
    // if the withdrawal envelope ever stops depending on outlet elevation.
    StratifiedReservoir reservoir;
    reservoir.build(/*surface=*/24.0, /*bed=*/6.0);

    const int deep    = reservoir.addPort("low-level", -21.0);
    const int shallow = reservoir.addPort("upper", -4.0);

    const int group = reservoir.addGroup({deep, shallow}, /*demand=*/60.0,
                                         /*target=*/15.0);
    ASSERT_GE(group, 0);

    reservoir.operate();

    const double deepTemperature = reservoir.solver.sources().withdrawalMean(
        reservoir.layered, reservoir.solver.structures().sourceIndex(deep), 60.0,
        reservoir.state.density, kGravity, 998.2,
        reservoir.state.speciesBlock(reservoir.temperature));

    const double shallowTemperature = reservoir.solver.sources().withdrawalMean(
        reservoir.layered, reservoir.solver.structures().sourceIndex(shallow), 60.0,
        reservoir.state.density, kGravity, 998.2,
        reservoir.state.speciesBlock(reservoir.temperature));

    EXPECT_LT(deepTemperature, 12.0)
        << "a port 21 m down in a 24-to-6 degree profile drew " << deepTemperature;
    EXPECT_GT(shallowTemperature, 18.0)
        << "a port 4 m down drew " << shallowTemperature;
    EXPECT_GT(shallowTemperature - deepTemperature, 6.0)
        << "the ports must offer a usable temperature range to blend across";
}

TEST(OutletGroups, TheBlendFollowsTheHeatConservationMixingLaw)
{
    // Q_w T_w + Q_c T_c = (Q_w + Q_c) T*, so the warm fraction is
    // (T* - T_c) / (T_w - T_c). Checked against the port temperatures the model
    // itself reports, so this tests the *split*, not the envelope.
    StratifiedReservoir reservoir;
    reservoir.build(24.0, 6.0);

    const int deep    = reservoir.addPort("low-level", -21.0);
    const int shallow = reservoir.addPort("upper", -4.0);

    constexpr double kDemand = 60.0;

    const int group = reservoir.addGroup({deep, shallow}, kDemand, /*target=*/15.0);
    ASSERT_GE(group, 0);

    reservoir.operate();

    const auto portTemperature = [&](int port) {
        return reservoir.solver.sources().withdrawalMean(
            reservoir.layered, reservoir.solver.structures().sourceIndex(port),
            kDemand, reservoir.state.density, kGravity, 998.2,
            reservoir.state.speciesBlock(reservoir.temperature));
    };

    const double cold = portTemperature(deep);
    const double warm = portTemperature(shallow);

    ASSERT_LT(cold, 15.0);
    ASSERT_GT(warm, 15.0) << "the target must be bracketed for this test to mean anything";

    const double expectedWarmFraction = (15.0 - cold) / (warm - cold);

    const double warmRelease =
        std::abs(reservoir.solver.structures().structure(shallow).scheduledDischarge);
    const double coldRelease =
        std::abs(reservoir.solver.structures().structure(deep).scheduledDischarge);

    EXPECT_NEAR(warmRelease, expectedWarmFraction * kDemand, 1.0e-6);
    EXPECT_NEAR(coldRelease, (1.0 - expectedWarmFraction) * kDemand, 1.0e-6);

    // Total release is what was asked for, and the blend hits the target.
    EXPECT_NEAR(warmRelease + coldRelease, kDemand, 1.0e-9);
    EXPECT_NEAR(reservoir.solver.outletGroups().deliveredRelease(group), kDemand, 1.0e-9);
    EXPECT_NEAR(reservoir.solver.outletGroups().achievedTemperature(group), 15.0, 1.0e-6);
}

TEST(OutletGroups, TheAchievedTemperatureTracksTheTargetAcrossTheAchievableRange)
{
    // Sweeping the target is what distinguishes a working controller from one
    // that happens to be right at a single set point.
    StratifiedReservoir reservoir;
    reservoir.build(24.0, 6.0);

    const int deep    = reservoir.addPort("low-level", -21.0);
    const int shallow = reservoir.addPort("upper", -4.0);

    const int group = reservoir.addGroup({deep, shallow}, 60.0, 15.0);
    ASSERT_GE(group, 0);

    std::string message;
    double      previous = -1.0e30;

    for (double target : {12.0, 14.0, 16.0, 18.0})
    {
        ASSERT_TRUE(reservoir.solver.outletGroups().setTarget(group, target, message))
            << message;

        reservoir.operate();

        const double achieved =
            reservoir.solver.outletGroups().achievedTemperature(group);

        EXPECT_NEAR(achieved, target, 1.0e-6)
            << "asked for " << target << " and scheduled a release of " << achieved;

        EXPECT_GT(achieved, previous) << "the response must be monotone in the target";
        previous = achieved;

        EXPECT_NEAR(reservoir.solver.outletGroups().deliveredRelease(group), 60.0,
                    1.0e-9);
    }
}

// ── Failing to hit the target, visibly ────────────────────────────────────

TEST(OutletGroups, ATargetColderThanTheReservoirGivesTheColdestPortAndMissesHonestly)
{
    // No combination of outlets can produce water the reservoir does not
    // contain. The model must open the coldest port and report an achieved
    // temperature above the target, not claim success.
    StratifiedReservoir reservoir;
    reservoir.build(24.0, 6.0);

    const int deep    = reservoir.addPort("low-level", -21.0);
    const int shallow = reservoir.addPort("upper", -4.0);

    const int group = reservoir.addGroup({deep, shallow}, 60.0, /*target=*/2.0);
    ASSERT_GE(group, 0);

    reservoir.operate();

    EXPECT_NEAR(
        std::abs(reservoir.solver.structures().structure(deep).scheduledDischarge),
        60.0, 1.0e-9)
        << "the coldest port should carry the whole release";
    EXPECT_NEAR(
        std::abs(reservoir.solver.structures().structure(shallow).scheduledDischarge),
        0.0, 1.0e-9);

    EXPECT_GT(reservoir.solver.outletGroups().achievedTemperature(group), 2.0)
        << "the model claimed to deliver water colder than anything in the lake";
}

TEST(OutletGroups, ATargetWarmerThanTheReservoirGivesTheWarmestPort)
{
    StratifiedReservoir reservoir;
    reservoir.build(24.0, 6.0);

    const int deep    = reservoir.addPort("low-level", -21.0);
    const int shallow = reservoir.addPort("upper", -4.0);

    const int group = reservoir.addGroup({deep, shallow}, 60.0, /*target=*/40.0);
    ASSERT_GE(group, 0);

    reservoir.operate();

    EXPECT_NEAR(
        std::abs(reservoir.solver.structures().structure(shallow).scheduledDischarge),
        60.0, 1.0e-9);
    EXPECT_LT(reservoir.solver.outletGroups().achievedTemperature(group), 40.0);
}

TEST(OutletGroups, ADemandBeyondCapacityIsReportedAsAShortfall)
{
    // The failure mode that would otherwise be diagnosed as a mass-balance bug
    // three modules away: a release that silently comes up short.
    StratifiedReservoir reservoir;
    reservoir.build(24.0, 6.0);

    // Deliberately tiny ports.
    const int deep    = reservoir.addPort("low-level", -21.0, /*area=*/0.01);
    const int shallow = reservoir.addPort("upper", -4.0, /*area=*/0.01);

    constexpr double kDemand = 5000.0;

    const int group = reservoir.addGroup({deep, shallow}, kDemand, 15.0);
    ASSERT_GE(group, 0);

    reservoir.operate();

    const double delivered = reservoir.solver.outletGroups().deliveredRelease(group);

    EXPECT_GT(delivered, 0.0);
    EXPECT_LT(delivered, kDemand)
        << "the ports cannot pass " << kDemand << " m³/s and the model should say so";

    // And what it does deliver is exactly the two capacities, not a guess.
    const double capacity =
        StructureSet::rate(reservoir.solver.structures().structure(deep), 0.0) +
        StructureSet::rate(reservoir.solver.structures().structure(shallow), 0.0);

    EXPECT_NEAR(delivered, capacity, 1.0e-9);
}

TEST(OutletGroups, StratificationIsWhatMakesBlendingPossibleAtAll)
{
    // The physical counterpart to the envelope theory: in a well-mixed column
    // every port draws the same water, so the ports offer no range and the
    // operator's lever is gone. This is why a dam that can manage release
    // temperature in August cannot in February, and a model that reported a
    // usable spread here would be wrong about the season.
    StratifiedReservoir reservoir;
    reservoir.build(/*surface=*/11.0, /*bed=*/11.0);

    const int deep    = reservoir.addPort("low-level", -21.0);
    const int shallow = reservoir.addPort("upper", -4.0);

    const int group = reservoir.addGroup({deep, shallow}, 60.0, /*target=*/18.0);
    ASSERT_GE(group, 0);

    reservoir.operate();

    const auto portTemperature = [&](int port) {
        return reservoir.solver.sources().withdrawalMean(
            reservoir.layered, reservoir.solver.structures().sourceIndex(port), 60.0,
            reservoir.state.density, kGravity, 998.2,
            reservoir.state.speciesBlock(reservoir.temperature));
    };

    EXPECT_NEAR(portTemperature(deep), portTemperature(shallow), 1.0e-6)
        << "an unstratified column cannot offer the ports different water";

    // The target is unreachable, and the achieved release says so plainly.
    EXPECT_NEAR(reservoir.solver.outletGroups().achievedTemperature(group), 11.0, 1.0e-6);
}

// ── End to end ────────────────────────────────────────────────────────────

TEST(OutletGroups, AGroupOperatesThroughAStepAndStaysMassConsistent)
{
    // The integration check: the group schedules, the ports rate against the
    // schedule, the subcycle enforces capacity, and the two modes still agree
    // about how much water left.
    StratifiedReservoir reservoir;
    reservoir.build(24.0, 6.0);

    const int deep    = reservoir.addPort("low-level", -21.0);
    const int shallow = reservoir.addPort("upper", -4.0);

    constexpr double kDemand = 40.0;

    const int group = reservoir.addGroup({deep, shallow}, kDemand, 15.0);
    ASSERT_GE(group, 0);

    const double before = reservoir.layered.totalVolume();

    std::string message;
    ASSERT_TRUE(reservoir.solver.advance(200.0, message)) << message;

    const double after = reservoir.layered.totalVolume();

    // Both ports are working.
    EXPECT_GT(std::abs(reservoir.solver.structures().discharge(deep)), 0.0);
    EXPECT_GT(std::abs(reservoir.solver.structures().discharge(shallow)), 0.0);

    // The volume that left is the release that was scheduled, to the accuracy
    // the free surface can be read back at.
    const double released = before - after;

    EXPECT_NEAR(released, kDemand * 200.0, 0.02 * kDemand * 200.0)
        << "scheduled " << kDemand << " m³/s for 200 s but the reservoir lost "
        << released << " m³";

    EXPECT_LT(reservoir.solver.diagnostics().surfaceFluxResidual, 1.0e-10);
}
