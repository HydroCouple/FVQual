// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_transfers.cpp
 * \brief Verification of internal transfers between columns (F4.3).
 *
 * A transfer is two sources that must stay welded together. The gates split
 * along that seam:
 *
 * - **Bookkeeping.** Volume moves from one column to another and none of it
 *   goes missing; constituent mass is conserved *exactly*, because the arriving
 *   composition is the mean through the same weights the withdrawal removes
 *   mass at. This is an identity, not an approximation, so it is asserted at
 *   round-off rather than against a tolerance.
 *
 * - **Physics.** The water that arrives is the water that left. A transfer from
 *   a deep intake delivers cold water and one from a shallow intake delivers
 *   warm water, in the same reservoir at the same rate — and that difference is
 *   the whole reason to model an intake elevation at all.
 */

#include "fvqual/forcing/transfers.h"
#include "fvqual/hydro/internalmode.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

using FVQual::Forcing::Placement;
using FVQual::Forcing::SourceSet;
using FVQual::Forcing::Transfer;
using FVQual::Forcing::TransferSet;
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
constexpr int    kColumns = 8;
constexpr int    kLayers  = 10;
constexpr double kLength  = 400.0;
constexpr double kWidth   = 60.0;
constexpr double kBed     = -20.0;

/*!
 * \brief A stratified basin used as both donor and recipient.
 */
struct Basin
{
    HorizontalMesh     horizontal;
    LayeredMesh        layered;
    SpeciesRegistry    registry;
    StateStore         state;
    InternalModeSolver solver;
    int                temperature = -1;
    int                tracer      = -1;

    void build(double surfaceTemperature, double bedTemperature,
               bool withTracer = false)
    {
        std::string message;

        ASSERT_TRUE(
            horizontal.build(Testing::channelStrip(kColumns, kLength, kWidth), message))
            << message;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = kBed;

        SigmaSpec spec;
        spec.layerCount        = kLayers;
        spec.minLayerThickness = 0.0;

        ASSERT_TRUE(layered.build(horizontal, spec, /*surface=*/0.0, message)) << message;

        temperature = registry.addTemperature(surfaceTemperature, message);
        ASSERT_GE(temperature, 0) << message;

        if (withTracer)
        {
            FVQual::State::SpeciesDefinition definition;
            definition.id           = "tracer";
            definition.units        = "g/m3";
            definition.initialValue = 0.0;

            tracer = registry.add(definition, message);
            ASSERT_GE(tracer, 0) << message;
        }

        state.allocate(layered, registry);

        for (int64_t c = 0; c < layered.columnCount(); ++c)
        {
            for (int k = 0; k < kLayers; ++k)
            {
                const int64_t cell     = layered.cell(c, k);
                const double  centre   = 0.5 * (layered.z(c, k) + layered.z(c, k + 1));
                const double  fraction = centre / kBed;

                state.speciesAt(temperature, cell) =
                    surfaceTemperature + fraction * (bedTemperature - surfaceTemperature);

                if (tracer >= 0)
                {
                    // A depth-varying tracer, so the transfer has something to
                    // carry that is not temperature and cannot ride the density
                    // machinery by accident.
                    state.speciesAt(tracer, cell) = 100.0 * fraction;
                }
            }
        }

        InternalModeOptions options;
        options.gravity   = kGravity;
        options.manningN  = 0.0;
        options.cflNumber = 0.5;

        ASSERT_TRUE(solver.initialize(layered, state, registry, options, message))
            << message;
    }

    [[nodiscard]] double totalMass(int species) const
    {
        return solver.totalSpeciesMass(species);
    }

    [[nodiscard]] double columnMeanTemperature(int64_t column) const
    {
        double weighted = 0.0;
        double volume   = 0.0;

        for (int k = 0; k < kLayers; ++k)
        {
            const int64_t cell = layered.cell(column, k);
            const double  v    = layered.cellVolume()[static_cast<size_t>(cell)];

            weighted += state.speciesAt(temperature, cell) * v;
            volume += v;
        }

        return weighted / volume;
    }

    int addTransfer(int64_t from, int64_t to, double rate, double intakeElevation)
    {
        Transfer transfer;
        transfer.id            = "pump";
        transfer.fromColumn    = from;
        transfer.toColumn      = to;
        transfer.discharge     = rate;
        transfer.fromElevation = intakeElevation;
        transfer.fromPlacement = Placement::SelectiveWithdrawal;
        transfer.toPlacement   = Placement::DensityMatched;

        std::string message;
        const int   index = solver.transfers().add(transfer, solver.sources(), message);

        EXPECT_GE(index, 0) << message;
        return index;
    }
};

} // namespace

// ── Configuration ─────────────────────────────────────────────────────────

TEST(Transfers, MalformedTransfersAreRefused)
{
    SourceSet   sources;
    TransferSet transfers;
    std::string message;

    Transfer selfLoop;
    selfLoop.fromColumn = 3;
    selfLoop.toColumn   = 3;
    selfLoop.discharge  = 5.0;
    EXPECT_LT(transfers.add(selfLoop, sources, message), 0);
    EXPECT_FALSE(message.empty());

    Transfer negativeColumn;
    negativeColumn.fromColumn = -1;
    negativeColumn.toColumn   = 2;
    EXPECT_LT(transfers.add(negativeColumn, sources, message), 0);

    Transfer backwards;
    backwards.fromColumn = 0;
    backwards.toColumn   = 4;
    backwards.discharge  = -10.0;
    EXPECT_LT(transfers.add(backwards, sources, message), 0);

    Transfer fine;
    fine.fromColumn = 0;
    fine.toColumn   = 4;
    fine.discharge  = 10.0;
    EXPECT_GE(transfers.add(fine, sources, message), 0) << message;

    // One transfer, two sources, and they must be equal and opposite.
    EXPECT_EQ(sources.count(), 2);
    ASSERT_TRUE(transfers.updateDischarges(sources, message)) << message;
    EXPECT_NEAR(sources.netDischarge(), 0.0, 1.0e-12)
        << "a transfer is internal and cannot be a net source";

    EXPECT_FALSE(transfers.setDischarge(0, -1.0, message));
}

// ── Bookkeeping ───────────────────────────────────────────────────────────

TEST(Transfers, VolumeMovesBetweenColumnsAndNoneIsLost)
{
    Basin basin;
    basin.build(22.0, 8.0);

    const int index = basin.addTransfer(/*from=*/1, /*to=*/6, /*rate=*/25.0,
                                        /*intake=*/-17.0);
    ASSERT_GE(index, 0);

    const double before = basin.layered.totalVolume();

    std::string message;
    ASSERT_TRUE(basin.solver.advance(150.0, message)) << message;

    // Total volume is untouched: the water went somewhere else in the model,
    // not out of it.
    EXPECT_NEAR(basin.layered.totalVolume(), before, 1.0e-6 * before);

    // And it went the right way.
    const auto &surface = basin.solver.external().surfaceElevation();

    EXPECT_LT(surface[1], surface[6])
        << "the donating column should sit below the receiving one";

    EXPECT_LT(basin.solver.diagnostics().surfaceFluxResidual, 1.0e-10);
}

TEST(Transfers, ConstituentMassIsConservedExactly)
{
    // The identity the design rests on. The withdrawal removes dt*Q*sum(w_k c_k)
    // and the inflow adds dt*Q*c_in with c_in the mean through those same
    // weights, so the two cancel term by term. Asserted at round-off, because a
    // tolerance here would hide precisely the bug worth catching.
    Basin basin;
    basin.build(22.0, 8.0, /*withTracer=*/true);

    const int index = basin.addTransfer(1, 6, 25.0, -17.0);
    ASSERT_GE(index, 0);

    const double heatBefore   = basin.totalMass(basin.temperature);
    const double tracerBefore = basin.totalMass(basin.tracer);

    ASSERT_GT(tracerBefore, 0.0) << "the tracer fixture must actually hold mass";

    std::string message;
    ASSERT_TRUE(basin.solver.advance(150.0, message)) << message;

    const double heatAfter   = basin.totalMass(basin.temperature);
    const double tracerAfter = basin.totalMass(basin.tracer);

    // Relative drift, at the level a sum of products can be expected to hold.
    EXPECT_LT(std::abs(heatAfter - heatBefore) / heatBefore, 1.0e-12)
        << "heat drifted by " << (heatAfter - heatBefore);

    EXPECT_LT(std::abs(tracerAfter - tracerBefore) / tracerBefore, 1.0e-12)
        << "tracer mass drifted by " << (tracerAfter - tracerBefore);
}

// ── Physics: the water that arrives is the water that left ────────────────

TEST(Transfers, IntakeElevationSetsWhatIsDelivered)
{
    // The reason an intake elevation is worth modeling. Same reservoir, same
    // rate, two different intakes — and the delivered water differs by most of
    // the reservoir's temperature range.
    Basin deep;
    deep.build(22.0, 8.0);
    ASSERT_GE(deep.addTransfer(1, 6, 25.0, /*intake=*/-18.0), 0);

    Basin shallow;
    shallow.build(22.0, 8.0);
    ASSERT_GE(shallow.addTransfer(1, 6, 25.0, /*intake=*/-2.0), 0);

    std::string message;
    ASSERT_TRUE(deep.solver.advance(60.0, message)) << message;
    ASSERT_TRUE(shallow.solver.advance(60.0, message)) << message;

    const double deepDelivered    = deep.solver.transfers().deliveredTemperature(0);
    const double shallowDelivered = shallow.solver.transfers().deliveredTemperature(0);

    EXPECT_LT(deepDelivered, 13.0)
        << "a deep intake in a 22-to-8 degree profile delivered " << deepDelivered;
    EXPECT_GT(shallowDelivered, 17.0)
        << "a shallow intake delivered " << shallowDelivered;

    EXPECT_GT(shallowDelivered - deepDelivered, 5.0)
        << "intake elevation barely mattered, which cannot be right in a "
           "stratified reservoir";
}

TEST(Transfers, TheDeliveredWaterIsWhatTheIntakeDraws)
{
    // Not a default, not the receiving column's own water: the donor's.
    Basin basin;
    basin.build(22.0, 8.0);

    const int index = basin.addTransfer(1, 6, 25.0, -18.0);
    ASSERT_GE(index, 0);

    std::string message;
    ASSERT_TRUE(basin.solver.advance(60.0, message)) << message;

    const double delivered = basin.solver.transfers().deliveredTemperature(index);

    // It matches the envelope-weighted mean of the donating column, not the
    // mean of the receiving one.
    const double donorMean     = basin.columnMeanTemperature(1);
    const double recipientMean = basin.columnMeanTemperature(6);

    EXPECT_LT(delivered, donorMean)
        << "a deep intake should draw colder than its own column's average";

    // The decisive comparison, stated as a property rather than a margin: what
    // arrives is hypolimnetic in character — nearer the bed temperature than
    // the surface one — and colder than the column it is arriving in.
    EXPECT_LT(std::abs(delivered - 8.0), std::abs(delivered - 22.0))
        << "delivered " << delivered
        << " degC, which is closer to the surface water than to the deep water "
           "the intake sits in";

    EXPECT_LT(delivered, recipientMean)
        << "delivered " << delivered << " into a column averaging " << recipientMean
        << ", which is the receiving column's own water rather than the donor's";

    // And it tracks the intake's weighted mean, which is what it is by
    // construction; the small gap is the drift since it was last recorded.
    EXPECT_NEAR(delivered,
                basin.solver.sources().weightedMean(
                    basin.layered, basin.solver.transfers().intakeSource(index),
                    basin.state.speciesBlock(basin.temperature)),
                0.5)
        << "the delivered temperature should track the intake's weighted mean";
}

TEST(Transfers, ADestratifyingPumpCoolsTheColumnItFeeds)
{
    // The end-to-end operational claim: lifting hypolimnetic water into another
    // column makes that column colder. A transfer that delivered a stated
    // composition instead of the donor's would not do this.
    Basin basin;
    basin.build(22.0, 8.0);

    const double recipientBefore = basin.columnMeanTemperature(6);

    ASSERT_GE(basin.addTransfer(/*from=*/1, /*to=*/6, /*rate=*/60.0,
                                /*intake=*/-18.0),
              0);

    std::string message;
    ASSERT_TRUE(basin.solver.advance(600.0, message)) << message;

    const double recipientAfter = basin.columnMeanTemperature(6);

    EXPECT_LT(recipientAfter, recipientBefore)
        << "pumping cold hypolimnetic water in warmed the column, from "
        << recipientBefore << " to " << recipientAfter;
}

TEST(Transfers, AnIdleTransferDoesNothing)
{
    // The control. Registered but not running, it must not perturb the model —
    // the failure mode where a paired source leaks because one half was written
    // and the other was not.
    Basin basin;
    basin.build(22.0, 8.0);

    ASSERT_GE(basin.addTransfer(1, 6, /*rate=*/0.0, -18.0), 0);

    const double before = basin.layered.totalVolume();
    const double heat   = basin.totalMass(basin.temperature);

    std::string message;
    ASSERT_TRUE(basin.solver.advance(300.0, message)) << message;

    EXPECT_NEAR(basin.layered.totalVolume(), before, 1.0e-9 * before);
    EXPECT_LT(std::abs(basin.totalMass(basin.temperature) - heat) / heat, 1.0e-12);
}
