// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_sources.cpp
 * \brief Verification of inflows and withdrawals (F3.4, Phase 4).
 *
 * Two families of property are at stake and they fail differently.
 *
 * *Bookkeeping* — volume and constituent budgets must close once water enters
 * and leaves. A source counted once in the surface equation and differently in
 * the layer equations produces fields that still look plausible while the
 * budget quietly drifts, which is precisely the failure a water-quality model
 * cannot tolerate.
 *
 * *Placement* — a plunging inflow must find its own density level, and a
 * withdrawal envelope must narrow as stratification strengthens. These are the
 * behaviors reservoir operations actually depend on, and getting the budget
 * right while placing the water in the wrong layer would be a model that
 * balances and lies.
 */

#include "fvqual/forcing/sources.h"
#include "fvqual/hydro/equationofstate.h"
#include "fvqual/hydro/internalmode.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

using FVQual::Forcing::Placement;
using FVQual::Forcing::Source;
using FVQual::Forcing::SourceSet;
using FVQual::Hydro::EquationOfState;
using FVQual::Hydro::InternalModeOptions;
using FVQual::Hydro::InternalModeSolver;
using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;
using FVQual::State::SpeciesDefinition;
using FVQual::State::SpeciesRegistry;
using FVQual::State::StateStore;
namespace Testing = FVQual::Testing;

namespace {

constexpr double kGravity = 9.80665;
constexpr double kReferenceDensity = 1000.0;

/*!
 * \brief A stratified channel model with sources attached.
 */
struct Model
{
    HorizontalMesh     horizontal;
    LayeredMesh        layered;
    SpeciesRegistry    registry;
    StateStore         state;
    InternalModeSolver solver;

    int temperature = -1;
    int tracer      = -1;

    template <typename TemperatureFunction>
    void build(int columns, int layers, double length, double depth,
               TemperatureFunction temp, bool withTracer = false)
    {
        std::string message;

        ASSERT_TRUE(horizontal.build(Testing::channelStrip(columns, length, 50.0), message))
            << message;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = -depth;

        SigmaSpec spec;
        spec.layerCount        = layers;
        spec.minLayerThickness = 0.0;

        ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

        ASSERT_GE(registry.addTemperature(10.0, message), 0) << message;
        temperature = registry.temperatureIndex();

        if (withTracer)
        {
            SpeciesDefinition dye;
            dye.id           = "DYE";
            dye.units        = "g/m3";
            dye.initialValue = 0.0;
            ASSERT_GE(registry.add(dye, message), 0) << message;
            tracer = registry.indexOf("DYE");
        }

        state.allocate(layered, registry);

        for (int64_t c = 0; c < layered.columnCount(); ++c)
        {
            for (int k = 0; k < layers; ++k)
            {
                const double z = 0.5 * (layered.z(c, k) + layered.z(c, k + 1));
                state.speciesAt(temperature, layered.cell(c, k)) = temp(z);
            }
        }

        InternalModeOptions options;
        options.gravity   = kGravity;
        options.manningN  = 0.0;
        options.cflNumber = 0.5;

        ASSERT_TRUE(solver.initialize(layered, state, registry, options, message))
            << message;
    }

    [[nodiscard]] double tracerMass() const
    {
        return tracer >= 0 ? state.totalMass(tracer, layered) : 0.0;
    }
};

/*!
 * \brief A sharp thermocline at -5 m: 25 °C above, 5 °C below.
 */
double twoLayerProfile(double z)
{
    return z > -5.0 ? 25.0 : 5.0;
}

/*!
 * \brief Builds a standalone layered mesh and density field for placement
 *        tests that need no time stepping.
 */
struct StaticColumn
{
    HorizontalMesh      horizontal;
    LayeredMesh         layered;
    std::vector<double> density;

    void build(int layers, double depth, double thermoclineDepth, double warm,
               double cold)
    {
        std::string message;

        ASSERT_TRUE(horizontal.build(Testing::channelStrip(3, 30.0, 10.0), message))
            << message;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = -depth;

        SigmaSpec spec;
        spec.layerCount        = layers;
        spec.minLayerThickness = 0.0;

        ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

        density.assign(static_cast<size_t>(layered.cellCount()), 0.0);

        for (int64_t c = 0; c < layered.columnCount(); ++c)
        {
            for (int k = 0; k < layers; ++k)
            {
                const double z = 0.5 * (layered.z(c, k) + layered.z(c, k + 1));
                density[static_cast<size_t>(layered.cell(c, k))] =
                    EquationOfState::density(z > -thermoclineDepth ? warm : cold);
            }
        }
    }

    /*!
     * \brief A column stratified in three bands: steep, gentle, steep.
     *
     * Gives a zone that develops freely (neither limit reaches a boundary) but
     * asymmetrically about a port placed off-center in the gentle band, which
     * is what eq. 4 is about.
     */
    void buildBanded(int layers, double depth, double gentleTop, double gentleBottom)
    {
        build(layers, depth, depth, 28.0, 4.0);

        for (int64_t c = 0; c < layered.columnCount(); ++c)
        {
            for (int k = 0; k < layers; ++k)
            {
                const double z = 0.5 * (layered.z(c, k) + layered.z(c, k + 1));

                double temperature;

                if (z > gentleTop)
                    temperature = 28.0 + (z - 0.0) / (gentleTop - 0.0) * (20.0 - 28.0);
                else if (z > gentleBottom)
                    temperature =
                        20.0 + (z - gentleTop) / (gentleBottom - gentleTop) * (18.0 - 20.0);
                else
                    temperature =
                        18.0 + (z - gentleBottom) / (-depth - gentleBottom) * (4.0 - 18.0);

                density[static_cast<size_t>(layered.cell(c, k))] =
                    EquationOfState::density(temperature);
            }
        }
    }

    /*!
     * \brief The same column but linearly stratified top to bottom.
     */
    void buildLinear(int layers, double depth, double warm, double cold)
    {
        build(layers, depth, depth, warm, cold);

        for (int64_t c = 0; c < layered.columnCount(); ++c)
        {
            for (int k = 0; k < layers; ++k)
            {
                const double z = 0.5 * (layered.z(c, k) + layered.z(c, k + 1));
                const double f = -z / depth;

                density[static_cast<size_t>(layered.cell(c, k))] =
                    EquationOfState::density(warm + f * (cold - warm));
            }
        }
    }
};

} // namespace

// ── Configuration ─────────────────────────────────────────────────────────

TEST(Sources, ContradictorySpecificationsAreRefused)
{
    SourceSet   sources;
    std::string message;

    // Selective withdrawal describes how water leaves, so it cannot be given a
    // positive discharge — the envelope theory has no meaning for an inflow.
    Source impossible;
    impossible.placement = Placement::SelectiveWithdrawal;
    impossible.discharge = 5.0;

    EXPECT_LT(sources.add(impossible, message), 0);
    EXPECT_FALSE(message.empty());

    Source backwards;
    backwards.placement = Placement::DensityMatched;
    backwards.discharge = -5.0;

    EXPECT_LT(sources.add(backwards, message), 0);

    Source valid;
    valid.placement = Placement::SelectiveWithdrawal;
    valid.discharge = -5.0;
    valid.elevation = -8.0;

    const int index = sources.add(valid, message);
    ASSERT_GE(index, 0) << message;

    // A schedule may close a gate, but must not reverse its physics.
    EXPECT_TRUE(sources.setDischarge(index, 0.0, message));
    EXPECT_FALSE(sources.setDischarge(index, 3.0, message));
}

TEST(Sources, DistributionSumsToTheColumnTotal)
{
    // The property the whole mass-consistency contract rests on: whatever is
    // spread over the cells must equal what the surface equation is told.
    StaticColumn column;
    column.build(12, 24.0, 6.0, 25.0, 5.0);

    SourceSet   sources;
    std::string message;

    Source inflow;
    inflow.column      = 0;
    inflow.discharge   = 12.0;
    inflow.placement   = Placement::DensityMatched;
    inflow.temperature = 14.0;
    ASSERT_GE(sources.add(inflow, message), 0) << message;

    Source outflow;
    outflow.column      = 2;
    outflow.discharge   = -7.0;
    outflow.placement   = Placement::SelectiveWithdrawal;
    outflow.elevation   = -10.0;
    ASSERT_GE(sources.add(outflow, message), 0) << message;

    Source spread;
    spread.column    = 1;
    spread.discharge = 3.0;
    spread.placement = Placement::Uniform;
    ASSERT_GE(sources.add(spread, message), 0) << message;

    std::vector<double> cells;
    std::vector<double> columns;

    sources.distribute(column.layered, column.density, kGravity, kReferenceDensity, cells,
                       columns);

    for (int64_t c = 0; c < column.layered.columnCount(); ++c)
    {
        double summed = 0.0;

        for (int k = 0; k < column.layered.layerCount(); ++k)
            summed += cells[static_cast<size_t>(column.layered.cell(c, k))];

        EXPECT_NEAR(summed, columns[static_cast<size_t>(c)], 1.0e-12) << "column " << c;
    }

    double total = 0.0;
    for (double value : cells)
        total += value;

    EXPECT_NEAR(total, sources.netDischarge(), 1.0e-12);
    EXPECT_NEAR(total, 8.0, 1.0e-12);
}

// ── Placement ─────────────────────────────────────────────────────────────

TEST(Sources, PlungingInflowFindsItsNeutralLevel)
{
    // The behavior that makes a tributary a three-dimensional problem: an
    // inflow colder than the epilimnion sinks past it rather than mixing in.
    StaticColumn column;
    column.build(20, 20.0, 5.0, 25.0, 5.0);

    SourceSet   sources;
    std::string message;

    const auto placedDepthFor = [&](double temperature) {
        SourceSet set;
        Source    inflow;
        inflow.column      = 0;
        inflow.discharge   = 10.0;
        inflow.placement   = Placement::DensityMatched;
        inflow.temperature = temperature;

        EXPECT_GE(set.add(inflow, message), 0) << message;

        std::vector<double> cells;
        std::vector<double> columns;
        set.distribute(column.layered, column.density, kGravity, kReferenceDensity, cells,
                       columns);

        return set.placedElevation()[0];
    };

    // Warmer than the whole reservoir: rides on the surface as an overflow.
    EXPECT_GT(placedDepthFor(30.0), -2.0);

    // Denser than everything: sinks to the bed as an underflow. Note the
    // temperature that achieves this is 4 C, not 0 — water is densest at 4 C,
    // so a 2 C inflow is *lighter* than a 5 C hypolimnion and correctly floats
    // on top of it. That is winter inverse stratification, and it is the whole
    // reason the equation of state carries the full polynomial rather than a
    // linear fit.
    EXPECT_LT(placedDepthFor(4.0), -18.0);
    EXPECT_GT(placedDepthFor(2.0), -8.0);

    // Between the two layers: settles at the thermocline as an interflow,
    // which is the case that matters and the one a depth-averaged or
    // surface-only treatment gets wrong.
    const double interflow = placedDepthFor(15.0);
    EXPECT_LT(interflow, -4.0);
    EXPECT_GT(interflow, -7.0);
}

TEST(Sources, WithdrawalEnvelopeNarrowsWithStratification)
{
    // Selective withdrawal only works as an operational lever because the
    // envelope is thin. Weak stratification widens it until the outlet draws
    // from most of the column — which is exactly why an operator can aim a
    // port in summer and not in winter.
    const auto envelopeFor = [&](double warm, double cold) {
        StaticColumn column;
        column.build(20, 20.0, 10.0, warm, cold);

        SourceSet   set;
        std::string message;

        Source outflow;
        outflow.column    = 1;
        outflow.discharge = -20.0;
        outflow.placement = Placement::SelectiveWithdrawal;
        outflow.elevation = -10.0;

        EXPECT_GE(set.add(outflow, message), 0) << message;

        std::vector<double> cells;
        std::vector<double> columns;
        set.distribute(column.layered, column.density, kGravity, kReferenceDensity, cells,
                       columns);

        return set.envelopeHalfHeight()[0];
    };

    const double strong = envelopeFor(28.0, 4.0);
    const double weak   = envelopeFor(12.0, 10.0);

    EXPECT_GT(strong, 0.0);
    EXPECT_LT(strong, weak) << "strong " << strong << " vs weak " << weak;
}

TEST(Sources, AThermoclineExcludesTheEpilimnionEntirely)
{
    // A step profile: well-mixed 28 degC above 10 m, well-mixed 4 degC below.
    // A low-level port in the hypolimnion cannot reach across the thermocline,
    // and inside the well-mixed layer there is no buoyancy to favor one depth
    // over another — so the draw is uniform there and exactly zero above.
    //
    // The uniformity is a property of the density-anomaly velocity profile the
    // reference formulation uses (Bohan and Grace 1973). A profile written in
    // distance from the port instead would peak at the port and taper through
    // water that is physically indistinguishable from it, which is the wrong
    // answer for a well-mixed layer.
    StaticColumn column;
    column.build(20, 20.0, 10.0, 28.0, 4.0);

    SourceSet   sources;
    std::string message;

    Source outflow;
    outflow.column    = 1;
    outflow.discharge = -15.0;
    outflow.placement = Placement::SelectiveWithdrawal;
    outflow.elevation = -15.0; // a low-level port

    ASSERT_GE(sources.add(outflow, message), 0) << message;

    std::vector<double> cells;
    std::vector<double> columns;
    sources.distribute(column.layered, column.density, kGravity, kReferenceDensity, cells,
                       columns);

    const auto drawnAt = [&](int k) {
        return std::abs(cells[static_cast<size_t>(column.layered.cell(1, k))]);
    };

    double epilimnion = 0.0;
    double hypolimnion = 0.0;

    for (int k = 0; k < column.layered.layerCount(); ++k)
    {
        const double centre =
            0.5 * (column.layered.z(1, k) + column.layered.z(1, k + 1));

        if (centre > -10.0)
            epilimnion += drawnAt(k);
        else
            hypolimnion += drawnAt(k);
    }

    EXPECT_NEAR(epilimnion, 0.0, 1.0e-12)
        << "the port drew " << epilimnion << " m3/s across a sharp thermocline";

    EXPECT_NEAR(hypolimnion, 15.0, 1.0e-9)
        << "the whole withdrawal must come from below the thermocline";

    // Uniform within the well-mixed layer.
    for (int k = 10; k < 20; ++k)
        EXPECT_NEAR(drawnAt(k), 1.5, 1.0e-9) << "at layer " << k;
}

TEST(Sources, AnAsymmetricZoneShiftsTheVelocityMaximumOffThePort)
{
    // Bohan and Grace (1973) eq. 4. The fastest approach velocity sits at the
    // port only when the withdrawal zone is symmetric about it; where the
    // stratification lets the zone grow further on one side than the other, the
    // maximum moves. The check computes eq. 4 from the model's *own* zone
    // limits and asks where the draw actually peaks, so it tests the relation
    // rather than a memorized elevation.
    // The gentle band is bounded above and below by steeply stratified water,
    // so the zone still develops freely; placing the port well off-center
    // *within* that band is what makes it lopsided.
    constexpr double kPort      = -6.0;
    constexpr double kDischarge = 8.0;

    StaticColumn column;
    column.buildBanded(80, 20.0, /*gentleTop=*/-4.0, /*gentleBottom=*/-16.0);

    SourceSet   sources;
    std::string message;

    Source outflow;
    outflow.column    = 1;
    outflow.discharge = -kDischarge;
    outflow.placement = Placement::SelectiveWithdrawal;
    outflow.elevation = kPort;

    ASSERT_GE(sources.add(outflow, message), 0) << message;

    std::vector<double> cells;
    std::vector<double> columns;
    sources.distribute(column.layered, column.density, kGravity, kReferenceDensity, cells,
                       columns);

    const LayeredMesh &mesh   = column.layered;
    const int          layers = mesh.layerCount();

    const auto centreOf = [&](int k) {
        return 0.5 * (mesh.z(1, k) + mesh.z(1, k + 1));
    };
    const auto drawnAt = [&](int k) {
        return std::abs(cells[static_cast<size_t>(mesh.cell(1, k))]);
    };

    int first = -1;
    int last  = -1;
    int peak  = -1;
    double peakDraw = 0.0;

    for (int k = 0; k < layers; ++k)
    {
        const double drawn = drawnAt(k);

        if (drawn > 0.0)
        {
            if (first < 0)
                first = k;
            last = k;
        }

        if (drawn > peakDraw)
        {
            peakDraw = drawn;
            peak     = k;
        }
    }

    ASSERT_GE(first, 1) << "the zone must not reach the surface for eq. 4 to apply";
    ASSERT_LE(last, layers - 2) << "the zone must not reach the bed for eq. 4 to apply";

    // The limit layers carry zero velocity, so the zone edges sit one outside
    // the first and last layers that draw.
    const double above = centreOf(first - 1) - kPort;
    const double below = kPort - centreOf(last + 1);

    ASSERT_GT(std::abs(above - below), 1.0)
        << "this fixture is meant to produce an asymmetric zone, but got above="
        << above << " below=" << below;

    // Eq. 4, computed here from the model's own limits.
    const double thickness = above + below;
    const double height =
        thickness * std::pow(std::sin(1.57 * below / thickness), 2.0);
    const double expected = kPort - below + height;

    int nearest = 0;
    double best = 1.0e30;

    for (int k = first - 1; k <= last + 1; ++k)
    {
        const double gap = std::abs(centreOf(k) - expected);

        if (gap < best)
        {
            best    = gap;
            nearest = k;
        }
    }

    EXPECT_EQ(peak, nearest)
        << "the draw peaked at " << centreOf(peak) << " m but eq. 4 puts the "
           "velocity maximum at " << expected << " m (layer centre "
        << centreOf(nearest) << " m)";

    // The discriminating part: the maximum is genuinely off the port. Without
    // eq. 4 the profile would be centered on the port layer and this would fail.
    int portLayer = 0;
    best = 1.0e30;

    for (int k = 0; k < layers; ++k)
    {
        const double gap = std::abs(centreOf(k) - kPort);

        if (gap < best)
        {
            best      = gap;
            portLayer = k;
        }
    }

    EXPECT_NE(peak, portLayer)
        << "the maximum stayed on the port layer, so the zone asymmetry was ignored";

    // And it moved toward the short side of the zone, as eq. 4 requires.
    if (below > above)
    {
        EXPECT_GT(centreOf(peak), kPort)
            << "the zone runs " << below << " m below the port and only " << above
            << " m above it, so the maximum should sit above the port";
    }
    else
    {
        EXPECT_LT(centreOf(peak), kPort)
            << "the zone runs " << above << " m above the port and only " << below
            << " m below it, so the maximum should sit below the port";
    }
}

TEST(Sources, TheZoneLimitSatisfiesThePublishedFroudeCondition)
{
    // The formulation check, and the one that would have caught the leading
    // coefficients being wrong. The withdrawal zone limit is where the
    // separation from the port catches up with the zone thickness the local
    // buoyancy frequency implies (Bohan and Grace 1973; Smith et al. 1987):
    //
    //     H  >=  ( c_bi Q / N(H) )^(1/3),     N(H) = sqrt( g dRho / (H rho) )
    //
    // Here the condition is evaluated independently, from the density field
    // and nothing else, and compared against where the model actually stopped
    // drawing. A test that only asked "does a deep port draw cold water" passes
    // whatever the constants are; this one does not.
    constexpr double kDischarge = 15.0;
    constexpr double kPort      = -12.0;

    StaticColumn column;
    column.buildLinear(20, 20.0, 28.0, 4.0);

    SourceSet   sources;
    std::string message;

    Source outflow;
    outflow.column    = 1;
    outflow.discharge = -kDischarge;
    outflow.placement = Placement::SelectiveWithdrawal;
    outflow.elevation = kPort;

    ASSERT_GE(sources.add(outflow, message), 0) << message;

    std::vector<double> cells;
    std::vector<double> columns;
    sources.distribute(column.layered, column.density, kGravity, kReferenceDensity, cells,
                       columns);

    const LayeredMesh &mesh   = column.layered;
    const int          layers = mesh.layerCount();

    const auto centreOf = [&](int k) {
        return 0.5 * (mesh.z(1, k) + mesh.z(1, k + 1));
    };
    const auto densityAt = [&](int k) {
        return column.density[static_cast<size_t>(mesh.cell(1, k))];
    };

    // The port's own layer, and the interference coefficient for its position.
    int outlet = 0;
    {
        double best = 1.0e30;
        for (int k = 0; k < layers; ++k)
        {
            const double gap = std::abs(centreOf(k) - kPort);
            if (gap < best) { best = gap; outlet = k; }
        }
    }

    const double surface  = mesh.z(1, 0);
    const double bed      = mesh.z(1, layers);
    const double position = (kPort - bed) / (surface - bed);
    const double cbi      = (position < 0.1 || position > 0.9) ? 2.0 : 1.0;

    ASSERT_DOUBLE_EQ(cbi, 1.0) << "this fixture is meant to be a free-developing zone";

    // Independently march to the upper limit.
    int expectedTop = 0;

    for (int k = outlet - 1; k >= 0; --k)
    {
        const double separation = centreOf(k) - kPort;

        if (separation <= 0.0)
            continue;

        const double contrast = std::abs(densityAt(k) - densityAt(outlet));

        if (contrast <= 0.0)
            continue;

        const double frequency =
            std::sqrt(kGravity * contrast / (separation * densityAt(outlet)));
        const double thickness = std::cbrt(cbi * kDischarge / frequency);

        if (separation >= thickness)
        {
            expectedTop = k;
            break;
        }
    }

    // Where the model actually stopped drawing. The limit layer itself carries
    // *zero* velocity by construction — the parabola is normalized by the
    // density contrast at the limit, so it vanishes exactly there — and the
    // reference implementation has the same property. The first layer with a
    // nonzero draw is therefore one step inside the limit.
    int actualTop = layers;

    for (int k = 0; k < layers; ++k)
    {
        if (std::abs(cells[static_cast<size_t>(mesh.cell(1, k))]) > 0.0)
        {
            actualTop = k;
            break;
        }
    }

    EXPECT_EQ(actualTop, expectedTop + 1)
        << "the draw begins at layer " << actualTop << " (centre "
        << centreOf(actualTop) << " m), which puts the zone limit at layer "
        << (actualTop - 1) << ", but the Froude condition puts it at layer "
        << expectedTop << " (centre " << centreOf(expectedTop) << " m)";

    // And the condition really does bracket the limit: satisfied at the limit,
    // not satisfied one layer inside it. Without this the equality above could
    // be matched by any monotone rule.
    const auto residual = [&](int k) {
        const double separation = centreOf(k) - kPort;
        const double contrast   = std::abs(densityAt(k) - densityAt(outlet));
        const double frequency =
            std::sqrt(kGravity * contrast / (separation * densityAt(outlet)));

        return separation - std::cbrt(cbi * kDischarge / frequency);
    };

    ASSERT_LT(expectedTop + 1, outlet) << "need a layer strictly inside the limit";

    EXPECT_GE(residual(expectedTop), 0.0) << "the limit does not satisfy the condition";
    EXPECT_LT(residual(expectedTop + 1), 0.0)
        << "the layer inside the limit already satisfied it, so the zone stopped short";
}

TEST(Sources, UnderLinearStratificationTheDrawPeaksAtThePort)
{
    // The counterpart. With a continuous density gradient the layers are no
    // longer interchangeable, and the parabolic profile does peak at the port
    // and fall away from it — the behavior the step-profile case cannot show.
    StaticColumn column;
    column.buildLinear(20, 20.0, /*warm=*/28.0, /*cold=*/4.0);

    SourceSet   sources;
    std::string message;

    Source outflow;
    outflow.column    = 1;
    outflow.discharge = -15.0;
    outflow.placement = Placement::SelectiveWithdrawal;
    outflow.elevation = -12.0;

    ASSERT_GE(sources.add(outflow, message), 0) << message;

    std::vector<double> cells;
    std::vector<double> columns;
    sources.distribute(column.layered, column.density, kGravity, kReferenceDensity, cells,
                       columns);

    int    strongest = 0;
    double peak      = 0.0;

    for (int k = 0; k < column.layered.layerCount(); ++k)
    {
        const double drawn =
            std::abs(cells[static_cast<size_t>(column.layered.cell(1, k))]);

        if (drawn > peak)
        {
            peak      = drawn;
            strongest = k;
        }
    }

    const double centre =
        0.5 * (column.layered.z(1, strongest) + column.layered.z(1, strongest + 1));

    EXPECT_NEAR(centre, -12.0, 1.5) << "withdrawal peaked at " << centre << " m";

    // And it falls away on both sides rather than sitting flat.
    const auto drawnAt = [&](int k) {
        return std::abs(cells[static_cast<size_t>(column.layered.cell(1, k))]);
    };

    EXPECT_LT(drawnAt(strongest - 3), peak) << "no fall-off above the port";
    EXPECT_LT(drawnAt(strongest + 3), peak) << "no fall-off below the port";
}

// ── Budgets under throughflow ─────────────────────────────────────────────

TEST(Sources, VolumeBudgetClosesUnderThroughflow)
{
    Model model;
    model.build(20, 12, 400.0, 20.0, twoLayerProfile);

    std::string message;

    Source inflow;
    inflow.column      = 0;
    inflow.discharge   = 25.0;
    inflow.placement   = Placement::DensityMatched;
    inflow.temperature = 12.0;
    ASSERT_GE(model.solver.sources().add(inflow, message), 0) << message;

    Source outflow;
    outflow.column    = 19;
    outflow.discharge = -25.0;
    outflow.placement = Placement::SelectiveWithdrawal;
    outflow.elevation = -12.0;
    ASSERT_GE(model.solver.sources().add(outflow, message), 0) << message;

    const double before = model.layered.totalVolume();

    ASSERT_TRUE(model.solver.advance(200.0, message)) << message;

    // Balanced throughflow: what comes in goes out, so the stored volume is
    // unchanged apart from the transient adjustment of the surface.
    const double drift = std::abs(model.layered.totalVolume() - before) / before;

    EXPECT_LT(drift, 1.0e-3) << "volume drifted by " << drift;

    // The mass-consistency gate must survive the sources — this is what would
    // break first if the column totals and cell shares disagreed.
    EXPECT_LT(model.solver.diagnostics().surfaceFluxResidual, 1.0e-10);
}

TEST(Sources, NetInflowRaisesTheSurfaceByTheRightVolume)
{
    Model model;
    model.build(16, 8, 320.0, 15.0, [](double) { return 10.0; });

    std::string message;

    Source inflow;
    inflow.column      = 0;
    inflow.discharge   = 40.0;
    inflow.placement   = Placement::Uniform;
    inflow.temperature = 10.0;
    ASSERT_GE(model.solver.sources().add(inflow, message), 0) << message;

    const double before   = model.layered.totalVolume();
    const double duration = 100.0;

    ASSERT_TRUE(model.solver.advance(duration, message)) << message;

    const double expected = before + 40.0 * duration;
    const double actual   = model.layered.totalVolume();

    EXPECT_NEAR(actual, expected, 1.0e-6 * expected)
        << "added " << (actual - before) << " m3, expected " << (expected - before);

    EXPECT_LT(model.solver.diagnostics().surfaceFluxResidual, 1.0e-10);
}

TEST(Sources, InflowDeliversItsConstituentLoad)
{
    // A tracer arriving with the inflow must appear in the domain at the rate
    // discharge x concentration — the elementary loading calculation every
    // water-quality answer is built on.
    Model model;
    model.build(16, 8, 320.0, 15.0, [](double) { return 10.0; }, /*withTracer=*/true);

    ASSERT_GE(model.tracer, 0);

    std::string message;

    Source inflow;
    inflow.column      = 0;
    inflow.discharge   = 20.0;
    inflow.placement   = Placement::Uniform;
    inflow.temperature = 10.0;
    inflow.concentration.assign(static_cast<size_t>(model.registry.count()), 0.0);
    inflow.concentration[static_cast<size_t>(model.tracer)] = 5.0;

    ASSERT_GE(model.solver.sources().add(inflow, message), 0) << message;

    EXPECT_NEAR(model.tracerMass(), 0.0, 1.0e-12);

    const double duration = 100.0;
    ASSERT_TRUE(model.solver.advance(duration, message)) << message;

    // 20 m3/s x 5 g/m3 x 100 s = 10000 g.
    const double expected = 20.0 * 5.0 * duration;

    EXPECT_NEAR(model.tracerMass(), expected, 0.05 * expected)
        << "delivered " << model.tracerMass() << " g, expected " << expected;
}

TEST(Sources, WithdrawalRemovesConstituentAtAmbientConcentration)
{
    // The complement: an outlet exports whatever is there, so a domain with a
    // uniform tracer and balanced throughflow of clean water must lose mass at
    // the rate the outlet carries it.
    Model model;
    model.build(16, 8, 320.0, 15.0, [](double) { return 10.0; }, /*withTracer=*/true);

    std::string message;

    for (int64_t i = 0; i < model.state.cellCount(); ++i)
        model.state.speciesAt(model.tracer, i) = 2.0;

    Source outflow;
    outflow.column    = 15;
    outflow.discharge = -30.0;
    outflow.placement = Placement::Uniform;
    ASSERT_GE(model.solver.sources().add(outflow, message), 0) << message;

    const double before = model.tracerMass();
    ASSERT_GT(before, 0.0);

    ASSERT_TRUE(model.solver.advance(60.0, message)) << message;

    const double removed = before - model.tracerMass();

    // 30 m3/s x 2 g/m3 x 60 s = 3600 g, approximately: the ambient
    // concentration near the outlet drifts a little as water leaves.
    EXPECT_GT(removed, 0.5 * 3600.0);
    EXPECT_LT(removed, 1.5 * 3600.0);
}
