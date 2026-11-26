// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_hypsography.cpp
 * \brief Depth-varying plan area in the layered mesh (F1.5).
 *
 * A prism is a poor model of a reservoir, and the numbers are not marginal.
 * DeGray's deepest section runs from 6.4 km wide at the surface to 75 m at the
 * bed — a ratio of 85 — and treating each section as a prism on its surface
 * area overstates the whole reservoir's volume by 2.07x. Almost all of that
 * error sits in the hypolimnion, which is exactly the water whose volume
 * governs how strongly a reservoir stratifies and how fast its oxygen runs out.
 *
 * Two things have to be true for the fix to be safe.
 *
 * **Nothing changes for a model that does not use it.** A column with no curve
 * must produce bit-for-bit the volumes it produced before, and the whole
 * existing physics suite must be unmoved. That is checked here and by the suite
 * as a whole.
 *
 * **The integral is exact, not sampled.** Cell volume is the integral of a
 * piecewise-linear area curve, and a cell that straddles a kink in that curve
 * is where a midpoint rule quietly stops being right — so the gates include a
 * kinked curve and compare against arithmetic done independently.
 */

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
#include <vector>

using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;
namespace Testing = FVQual::Testing;

namespace {

/*!
 * \brief A one-column basin, so the arithmetic is checkable by hand.
 */
struct Column
{
    HorizontalMesh horizontal;
    LayeredMesh    layered;

    double planArea = 0.0;

    void build(int layers, double bed, double surface)
    {
        std::string message;

        ASSERT_TRUE(horizontal.build(Testing::channelStrip(1, 100.0, 50.0), message))
            << message;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = bed;

        SigmaSpec spec;
        spec.layerCount        = layers;
        spec.minLayerThickness = 0.0;

        ASSERT_TRUE(layered.build(horizontal, spec, surface, message)) << message;

        planArea = horizontal.faceArea()[0];
    }

    [[nodiscard]] bool applyCurve(const std::vector<double> &elevation,
                                  const std::vector<double> &area,
                                  std::string &message)
    {
        if (!layered.setHypsography(0, elevation, area, message))
            return false;

        const std::vector<double> surface(1, layered.z(0, 0));

        return layered.updateSurface(surface, message);
    }

    [[nodiscard]] double volume() const { return layered.totalVolume(); }
};

} // namespace

// ── The default is unchanged ──────────────────────────────────────────────

TEST(Hypsography, WithoutACurveTheMeshIsExactlyAPrism)
{
    Column column;
    column.build(/*layers=*/8, /*bed=*/-20.0, /*surface=*/0.0);

    EXPECT_FALSE(column.layered.hasHypsography());

    const double expected = column.planArea * 20.0;

    EXPECT_NEAR(column.volume(), expected, 1.0e-9 * expected);

    for (int k = 0; k < 8; ++k)
    {
        const int64_t cell = column.layered.cell(0, k);

        EXPECT_EQ(column.layered.cellVolume()[static_cast<size_t>(cell)],
                  column.layered.thickness()[static_cast<size_t>(cell)] * column.planArea)
            << "layer " << k;
    }

    // Interface areas fall back to the constant plan area too.
    for (int k = 0; k <= 8; ++k)
    {
        EXPECT_EQ(column.layered.interfaceArea()[static_cast<size_t>(
                      column.layered.interfaceSlot(0, k))],
                  column.planArea);
    }
}

TEST(Hypsography, AConstantCurveReproducesThePrismBitForBit)
{
    // The strongest possible statement that the new path degenerates correctly:
    // not "close to", but the same doubles.
    Column reference;
    reference.build(/*layers=*/6, /*bed=*/-15.0, /*surface=*/0.0);

    Column curved;
    curved.build(/*layers=*/6, /*bed=*/-15.0, /*surface=*/0.0);

    std::string message;
    ASSERT_TRUE(curved.applyCurve({-15.0, 0.0}, {curved.planArea, curved.planArea},
                                  message))
        << message;

    ASSERT_TRUE(curved.layered.hasHypsography());

    for (int k = 0; k < 6; ++k)
    {
        const size_t cell = static_cast<size_t>(reference.layered.cell(0, k));

        EXPECT_DOUBLE_EQ(curved.layered.cellVolume()[cell],
                         reference.layered.cellVolume()[cell])
            << "layer " << k;
    }
}

// ── The integral is exact ─────────────────────────────────────────────────

TEST(Hypsography, ALinearTaperMatchesItsAnalyticVolume)
{
    // A(z) = a + b (z - bed). The exact volume between two elevations is the
    // trapezoid, and the implementation should reproduce it to round-off rather
    // than approximately.
    constexpr double bed     = -20.0;
    constexpr double surface = 0.0;
    constexpr double bedArea = 1000.0;
    constexpr double topArea = 9000.0;

    Column column;
    column.build(/*layers=*/10, bed, surface);

    std::string message;
    ASSERT_TRUE(column.applyCurve({bed, surface}, {bedArea, topArea}, message))
        << message;

    // Whole column: mean area times depth, because the taper is linear.
    const double expected = 0.5 * (bedArea + topArea) * (surface - bed);

    EXPECT_NEAR(column.volume(), expected, 1.0e-9 * expected);

    // And each layer individually.
    const double slope = (topArea - bedArea) / (surface - bed);

    for (int k = 0; k < 10; ++k)
    {
        const double top    = column.layered.z(0, k);
        const double bottom = column.layered.z(0, k + 1);

        const double areaTop    = bedArea + slope * (top - bed);
        const double areaBottom = bedArea + slope * (bottom - bed);

        const double cellExpected = 0.5 * (areaTop + areaBottom) * (top - bottom);

        const size_t cell = static_cast<size_t>(column.layered.cell(0, k));

        EXPECT_NEAR(column.layered.cellVolume()[cell], cellExpected,
                    1.0e-9 * cellExpected)
            << "layer " << k;
    }

    // The taper actually bites: the bottom layer holds far less than the top.
    const size_t top = static_cast<size_t>(column.layered.cell(0, 0));
    const size_t bot = static_cast<size_t>(column.layered.cell(0, 9));

    EXPECT_GT(column.layered.cellVolume()[top],
              3.0 * column.layered.cellVolume()[bot]);
}

TEST(Hypsography, ACellStraddlingAKinkIsStillExact)
{
    // The case that separates an exact integral from a sampled one. The curve
    // bends sharply at mid-depth; a cell spanning the bend gets the right
    // answer only if the integration splits there.
    constexpr double bed     = -10.0;
    constexpr double surface = 0.0;
    constexpr double kink    = -5.0;

    // Nearly flat below the kink, steep above it — a reservoir with a narrow
    // gorge opening into a broad pool.
    const std::vector<double> elevation{bed, kink, surface};
    const std::vector<double> area{100.0, 200.0, 4200.0};

    // Two layers, so the *upper* cell spans [-5, 0] and the lower [-10, -5]:
    // each lands exactly on the kink. Use three layers instead so the middle
    // one straddles it.
    Column column;
    column.build(/*layers=*/3, bed, surface);

    std::string message;
    ASSERT_TRUE(column.applyCurve(elevation, area, message)) << message;

    // Independent arithmetic, splitting at the kink by hand.
    const auto areaAt = [&](double z) {
        if (z <= kink)
            return 100.0 + (200.0 - 100.0) * (z - bed) / (kink - bed);

        return 200.0 + (4200.0 - 200.0) * (z - kink) / (surface - kink);
    };

    const auto exact = [&](double lower, double upper) {
        if (lower >= kink || upper <= kink)
            return 0.5 * (areaAt(lower) + areaAt(upper)) * (upper - lower);

        return 0.5 * (areaAt(lower) + areaAt(kink)) * (kink - lower) +
               0.5 * (areaAt(kink) + areaAt(upper)) * (upper - kink);
    };

    double total = 0.0;

    for (int k = 0; k < 3; ++k)
    {
        const double top    = column.layered.z(0, k);
        const double bottom = column.layered.z(0, k + 1);

        const double cellExpected = exact(bottom, top);

        const size_t cell = static_cast<size_t>(column.layered.cell(0, k));

        EXPECT_NEAR(column.layered.cellVolume()[cell], cellExpected,
                    1.0e-9 * std::max(cellExpected, 1.0))
            << "layer " << k << " spanning " << bottom << " to " << top;

        total += cellExpected;
    }

    EXPECT_NEAR(column.volume(), total, 1.0e-9 * total);

    // The middle layer really does straddle the kink, or this gate is testing
    // nothing it claims to.
    EXPECT_LT(column.layered.z(0, 2), kink);
    EXPECT_GT(column.layered.z(0, 1), kink);
}

TEST(Hypsography, InterfaceAreasAreTheCurveAtThatElevation)
{
    constexpr double bed     = -12.0;
    constexpr double surface = 0.0;

    Column column;
    column.build(/*layers=*/6, bed, surface);

    std::string message;
    ASSERT_TRUE(column.applyCurve({bed, surface}, {500.0, 6500.0}, message)) << message;

    const double slope = (6500.0 - 500.0) / (surface - bed);

    for (int k = 0; k <= 6; ++k)
    {
        const double z        = column.layered.z(0, k);
        const double expected = 500.0 + slope * (z - bed);

        EXPECT_NEAR(column.layered.interfaceArea()[static_cast<size_t>(
                        column.layered.interfaceSlot(0, k))],
                    expected, 1.0e-9 * expected)
            << "interface " << k;
    }
}

// ── The storage curve ─────────────────────────────────────────────────────

TEST(Hypsography, TheVolumeElevationCurveFollowsTheHypsography)
{
    // What a reservoir operator actually reads off: how much water is in the
    // pool at a given level. It is the integral of the area curve, and getting
    // it right is the whole reason for this feature.
    constexpr double bed     = -20.0;
    constexpr double bedArea = 400.0;
    constexpr double topArea = 8400.0;
    constexpr double slope   = (topArea - bedArea) / 20.0;

    Column column;
    column.build(/*layers=*/10, bed, /*surface=*/0.0);

    std::string message;
    ASSERT_TRUE(column.applyCurve({bed, 0.0}, {bedArea, topArea}, message)) << message;

    for (double level : {-15.0, -10.0, -5.0, -1.0})
    {
        const std::vector<double> surface(1, level);

        ASSERT_TRUE(column.layered.updateSurface(surface, message)) << message;

        const double depth    = level - bed;
        const double areaHere = bedArea + slope * depth;
        const double expected = 0.5 * (bedArea + areaHere) * depth;

        EXPECT_NEAR(column.volume(), expected, 1.0e-9 * expected)
            << "at elevation " << level;
    }

    // Storage falls away faster than linearly as the pool drops — which is the
    // behaviour a prism cannot reproduce and the reason this matters.
    const std::vector<double> half(1, -10.0);
    ASSERT_TRUE(column.layered.updateSurface(half, message)) << message;
    const double atHalfDepth = column.volume();

    const std::vector<double> full(1, 0.0);
    ASSERT_TRUE(column.layered.updateSurface(full, message)) << message;
    const double atFull = column.volume();

    EXPECT_LT(atHalfDepth, 0.35 * atFull)
        << "half the depth held " << atHalfDepth / atFull << " of the volume";
}

TEST(Hypsography, ASteepTaperConcentratesVolumeNearTheSurface)
{
    // DeGray's shape, in miniature: 85:1 from surface to bed. The point is that
    // the hypolimnion is a small fraction of the reservoir, which is precisely
    // what a prism gets wrong and what governs stratification and oxygen.
    constexpr double bed     = -26.0;
    constexpr double bedArea = 75.0;
    constexpr double topArea = 6418.0;

    Column column;
    column.build(/*layers=*/26, bed, /*surface=*/0.0);

    std::string message;
    ASSERT_TRUE(column.applyCurve({bed, 0.0}, {bedArea, topArea}, message)) << message;

    const double prism = topArea * 26.0;

    // The prism overstates by close to a factor of two, as it does for DeGray.
    EXPECT_GT(prism / column.volume(), 1.8);
    EXPECT_LT(prism / column.volume(), 2.1);

    // The bottom half of the column holds a small share of the water.
    double lower = 0.0;

    for (int k = 13; k < 26; ++k)
        lower += column.layered.cellVolume()[static_cast<size_t>(column.layered.cell(0, k))];

    EXPECT_LT(lower / column.volume(), 0.3)
        << "the lower half held " << lower / column.volume();
}

// ── Refusals ──────────────────────────────────────────────────────────────

TEST(Hypsography, MalformedCurvesAreRefused)
{
    Column column;
    column.build(/*layers=*/4, /*bed=*/-10.0, /*surface=*/0.0);

    std::string message;

    // Mismatched lengths.
    EXPECT_FALSE(column.layered.setHypsography(0, {{-10.0, 0.0}}, {{100.0}}, message));

    // Too few samples to interpolate between.
    EXPECT_FALSE(column.layered.setHypsography(0, {{-10.0}}, {{100.0}}, message));

    // Elevations must increase.
    EXPECT_FALSE(
        column.layered.setHypsography(0, {{0.0, -10.0}}, {{100.0, 200.0}}, message));
    EXPECT_FALSE(
        column.layered.setHypsography(0, {{-10.0, -10.0}}, {{100.0, 200.0}}, message));

    // A body that widens downward would make the storage curve non-monotone, so
    // a free surface would no longer follow uniquely from a volume.
    EXPECT_FALSE(
        column.layered.setHypsography(0, {{-10.0, 0.0}}, {{500.0, 100.0}}, message));

    // Negative area.
    EXPECT_FALSE(
        column.layered.setHypsography(0, {{-10.0, 0.0}}, {{-1.0, 100.0}}, message));

    // Out of range.
    EXPECT_FALSE(
        column.layered.setHypsography(7, {{-10.0, 0.0}}, {{100.0, 200.0}}, message));

    // A curve set before build() has no mesh to index into.
    LayeredMesh fresh;
    EXPECT_FALSE(fresh.setHypsography(0, {{-10.0, 0.0}}, {{100.0, 200.0}}, message));
}

TEST(Hypsography, TheCurveIsClampedOutsideItsSamples)
{
    // A survey stops somewhere. Extrapolating a linear fit below the deepest
    // sample runs the area negative; above the highest it invents storage the
    // survey never saw.
    Column column;
    column.build(/*layers=*/4, /*bed=*/-10.0, /*surface=*/0.0);

    std::string message;
    ASSERT_TRUE(column.applyCurve({-6.0, -4.0}, {200.0, 600.0}, message)) << message;

    EXPECT_NEAR(column.layered.planArea(0, -20.0), 200.0, 1.0e-12);
    EXPECT_NEAR(column.layered.planArea(0, -6.0), 200.0, 1.0e-12);
    EXPECT_NEAR(column.layered.planArea(0, -5.0), 400.0, 1.0e-12);
    EXPECT_NEAR(column.layered.planArea(0, -4.0), 600.0, 1.0e-12);
    EXPECT_NEAR(column.layered.planArea(0, 50.0), 600.0, 1.0e-12);

    // Nothing is negative anywhere, which is what clamping buys.
    for (double z = -30.0; z <= 30.0; z += 0.5)
        EXPECT_GT(column.layered.planArea(0, z), 0.0) << "at " << z;
}

TEST(Hypsography, ModeConsistencyHoldsOnAHypsographicMesh)
{
    // The gate stage 2 exists to satisfy. The vertical velocity is integrated up
    // from a no-flux bed, so its value at the free surface vanishes only if the
    // internal mode's cell volumes and the external mode's stage agree about how
    // much water a column holds. With a storage curve those are two different
    // calculations of the same thing, and they have to match.
    HorizontalMesh horizontal;
    LayeredMesh    layered;

    std::string message;

    ASSERT_TRUE(horizontal.build(Testing::channelStrip(6, 600.0, 40.0), message))
        << message;

    for (int64_t f = 0; f < horizontal.faceCount(); ++f)
        horizontal.faceBed()[static_cast<size_t>(f)] = -20.0;

    SigmaSpec spec;
    spec.layerCount        = 6;
    spec.minLayerThickness = 0.0;

    ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

    // A steep taper, so the curve is doing real work rather than sitting near
    // the constant it replaced.
    for (int64_t c = 0; c < layered.columnCount(); ++c)
    {
        const double top = horizontal.faceArea()[static_cast<size_t>(c)];

        ASSERT_TRUE(layered.setHypsography(c, {{-20.0, 5.0}}, {{0.15 * top, top}},
                                           message))
            << message;
    }

    const std::vector<double> surface(static_cast<size_t>(layered.columnCount()), 0.0);
    ASSERT_TRUE(layered.updateSurface(surface, message)) << message;

    FVQual::State::SpeciesRegistry registry;
    ASSERT_GE(registry.addTemperature(10.0, message), 0) << message;

    FVQual::State::StateStore state;
    state.allocate(layered, registry);

    for (int64_t cell = 0; cell < layered.cellCount(); ++cell)
        state.speciesAt(registry.temperatureIndex(), cell) = 10.0;

    FVQual::Hydro::InternalModeSolver  solver;
    FVQual::Hydro::InternalModeOptions options;
    options.manningN = 0.0;

    ASSERT_TRUE(solver.initialize(layered, state, registry, options, message)) << message;

    // Push water in at one end so the surface actually moves and the storage
    // curve is exercised, rather than sitting at rest where any formulation
    // would agree.
    FVQual::Forcing::Source inflow;
    inflow.column    = 0;
    inflow.discharge = 40.0;
    inflow.placement = FVQual::Forcing::Placement::Uniform;

    ASSERT_GE(solver.sources().add(inflow, message), 0) << message;

    double worst = 0.0;

    for (int i = 0; i < 20; ++i)
    {
        ASSERT_TRUE(solver.step(5.0, message)) << message;

        worst = std::max(worst, std::abs(solver.diagnostics().surfaceFluxResidual));
    }

    EXPECT_LT(worst, 1.0e-10) << "surface flux residual reached " << worst;

    // The surface really did move, or the gate is trivial.
    EXPECT_GT(layered.surfaceElevation()[0], 1.0e-3);
}

TEST(Hypsography, StillWaterStaysStillOnAHypsographicMesh)
{
    // Well-balancedness is the property the whole barotropic solver is built
    // around, and a storage curve is a change to the continuity update — so it
    // is exactly the sort of change that breaks it. A lake at rest over sloping
    // bathymetry must still not develop currents.
    HorizontalMesh horizontal;
    LayeredMesh    layered;

    std::string message;

    ASSERT_TRUE(horizontal.build(Testing::channelStrip(8, 800.0, 40.0), message))
        << message;

    // Sloping bed, so the bed-slope source term is genuinely active.
    for (int64_t f = 0; f < horizontal.faceCount(); ++f)
    {
        const double x = horizontal.faceCentroidX()[static_cast<size_t>(f)];

        horizontal.faceBed()[static_cast<size_t>(f)] = -30.0 + 0.02 * x;
    }

    SigmaSpec spec;
    spec.layerCount        = 5;
    spec.minLayerThickness = 0.0;

    ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

    for (int64_t c = 0; c < layered.columnCount(); ++c)
    {
        const double top = horizontal.faceArea()[static_cast<size_t>(c)];
        const double bed = horizontal.faceBed()[static_cast<size_t>(c)];

        ASSERT_TRUE(layered.setHypsography(c, {{bed, 5.0}}, {{0.2 * top, top}}, message))
            << message;
    }

    const std::vector<double> surface(static_cast<size_t>(layered.columnCount()), 0.0);
    ASSERT_TRUE(layered.updateSurface(surface, message)) << message;

    FVQual::State::SpeciesRegistry registry;
    ASSERT_GE(registry.addTemperature(10.0, message), 0) << message;

    FVQual::State::StateStore state;
    state.allocate(layered, registry);

    for (int64_t cell = 0; cell < layered.cellCount(); ++cell)
        state.speciesAt(registry.temperatureIndex(), cell) = 10.0;

    FVQual::Hydro::InternalModeSolver  solver;
    FVQual::Hydro::InternalModeOptions options;
    options.manningN = 0.0;

    ASSERT_TRUE(solver.initialize(layered, state, registry, options, message)) << message;

    for (int i = 0; i < 20; ++i)
        ASSERT_TRUE(solver.step(5.0, message)) << message;

    // The surface has not moved and nothing is flowing.
    for (int64_t c = 0; c < layered.columnCount(); ++c)
    {
        EXPECT_NEAR(layered.surfaceElevation()[static_cast<size_t>(c)], 0.0, 1.0e-12)
            << "column " << c << " surface drifted";
    }

    EXPECT_LT(solver.diagnostics().maxHorizontalVelocity, 1.0e-12)
        << "still water developed a current of "
        << solver.diagnostics().maxHorizontalVelocity << " m/s";
}

TEST(Hypsography, TheStageFollowsTheStorageCurveWhenWaterIsAdded)
{
    // The inverse in action. A known volume goes in; the stage that comes out
    // must be the one the curve implies, not the one a constant area would.
    // For a linear taper A(z) = a + s (z - bed), V = a d + (s/2) d^2, so
    // d = (-a + sqrt(a^2 + 2 s V)) / s.
    HorizontalMesh horizontal;
    LayeredMesh    layered;

    std::string message;

    ASSERT_TRUE(horizontal.build(Testing::channelStrip(1, 100.0, 50.0), message))
        << message;

    horizontal.faceBed()[0] = -20.0;

    SigmaSpec spec;
    spec.layerCount        = 4;
    spec.minLayerThickness = 0.0;

    ASSERT_TRUE(layered.build(horizontal, spec, -10.0, message)) << message;

    const double top     = horizontal.faceArea()[0];
    const double bedArea = 0.25 * top;

    ASSERT_TRUE(layered.setHypsography(0, {{-20.0, 0.0}}, {{bedArea, top}}, message))
        << message;

    const std::vector<double> start(1, -10.0);
    ASSERT_TRUE(layered.updateSurface(start, message)) << message;

    FVQual::State::SpeciesRegistry registry;
    ASSERT_GE(registry.addTemperature(10.0, message), 0) << message;

    FVQual::State::StateStore state;
    state.allocate(layered, registry);
    state.speciesAt(registry.temperatureIndex(), 0) = 10.0;

    for (int64_t cell = 0; cell < layered.cellCount(); ++cell)
        state.speciesAt(registry.temperatureIndex(), cell) = 10.0;

    FVQual::Hydro::InternalModeSolver  solver;
    FVQual::Hydro::InternalModeOptions options;
    options.manningN = 0.0;

    ASSERT_TRUE(solver.initialize(layered, state, registry, options, message)) << message;

    const double slope    = (top - bedArea) / 20.0;
    const double before   = 10.0; // depth above the bed at the start
    const double vBefore  = bedArea * before + 0.5 * slope * before * before;

    // A known inflow for a known time.
    constexpr double discharge = 50.0;
    constexpr double dt        = 5.0;
    constexpr int    steps     = 20;

    FVQual::Forcing::Source inflow;
    inflow.column    = 0;
    inflow.discharge = discharge;
    inflow.placement = FVQual::Forcing::Placement::Uniform;

    ASSERT_GE(solver.sources().add(inflow, message), 0) << message;

    for (int i = 0; i < steps; ++i)
        ASSERT_TRUE(solver.step(dt, message)) << message;

    const double vAfter = vBefore + discharge * dt * steps;

    const double depth =
        (-bedArea + std::sqrt(bedArea * bedArea + 2.0 * slope * vAfter)) / slope;

    EXPECT_NEAR(layered.surfaceElevation()[0], -20.0 + depth, 1.0e-8)
        << "stage " << layered.surfaceElevation()[0] << " against " << -20.0 + depth;

    // A constant-area answer would be visibly different, or this gate would pass
    // on the old formulation too.
    const double naive = -10.0 + discharge * dt * steps / top;

    EXPECT_GT(std::abs((-20.0 + depth) - naive), 0.05)
        << "the curve made no difference; the gate is vacuous";
}

// Independent two-box oracle: volumes 160/80 m3, contact area 60 m2,
// spacing 2 m. Conservative backward Euler diffusion has conductance 3 m3/s.
TEST(Hypsography, VerticalDiffusionUsesInterfaceAreaAndUnequalCellVolumes) {
    Column column; ASSERT_NO_FATAL_FAILURE(column.build(2,-4,0));
    std::string message;
    ASSERT_TRUE(column.applyCurve({-4,0},{20,100},message)) << message;
    FVQual::State::SpeciesRegistry registry;
    ASSERT_GE(registry.addTemperature(10,message),0);
    FVQual::State::SpeciesDefinition tracer; tracer.id="tracer"; tracer.units="kg/m3";
    ASSERT_GE(registry.add(tracer,message),0);
    FVQual::State::StateStore state; state.allocate(column.layered,registry);
    state.speciesAt(0,0)=20; state.speciesAt(0,1)=10;
    state.speciesAt(1,0)=0; state.speciesAt(1,1)=3;
    FVQual::Hydro::InternalModeOptions options; options.neutralViscosity=0;
    options.backgroundDiffusivity=.1;
    options.auditHeatTransport=true;
    FVQual::Hydro::InternalModeSolver solver;
    ASSERT_TRUE(solver.initialize(column.layered,state,registry,options,message));
    ASSERT_TRUE(solver.step(10,message)) << message;
    EXPECT_NEAR(state.speciesAt(1,0),.36,1e-12);
    EXPECT_NEAR(state.speciesAt(1,1),2.28,1e-12);
    EXPECT_NEAR(solver.totalSpeciesMass(1),240,1e-10);
    EXPECT_NEAR(solver.totalSpeciesMass(0),4000,1e-9);
    EXPECT_NEAR(solver.scalarBudgetResidual(1),0,1e-10);
    EXPECT_NEAR(solver.scalarBudgetResidual(0),0,1e-9);
    const auto &audit=solver.heatTransportAudit();
    const double capacity=1000*4186;
    // Heat conductance is -30*(T_top-T_bottom) for this ten-second solve.
    EXPECT_NEAR(audit.upwardDiffusion[1]/capacity,-192,1e-9);
    EXPECT_NEAR(audit.verticalDiffusion[0]/capacity,-192,1e-9);
    EXPECT_NEAR(audit.verticalDiffusion[1]/capacity,192,1e-9);
    EXPECT_NEAR(audit.horizontalAdvection[0],0,1e-9);
    auto checkpoint=solver.checkpoint(message); ASSERT_TRUE(checkpoint) << message;
    const auto saved=audit.verticalDiffusion;
    ASSERT_TRUE(solver.step(10,message)) << message;
    const auto replay=solver.heatTransportAudit().verticalDiffusion;
    ASSERT_TRUE(solver.restoreCheckpoint(*checkpoint,message)) << message;
    EXPECT_EQ(solver.heatTransportAudit().verticalDiffusion,saved);
    ASSERT_TRUE(solver.step(10,message)) << message;
    EXPECT_EQ(solver.heatTransportAudit().verticalDiffusion,replay);
}
TEST(Hypsography, SurfaceAndPenetratingSolarEnergyUseActualSurfaceAreaAndStorage) {
    for (bool weather : {false,true}) {
        Column column; ASSERT_NO_FATAL_FAILURE(column.build(2,-4,0));
        std::string message;
        // The wet surface is 100 m2; the polygon's prism is 5000 m2.
        ASSERT_TRUE(column.applyCurve({-4,0},{20,100},message));
        FVQual::State::SpeciesRegistry registry;
        ASSERT_GE(registry.addTemperature(10,message),0);
        FVQual::State::StateStore state; state.allocate(column.layered,registry);
        FVQual::Hydro::InternalModeOptions options;
        options.neutralViscosity=0; options.backgroundDiffusivity=0;
        options.surfaceHeatFlux=weather ? 0 : 100;
        options.auditHeatTransport=true;
        FVQual::Hydro::InternalModeSolver solver;
        ASSERT_TRUE(solver.initialize(column.layered,state,registry,options,message));
        double watts=10000;
        if (weather) {
            solver.enableMeteorology(.01); // rays reach both sloping banks and bed
            auto &met=solver.meteorology();
            met.airTemperature.setUniform(10); met.dewPoint.setUniform(10);
            met.windSpeed.setUniform(0); met.shortwave.setUniform(1000);
            met.base.dragCoefficient=0;
            watts=100*FVQual::Hydro::surfaceHeatFlux(met.at(0),10).total();
        }
        ASSERT_TRUE(solver.step(1,message)) << message;
        EXPECT_NEAR((solver.totalSpeciesMass(0)-2400)*1000*4186,watts,1e-5);
        EXPECT_NEAR(solver.scalarBudgetResidual(0)*1000*4186,0,1e-5);
        const auto &audit=solver.heatTransportAudit();
        EXPECT_NEAR(audit.surface[0]+audit.surface[1],watts,1e-5);
        for (int k=0;k<2;++k)
            EXPECT_NEAR(audit.initial[k]+audit.surface[k],
                state.speciesAt(0,k)*column.layered.cellVolume()[k]*1000*4186,1e-5);
        if (!weather) {
            EXPECT_NEAR(solver.budget().surfaceHeat[5],10000,1e-9);
            EXPECT_NEAR(state.speciesAt(0,0)-10,10000./(160*1000*4186),1e-12);
        } else {
            EXPECT_NEAR(solver.budget().surfaceHeat[0],100*1000*.94,1e-8);
        }
    }
}

TEST(Hypsography, SettlingAccountsForProjectedBanksAndTheBottom) {
    for (bool receiver : {false,true}) {
        Column column; ASSERT_NO_FATAL_FAILURE(column.build(2,-4,0));
        std::string message;
        ASSERT_TRUE(column.applyCurve({-4,0},{20,100},message));
        FVQual::State::SpeciesRegistry registry;
        ASSERT_GE(registry.addTemperature(10,message),0);
        FVQual::State::SpeciesDefinition particle; particle.id="particles"; particle.units="kg/m3";
        particle.initialValue=1;
        ASSERT_GE(registry.add(particle,message),0);
        FVQual::State::SpeciesDefinition sediment; sediment.id="sediment"; sediment.units="kg/m2";
        sediment.kind=FVQual::State::SpeciesKind::MsxBenthic;
        ASSERT_GE(registry.add(sediment,message),0);
        FVQual::State::StateStore state; state.allocate(column.layered,registry);
        FVQual::Hydro::InternalModeOptions options;
        options.neutralViscosity=0; options.backgroundDiffusivity=0;
        FVQual::Hydro::InternalModeSolver solver;
        ASSERT_TRUE(solver.initialize(column.layered,state,registry,options,message));
        ASSERT_TRUE(solver.settling().setVelocity(1,1,registry,message));
        if (receiver) ASSERT_TRUE(solver.settling().setBenthicReceiver(1,2,registry,message));
        ASSERT_TRUE(solver.step(1,message));
        // Independent two-box solution: V=(160,80), interfaces=(100,60,20).
        EXPECT_NEAR(state.speciesAt(1,0),8./13,1e-12);
        EXPECT_NEAR(state.speciesAt(1,1),76./91,1e-12);
        const double deposited=6800./91;
        EXPECT_NEAR(solver.settling().deliveredToBed(1),deposited,1e-10);
        EXPECT_NEAR(solver.totalSpeciesMass(1)+(receiver ? solver.totalSpeciesMass(2)
            : solver.settling().lostToBed(1)),240,1e-10);
        EXPECT_NEAR(solver.scalarBudgetResidual(1),0,1e-10);
        EXPECT_NEAR(solver.scalarBudgetResidual(2),0,1e-10);
        if (receiver) EXPECT_NEAR(solver.totalSpeciesMass(2),deposited,1e-10);
    }
}

TEST(Hypsography, ScopedWaterReactionUsesActualSurfaceAndBottomContact) {
    Column column; ASSERT_NO_FATAL_FAILURE(column.build(2,-4,0));
    std::string message;
    ASSERT_TRUE(column.applyCurve({-4,0},{20,100},message));
    FVQual::State::SpeciesRegistry registry;
    ASSERT_GE(registry.addTemperature(10,message),0);
    FVQual::State::SpeciesDefinition tracer; tracer.id="tracer"; tracer.units="kg/m3"; tracer.initialValue=1;
    ASSERT_GE(registry.add(tracer,message),0);
    FVQual::State::StateStore state; state.allocate(column.layered,registry);
    FVQual::Hydro::InternalModeOptions options;
    options.neutralViscosity=0; options.backgroundDiffusivity=0;
    FVQual::Hydro::InternalModeSolver solver;
    ASSERT_TRUE(solver.initialize(column.layered,state,registry,options,message));
    ASSERT_TRUE(solver.reactions().addRate("tracer","0.1",FVQual::Kinetics::RateScope::Surface,message));
    ASSERT_TRUE(solver.reactions().addRate("tracer","0 - 0.1",FVQual::Kinetics::RateScope::Bed,message));
    ASSERT_TRUE(solver.reactions().compile(registry,message));
    ASSERT_TRUE(solver.step(10,message));
    EXPECT_NEAR(state.speciesAt(1,0),1.625,1e-10); // 100 kg into 160 m3
    EXPECT_NEAR(state.speciesAt(1,1),.75,1e-10);   // 20 kg out of 80 m3
    EXPECT_NEAR(solver.totalSpeciesMass(1),320,1e-8);
    EXPECT_NEAR(solver.budget().reactions[1],80,1e-8);
}

TEST(Hypsography, BenthicInventoryRatesRetainArealUnitsAndActOnlyAtBed) {
    Column column; ASSERT_NO_FATAL_FAILURE(column.build(2,-4,0));
    std::string message;
    ASSERT_TRUE(column.applyCurve({-4,0},{20,100},message));
    FVQual::State::SpeciesRegistry registry;
    ASSERT_GE(registry.addTemperature(10,message),0);
    FVQual::State::SpeciesDefinition sediment; sediment.id="sediment"; sediment.units="kg/m2";
    sediment.kind=FVQual::State::SpeciesKind::MsxBenthic;
    ASSERT_GE(registry.add(sediment,message),0);
    FVQual::State::StateStore state; state.allocate(column.layered,registry);
    FVQual::Hydro::InternalModeOptions options;
    options.neutralViscosity=0; options.backgroundDiffusivity=0;
    FVQual::Hydro::InternalModeSolver solver;
    ASSERT_TRUE(solver.initialize(column.layered,state,registry,options,message));
    ASSERT_TRUE(solver.reactions().addRate("sediment","0.1",FVQual::Kinetics::RateScope::Bed,message));
    ASSERT_TRUE(solver.reactions().compile(registry,message));
    ASSERT_TRUE(solver.step(10,message));
    EXPECT_EQ(state.speciesAt(1,0),0);
    EXPECT_NEAR(state.speciesAt(1,1),1,1e-10);
    EXPECT_NEAR(solver.totalSpeciesMass(1),column.planArea,1e-8);
}
