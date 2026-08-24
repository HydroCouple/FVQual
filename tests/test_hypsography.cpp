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

TEST(Hypsography, TheInternalModeRefusesAHypsographicMeshForNow)
{
    // Stage 1 makes the mesh's volumes right. It does not yet make the two
    // modes agree about them: the external solver still evolves depth on a
    // constant plan area, which is the form its well-balancedness rests on.
    // Running anyway would give plausible output with a silently broken mass
    // budget, so it is refused until stage 2 lands.
    HorizontalMesh horizontal;
    LayeredMesh    layered;

    std::string message;

    ASSERT_TRUE(horizontal.build(Testing::channelStrip(3, 300.0, 40.0), message))
        << message;

    for (int64_t f = 0; f < horizontal.faceCount(); ++f)
        horizontal.faceBed()[static_cast<size_t>(f)] = -10.0;

    SigmaSpec spec;
    spec.layerCount        = 4;
    spec.minLayerThickness = 0.0;

    ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

    FVQual::State::SpeciesRegistry registry;
    ASSERT_GE(registry.addTemperature(10.0, message), 0) << message;

    FVQual::State::StateStore state;
    state.allocate(layered, registry);

    FVQual::Hydro::InternalModeSolver solver;
    FVQual::Hydro::InternalModeOptions options;

    // Without a curve it initializes, as it always has.
    ASSERT_TRUE(solver.initialize(layered, state, registry, options, message))
        << message;

    // With one, it refuses and says why.
    ASSERT_TRUE(layered.setHypsography(0, {{-10.0, 0.0}}, {{100.0, 900.0}}, message))
        << message;

    const std::vector<double> surface(static_cast<size_t>(layered.columnCount()), 0.0);
    ASSERT_TRUE(layered.updateSurface(surface, message)) << message;

    FVQual::Hydro::InternalModeSolver second;

    EXPECT_FALSE(second.initialize(layered, state, registry, options, message));
    EXPECT_NE(message.find("depth-varying plan area"), std::string::npos) << message;
}
