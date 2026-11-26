// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_verticalcoordinate.cpp
 * \brief Gates for SigmaVerticalCoordinate against HydroCouple's
 *        IVerticalCoordinate contract (F8.2e).
 *
 * Written from the coupling partner's side: every gate reaches the coordinate
 * through the `IVerticalCoordinate` base, never through `LayeredMesh`, so
 * anything asserted here is something a partner that knows nothing about sigma
 * can actually rely on.
 */

#include "fvqual/component/sigmaverticalcoordinate.h"
#include "fvqual/mesh/horizontalmesh.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/mesh/sigmaspec.h"

#include "testmeshes.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using FVQual::Component::SigmaVerticalCoordinate;
using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;
using HydroCouple::Spatial::IVerticalCoordinate;
using HydroCouple::Spatial::VerticalCoordinateKind;

namespace {

constexpr int    kColumns = 6;
constexpr int    kLayers  = 10;
constexpr double kStage   = 2.0;

//! A basin that deepens along its length, so no two columns share a bed and an
//! implementation that ignored the cell index could not pass.
struct Basin
{
    HorizontalMesh horizontal;
    LayeredMesh    layered;

    [[nodiscard]] double bed(int64_t column) const
    {
        return horizontal.faceBed()[static_cast<size_t>(column)];
    }
};

void makeBasin(Basin &basin)
{
    std::string message;

    ASSERT_TRUE(basin.horizontal.build(
        FVQual::Testing::basinGrid(kColumns, 1, 1000.0, 100.0), message))
        << message;

    for (int64_t f = 0; f < basin.horizontal.faceCount(); ++f)
    {
        basin.horizontal.faceBed()[static_cast<size_t>(f)] =
            -5.0 - 0.5 * static_cast<double>(f);
    }

    SigmaSpec spec;
    spec.layerCount        = kLayers;
    spec.minLayerThickness = 0.0;

    ASSERT_TRUE(basin.layered.build(basin.horizontal, spec, kStage, message)) << message;
}

} // namespace

TEST(VerticalCoordinate, ReportsItselfAsSigmaAndAsTimeVarying)
{
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(makeBasin(basin));

    SigmaVerticalCoordinate  coordinate(basin.layered);
    const IVerticalCoordinate &vertical = coordinate;

    EXPECT_EQ(vertical.kind(), VerticalCoordinateKind::Sigma);
    EXPECT_EQ(vertical.layerCount(), kLayers);

    // Not a formality. A partner that caches elevations across an exchange is
    // correct for a z-level model and wrong for this one, and this flag is the
    // only thing that tells it which it is talking to.
    EXPECT_TRUE(vertical.isTimeVarying());
}

TEST(VerticalCoordinate, InterfacesRunFromTheFreeSurfaceDownToTheBed)
{
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(makeBasin(basin));

    SigmaVerticalCoordinate  coordinate(basin.layered);
    const IVerticalCoordinate &vertical = coordinate;

    for (int64_t c = 0; c < basin.layered.columnCount(); ++c)
    {
        EXPECT_NEAR(vertical.interfaceElevation(c, 0), kStage, 1.0e-12) << "column " << c;
        EXPECT_NEAR(vertical.interfaceElevation(c, kLayers), basin.bed(c), 1.0e-12)
            << "column " << c;

        double previous = vertical.interfaceElevation(c, 0);
        for (int k = 1; k <= kLayers; ++k)
        {
            const double current = vertical.interfaceElevation(c, k);
            ASSERT_LT(current, previous) << "column " << c << " interface " << k;
            previous = current;
        }
    }
}

TEST(VerticalCoordinate, ReportedElevationsFollowAFreeSurfaceThatMoves)
{
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(makeBasin(basin));

    SigmaVerticalCoordinate  coordinate(basin.layered);
    const IVerticalCoordinate &vertical = coordinate;

    constexpr int64_t kColumn = 3;
    constexpr int     kMiddle = 5;

    const double before = vertical.interfaceElevation(kColumn, kMiddle);

    const std::vector<double> raised(
        static_cast<size_t>(basin.layered.columnCount()), kStage + 1.0);

    std::string message;
    ASSERT_TRUE(basin.layered.updateSurface(raised, message)) << message;

    // The whole reason this is an interface rather than a cached array. A
    // mid-column sigma interface must rise when the stage does.
    EXPECT_GT(vertical.interfaceElevation(kColumn, kMiddle), before);
    EXPECT_NEAR(vertical.interfaceElevation(kColumn, 0), kStage + 1.0, 1.0e-12);

    // And the bed must not have moved with it.
    EXPECT_NEAR(vertical.interfaceElevation(kColumn, kLayers), basin.bed(kColumn),
                1.0e-12);
}

TEST(VerticalCoordinate, BulkAccessorAgreesWithThePerInterfaceOne)
{
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(makeBasin(basin));

    SigmaVerticalCoordinate  coordinate(basin.layered);
    const IVerticalCoordinate &vertical = coordinate;

    // Two ways to ask one question is two chances to answer differently, and a
    // remap driven by the bulk path would inherit the discrepancy in silence.
    const auto bulk = vertical.interfaceElevations();

    for (int64_t c = 0; c < basin.layered.columnCount(); ++c)
    {
        for (int k = 0; k <= kLayers; ++k)
        {
            ASSERT_DOUBLE_EQ(bulk[static_cast<size_t>(basin.layered.interfaceSlot(c, k))],
                             vertical.interfaceElevation(c, k))
                << "column " << c << " interface " << k;
        }
    }
}

TEST(VerticalCoordinate, ColumnsWithDifferentBedsReportDifferentElevations)
{
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(makeBasin(basin));

    SigmaVerticalCoordinate  coordinate(basin.layered);
    const IVerticalCoordinate &vertical = coordinate;

    // Guards against an implementation that ignores the cell index, which would
    // pass every gate above if every column happened to share a bed.
    EXPECT_LT(vertical.interfaceElevation(kColumns - 1, kLayers / 2),
              vertical.interfaceElevation(0, kLayers / 2));
}

TEST(VerticalCoordinate, AnOutOfRangeIndexIsUnmistakableRatherThanPlausible)
{
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(makeBasin(basin));

    SigmaVerticalCoordinate  coordinate(basin.layered);
    const IVerticalCoordinate &vertical = coordinate;

    // A partner indexing off the end must see NaN, not a neighbouring column's
    // elevation — which would look entirely reasonable plotted in a profile.
    EXPECT_TRUE(std::isnan(vertical.interfaceElevation(-1, 0)));
    EXPECT_TRUE(std::isnan(
        vertical.interfaceElevation(basin.layered.columnCount(), 0)));
    EXPECT_TRUE(std::isnan(vertical.interfaceElevation(0, -1)));
    EXPECT_TRUE(std::isnan(vertical.interfaceElevation(0, kLayers + 1)));


}

TEST(VerticalCoordinate, ElevationsAreTheMeshsOwnAndNotARederivationFromSigma)
{
    // The design decision this pins. It would be natural to answer
    // `interfaceElevation` by recomputing `eta + sigma * (bed - eta)` from the
    // published formula terms, and under a uniform grid that is exactly right,
    // which is why every other gate here passes either way.
    //
    // It stops being right the moment the generator does anything else. With a
    // minimum layer thickness enforced, `SigmaGridGenerator` lifts thin layers
    // to the floor and takes the surplus proportionally from the thick ones, so
    // the naive formula and the mesh the model actually integrates over part
    // company — quietly, and only in the shallow columns. This gate holds the
    // coordinate to the real mesh.
    //
    // Stretched rather than uniform layering, deliberately: under uniform sigma
    // every layer has the same thickness, so the floor either binds on all of
    // them at once — where the generator falls back to equal layers and the
    // naive formula is right again — or on none. The redistribution only has
    // something to redistribute when the layers differ. An earlier version of
    // this gate used uniform layering and measured a disagreement of 3.6e-15,
    // and it failed loudly rather than passing on a grid where the mechanism it
    // claims to test was never engaged.
    std::string message;

    HorizontalMesh horizontal;
    ASSERT_TRUE(horizontal.build(FVQual::Testing::basinGrid(kColumns, 1, 1000.0, 100.0),
                                 message))
        << message;

    // A shoaling arm: shallow enough at one end that the surface-clustered
    // layers there fall under the floor.
    for (int64_t f = 0; f < horizontal.faceCount(); ++f)
    {
        horizontal.faceBed()[static_cast<size_t>(f)] =
            -2.0 - 3.0 * static_cast<double>(f);
    }

    SigmaSpec spec;
    spec.layerCount        = kLayers;
    spec.distribution      = SigmaSpec::Distribution::Stretched;
    spec.theta             = 6.0;
    spec.b                 = 0.0;
    spec.minLayerThickness = 0.1;

    LayeredMesh layered;
    ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

    SigmaVerticalCoordinate  coordinate(layered);
    const IVerticalCoordinate &vertical = coordinate;

    const std::vector<double> &sigma = layered.sigmaInterfaces();

    double largestDisagreement = 0.0;

    for (int64_t c = 0; c < layered.columnCount(); ++c)
    {
        const double bed = horizontal.faceBed()[static_cast<size_t>(c)];
        const double eta = layered.surfaceElevation()[static_cast<size_t>(c)];

        for (int k = 0; k <= kLayers; ++k)
        {
            const double reported = vertical.interfaceElevation(c, k);

            // Agrees with the mesh, exactly.
            ASSERT_DOUBLE_EQ(reported, layered.z(c, k))
                << "column " << c << " interface " << k;

            const double naive = eta + sigma[static_cast<size_t>(k)] * (bed - eta);
            largestDisagreement =
                std::max(largestDisagreement, std::abs(reported - naive));
        }
    }

    // And the two really do differ, so the assertion above is load-bearing
    // rather than a tautology restated on a grid where nothing can go wrong.
    EXPECT_GT(largestDisagreement, 1.0e-6)
        << "the naive sigma formula matched the mesh everywhere, so this gate "
           "proves nothing as written";
}

TEST(VerticalCoordinate, LayerThicknessesSumToTheWaterColumnDepth)
{
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(makeBasin(basin));

    SigmaVerticalCoordinate  coordinate(basin.layered);
    const IVerticalCoordinate &vertical = coordinate;

    // What a conservative remap rests on: the layers must tile the column with
    // no gap and no overlap, or mass changes as it crosses the boundary.
    for (int64_t c = 0; c < basin.layered.columnCount(); ++c)
    {
        double sum = 0.0;
        for (int k = 0; k < kLayers; ++k)
        {
            sum += vertical.interfaceElevation(c, k)
                   - vertical.interfaceElevation(c, k + 1);
        }

        EXPECT_NEAR(sum, kStage - basin.bed(c), 1.0e-10) << "column " << c;
    }
}

TEST(VerticalCoordinate, EpochAndBulkViewsFollowGeometryAndRebinding)
{
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(makeBasin(basin));
    SigmaVerticalCoordinate coordinate(basin.layered);
    EXPECT_EQ(coordinate.columnCount(), basin.layered.columnCount());
    EXPECT_EQ(coordinate.interfaceElevations().size(),
              static_cast<size_t>(basin.layered.interfaceSlotCount()));
    auto epoch = coordinate.geometryEpoch();
    std::vector<double> surface(static_cast<size_t>(basin.layered.columnCount()), 1.0);
    std::string message;
    ASSERT_TRUE(basin.layered.updateSurface(surface, message)) << message;
    EXPECT_GT(coordinate.geometryEpoch(), epoch);
    EXPECT_DOUBLE_EQ(coordinate.interfaceElevations().front(), 1.0);
    epoch = coordinate.geometryEpoch();
    ASSERT_FALSE(basin.layered.updateSurface({}, message));
    EXPECT_EQ(coordinate.geometryEpoch(), epoch);
    FVQual::Mesh::LayeredMesh replacement;
    coordinate.refresh(replacement);
    EXPECT_GT(coordinate.geometryEpoch(), epoch);
    EXPECT_EQ(coordinate.columnCount(), 0);
    EXPECT_TRUE(coordinate.interfaceElevations().empty());
}
