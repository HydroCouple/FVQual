// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/mesh/layeredmesh.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <numeric>

using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;
namespace Testing = FVQual::Testing;

namespace {

/*!
 * \brief A flat-bedded quad strip at a uniform bed elevation.
 */
HorizontalMesh flatStrip(int columns, double bedElevation)
{
    auto definition = Testing::quadStrip(columns);
    definition.nodeZ.assign(definition.nodeX.size(), bedElevation);

    HorizontalMesh mesh;
    std::string    message;
    EXPECT_TRUE(mesh.build(definition, message)) << message;
    return mesh;
}

SigmaSpec uniformSpec(int layers)
{
    SigmaSpec spec;
    spec.distribution      = SigmaSpec::Distribution::Uniform;
    spec.layerCount        = layers;
    spec.minLayerThickness = 0.0;
    return spec;
}

} // namespace

TEST(LayeredMesh, DimensionsAndIndexing)
{
    const HorizontalMesh horizontal = flatStrip(3, -10.0);

    LayeredMesh mesh;
    std::string message;

    ASSERT_TRUE(mesh.build(horizontal, uniformSpec(5), 0.0, message)) << message;

    EXPECT_EQ(mesh.columnCount(), 3);
    EXPECT_EQ(mesh.layerCount(), 5);
    EXPECT_EQ(mesh.cellCount(), 15);
    EXPECT_EQ(mesh.interfaceCount(), 6);
    EXPECT_EQ(mesh.interfaceSlotCount(), 18);

    // Cells of one column must be contiguous: that is what makes the vertical
    // implicit solves cache-friendly.
    for (int64_t c = 0; c < mesh.columnCount(); ++c)
    {
        for (int k = 0; k + 1 < mesh.layerCount(); ++k)
            EXPECT_EQ(mesh.cell(c, k) + 1, mesh.cell(c, k + 1));
    }

    EXPECT_EQ(mesh.cell(0, 0), 0);
    EXPECT_EQ(mesh.cell(1, 0), 5);
    EXPECT_EQ(mesh.cell(2, 4), 14);
}

TEST(LayeredMesh, InterfacesDescendFromSurfaceToBed)
{
    const HorizontalMesh horizontal = flatStrip(2, -25.0);

    LayeredMesh mesh;
    std::string message;

    ASSERT_TRUE(mesh.build(horizontal, uniformSpec(8), 3.0, message)) << message;

    for (int64_t c = 0; c < mesh.columnCount(); ++c)
    {
        EXPECT_NEAR(mesh.z(c, 0), 3.0, 1e-12) << "column " << c << " surface";
        EXPECT_NEAR(mesh.z(c, mesh.layerCount()), -25.0, 1e-9)
            << "column " << c << " bed";

        for (int k = 0; k < mesh.layerCount(); ++k)
            EXPECT_GT(mesh.z(c, k), mesh.z(c, k + 1)) << "column " << c << " layer " << k;
    }
}

TEST(LayeredMesh, ThicknessesSumToColumnDepth)
{
    const HorizontalMesh horizontal = flatStrip(4, -12.0);

    LayeredMesh mesh;
    std::string message;

    ASSERT_TRUE(mesh.build(horizontal, uniformSpec(7), 1.5, message)) << message;

    for (int64_t c = 0; c < mesh.columnCount(); ++c)
    {
        double sum = 0.0;
        for (int k = 0; k < mesh.layerCount(); ++k)
            sum += mesh.thickness()[static_cast<size_t>(mesh.cell(c, k))];

        EXPECT_NEAR(sum, mesh.columnDepth()[static_cast<size_t>(c)], 1e-9)
            << "column " << c;
        EXPECT_NEAR(mesh.columnDepth()[static_cast<size_t>(c)], 13.5, 1e-12);
    }
}

TEST(LayeredMesh, TotalVolumeMatchesAreaTimesDepth)
{
    const HorizontalMesh horizontal = flatStrip(6, -20.0);

    LayeredMesh mesh;
    std::string message;

    ASSERT_TRUE(mesh.build(horizontal, uniformSpec(10), 0.0, message)) << message;

    // Six unit-area columns, each 20 m deep.
    EXPECT_NEAR(mesh.totalVolume(), 6.0 * 20.0, 1e-9);
    EXPECT_NEAR(mesh.totalVolume(), horizontal.totalArea() * 20.0, 1e-9);
}

TEST(LayeredMesh, RaisingTheSurfaceStretchesLayersNotCounts)
{
    const HorizontalMesh horizontal = flatStrip(3, -10.0);

    LayeredMesh mesh;
    std::string message;

    ASSERT_TRUE(mesh.build(horizontal, uniformSpec(5), 0.0, message)) << message;

    const double before = mesh.totalVolume();
    const int    layers = mesh.layerCount();

    const std::vector<double> raised(3, 2.0);
    ASSERT_TRUE(mesh.updateSurface(raised, message)) << message;

    // Terrain-following layers stretch with the free surface: the layer count
    // never changes, which is the whole point of sigma over z-layers.
    EXPECT_EQ(mesh.layerCount(), layers);
    EXPECT_EQ(mesh.cellCount(), 15);
    EXPECT_NEAR(mesh.totalVolume(), before * 12.0 / 10.0, 1e-9);
}

TEST(LayeredMesh, SurfaceBelowBedIsRejected)
{
    const HorizontalMesh horizontal = flatStrip(2, -5.0);

    LayeredMesh mesh;
    std::string message;

    EXPECT_FALSE(mesh.build(horizontal, uniformSpec(4), -6.0, message));
    EXPECT_FALSE(message.empty());
}

TEST(LayeredMesh, MinimumThicknessIsEnforcedInShallowColumns)
{
    auto definition = Testing::quadStrip(2);
    definition.nodeZ.assign(definition.nodeX.size(), -0.20);

    HorizontalMesh horizontal;
    std::string    message;
    ASSERT_TRUE(horizontal.build(definition, message)) << message;

    SigmaSpec spec;
    spec.distribution      = SigmaSpec::Distribution::Stretched;
    spec.layerCount        = 4;
    spec.theta             = 6.0;
    spec.minLayerThickness = 0.04;

    LayeredMesh mesh;
    ASSERT_TRUE(mesh.build(horizontal, spec, 0.0, message)) << message;

    for (int64_t c = 0; c < mesh.columnCount(); ++c)
    {
        for (int k = 0; k < mesh.layerCount(); ++k)
        {
            const double dz = mesh.thickness()[static_cast<size_t>(mesh.cell(c, k))];
            EXPECT_GE(dz, spec.minLayerThickness - 1e-12)
                << "column " << c << " layer " << k;
        }
    }
}

TEST(LayeredMesh, HaneyReportPassesOnFlatBathymetry)
{
    const HorizontalMesh horizontal = flatStrip(5, -30.0);

    LayeredMesh mesh;
    std::string message;

    ASSERT_TRUE(mesh.build(horizontal, uniformSpec(6), 0.0, message)) << message;

    const auto report = mesh.haneyReport(0.2);

    EXPECT_TRUE(report.passed()) << report.summary();
    EXPECT_NEAR(report.maxRatio, 0.0, 1e-12);
    EXPECT_EQ(report.exceedances, 0);
}

TEST(LayeredMesh, HaneyReportFlagsSteepBathymetry)
{
    auto definition = Testing::quadStrip(2);
    definition.nodeZ.assign(definition.nodeX.size(), -10.0);

    HorizontalMesh horizontal;
    std::string    message;
    ASSERT_TRUE(horizontal.build(definition, message)) << message;

    // A cliff: one column 2 m deep beside one 40 m deep gives
    // r = 38/42 = 0.905, far beyond the 0.2 sigma-coordinate guidance.
    horizontal.faceBed()[0] = -2.0;
    horizontal.faceBed()[1] = -40.0;

    LayeredMesh mesh;
    ASSERT_TRUE(mesh.build(horizontal, uniformSpec(6), 0.0, message)) << message;

    const auto report = mesh.haneyReport(0.2);

    EXPECT_FALSE(report.passed());
    EXPECT_EQ(report.exceedances, 1);
    EXPECT_NEAR(report.maxRatio, 38.0 / 42.0, 1e-12);
    EXPECT_NE(report.summary().find("hybrid z-sigma"), std::string::npos)
        << report.summary();
}
