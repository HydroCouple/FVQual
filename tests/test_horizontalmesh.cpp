// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/mesh/horizontalmesh.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::kNoFace;
namespace Testing = FVQual::Testing;

TEST(HorizontalMesh, SingleTriangleGeometry)
{
    HorizontalMesh mesh;
    std::string    message;

    ASSERT_TRUE(mesh.build(Testing::unitTriangle(), message)) << message;

    EXPECT_EQ(mesh.faceCount(), 1);
    EXPECT_EQ(mesh.nodeCount(), 3);
    EXPECT_EQ(mesh.edgeCount(), 3);

    EXPECT_DOUBLE_EQ(mesh.faceArea()[0], 0.5);
    EXPECT_DOUBLE_EQ(mesh.faceCentroidX()[0], 1.0 / 3.0);
    EXPECT_DOUBLE_EQ(mesh.faceCentroidY()[0], 1.0 / 3.0);

    for (int64_t e = 0; e < mesh.edgeCount(); ++e)
        EXPECT_TRUE(mesh.isBoundaryEdge(e)) << "edge " << e;
}

TEST(HorizontalMesh, ClockwiseWindingIsNormalized)
{
    // The same triangle wound the wrong way must produce the same positive
    // area — winding is normalized so edge normals have a consistent sign.
    auto definition       = Testing::unitTriangle();
    definition.faceNodes  = {2, 1, 0};

    HorizontalMesh mesh;
    std::string    message;

    ASSERT_TRUE(mesh.build(definition, message)) << message;
    EXPECT_DOUBLE_EQ(mesh.faceArea()[0], 0.5);
}

TEST(HorizontalMesh, SharedEdgeIsCountedOnce)
{
    HorizontalMesh mesh;
    std::string    message;

    ASSERT_TRUE(mesh.build(Testing::twoTriangleSquare(), message)) << message;

    EXPECT_EQ(mesh.faceCount(), 2);
    EXPECT_EQ(mesh.edgeCount(), 5);
    EXPECT_DOUBLE_EQ(mesh.totalArea(), 1.0);

    int64_t interior = 0;
    for (int64_t e = 0; e < mesh.edgeCount(); ++e)
    {
        if (!mesh.isBoundaryEdge(e))
        {
            ++interior;
            EXPECT_NE(mesh.edgeLeftFace()[static_cast<size_t>(e)],
                      mesh.edgeRightFace()[static_cast<size_t>(e)]);
        }
    }

    EXPECT_EQ(interior, 1);
}

TEST(HorizontalMesh, EdgeNormalsAreUnitAndPointLeftToRight)
{
    HorizontalMesh mesh;
    std::string    message;

    ASSERT_TRUE(mesh.build(Testing::twoTriangleSquare(), message)) << message;

    for (int64_t e = 0; e < mesh.edgeCount(); ++e)
    {
        const size_t i  = static_cast<size_t>(e);
        const double nx = mesh.edgeNormalX()[i];
        const double ny = mesh.edgeNormalY()[i];

        EXPECT_NEAR(std::hypot(nx, ny), 1.0, 1e-12) << "edge " << e;

        // The normal must have a positive component along the vector from the
        // left face's centroid to the edge midpoint: it points outward.
        const size_t left = static_cast<size_t>(mesh.edgeLeftFace()[i]);
        const double dx   = mesh.edgeMidX()[i] - mesh.faceCentroidX()[left];
        const double dy   = mesh.edgeMidY()[i] - mesh.faceCentroidY()[left];

        EXPECT_GT(nx * dx + ny * dy, 0.0) << "edge " << e << " normal points inward";
    }
}

TEST(HorizontalMesh, MixedTriQuadMesh)
{
    HorizontalMesh mesh;
    std::string    message;

    ASSERT_TRUE(mesh.build(Testing::quadPlusTriangle(), message)) << message;

    EXPECT_EQ(mesh.faceCount(), 2);
    EXPECT_EQ(mesh.faceNodeCount(0), 4);
    EXPECT_EQ(mesh.faceNodeCount(1), 3);

    EXPECT_DOUBLE_EQ(mesh.faceArea()[0], 2.0);
    EXPECT_DOUBLE_EQ(mesh.faceArea()[1], 0.5);
    EXPECT_DOUBLE_EQ(mesh.totalArea(), 2.5);

    // Quad has 4 edges, triangle 3, one shared: 6 unique.
    EXPECT_EQ(mesh.edgeCount(), 6);
}

TEST(HorizontalMesh, EdgeSumPerFaceMatchesNodeCount)
{
    HorizontalMesh mesh;
    std::string    message;

    ASSERT_TRUE(mesh.build(Testing::quadStrip(4), message)) << message;

    for (int64_t f = 0; f < mesh.faceCount(); ++f)
    {
        const int64_t begin = mesh.faceEdgeOffsets()[static_cast<size_t>(f)];
        const int64_t end   = mesh.faceEdgeOffsets()[static_cast<size_t>(f) + 1];

        EXPECT_EQ(end - begin, mesh.faceNodeCount(f)) << "face " << f;
    }
}

TEST(HorizontalMesh, QuadStripHasExpectedTopology)
{
    HorizontalMesh mesh;
    std::string    message;

    ASSERT_TRUE(mesh.build(Testing::quadStrip(5), message)) << message;

    EXPECT_EQ(mesh.faceCount(), 5);
    EXPECT_DOUBLE_EQ(mesh.totalArea(), 5.0);

    int64_t interior = 0;
    for (int64_t e = 0; e < mesh.edgeCount(); ++e)
        interior += mesh.isBoundaryEdge(e) ? 0 : 1;

    EXPECT_EQ(interior, 4);
    EXPECT_DOUBLE_EQ(mesh.minimumCentroidDistance(), 1.0);
}

TEST(HorizontalMesh, NonManifoldMeshIsRejected)
{
    // Three triangles hinged on the same edge 0-1.
    HydroCouple::SDK::IO::MeshDefinition definition;
    definition.nodeX           = {0.0, 1.0, 0.5, 0.5, 2.0};
    definition.nodeY           = {0.0, 0.0, 1.0, -1.0, 0.5};
    definition.faceNodeOffsets = {0, 3, 6, 9};
    definition.faceNodes       = {0, 1, 2, 1, 0, 3, 0, 1, 4};

    HorizontalMesh mesh;
    std::string    message;

    EXPECT_FALSE(mesh.build(definition, message));
    EXPECT_NE(message.find("non-manifold"), std::string::npos) << message;
}

TEST(HorizontalMesh, DegenerateFaceIsRejected)
{
    HydroCouple::SDK::IO::MeshDefinition definition;
    definition.nodeX           = {0.0, 1.0, 2.0};
    definition.nodeY           = {0.0, 0.0, 0.0}; // collinear
    definition.faceNodeOffsets = {0, 3};
    definition.faceNodes       = {0, 1, 2};

    HorizontalMesh mesh;
    std::string    message;

    EXPECT_FALSE(mesh.build(definition, message));
    EXPECT_FALSE(message.empty());
}

TEST(HorizontalMesh, NodeElevationsAverageOntoFaces)
{
    auto definition   = Testing::unitTriangle();
    definition.nodeZ  = {-10.0, -20.0, -30.0};

    HorizontalMesh mesh;
    std::string    message;

    ASSERT_TRUE(mesh.build(definition, message)) << message;
    EXPECT_DOUBLE_EQ(mesh.faceBed()[0], -20.0);
}
