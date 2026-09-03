// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
//
// Gates for the mesh argument: a supplied mesh replaces the rectangular
// basin (its faces become the columns, its node z the bed), the same mesh
// arrives identically by UGRID file and by adoption, and the default
// rectangle still stands when nobody supplies one.

#include "fvqual/component/fvqualcomponent.h"

#include "hydrocouplesdk/spatial/polyhedralsurfaceargument.h"

#include <gtest/gtest.h>

#include <filesystem>
#include <string>
#include <vector>

using FVQual::Component::FVQualComponent;
using ComponentStatus = HydroCouple::IModelComponent::ComponentStatus;
using HydroCouple::SDK::IO::MeshDefinition;
using HydroCouple::SDK::Spatial::PolyhedralSurfaceArgument;

namespace {

//! Three quads in a row, with a bed that deepens along the run.
MeshDefinition slopingChannel()
{
    MeshDefinition mesh;
    mesh.meshName = "channel";
    // Two rows of four nodes; z drops from -5 at the head to -8 at the toe.
    for (int j = 0; j <= 1; ++j)
        for (int i = 0; i <= 3; ++i)
        {
            mesh.nodeX.push_back(100.0 * i);
            mesh.nodeY.push_back(50.0 * j);
            mesh.nodeZ.push_back(-5.0 - 1.0 * i);
        }
    mesh.faceNodeOffsets.push_back(0);
    const auto node = [](int i, int j) { return j * 4 + i; };
    for (int i = 0; i < 3; ++i)
    {
        mesh.faceNodes.push_back(node(i, 0));
        mesh.faceNodes.push_back(node(i + 1, 0));
        mesh.faceNodes.push_back(node(i + 1, 1));
        mesh.faceNodes.push_back(node(i, 1));
        mesh.faceNodeOffsets.push_back(
            static_cast<int64_t>(mesh.faceNodes.size()));
    }
    return mesh;
}

PolyhedralSurfaceArgument *meshArgumentOf(FVQualComponent &component)
{
    for (HydroCouple::IArgument *argument : component.arguments())
    {
        if (argument->id() == "mesh")
            return dynamic_cast<PolyhedralSurfaceArgument *>(argument);
    }
    return nullptr;
}

} // namespace

TEST(FVQualMeshArgument, ASuppliedMeshBecomesTheColumnsAndTheBed)
{
    FVQualComponent component("fvqual.mesharg");
    std::string message;
    // arguments() is valid in Created — the lazy-create contract.
    PolyhedralSurfaceArgument *meshArgument = meshArgumentOf(component);
    ASSERT_NE(meshArgument, nullptr);
    ASSERT_TRUE(meshArgument->setMesh(slopingChannel(), message)) << message;

    component.initialize();
    ASSERT_EQ(component.status(), ComponentStatus::Initialized);

    // Its three faces are the three columns, so the met input sizes to
    // them.
    ASSERT_EQ(component.inputs().size(), 1u);
    EXPECT_EQ(component.inputs()[0]->shape()[0], 3);

    component.prepare();
    ASSERT_NE(component.status(), ComponentStatus::Failed);

    // The bed is the mesh's own: per-face mean of node z, deepening
    // along the run — not the arguments' uniform -depth.
    const auto &horizontal = component.mesh().horizontal();
    ASSERT_EQ(horizontal.faceCount(), 3);
    EXPECT_DOUBLE_EQ(horizontal.faceBed()[0], -5.5);
    EXPECT_DOUBLE_EQ(horizontal.faceBed()[1], -6.5);
    EXPECT_DOUBLE_EQ(horizontal.faceBed()[2], -7.5);
}

TEST(FVQualMeshArgument, NobodySuppliesOneAndTheRectangleStillStands)
{
    FVQualComponent component("fvqual.defaultmesh");
    component.initialize();
    ASSERT_EQ(component.status(), ComponentStatus::Initialized);
    component.prepare();
    ASSERT_NE(component.status(), ComponentStatus::Failed);

    const auto &horizontal = component.mesh().horizontal();
    EXPECT_EQ(horizontal.faceCount(), 10); // the default 10-column basin
    EXPECT_DOUBLE_EQ(horizontal.faceBed()[0], -20.0); // uniform -depth
}

#ifdef USE_NETCDF
#include "hydrocouplesdk/io/netcdfugridwriter.h"

TEST(FVQualMeshArgument, TheSameMeshByFileAndByAdoptionAgree)
{
    const std::filesystem::path path =
        std::filesystem::path(FVQual_TEST_OUTPUT_DIR) /
        "mesh_argument_channel.nc";
    {
        HydroCouple::SDK::IO::NetCDFUGRIDWriter writer(path,
                                                       slopingChannel());
        ASSERT_EQ(writer.initialize(), 0) << writer.errorMessage();
        ASSERT_EQ(writer.finalize(), 0) << writer.errorMessage();
    }

    FVQualComponent byFile("fvqual.byfile");
    std::string message;
    ASSERT_TRUE(meshArgumentOf(byFile)->initialize(
        path.string(), HydroCouple::IArgument::ArgumentInputType::File,
        message))
        << message;
    byFile.initialize();
    ASSERT_EQ(byFile.status(), ComponentStatus::Initialized);
    byFile.prepare();
    ASSERT_NE(byFile.status(), ComponentStatus::Failed);

    FVQualComponent adopted("fvqual.adopted");
    ASSERT_TRUE(meshArgumentOf(adopted)->setMesh(slopingChannel(), message))
        << message;
    adopted.initialize();
    adopted.prepare();
    ASSERT_NE(adopted.status(), ComponentStatus::Failed);

    const auto &a = byFile.mesh().horizontal();
    const auto &b = adopted.mesh().horizontal();
    ASSERT_EQ(a.faceCount(), b.faceCount());
    for (int64_t f = 0; f < a.faceCount(); ++f)
        EXPECT_DOUBLE_EQ(a.faceBed()[static_cast<size_t>(f)],
                         b.faceBed()[static_cast<size_t>(f)]);
}
#endif
