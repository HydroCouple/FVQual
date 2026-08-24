// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_ugrid_roundtrip.cpp
 * \brief Verifies that what LayeredUgridWriter writes, UgridReader reads back
 *        — the F1.2/F1.4 acceptance check.
 *
 * Artifacts land in tests/output/ so a failed run can be opened in xugrid,
 * ncdump, or QGIS rather than vanishing with a temporary directory.
 */

#include "fvqual/io/layeredugridwriter.h"
#include "fvqual/io/ugridreader.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/state/statestore.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <filesystem>

using FVQual::IO::LayeredUgridMetadata;
using FVQual::IO::LayeredUgridWriter;
using FVQual::IO::UgridReader;
using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;
using FVQual::State::SpeciesRegistry;
using FVQual::State::StateStore;
namespace SDKIO   = HydroCouple::SDK::IO;
namespace Testing = FVQual::Testing;

namespace {

std::filesystem::path outputPath(const std::string &filename)
{
    const std::filesystem::path directory{FVQual_TEST_OUTPUT_DIR};
    std::error_code             ec;
    std::filesystem::create_directories(directory, ec);
    return directory / filename;
}

} // namespace

TEST(LayeredUgridWriter, ReportsAvailabilityConsistently)
{
#ifdef FVQual_HAVE_NETCDF
    EXPECT_TRUE(LayeredUgridWriter::available());
    EXPECT_TRUE(UgridReader::available());
#else
    EXPECT_FALSE(LayeredUgridWriter::available());
    EXPECT_FALSE(UgridReader::available());
#endif
}

TEST(LayeredUgridWriter, WithoutNetCDFTheFailureIsExplained)
{
#ifdef FVQual_HAVE_NETCDF
    GTEST_SKIP() << "build has NetCDF enabled";
#else
    LayeredUgridMetadata vertical;
    vertical.sigmaInterfaces = {0.0, 0.5, 1.0};
    vertical.bedElevation    = {-10.0};

    LayeredUgridWriter writer(outputPath("unused.nc"), Testing::unitTriangle(), vertical);

    EXPECT_NE(writer.initialize(), 0);
    EXPECT_NE(writer.errorMessage().find("NetCDF"), std::string::npos)
        << writer.errorMessage();

    SDKIO::MeshDefinition mesh;
    std::string           message;
    EXPECT_FALSE(UgridReader::read(outputPath("unused.nc"), mesh, message));
    EXPECT_NE(message.find("NetCDF"), std::string::npos) << message;
#endif
}

#ifdef FVQual_HAVE_NETCDF

TEST(LayeredUgridWriter, RejectsAMissingSurfaceField)
{
    LayeredUgridMetadata vertical;
    vertical.sigmaInterfaces = {0.0, 0.5, 1.0};
    vertical.bedElevation    = {-10.0, -10.0};

    LayeredUgridWriter writer(outputPath("no_surface.nc"), Testing::twoTriangleSquare(),
                              vertical);
    writer.addLayeredField("temperature", "degC");

    // The CF sigma formula needs eta; refusing early beats writing a file
    // whose vertical coordinate cannot be reconstructed.
    EXPECT_NE(writer.initialize(), 0);
    EXPECT_NE(writer.errorMessage().find("free-surface"), std::string::npos)
        << writer.errorMessage();
}

TEST(LayeredUgridWriter, RoundTripsMeshAndFields)
{
    // ── Build a small layered model state ─────────────────────────────────
    auto definition = Testing::quadPlusTriangle();
    definition.nodeZ.assign(definition.nodeX.size(), -15.0);

    HorizontalMesh horizontal;
    std::string    message;
    ASSERT_TRUE(horizontal.build(definition, message)) << message;

    SigmaSpec spec;
    spec.layerCount        = 4;
    spec.minLayerThickness = 0.0;

    LayeredMesh layered;
    ASSERT_TRUE(layered.build(horizontal, spec, 0.5, message)) << message;

    SpeciesRegistry registry;
    const int       temperature = registry.addTemperature(14.0, message);
    ASSERT_GE(temperature, 0) << message;

    StateStore state;
    state.allocate(layered, registry);

    // ── Write two time levels ─────────────────────────────────────────────
    const std::filesystem::path path = outputPath("layered_roundtrip.nc");

    LayeredUgridMetadata vertical;
    vertical.sigmaInterfaces = layered.sigmaInterfaces();
    vertical.bedElevation    = horizontal.faceBed();

    {
        LayeredUgridWriter writer(path, definition, vertical);
        writer.addColumnField("eta", "m", "sea_surface_height_above_geoid");
        writer.addLayeredField("temperature", "degC", "sea_water_temperature");

        ASSERT_EQ(writer.initialize(), 0) << writer.errorMessage();

        for (int step = 0; step < 2; ++step)
        {
            SDKIO::Snapshot snapshot(static_cast<double>(step));
            snapshot.addField(state.columnSlice("eta", state.eta));
            snapshot.addField(
                state.layeredSlice("temperature", state.speciesBlock(temperature)));

            ASSERT_EQ(writer.write(snapshot), 0) << writer.errorMessage();
        }

        ASSERT_EQ(writer.finalize(), 0) << writer.errorMessage();

        const auto catalog = writer.catalog();
        ASSERT_EQ(catalog.size(), 2u);
        EXPECT_EQ(catalog[1].itemId, "temperature");
        ASSERT_EQ(catalog[1].shape.size(), 3u);
        EXPECT_EQ(catalog[1].shape[0], 2);  // two time levels
        EXPECT_EQ(catalog[1].shape[2], 4);  // four layers
        EXPECT_EQ(catalog[1].format, "netcdf_ugrid");
    }

    ASSERT_TRUE(std::filesystem::exists(path));

    // ── Read the mesh back ────────────────────────────────────────────────
    SDKIO::MeshDefinition readBack;
    ASSERT_TRUE(UgridReader::read(path, readBack, message)) << message;

    EXPECT_EQ(readBack.nodeCount(), definition.nodeCount());
    EXPECT_EQ(readBack.faceCount(), definition.faceCount());
    EXPECT_EQ(readBack.faceNodes, definition.faceNodes);
    EXPECT_EQ(readBack.faceNodeOffsets, definition.faceNodeOffsets);

    for (int64_t n = 0; n < readBack.nodeCount(); ++n)
    {
        EXPECT_DOUBLE_EQ(readBack.nodeX[static_cast<size_t>(n)],
                         definition.nodeX[static_cast<size_t>(n)]);
        EXPECT_DOUBLE_EQ(readBack.nodeY[static_cast<size_t>(n)],
                         definition.nodeY[static_cast<size_t>(n)]);
    }

    // The mixed tri-quad case is the one that exercises fill-value padding:
    // the quad keeps 4 nodes and the triangle 3, not 4 with a stray entry.
    EXPECT_EQ(readBack.faceNodeOffsets[1] - readBack.faceNodeOffsets[0], 4);
    EXPECT_EQ(readBack.faceNodeOffsets[2] - readBack.faceNodeOffsets[1], 3);

    // Bathymetry must survive the trip, or the file is a results archive
    // rather than a model input and FVQual cannot re-read its own output.
    ASSERT_EQ(readBack.nodeZ.size(), definition.nodeZ.size());
    for (size_t n = 0; n < readBack.nodeZ.size(); ++n)
        EXPECT_DOUBLE_EQ(readBack.nodeZ[n], definition.nodeZ[n]);

    // A rebuilt mesh must be usable, not merely readable: the full pipeline
    // has to run again off the file alone.
    HorizontalMesh rebuilt;
    ASSERT_TRUE(rebuilt.build(readBack, message)) << message;
    EXPECT_DOUBLE_EQ(rebuilt.totalArea(), horizontal.totalArea());

    LayeredMesh relayered;
    ASSERT_TRUE(relayered.build(rebuilt, spec, 0.5, message)) << message;
    EXPECT_EQ(relayered.cellCount(), layered.cellCount());
    EXPECT_NEAR(relayered.totalVolume(), layered.totalVolume(), 1e-9);
}

TEST(LayeredUgridWriter, RejectsMisSizedSnapshotFields)
{
    auto definition = Testing::twoTriangleSquare();
    definition.nodeZ.assign(definition.nodeX.size(), -5.0);

    LayeredUgridMetadata vertical;
    vertical.sigmaInterfaces = {0.0, 0.5, 1.0};
    vertical.bedElevation    = {-5.0, -5.0};

    LayeredUgridWriter writer(outputPath("mis_sized.nc"), definition, vertical);
    writer.addColumnField("eta", "m");

    ASSERT_EQ(writer.initialize(), 0) << writer.errorMessage();

    SDKIO::FieldSlice wrong;
    wrong.itemId = "eta";
    wrong.kind   = HydroCouple::DataKind::Float64;
    wrong.shape  = {5}; // mesh has 2 faces
    wrong.data.resize(5 * sizeof(double));

    SDKIO::Snapshot snapshot(0.0);
    snapshot.addField(std::move(wrong));

    EXPECT_NE(writer.write(snapshot), 0);
    EXPECT_NE(writer.errorMessage().find("expected"), std::string::npos)
        << writer.errorMessage();

    EXPECT_EQ(writer.finalize(), 0);
}

#endif // FVQual_HAVE_NETCDF
