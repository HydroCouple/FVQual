// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_w2bathymetry.cpp
 * \brief Reading CE-QUAL-W2 bathymetry into a FVQual mesh (F8.1).
 *
 * Two kinds of gate, and the second is the one that matters.
 *
 * **Synthetic**, where the answer is arithmetic: a rectangular basin whose
 * volume and area can be written down, and every malformed file refused with a
 * line number.
 *
 * **Real**, against DeGray's own bathymetry — which ships with CE-QUAL-W2 and
 * is in this repository. A reader that parses a hand-built file and falls over
 * on a real one has verified nothing, and the reservoir supplies its own
 * reference values: with `EBOT = 66.35 m` from the control file's location card,
 * the grid holds 7.29e8 m³ over 5.09e7 m² at the initial surface of 123.80 m.
 * Those match DeGray's published surface area and storage, which is how the
 * datum interpretation was checked in the first place — the file does not
 * contain elevations, and getting them wrong is the single easiest way to
 * produce a reservoir of the wrong size.
 *
 * The mesh gate is the one with teeth: the volume FVQual's own hypsographic
 * mesh computes must equal the volume W2's cell arithmetic gives, because those
 * are two independent calculations that only agree if the storage curve was
 * sampled correctly.
 */

#include "fvqual/io/w2bathymetry.h"
#include "fvqual/mesh/horizontalmesh.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/mesh/sigmaspec.h"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>

using FVQual::IO::W2Bathymetry;
using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;

namespace {

//! Path to DeGray's bathymetry as it ships with CE-QUAL-W2 V5.
const char *kDeGray =
    "/sessions/festive-sharp-feynman/mnt/v5/Examples/DeGray/InputFiles/bth1.csv";

//! EBOT for DeGray's water body, from the control file's location card.
constexpr double kDeGrayBottom = 66.35;

[[nodiscard]] bool haveDeGray()
{
    std::ifstream probe(kDeGray);
    return probe.good();
}

/*!
 * \brief A rectangular two-segment basin whose volume is trivial arithmetic.
 *
 * Three layers of 2 m, constant 100 m width, 500 m segments, padded with a
 * zero-width segment at each end as W2 does.
 */
std::string rectangularBasin()
{
    return
        "$\n"
        ",1,2,3,4,\n"
        "DLX,500.0,500.0,500.0,500.0,\n"
        "ELWS,106.0,106.0,106.0,106.0,\n"
        "PHI0,0.0,0.0,0.0,0.0,\n"
        "FRIC,0.03,0.03,0.03,0.03,\n"
        "LAYERH,,,,,K\n"
        "2.0,0.0,100.0,100.0,0.0,1\n"
        "2.0,0.0,100.0,100.0,0.0,2\n"
        "2.0,0.0,100.0,100.0,0.0,3\n";
}

[[nodiscard]] bool parseText(W2Bathymetry &bathymetry, const std::string &text,
                             double bottom, std::string &message)
{
    std::istringstream input(text);

    return bathymetry.parse(input, "<test>", bottom, message);
}

} // namespace

// ── Synthetic ─────────────────────────────────────────────────────────────

TEST(W2Bathymetry, ARectangularBasinReadsWithTheRightDimensions)
{
    W2Bathymetry bathymetry;
    std::string  message;

    ASSERT_TRUE(parseText(bathymetry, rectangularBasin(), 100.0, message)) << message;

    EXPECT_EQ(bathymetry.segmentCount(), 4);
    EXPECT_EQ(bathymetry.layerCount(), 3);

    // The zero-width padding segments are not reservoir.
    ASSERT_EQ(bathymetry.activeSegments().size(), 2u);
    EXPECT_EQ(bathymetry.activeSegments()[0], 1);
    EXPECT_EQ(bathymetry.activeSegments()[1], 2);

    EXPECT_NEAR(bathymetry.segmentLength()[1], 500.0, 1.0e-12);
    EXPECT_NEAR(bathymetry.initialSurface()[1], 106.0, 1.0e-12);

    // Layer 0 is the top of the grid: three 2 m layers over a bed at 100 m puts
    // its bottom at 104 m and the top of the grid at 106 m.
    EXPECT_NEAR(bathymetry.layerBottom(2), 100.0, 1.0e-12);
    EXPECT_NEAR(bathymetry.layerBottom(1), 102.0, 1.0e-12);
    EXPECT_NEAR(bathymetry.layerBottom(0), 104.0, 1.0e-12);
}

TEST(W2Bathymetry, VolumeAndAreaAreTheArithmeticTheyShouldBe)
{
    W2Bathymetry bathymetry;
    std::string  message;

    ASSERT_TRUE(parseText(bathymetry, rectangularBasin(), 100.0, message)) << message;

    // Full: two segments, 100 m x 500 m x 6 m.
    EXPECT_NEAR(bathymetry.volumeAt(106.0), 2 * 100.0 * 500.0 * 6.0, 1.0e-6);

    // Drawn down to mid-grid: 3 m of water.
    EXPECT_NEAR(bathymetry.volumeAt(103.0), 2 * 100.0 * 500.0 * 3.0, 1.0e-6);

    // Empty.
    EXPECT_NEAR(bathymetry.volumeAt(100.0), 0.0, 1.0e-9);

    // Area is width x length, constant here at every stage.
    EXPECT_NEAR(bathymetry.surfaceAreaAt(103.0), 2 * 100.0 * 500.0, 1.0e-6);
}

TEST(W2Bathymetry, TheMeshReproducesTheBathymetrysOwnVolume)
{
    // The gate with teeth, in miniature: FVQual's hypsographic mesh and W2's
    // cell arithmetic are two independent calculations of the same water.
    W2Bathymetry bathymetry;
    std::string  message;

    ASSERT_TRUE(parseText(bathymetry, rectangularBasin(), 100.0, message)) << message;

    HydroCouple::SDK::IO::MeshDefinition definition;
    std::vector<double>                  bed;

    ASSERT_TRUE(bathymetry.toMeshDefinition(definition, bed, message)) << message;

    HorizontalMesh horizontal;
    ASSERT_TRUE(horizontal.build(definition, message)) << message;

    ASSERT_EQ(horizontal.faceCount(), 2);

    for (int64_t f = 0; f < horizontal.faceCount(); ++f)
        horizontal.faceBed()[static_cast<size_t>(f)] = bed[static_cast<size_t>(f)];

    // Plan geometry: 100 m wide by 500 m long.
    EXPECT_NEAR(horizontal.faceArea()[0], 100.0 * 500.0, 1.0e-6);

    SigmaSpec spec;
    spec.layerCount        = 6;
    spec.minLayerThickness = 0.0;

    LayeredMesh layered;
    ASSERT_TRUE(layered.build(horizontal, spec, 106.0, message)) << message;

    ASSERT_TRUE(bathymetry.applyHypsography(layered, message)) << message;

    const std::vector<double> surface(2, 106.0);
    ASSERT_TRUE(layered.updateSurface(surface, message)) << message;

    EXPECT_NEAR(layered.totalVolume(), bathymetry.volumeAt(106.0),
                1.0e-6 * bathymetry.volumeAt(106.0));
}

TEST(W2Bathymetry, MalformedFilesAreRefusedWithALineNumber)
{
    struct Case
    {
        const char *text;
        const char *needle;
    };

    const Case cases[] = {
        // A data row before the segment-index row.
        {"$\nDLX,500.0,\n", "before the segment-index row"},
        // Wrong number of values.
        {"$\n,1,2,\nDLX,500.0,\n", "expected 2"},
        // Not a number.
        {"$\n,1,2,\nDLX,500.0,wide,\n", "not a number"},
        // An unknown row label.
        {"$\n,1,2,\nWHAT,1.0,2.0,\n", "unknown row"},
        // A negative width.
        {"$\n,1,2,\nDLX,500.0,500.0,\nLAYERH,,,K\n2.0,-5.0,100.0,1\n", "negative"},
        // A non-positive layer thickness.
        {"$\n,1,2,\nDLX,500.0,500.0,\nLAYERH,,,K\n0.0,50.0,100.0,1\n", "positive"},
    };

    for (const Case &c : cases)
    {
        W2Bathymetry bathymetry;
        std::string  message;

        EXPECT_FALSE(parseText(bathymetry, c.text, 100.0, message))
            << "accepted malformed input: " << c.text;

        EXPECT_NE(message.find("<test>:"), std::string::npos)
            << "no line number in: " << message;

        EXPECT_NE(message.find(c.needle), std::string::npos)
            << "message did not explain the problem: " << message;
    }

    // A file with no water anywhere.
    W2Bathymetry empty;
    std::string  message;

    EXPECT_FALSE(parseText(empty,
                           "$\n,1,2,\nDLX,500.0,500.0,\nLAYERH,,,K\n2.0,0.0,0.0,1\n",
                           100.0, message));
    EXPECT_NE(message.find("zero width"), std::string::npos) << message;
}

// ── DeGray, the real thing ────────────────────────────────────────────────

TEST(W2Bathymetry, DeGrayReadsWithItsPublishedDimensions)
{
    if (!haveDeGray())
        GTEST_SKIP() << "DeGray bathymetry not available at " << kDeGray;

    W2Bathymetry bathymetry;
    std::string  message;

    ASSERT_TRUE(bathymetry.read(kDeGray, kDeGrayBottom, message)) << message;

    // IMX and KMX from the control file.
    EXPECT_EQ(bathymetry.segmentCount(), 32);
    EXPECT_EQ(bathymetry.layerCount(), 36);

    // Branch 1 runs segments 2..31 in W2's one-based numbering, so the padded
    // ends are excluded and 30 segments carry water.
    EXPECT_EQ(bathymetry.activeSegments().size(), 30u);

    EXPECT_NEAR(bathymetry.segmentLength()[1], 1000.0, 1.0e-9);
    EXPECT_NEAR(bathymetry.initialSurface()[1], 123.8, 1.0e-9);
    EXPECT_NEAR(bathymetry.layerHeight()[0], 2.0, 1.0e-9);

    // 36 two-metre layers on a bed at 66.35 m.
    EXPECT_NEAR(bathymetry.layerBottom(35), kDeGrayBottom, 1.0e-9);
    EXPECT_NEAR(bathymetry.layerBottom(0), kDeGrayBottom + 70.0, 1.0e-9);
}

TEST(W2Bathymetry, DeGrayHoldsTheWaterItIsKnownToHold)
{
    if (!haveDeGray())
        GTEST_SKIP() << "DeGray bathymetry not available";

    W2Bathymetry bathymetry;
    std::string  message;

    ASSERT_TRUE(bathymetry.read(kDeGray, kDeGrayBottom, message)) << message;

    const double volume = bathymetry.volumeAt(123.8);
    const double area   = bathymetry.surfaceAreaAt(123.8);

    // Against the reservoir's published figures — roughly 5.4e7 m² of surface
    // and of order 7e8 m³ of storage. Loose bounds on purpose: the point is to
    // catch a datum off by a layer or a grid read upside down, both of which
    // would be wrong by a large factor, not to pin a number the source does not
    // define to three digits.
    EXPECT_GT(volume, 5.0e8) << "volume " << volume;
    EXPECT_LT(volume, 1.0e9) << "volume " << volume;

    EXPECT_GT(area, 3.0e7) << "area " << area;
    EXPECT_LT(area, 7.0e7) << "area " << area;

    // Mean depth of a large impoundment: metres, not centimetres or hundreds.
    EXPECT_GT(volume / area, 5.0);
    EXPECT_LT(volume / area, 40.0);

    // Storage falls away as the pool drops, and faster than linearly — the
    // hypsography a prism cannot represent.
    const double drawn = bathymetry.volumeAt(113.8); // ten metres down

    EXPECT_LT(drawn, volume);
    EXPECT_LT(drawn / volume, 0.6) << "ten metres of drawdown left " << drawn / volume;
}

TEST(W2Bathymetry, TheDeGrayMeshReproducesTheBathymetrysVolume)
{
    // The gate this whole feature exists to pass. FVQual builds a mesh, gives
    // every column a storage curve, and integrates it; W2's own cell arithmetic
    // says what the answer should be. They are independent calculations, and
    // they agree only if the curve was sampled at layer midpoints — sampling at
    // boundaries averages adjacent widths and loses storage wherever the
    // reservoir narrows, which for DeGray is everywhere.
    if (!haveDeGray())
        GTEST_SKIP() << "DeGray bathymetry not available";

    W2Bathymetry bathymetry;
    std::string  message;

    ASSERT_TRUE(bathymetry.read(kDeGray, kDeGrayBottom, message)) << message;

    HydroCouple::SDK::IO::MeshDefinition definition;
    std::vector<double>                  bed;

    ASSERT_TRUE(bathymetry.toMeshDefinition(definition, bed, message)) << message;

    HorizontalMesh horizontal;
    ASSERT_TRUE(horizontal.build(definition, message)) << message;

    ASSERT_EQ(horizontal.faceCount(), 30);

    for (int64_t f = 0; f < horizontal.faceCount(); ++f)
        horizontal.faceBed()[static_cast<size_t>(f)] = bed[static_cast<size_t>(f)];

    // Enough layers to resolve a 57 m column; sigma layers are not W2 layers, so
    // the two grids deliberately do not line up — which is what makes the volume
    // agreement a statement about the storage curve rather than about indexing.
    SigmaSpec spec;
    spec.layerCount        = 30;
    spec.minLayerThickness = 0.0;

    LayeredMesh layered;
    ASSERT_TRUE(layered.build(horizontal, spec, 123.8, message)) << message;

    ASSERT_TRUE(bathymetry.applyHypsography(layered, message)) << message;

    const std::vector<double> surface(static_cast<size_t>(layered.columnCount()), 123.8);
    ASSERT_TRUE(layered.updateSurface(surface, message)) << message;

    const double reference = bathymetry.volumeAt(123.8);
    const double modelled  = layered.totalVolume();

    // Measured agreement is 0.03% — 7.28973e8 against 7.29223e8 — so the gate
    // is set at 0.2%, tight enough to catch a real error and loose enough to
    // survive a change to the sigma layering. A tolerance the answer beats by
    // sixty times would not be a gate at all. For scale, the prism mesh below
    // is out by 140%.
    EXPECT_NEAR(modelled, reference, 0.002 * reference)
        << "mesh " << modelled << " m³ against bathymetry " << reference << " m³";

    // And the mesh is emphatically not a prism.
    double prism = 0.0;

    for (int64_t c = 0; c < layered.columnCount(); ++c)
    {
        prism += horizontal.faceArea()[static_cast<size_t>(c)] *
                 layered.columnDepth()[static_cast<size_t>(c)];
    }

    EXPECT_GT(prism / modelled, 2.0)
        << "the storage curve made little difference; is it being applied? "
        << "prism " << prism << " against " << modelled;
}

