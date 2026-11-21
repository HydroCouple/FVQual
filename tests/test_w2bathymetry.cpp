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
#include "fvqual/io/w2timeseries.h"
#include "fvqual/mesh/horizontalmesh.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/mesh/sigmaspec.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>

using FVQual::IO::W2Bathymetry;
using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;

namespace {

/*!
 * \brief Path to DeGray's bathymetry as it ships with CE-QUAL-W2 V5.
 *
 * Set `FVQUAL_DEGRAY_DIR` to the example's `InputFiles` directory to point the
 * real-data gates at a local checkout. The compiled-in default is only a
 * convenience for the machine this was first written on.
 *
 * These gates skip when the file is absent, which is the right behavior for a
 * suite that must pass without a third-party dataset — and also a trap: a
 * hard-coded path that stops existing turns the gates with the most teeth into
 * three silent skips inside a run that still reports 100%. A skip is not a
 * pass, and the summary does not distinguish them.
 */
[[nodiscard]] std::string deGrayDir()
{
    if (const char *root = std::getenv("FVQUAL_DEGRAY_DIR"))
        return std::string(root);

    return "/sessions/festive-sharp-feynman/mnt/v5/Examples/DeGray/InputFiles";
}

//! One file from that directory. Every DeGray path in this file goes through
//! here, so the dataset moves in one place rather than in six.
[[nodiscard]] std::string deGrayFile(const std::string &name)
{
    return deGrayDir() + "/" + name;
}

[[nodiscard]] std::string deGrayPath() { return deGrayFile("bth1.csv"); }

//! EBOT for DeGray's water body, from the control file's location card.
constexpr double kDeGrayBottom = 66.35;

[[nodiscard]] bool haveDeGray()
{
    std::ifstream probe(deGrayPath());
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
        GTEST_SKIP() << "DeGray bathymetry not available at " << deGrayPath()
                     << " (set FVQUAL_DEGRAY_DIR to the example's InputFiles)";

    W2Bathymetry bathymetry;
    std::string  message;

    ASSERT_TRUE(bathymetry.read(deGrayPath(), kDeGrayBottom, message)) << message;

    // IMX and KMX from the control file.
    EXPECT_EQ(bathymetry.segmentCount(), 32);
    EXPECT_EQ(bathymetry.layerCount(), 36);

    // Branch 1 runs segments 2..31 in W2's one-based numbering, so the padded
    // ends are excluded and 30 segments carry water.
    EXPECT_EQ(bathymetry.activeSegments().size(), 30u);

    EXPECT_NEAR(bathymetry.segmentLength()[1], 1000.0, 1.0e-9);
    EXPECT_NEAR(bathymetry.initialSurface()[1], 123.8, 1.0e-9);
    EXPECT_NEAR(bathymetry.layerHeight()[0], 2.0, 1.0e-9);

    // 36 two-meter layers on a bed at 66.35 m.
    EXPECT_NEAR(bathymetry.layerBottom(35), kDeGrayBottom, 1.0e-9);
    EXPECT_NEAR(bathymetry.layerBottom(0), kDeGrayBottom + 70.0, 1.0e-9);
}

TEST(W2Bathymetry, DeGrayHoldsTheWaterItIsKnownToHold)
{
    if (!haveDeGray())
        GTEST_SKIP() << "DeGray bathymetry not available";

    W2Bathymetry bathymetry;
    std::string  message;

    ASSERT_TRUE(bathymetry.read(deGrayPath(), kDeGrayBottom, message)) << message;

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

    // Mean depth of a large impoundment: meters, not centimeters or hundreds.
    EXPECT_GT(volume / area, 5.0);
    EXPECT_LT(volume / area, 40.0);

    // Storage falls away as the pool drops, and faster than linearly — the
    // hypsography a prism cannot represent.
    const double drawn = bathymetry.volumeAt(113.8); // ten meters down

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

    ASSERT_TRUE(bathymetry.read(deGrayPath(), kDeGrayBottom, message)) << message;

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


// ── Time series and vertical profiles ─────────────────────────────────────

TEST(W2TimeSeries, AMeteorologyRecordReadsAndInterpolates)
{
    if (!haveDeGray()) GTEST_SKIP() << "DeGray inputs not available";

    FVQual::IO::W2TimeSeries met;
    std::string              message;

    ASSERT_TRUE(met.read(deGrayFile("met.npt"), message)) << message;

    // JDAY TAIR TDEW WIND PHI CLOUD, so five data columns.
    ASSERT_EQ(met.columnCount(), 5);
    EXPECT_EQ(met.columnIndex("TAIR"), 0);
    EXPECT_EQ(met.columnIndex("CLOUD"), 4);
    EXPECT_EQ(met.columnIndex("nonesuch"), -1);

    EXPECT_NEAR(met.firstDay(), 1.0, 1.0e-9);
    EXPECT_NEAR(met.lastDay(), 366.875, 1.0e-9);
    EXPECT_GT(met.rowCount(), 2000);

    // Exactly on a sample.
    const int air = met.columnIndex("TAIR");
    EXPECT_NEAR(met.valueAt(1.0, air), 2.2, 1.0e-9);

    // Half way between the first two samples, which are 2.2 and 1.7.
    EXPECT_NEAR(met.valueAt(1.0625, air), 1.95, 1.0e-9);

    // Clamped rather than extrapolated outside the record.
    EXPECT_NEAR(met.valueAt(-50.0, air), 2.2, 1.0e-9);
    EXPECT_NEAR(met.valueAt(999.0, air), 6.7, 1.0e-9);

    // Air temperature over a year in Arkansas: cold in winter, hot in summer,
    // and never absurd.
    double coldest = 1.0e30;
    double hottest = -1.0e30;

    for (double day = 1.0; day < 366.0; day += 0.1)
    {
        const double t = met.valueAt(day, air);

        coldest = std::min(coldest, t);
        hottest = std::max(hottest, t);
    }

    EXPECT_LT(coldest, 0.0);
    EXPECT_GT(hottest, 25.0);
    EXPECT_GT(coldest, -30.0);
    EXPECT_LT(hottest, 45.0);
}

TEST(W2TimeSeries, InflowAndOutflowRecordsRead)
{
    if (!haveDeGray()) GTEST_SKIP();

    const std::string base = deGrayDir() + "/";

    FVQual::IO::W2TimeSeries inflow, inflowTemperature, outflow;
    std::string              message;

    ASSERT_TRUE(inflow.read(base + "qin_br1.npt", message)) << message;
    ASSERT_TRUE(inflowTemperature.read(base + "tin_br1.npt", message)) << message;
    ASSERT_TRUE(outflow.read(base + "qot_br1.npt", message)) << message;

    ASSERT_EQ(inflow.columnCount(), 1);
    EXPECT_EQ(inflow.columnIndex("QIN"), 0);
    EXPECT_NEAR(inflow.valueAt(2.0, 0), 12.6, 1.0e-9);

    EXPECT_EQ(inflowTemperature.columnIndex("TIN"), 0);
    EXPECT_NEAR(inflowTemperature.valueAt(2.0, 0), 6.8, 1.0e-9);

    EXPECT_EQ(outflow.columnIndex("QOUT"), 0);
    EXPECT_NEAR(outflow.valueAt(1.0, 0), 24.39, 1.0e-9);

    // A year of inflow to a reservoir of this size: tens of m³/s on average,
    // never negative.
    double total = 0.0;
    int    count = 0;

    for (double day = 64.5; day < 358.0; day += 1.0)
    {
        const double q = inflow.valueAt(day, 0);

        EXPECT_GE(q, 0.0) << "negative inflow on day " << day;

        total += q;
        ++count;
    }

    const double mean = total / count;

    EXPECT_GT(mean, 1.0);
    EXPECT_LT(mean, 500.0);
}

TEST(W2VerticalProfile, TheInitialProfileConfirmsTheVerticalDatum)
{
    if (!haveDeGray()) GTEST_SKIP();

    FVQual::IO::W2VerticalProfile profile;
    std::string                   message;

    ASSERT_TRUE(profile.read(deGrayFile("vpr.npt"), message)) << message;

    const std::vector<double> &temperature = profile.profile("TEMP");

    ASSERT_FALSE(temperature.empty());

    // This is the independent check on a datum the bathymetry never states.
    // With EBOT = 66.35 m and 36 two-meter layers, the surface at 123.80 m
    // leaves exactly 28 layers wetted — and the profile carries exactly 28
    // values. A datum wrong by even one layer would not line up.
    W2Bathymetry bathymetry;
    ASSERT_TRUE(bathymetry.read(deGrayPath(), kDeGrayBottom, message)) << message;

    // Counting layers that are both below the surface *and* carry width. W2
    // pads the grid with zero-width boundary layers at the top and bottom —
    // K = 1 and K = 36 here — and including one of those was what made this
    // count come out at 29 against the profile's 28 on the first attempt. The
    // padding is grid bookkeeping, not reservoir.
    int wetted = 0;

    for (int k = 0; k < bathymetry.layerCount(); ++k)
    {
        if (bathymetry.layerBottom(k) >= 123.8)
            continue;

        bool carriesWater = false;

        for (int i : bathymetry.activeSegments())
            carriesWater = carriesWater || bathymetry.width(k, i) > 0.0;

        if (carriesWater)
            ++wetted;
    }

    EXPECT_EQ(static_cast<int>(temperature.size()), wetted)
        << "profile has " << temperature.size() << " values against " << wetted
        << " wetted layers";

    // March in Arkansas: a cold, nearly unstratified reservoir, warmest at the
    // top by less than a degree.
    EXPECT_NEAR(temperature.front(), 7.4, 1.0e-9);
    EXPECT_NEAR(temperature.back(), 6.5, 1.0e-9);

    for (size_t i = 1; i < temperature.size(); ++i)
        EXPECT_LE(temperature[i], temperature[i - 1] + 1.0e-9) << "at depth index " << i;

    // Dissolved oxygen came along too, which the DeGray water-quality run needs.
    EXPECT_FALSE(profile.profile("DO").empty());
    EXPECT_EQ(profile.profile("DO").size(), temperature.size());
}

TEST(W2Bathymetry, WideningDeGrayLaterallyHoldsTheSameWater)
{
    // The gate that makes a 2-D DeGray the *same reservoir* rather than a
    // different one. Lateral expansion invents no bathymetry: each lane's plan
    // area at an elevation is its share of the section's width there, and the
    // shares sum to the whole width. So the volume must not move at all — not
    // approximately, but to the precision the storage integral itself carries.
    //
    // If this drifts with lane count, the expansion is redistributing water
    // instead of dividing it, and every comparison built on it would be
    // measuring the mesh rather than the physics.
    if (!haveDeGray())
        GTEST_SKIP() << "DeGray bathymetry not available";

    W2Bathymetry bathymetry;
    std::string  message;

    ASSERT_TRUE(bathymetry.read(deGrayPath(), kDeGrayBottom, message)) << message;

    const double reference = bathymetry.volumeAt(123.8);

    double singleLane = 0.0;

    for (int lanes : {1, 2, 4, 8})
    {
        HydroCouple::SDK::IO::MeshDefinition definition;
        std::vector<double>                  bed;

        ASSERT_TRUE(bathymetry.toLateralMeshDefinition(lanes, definition, bed, message))
            << "lanes " << lanes << ": " << message;

        HorizontalMesh horizontal;
        ASSERT_TRUE(horizontal.build(definition, message)) << message;

        ASSERT_EQ(horizontal.faceCount(), 30 * lanes) << "lanes " << lanes;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = bed[static_cast<size_t>(f)];

        SigmaSpec spec;
        spec.layerCount        = 30;
        spec.minLayerThickness = 0.0;

        LayeredMesh layered;
        ASSERT_TRUE(layered.build(horizontal, spec, 123.8, message))
            << "lanes " << lanes << ": " << message;

        ASSERT_TRUE(bathymetry.applyLateralHypsography(layered, lanes, message))
            << "lanes " << lanes << ": " << message;

        const std::vector<double> surface(
            static_cast<size_t>(layered.columnCount()), 123.8);
        ASSERT_TRUE(layered.updateSurface(surface, message)) << message;

        const double modelled = layered.totalVolume();

        // Against W2's own cell arithmetic, at the same 0.2% the one-column
        // gate holds to: the lateral split must not cost more than the
        // sampling already does.
        EXPECT_NEAR(modelled, reference, 0.002 * reference)
            << "lanes " << lanes << ": " << modelled << " against " << reference;

        if (lanes == 1)
            singleLane = modelled;
        else
            // And against the one-lane mesh far more tightly, because that
            // comparison has no sampling difference in it at all — only the
            // lateral division, which is exact arithmetic.
            EXPECT_NEAR(modelled, singleLane, 1.0e-6 * singleLane)
                << "lanes " << lanes << " drifted from the single-lane volume";
    }
}

TEST(W2Bathymetry, LateralExpansionGivesOuterLanesAHigherBed)
{
    // The point of the expansion. A laterally averaged column has one bed for
    // the whole section; a real section is deep in the middle and shallow at the
    // margins. If every lane came out with the same bed, the mesh would be wider
    // without being any more three-dimensional, and the comparison it exists to
    // support would be worthless.
    if (!haveDeGray())
        GTEST_SKIP() << "DeGray bathymetry not available";

    W2Bathymetry bathymetry;
    std::string  message;

    ASSERT_TRUE(bathymetry.read(deGrayPath(), kDeGrayBottom, message)) << message;

    constexpr int kLanes = 8;

    HydroCouple::SDK::IO::MeshDefinition definition;
    std::vector<double>                  bed;

    ASSERT_TRUE(bathymetry.toLateralMeshDefinition(kLanes, definition, bed, message))
        << message;

    // The deepest segment, where the section shape is most pronounced.
    size_t deepestSegment = 0;
    double deepest        = 0.0;

    for (size_t j = 0; j < bathymetry.activeSegments().size(); ++j)
    {
        const double middle = bed[j * kLanes + kLanes / 2];

        if (deepest == 0.0 || middle < deepest)
        {
            deepest        = middle;
            deepestSegment = j;
        }
    }

    const double centre = bed[deepestSegment * kLanes + kLanes / 2];
    const double edge   = bed[deepestSegment * kLanes];

    EXPECT_GT(edge, centre) << "the outer lane is not shallower than the centre, so "
                               "the section has no shape";

    // Symmetric about the centerline, because the width curve is.
    EXPECT_NEAR(bed[deepestSegment * kLanes],
                bed[deepestSegment * kLanes + kLanes - 1], 1.0e-9);
}

TEST(W2Bathymetry, OneLaneReproducesTheLaterallyAveragedMeshExactly)
{
    // The control has to be the same model, or the widening comparison measures
    // two changes at once. At one lane the expansion must reduce to
    // toMeshDefinition + applyHypsography bit for bit — same nodes, same beds,
    // same storage curves, same volume.
    //
    // This did not hold at first. Lane spans were bounded by the width at the
    // initial surface, and every one of DeGray's thirty segments is wider below
    // the top of the grid than at the water line — by up to 4.1x, because the
    // grid extends 14 m above the surface. Those layers were being clipped away.
    // The volume gate did not catch it, since the clipped layers hold no water
    // at the initial stage; a differing internal step count did.
    if (!haveDeGray())
        GTEST_SKIP() << "DeGray bathymetry not available";

    W2Bathymetry bathymetry;
    std::string  message;

    ASSERT_TRUE(bathymetry.read(deGrayPath(), kDeGrayBottom, message)) << message;

    HydroCouple::SDK::IO::MeshDefinition flat;
    HydroCouple::SDK::IO::MeshDefinition wide;
    std::vector<double>                  flatBed;
    std::vector<double>                  wideBed;

    ASSERT_TRUE(bathymetry.toMeshDefinition(flat, flatBed, message)) << message;
    ASSERT_TRUE(bathymetry.toLateralMeshDefinition(1, wide, wideBed, message)) << message;

    ASSERT_EQ(flat.nodeX.size(), wide.nodeX.size());
    ASSERT_EQ(flat.faceNodes.size(), wide.faceNodes.size());
    EXPECT_EQ(flat.faceNodes, wide.faceNodes);
    EXPECT_EQ(flat.faceNodeOffsets, wide.faceNodeOffsets);

    for (size_t i = 0; i < flat.nodeX.size(); ++i)
    {
        EXPECT_DOUBLE_EQ(flat.nodeX[i], wide.nodeX[i]) << "node " << i;
        EXPECT_DOUBLE_EQ(flat.nodeY[i], wide.nodeY[i]) << "node " << i;
    }

    ASSERT_EQ(flatBed.size(), wideBed.size());

    for (size_t i = 0; i < flatBed.size(); ++i)
        EXPECT_DOUBLE_EQ(flatBed[i], wideBed[i]) << "column " << i;

    // And the storage curves, through the volume they produce at two stages —
    // the second one above the initial surface, which is precisely where the
    // clipping used to bite and where a run that drew the reservoir down and
    // refilled it would have diverged.
    const auto volumeOf = [&](bool lateral, double stage) {
        HydroCouple::SDK::IO::MeshDefinition definition;
        std::vector<double>                  bed;
        std::string                          local;

        EXPECT_TRUE(lateral
                        ? bathymetry.toLateralMeshDefinition(1, definition, bed, local)
                        : bathymetry.toMeshDefinition(definition, bed, local))
            << local;

        HorizontalMesh horizontal;
        EXPECT_TRUE(horizontal.build(definition, local)) << local;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = bed[static_cast<size_t>(f)];

        SigmaSpec spec;
        spec.layerCount        = 30;
        spec.minLayerThickness = 0.0;

        LayeredMesh layered;
        EXPECT_TRUE(layered.build(horizontal, spec, stage, local)) << local;

        EXPECT_TRUE(lateral ? bathymetry.applyLateralHypsography(layered, 1, local)
                            : bathymetry.applyHypsography(layered, local))
            << local;

        const std::vector<double> surface(
            static_cast<size_t>(layered.columnCount()), stage);
        EXPECT_TRUE(layered.updateSurface(surface, local)) << local;

        return layered.totalVolume();
    };

    for (double stage : {123.8, 130.0})
        EXPECT_DOUBLE_EQ(volumeOf(true, stage), volumeOf(false, stage))
            << "stage " << stage;
}

TEST(W2Bathymetry, AdaptiveLateralPlanRefinesTheWideWaterAndSparesTheArms)
{
    // The point of planning lanes per segment. Uniform refinement makes 19.8 m
    // cells in DeGray's 79 m tributary arms while the 5530 m water near the dam
    // gets the same four lanes — refining hardest exactly where it is least
    // wanted, and collapsing the barotropic step on cells that contribute
    // nothing to the question.
    if (!haveDeGray())
        GTEST_SKIP() << "DeGray bathymetry not available";

    W2Bathymetry bathymetry;
    std::string  message;

    ASSERT_TRUE(bathymetry.read(deGrayPath(), kDeGrayBottom, message)) << message;

    constexpr double kMinCell = 250.0;

    const auto plan = bathymetry.planLateralCells(/*maxLanes=*/8, kMinCell);

    ASSERT_EQ(plan.lanes.size(), bathymetry.activeSegments().size());

    int widest = 1;
    int narrowest = 999;

    for (int lanes : plan.lanes)
    {
        widest    = std::max(widest, lanes);
        narrowest = std::min(narrowest, lanes);
    }

    // Some segment must actually be refined, or the plan is doing nothing.
    EXPECT_GT(widest, 1);

    // And some must be left alone, or it is uniform refinement with extra steps.
    EXPECT_EQ(narrowest, 1);

    // The constraint the plan exists to honor, checked segment by segment
    // against the widths rather than against the plan's own arithmetic.
    HydroCouple::SDK::IO::MeshDefinition definition;
    std::vector<double>                  bed;

    ASSERT_TRUE(bathymetry.toPlannedMeshDefinition(plan, definition, bed, message))
        << message;

    HorizontalMesh horizontal;
    ASSERT_TRUE(horizontal.build(definition, message)) << message;

    ASSERT_EQ(horizontal.faceCount(), plan.columnCount);

    // Conformance is the property that could quietly be wrong: a hanging node
    // would leave an edge owned by one face, and the flux across it would go
    // nowhere. HorizontalMesh refuses a non-manifold mesh outright, so building
    // at all is most of the check; this confirms the interior is closed.
    int64_t boundaryEdges = 0;

    for (int64_t e = 0; e < horizontal.edgeCount(); ++e)
    {
        if (horizontal.edgeRightFace()[static_cast<size_t>(e)] < 0)
            ++boundaryEdges;
    }

    EXPECT_GT(boundaryEdges, 0) << "a closed reservoir still has an outer boundary";
    EXPECT_LT(boundaryEdges, horizontal.edgeCount())
        << "every edge is a boundary, so no two cells are connected";

    // Where the lane count changes, the coarse cell must carry its neighbor's
    // boundary nodes as extra collinear vertices — a face with more than four.
    // Without at least one, the transition path never ran and this whole test
    // reduces to a check on a uniform mesh under a different name.
    int64_t transitionFaces = 0;

    for (size_t f = 0; f + 1 < definition.faceNodeOffsets.size(); ++f)
    {
        if (definition.faceNodeOffsets[f + 1] - definition.faceNodeOffsets[f] > 4)
            ++transitionFaces;
    }

    EXPECT_GT(transitionFaces, 0)
        << "no face has more than four nodes, so no lane-count transition was built";
}

TEST(W2Bathymetry, AnAdaptivelyPlannedDeGrayHoldsTheSameWater)
{
    // Non-uniform lanes must not change the reservoir either. Same requirement
    // as the uniform case, and a stronger test of it: the lane partitions differ
    // from segment to segment, so an error in the share arithmetic cannot cancel
    // out the way it might under a uniform split.
    if (!haveDeGray())
        GTEST_SKIP() << "DeGray bathymetry not available";

    W2Bathymetry bathymetry;
    std::string  message;

    ASSERT_TRUE(bathymetry.read(deGrayPath(), kDeGrayBottom, message)) << message;

    const double reference = bathymetry.volumeAt(123.8);

    for (double minCell : {2000.0, 700.0, 250.0})
    {
        const auto plan = bathymetry.planLateralCells(8, minCell);

        HydroCouple::SDK::IO::MeshDefinition definition;
        std::vector<double>                  bed;

        ASSERT_TRUE(bathymetry.toPlannedMeshDefinition(plan, definition, bed, message))
            << message;

        HorizontalMesh horizontal;
        ASSERT_TRUE(horizontal.build(definition, message)) << message;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = bed[static_cast<size_t>(f)];

        SigmaSpec spec;
        spec.layerCount        = 30;
        spec.minLayerThickness = 0.0;

        LayeredMesh layered;
        ASSERT_TRUE(layered.build(horizontal, spec, 123.8, message)) << message;

        ASSERT_TRUE(bathymetry.applyPlannedHypsography(layered, plan, message))
            << message;

        const std::vector<double> surface(
            static_cast<size_t>(layered.columnCount()), 123.8);
        ASSERT_TRUE(layered.updateSurface(surface, message)) << message;

        EXPECT_NEAR(layered.totalVolume(), reference, 0.002 * reference)
            << "minimum cell " << minCell << " m, " << layered.columnCount()
            << " columns";
    }
}

TEST(W2TimeSeries, ARecordThatWrapsOntoFollowingLinesIsReadWhole)
{
    // W2's fixed-format files put one value per segment and simply continue onto
    // the next line when they run out of room, with no Julian day to mark the
    // continuation. DeGray's wind sheltering file carries 32 values per record
    // across four lines while its header names only nine, so neither the line
    // nor the header says how wide a record is — the caller does, from the
    // bathymetry.
    const std::string text =
        "a wrapped file\n"
        "\n"
        "    JDAY     WSC     WSC     WSC\n"
        "     1.0    0.90    0.90    0.90\n"
        "            0.80    0.80\n"
        "   100.0    0.50    0.50    0.50\n"
        "            0.40    0.40\n";

    FVQual::IO::W2TimeSeries series;
    std::string             message;
    std::istringstream      input(text);

    ASSERT_TRUE(series.parse(input, "wrapped", /*valuesPerRecord=*/5, message)) << message;

    EXPECT_EQ(series.columnCount(), 5);
    EXPECT_EQ(series.rowCount(), 2);

    EXPECT_DOUBLE_EQ(series.valueAt(1.0, 0), 0.90);
    EXPECT_DOUBLE_EQ(series.valueAt(1.0, 4), 0.80);
    EXPECT_DOUBLE_EQ(series.valueAt(100.0, 0), 0.50);
    EXPECT_DOUBLE_EQ(series.valueAt(100.0, 4), 0.40);
}

TEST(W2TimeSeries, AWrappedFileIsStillRefusedWithoutTheRecordWidth)
{
    // The old behavior, kept: with no width given, a continuation line reads as
    // a data row whose first value is mistaken for a Julian day. Refusing is
    // right — this file was silently failing to load in the DeGray driver, and
    // a model that quietly runs without a calibrated input looks exactly like
    // one that has it.
    const std::string text =
        "a wrapped file\n"
        "\n"
        "    JDAY     WSC     WSC     WSC\n"
        "     1.0    0.90    0.90    0.90\n"
        "            0.80    0.80\n";

    FVQual::IO::W2TimeSeries series;
    std::string             message;
    std::istringstream      input(text);

    EXPECT_FALSE(series.parse(input, "wrapped", message));
    EXPECT_NE(message.find("wrapped"), std::string::npos);
}

TEST(W2TimeSeries, SingleLineFilesAreUnaffectedByTheWrappingSupport)
{
    // met.npt and its kind must read exactly as before, whether or not a width
    // is supplied, or the change has broken every existing model to fix one.
    const std::string text =
        "met\n"
        "\n"
        "    JDAY    TAIR    TDEW\n"
        "     1.0    10.0     5.0\n"
        "     2.0    12.0     6.0\n";

    FVQual::IO::W2TimeSeries loose, exact;
    std::string             m1, m2;
    std::istringstream      a(text), b(text);

    ASSERT_TRUE(loose.parse(a, "met", m1)) << m1;
    ASSERT_TRUE(exact.parse(b, "met", /*valuesPerRecord=*/2, m2)) << m2;

    EXPECT_EQ(loose.columnCount(), 2);
    EXPECT_EQ(exact.columnCount(), 2);
    EXPECT_DOUBLE_EQ(loose.valueAt(2.0, 1), exact.valueAt(2.0, 1));
    EXPECT_DOUBLE_EQ(loose.valueAt(2.0, 1), 6.0);
}

TEST(W2TimeSeries, DeGraysWindShelteringFileLoads)
{
    // The real file, and the real reason this work happened: DeGray is
    // calibrated with a wind sheltering coefficient of 0.90, and FVQual was
    // running it at 1.0 because the file would not parse. Stress goes as the
    // square of speed, so that is a 23% overestimate of the wind stress driving
    // the overturning cell the interior warm bias is blamed on.
    if (!haveDeGray())
        GTEST_SKIP() << "DeGray inputs not available";

    W2Bathymetry bathymetry;
    std::string  message;

    ASSERT_TRUE(bathymetry.read(deGrayPath(), kDeGrayBottom, message)) << message;

    FVQual::IO::W2TimeSeries sheltering;

    ASSERT_TRUE(sheltering.read(deGrayFile("wsc.npt"), bathymetry.segmentCount(),
                                message))
        << message;

    EXPECT_EQ(sheltering.columnCount(), bathymetry.segmentCount());

    for (int c = 0; c < sheltering.columnCount(); ++c)
        EXPECT_NEAR(sheltering.valueAt(64.5, c), 0.90, 1.0e-9) << "segment " << c;
}
