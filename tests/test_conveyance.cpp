// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_conveyance.cpp
 * \brief Gates for the hypsographic conveyance depth (F1.5a, task 83).
 *
 * The derivation claims to have no free parameter: the shape of the storage
 * curve alone decides the answer. That is only worth anything if the two shapes
 * with a closed form come out right, so those are the first two gates, and they
 * are checked against arithmetic done outside the code rather than against the
 * function's own output.
 */

#include "fvqual/hydro/barotropicsolver.h"
#include "fvqual/mesh/horizontalmesh.h"

#include "testmeshes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

using FVQual::Hydro::BarotropicOptions;
using FVQual::Hydro::BarotropicSolver;
using FVQual::Mesh::HorizontalMesh;

namespace {

constexpr int    kColumns = 4;
constexpr double kBed     = -10.0;
constexpr double kStage   = 0.0;
constexpr double kDepth   = kStage - kBed;

struct Basin
{
    HorizontalMesh  mesh;
    BarotropicSolver solver;
};

void makeBasin(Basin &basin)
{
    std::string message;

    ASSERT_TRUE(basin.mesh.build(FVQual::Testing::basinGrid(kColumns, 1, 800.0, 200.0),
                                 message))
        << message;

    for (int64_t f = 0; f < basin.mesh.faceCount(); ++f)
        basin.mesh.faceBed()[static_cast<size_t>(f)] = kBed;

    BarotropicOptions options;
    ASSERT_TRUE(basin.solver.initialize(basin.mesh, options, message)) << message;
}

//! The smallest area the curve carries, as a fraction of the surface area.
//! setHypsography rejects a zero area, so a mathematically triangular column is
//! not representable and the shapes below are all trapezoidal by this much.
constexpr double kBedFraction = 0.02;

/*!
 * \brief Closed form for the conveyance depth of area(f) = A[eps + (1-eps) f^p].
 *
 * With z = bed + f H and depth (H - z) = H(1 - f),
 *
 *   INTEGRAL (H-z)^{5/3} dA = A H^{5/3} [ eps + (1-eps) p B(p, 8/3) ]
 *
 * where B is the Euler beta function, so
 *
 *   h_c = H [ eps + (1-eps) p B(p, 8/3) ]^{3/5}.
 *
 * Computed here from tgamma — a different route entirely from the piecewise
 * quadrature under test, which is the point of checking against it.
 */
double expectedConveyanceDepth(double power)
{
    const double beta = std::tgamma(power) * std::tgamma(8.0 / 3.0) /
                        std::tgamma(power + 8.0 / 3.0);

    const double shaped = kBedFraction + (1.0 - kBedFraction) * power * beta;

    return kDepth * std::pow(shaped, 3.0 / 5.0);
}

//! Gives every column a curve whose area grows as (z - bed)^power.
void setShapedCurve(Basin &basin, double power)
{
    std::string message;

    for (int64_t c = 0; c < basin.mesh.faceCount(); ++c)
    {
        const size_t i    = static_cast<size_t>(c);
        const double full = basin.mesh.faceArea()[i];

        std::vector<double> elevation;
        std::vector<double> area;

        constexpr int kSamples = 400; // Fine, so the shape is resolved.

        for (int s = 0; s < kSamples; ++s)
        {
            const double f = static_cast<double>(s) / (kSamples - 1);

            elevation.push_back(kBed + f * kDepth);
            area.push_back(full * (kBedFraction +
                                   (1.0 - kBedFraction) * std::pow(f, power)));
        }

        ASSERT_TRUE(basin.solver.setHypsography(c, elevation, area, message)) << message;
    }
}

double areaOf(const Basin &basin, int64_t column)
{
    return basin.mesh.faceArea()[static_cast<size_t>(column)];
}

} // namespace

TEST(Conveyance, APrismaticColumnConveysAtItsFullDepth)
{
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(makeBasin(basin));

    // No curve at all: the fallback path. This is the one that decides whether
    // every existing model is untouched, so it is checked exactly.
    for (int64_t c = 0; c < basin.mesh.faceCount(); ++c)
        EXPECT_DOUBLE_EQ(basin.solver.conveyanceDepth(c, kStage, areaOf(basin, c), kBed),
                         kDepth);
}

TEST(Conveyance, AColumnWhoseAreaIsConstantWithDepthAlsoConveysAtFullDepth)
{
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(makeBasin(basin));

    // A curve that says "the same area at every elevation" is a prism described
    // the long way, and must give the same answer as no curve at all. Without
    // the bed-area term in the integral it would not.
    setShapedCurve(basin, 0.0);

    for (int64_t c = 0; c < basin.mesh.faceCount(); ++c)
        EXPECT_NEAR(basin.solver.conveyanceDepth(c, kStage, areaOf(basin, c), kBed),
                    kDepth, 1.0e-9);
}

TEST(Conveyance, ATriangularColumnMatchesItsClosedForm)
{
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(makeBasin(basin));

    // Area growing linearly with elevation. For the pure triangle dA/dz is
    // constant and the integral is A H^{5/3} * 3/8, giving h_c = (3/8)^{3/5} H;
    // the small bed area the curve is obliged to carry shifts that slightly,
    // and expectedConveyanceDepth() accounts for it through the beta function.
    setShapedCurve(basin, 1.0);

    const double expected = expectedConveyanceDepth(1.0);

    // Sanity on the closed form itself: it must stay near the clean triangle.
    ASSERT_NEAR(expected, std::pow(3.0 / 8.0, 3.0 / 5.0) * kDepth, 0.05 * kDepth);

    for (int64_t c = 0; c < basin.mesh.faceCount(); ++c)
        EXPECT_NEAR(basin.solver.conveyanceDepth(c, kStage, areaOf(basin, c), kBed),
                    expected, 0.02 * kDepth)
            << "column " << c;

    // And that number must sit between the storage depth and the full depth.
    // A triangular column's storage depth is half its depth; if the conveyance
    // depth came out below that, it would be claiming the shallow margins carry
    // more than their share, which is the error this whole change exists to fix.
    EXPECT_GT(expected, 0.5 * kDepth);
    EXPECT_LT(expected, kDepth);
}

TEST(Conveyance, ANarrowingColumnConveysDeeperThanItsStorageAverage)
{
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(makeBasin(basin));

    // The DeGray shape: wide at the surface, narrow at the bed.
    setShapedCurve(basin, 2.0);

    for (int64_t c = 0; c < basin.mesh.faceCount(); ++c)
    {
        const double conveyance =
            basin.solver.conveyanceDepth(c, kStage, areaOf(basin, c), kBed);

        EXPECT_NEAR(conveyance, expectedConveyanceDepth(2.0), 0.02 * kDepth)
            << "column " << c;

        // Storage depth for area ~ (z-bed)^2 is about H/3. The conveyance depth
        // must exceed it: that gap is the entire correction this change makes,
        // and if it vanished the change would be doing nothing.
        const double storageDepth = kDepth / 3.0;

        EXPECT_GT(conveyance, storageDepth) << "column " << c;
        EXPECT_LT(conveyance, kDepth) << "column " << c;
    }
}

TEST(Conveyance, ConveyanceDepthRisesWithStage)
{
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(makeBasin(basin));
    setShapedCurve(basin, 2.0);

    double previous = 0.0;

    for (int step = 1; step <= 10; ++step)
    {
        const double eta = kBed + 0.1 * static_cast<double>(step) * kDepth;
        const double conveyance =
            basin.solver.conveyanceDepth(0, eta, areaOf(basin, 0), kBed);

        ASSERT_GT(conveyance, previous) << "stage " << eta;
        ASSERT_LE(conveyance, eta - kBed + 1.0e-12) << "stage " << eta;

        previous = conveyance;
    }
}

TEST(Conveyance, AShelvedColumnFeelsMoreFrictionThanAPrismaticOne)
{
    // F1.5b. The whole suite passed unchanged when the friction term was
    // switched to the conveyance depth, which means nothing existing exercises
    // friction on a hypsographic column. This is that gate.
    //
    // Two basins alike in every respect except the storage curve, driven by the
    // same wind for the same time. The curved one conveys at h_c rather than at
    // its maximum depth, so Manning bites harder and it ends up slower.
    const auto peakSpeedAfterWind = [](bool withCurve) {
        Basin basin;
        makeBasin(basin);

        if (withCurve)
            setShapedCurve(basin, 2.0);

        std::string message;

        const std::vector<double> flat(static_cast<size_t>(kColumns), kStage);
        EXPECT_TRUE(basin.solver.setSurface(flat, message)) << message;

        const std::vector<double> stressX(static_cast<size_t>(kColumns), 1.5);
        const std::vector<double> stressY(static_cast<size_t>(kColumns), 0.0);
        EXPECT_TRUE(basin.solver.setWindStress(stressX, stressY, message)) << message;

        double peak = 0.0;

        for (int i = 0; i < 200; ++i)
        {
            if (!basin.solver.step(1.0, message))
                break;

            for (int64_t c = 0; c < kColumns; ++c)
            {
                const size_t k = static_cast<size_t>(c);
                peak = std::max(peak, std::hypot(basin.solver.velocityX()[k],
                                                 basin.solver.velocityY()[k]));
            }
        }

        return peak;
    };

    const double prismatic = peakSpeedAfterWind(false);
    const double shelved   = peakSpeedAfterWind(true);

    ASSERT_GT(prismatic, 0.0) << "the prismatic basin never moved, so there is "
                                 "nothing to compare";

    // The direction is the claim. The magnitude is not gated tightly because
    // the curve changes storage as well as conveyance, so the two runs do not
    // hold the same water and the ratio is not purely a friction ratio.
    EXPECT_LT(shelved, prismatic)
        << "prismatic " << prismatic << " m/s, shelved " << shelved << " m/s";
}

TEST(Conveyance, ADryOrInvertedColumnConveysNothing)
{
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(makeBasin(basin));
    setShapedCurve(basin, 2.0);

    EXPECT_DOUBLE_EQ(basin.solver.conveyanceDepth(0, kBed, areaOf(basin, 0), kBed), 0.0);
    EXPECT_DOUBLE_EQ(
        basin.solver.conveyanceDepth(0, kBed - 1.0, areaOf(basin, 0), kBed), 0.0);
}
