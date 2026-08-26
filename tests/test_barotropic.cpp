// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_barotropic.cpp
 * \brief The Phase 2 verification gate for the external mode.
 *
 * These are not smoke tests. Each one pins a property that a shallow-water
 * solver can lose silently and expensively:
 *
 * - well-balancedness (still water over sloping bathymetry stays still),
 * - mass conservation to machine precision,
 * - the correct wave speed, measured as a seiche period against theory,
 * - the wind-stress magnitude, measured as a steady surface slope,
 * - shock capturing, against Ritter's dam-break solution,
 * - positivity at a wetting front.
 */

#include "fvqual/hydro/barotropicsolver.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <numeric>
#include <vector>

using FVQual::Hydro::BarotropicOptions;
using FVQual::Hydro::BarotropicSolver;
using FVQual::Mesh::HorizontalMesh;
namespace Testing = FVQual::Testing;

namespace {

constexpr double kGravity = 9.80665;

/*!
 * \brief Builds a channel mesh with a caller-supplied bed profile.
 * \param bedOfX Bed elevation as a function of cell-centre x.
 */
template <typename BedFunction>
HorizontalMesh channel(int columns, double length, double width, BedFunction bedOfX)
{
    HorizontalMesh mesh;
    std::string    message;

    EXPECT_TRUE(mesh.build(Testing::channelStrip(columns, length, width), message))
        << message;

    for (int64_t f = 0; f < mesh.faceCount(); ++f)
        mesh.faceBed()[static_cast<size_t>(f)] = bedOfX(mesh.faceCentroidX()[static_cast<size_t>(f)]);

    return mesh;
}

BarotropicOptions frictionlessOptions()
{
    BarotropicOptions options;
    options.gravity         = kGravity;
    options.frictionEnabled = false;
    options.cflNumber       = 0.45;
    return options;
}

double maximumSpeed(const BarotropicSolver &solver)
{
    double peak = 0.0;

    for (int64_t c = 0; c < solver.columnCount(); ++c)
    {
        peak = std::max(peak, std::hypot(solver.velocityX()[static_cast<size_t>(c)],
                                         solver.velocityY()[static_cast<size_t>(c)]));
    }

    return peak;
}

} // namespace

// ── Well-balancedness ─────────────────────────────────────────────────────

TEST(Barotropic, LakeAtRestOverIrregularBedIsExact)
{
    // A bumpy bed under flat water. Without hydrostatic reconstruction this
    // is where a shallow-water solver invents currents out of truncation
    // error; with it, the balance is exact to round-off.
    const HorizontalMesh mesh = channel(60, 600.0, 20.0, [](double x) {
        return -20.0 + 6.0 * std::sin(x / 90.0) + 2.5 * std::cos(x / 31.0);
    });

    BarotropicSolver solver;
    std::string      message;

    ASSERT_TRUE(solver.initialize(mesh, frictionlessOptions(), message)) << message;
    solver.setUniformSurface(0.0);

    const std::vector<double> initialDepth  = solver.depth();
    const double              initialVolume = solver.totalVolume();

    for (int i = 0; i < 200; ++i)
        ASSERT_TRUE(solver.step(0.05, message)) << message;

    EXPECT_LE(maximumSpeed(solver), 1.0e-13);

    for (int64_t c = 0; c < solver.columnCount(); ++c)
    {
        const size_t i = static_cast<size_t>(c);
        EXPECT_NEAR(solver.depth()[i], initialDepth[i], 1.0e-12) << "column " << c;
        EXPECT_NEAR(solver.surfaceElevation()[i], 0.0, 1.0e-12) << "column " << c;
    }

    EXPECT_NEAR(solver.totalVolume(), initialVolume, 1.0e-9);
}

TEST(Barotropic, LakeAtRestWithEmergentBedIsExact)
{
    // The same test with an island: part of the bed rises above the water
    // line, so the reconstruction has to handle a wet/dry interface at rest
    // without leaking momentum across it.
    const HorizontalMesh mesh = channel(80, 800.0, 20.0, [](double x) {
        return x > 300.0 && x < 500.0 ? 2.0 : -15.0;
    });

    BarotropicSolver solver;
    std::string      message;

    ASSERT_TRUE(solver.initialize(mesh, frictionlessOptions(), message)) << message;
    solver.setUniformSurface(0.0);

    int64_t dryColumns = 0;
    for (int64_t c = 0; c < solver.columnCount(); ++c)
        dryColumns += solver.isWet(c) ? 0 : 1;

    ASSERT_GT(dryColumns, 0) << "the fixture must actually contain dry columns";

    const std::vector<double> initialDepth = solver.depth();

    for (int i = 0; i < 200; ++i)
        ASSERT_TRUE(solver.step(0.05, message)) << message;

    EXPECT_LE(maximumSpeed(solver), 1.0e-13);

    for (int64_t c = 0; c < solver.columnCount(); ++c)
        EXPECT_NEAR(solver.depth()[static_cast<size_t>(c)],
                    initialDepth[static_cast<size_t>(c)], 1.0e-12)
            << "column " << c;
}

// ── Conservation ──────────────────────────────────────────────────────────

TEST(Barotropic, MassIsConservedInAClosedBasin)
{
    const HorizontalMesh mesh =
        channel(100, 1000.0, 20.0, [](double x) { return -10.0 + 3.0 * std::sin(x / 120.0); });

    BarotropicSolver solver;
    std::string      message;

    ASSERT_TRUE(solver.initialize(mesh, frictionlessOptions(), message)) << message;

    // A displaced surface that will slosh for the whole run.
    std::vector<double> surface(static_cast<size_t>(mesh.faceCount()));
    for (int64_t c = 0; c < mesh.faceCount(); ++c)
    {
        surface[static_cast<size_t>(c)] =
            0.30 * std::cos(3.14159265358979 * mesh.faceCentroidX()[static_cast<size_t>(c)] / 1000.0);
    }

    ASSERT_TRUE(solver.setSurface(surface, message)) << message;

    const double initialVolume = solver.totalVolume();
    ASSERT_GT(initialVolume, 0.0);

    ASSERT_TRUE(solver.advance(400.0, message)) << message;
    EXPECT_GT(solver.diagnostics().steps, 100);

    // The Phase 2 gate: relative mass closure below 1e-10.
    const double relativeError =
        std::abs(solver.totalVolume() - initialVolume) / initialVolume;

    EXPECT_LT(relativeError, 1.0e-10) << "volume drift " << relativeError;

    // Sloshing must actually have happened, or conservation is trivial.
    EXPECT_GT(maximumSpeed(solver), 1.0e-3);
}

// ── Wave speed ────────────────────────────────────────────────────────────

TEST(Barotropic, SeichePeriodMatchesLinearTheory)
{
    // Fundamental seiche in a closed rectangular basin:
    //     T = 2 L / sqrt(g H)
    // Getting this right means the pressure gradient and the wave speed are
    // both right; an error in either shows up as a detuned period.
    constexpr int    kColumns = 200;
    constexpr double kLength  = 1000.0;
    constexpr double kDepth   = 10.0;
    constexpr double kAmplitude = 0.02;

    const double expectedPeriod = 2.0 * kLength / std::sqrt(kGravity * kDepth);

    const HorizontalMesh mesh =
        channel(kColumns, kLength, 20.0, [](double) { return -kDepth; });

    BarotropicSolver solver;
    std::string      message;

    ASSERT_TRUE(solver.initialize(mesh, frictionlessOptions(), message)) << message;

    std::vector<double> surface(static_cast<size_t>(kColumns));
    for (int64_t c = 0; c < mesh.faceCount(); ++c)
    {
        const double x = mesh.faceCentroidX()[static_cast<size_t>(c)];
        surface[static_cast<size_t>(c)] =
            kAmplitude * std::cos(3.14159265358979 * x / kLength);
    }

    ASSERT_TRUE(solver.setSurface(surface, message)) << message;

    // Track the surface at the closed end and time its zero crossings.
    std::vector<double> crossings;
    double              previous = solver.surfaceElevation()[0];
    double              previousTime = solver.time();

    constexpr double kSampleInterval = 0.25;
    const double     endTime         = 3.5 * expectedPeriod;

    while (solver.time() < endTime)
    {
        ASSERT_TRUE(solver.advance(kSampleInterval, message)) << message;

        const double current = solver.surfaceElevation()[0];

        if ((previous > 0.0 && current <= 0.0) || (previous < 0.0 && current >= 0.0))
        {
            // Linear interpolation to the crossing instant.
            const double fraction = previous / (previous - current);
            crossings.push_back(previousTime + fraction * kSampleInterval);
        }

        previous     = current;
        previousTime = solver.time();
    }

    ASSERT_GE(crossings.size(), 4u) << "too few zero crossings to time the period";

    double intervalSum = 0.0;
    for (size_t i = 1; i < crossings.size(); ++i)
        intervalSum += crossings[i] - crossings[i - 1];

    const double measuredPeriod =
        2.0 * intervalSum / static_cast<double>(crossings.size() - 1);

    EXPECT_NEAR(measuredPeriod, expectedPeriod, 0.05 * expectedPeriod)
        << "measured " << measuredPeriod << " s vs analytic " << expectedPeriod << " s";
}

// ── Wind forcing ──────────────────────────────────────────────────────────

/*!
 * \brief Wind setup reached from rest, over a bed that shoals.
 *
 * The existing gate below starts the basin *at* the analytic tilt and checks
 * that wind holds it. That tests the equilibrium and says nothing about the
 * approach to it, which is where DeGray fails: under the year's strongest wind
 * a shallow column drains to its bed, and which resolutions survive reshuffles
 * when the formulation changes. An instability in the transient is invisible to
 * a test that begins at the answer.
 *
 * So: start flat, switch the wind on, and require the basin to arrive at the
 * analytic tilt and stay bounded on the way. Over a sloping bed, because
 * `d(eta)/dx = tau/(rho g H)` steepens as `H` falls, and the shallow end is
 * exactly where the real model comes apart.
 */
TEST(Barotropic, WindSetupIsReachedFromRestOverAShoalingBed)
{
    constexpr int    kColumns = 60;
    constexpr double kLength  = 6000.0;
    constexpr double kDeep    = 30.0;
    constexpr double kShallow = 4.0;
    constexpr double kStress  = 0.25; // Pa, near DeGray's peak wind.

    BarotropicOptions options;
    options.gravity         = kGravity;
    options.windStressX     = kStress;
    options.frictionEnabled = true;

    // Linear shoal from deep at x = 0 to shallow at x = L.
    const HorizontalMesh mesh =
        channel(kColumns, kLength, 200.0, [](double x) {
            const double f = x / kLength;
            return -(kDeep + (kShallow - kDeep) * f);
        });

    BarotropicSolver solver;
    std::string      message;

    ASSERT_TRUE(solver.initialize(mesh, options, message)) << message;

    // Flat surface, at rest. Everything that follows is the wind's doing.
    const std::vector<double> flat(static_cast<size_t>(kColumns), 0.0);
    ASSERT_TRUE(solver.setSurface(flat, message)) << message;

    // Long enough for the basin to fill and settle: the gravity-wave traverse
    // is ~L/sqrt(gH) ≈ 350 s, so several thousand seconds is many crossings.
    for (int i = 0; i < 40; ++i)
    {
        ASSERT_TRUE(solver.advance(300.0, message)) << message;

        // Bounded on the way, not merely at the end. A run that blows up and
        // recovers is not a run that worked.
        for (int64_t c = 0; c < solver.columnCount(); ++c)
        {
            const size_t i2 = static_cast<size_t>(c);

            const double depth = solver.surfaceElevation()[i2] -
                                 mesh.faceBed()[i2];

            ASSERT_GT(depth, 0.1)
                << "column " << c << " drained to " << depth
                << " m at sweep " << i << " — the shallow end is where DeGray fails";
        }
    }

    // The analytic setup, integrated across the shoal:
    //   d(eta)/dx = tau/(rho g H(x))  =>  Delta eta = tau L ln(Hd/Hs) /
    //                                     (rho g (Hd - Hs))
    const double expected = kStress * kLength * std::log(kDeep / kShallow) /
                            (options.waterDensity * kGravity * (kDeep - kShallow));

    const size_t last = static_cast<size_t>(solver.columnCount() - 1);

    const double measured =
        solver.surfaceElevation()[last] - solver.surfaceElevation()[0];

    EXPECT_GT(measured, 0.0) << "the surface did not tilt downwind";

    EXPECT_NEAR(measured, expected, 0.35 * expected)
        << "setup " << measured << " m against an analytic " << expected << " m";

    // And the basin is closed, so at steady state it must be near rest.
    EXPECT_LT(maximumSpeed(solver), 0.10)
        << "a closed basin under steady wind should settle, not circulate";
}

/*!
 * \brief The same shoal, but with a storage curve, under an impulsive gust.
 *
 * DeGray fails at day 96.75 and nowhere else. The met record there reads
 * 3.4, **12.1**, 2.3, 0.0 m/s at three-hour spacing — one isolated gust, a
 * 12.7-fold jump in stress and back down within nine hours.
 *
 * The gate above passes over a shoaling bed, so a sloping bed under steady wind
 * is not the problem. Two things separate it from DeGray: the storage curve, and
 * the impulse. This adds both. A reservoir that narrows with depth holds far
 * less than its surface area suggests, so the same outflux empties it faster,
 * and an impulsive stress is what discovers that.
 */
TEST(Barotropic, AnImpulsiveGustOverATaperingBasinDoesNotEmptyIt)
{
    constexpr int    kColumns = 40;
    constexpr double kLength  = 8000.0;
    constexpr double kDeep    = 25.0;
    constexpr double kShallow = 4.0;
    constexpr double kGust    = 0.25; // Pa, DeGray's peak.

    BarotropicOptions options;
    options.gravity         = kGravity;
    options.frictionEnabled = true;

    const HorizontalMesh mesh =
        channel(kColumns, kLength, 300.0, [](double x) {
            const double f = x / kLength;
            return -(kDeep + (kShallow - kDeep) * f);
        });

    BarotropicSolver solver;
    std::string      message;

    ASSERT_TRUE(solver.initialize(mesh, options, message)) << message;

    // A storage curve that tapers hard with depth: full plan area at the
    // surface falling to a tenth of it at the bed, which is the shape DeGray
    // has and a prismatic channel does not.
    for (int64_t c = 0; c < mesh.faceCount(); ++c)
    {
        const size_t i   = static_cast<size_t>(c);
        const double bed = mesh.faceBed()[i];
        const double top = 1.0;

        std::vector<double> elevation;
        std::vector<double> planArea;

        constexpr int kSamples = 12;

        for (int s = 0; s < kSamples; ++s)
        {
            const double f = static_cast<double>(s) / (kSamples - 1);
            const double z = bed + f * (top - bed);

            elevation.push_back(z);
            planArea.push_back(mesh.faceArea()[i] * (0.1 + 0.9 * f));
        }

        ASSERT_TRUE(solver.setHypsography(c, elevation, planArea, message)) << message;
    }

    const std::vector<double> flat(static_cast<size_t>(kColumns), 0.0);
    ASSERT_TRUE(solver.setSurface(flat, message)) << message;

    // Driven through setWindStress, which is the path the internal mode uses,
    // so this exercises what DeGray exercises rather than a parallel one.
    const auto sweep = [&](double stress, double seconds) {
        const std::vector<double> taux(static_cast<size_t>(kColumns), stress);
        const std::vector<double> tauy(static_cast<size_t>(kColumns), 0.0);

        ASSERT_TRUE(solver.setWindStress(taux, tauy, message)) << message;
        ASSERT_TRUE(solver.advance(seconds, message)) << message;

        for (int64_t c = 0; c < solver.columnCount(); ++c)
        {
            const size_t i = static_cast<size_t>(c);

            const double depth = solver.surfaceElevation()[i] - mesh.faceBed()[i];

            ASSERT_GT(depth, 0.1)
                << "column " << c << " drained to " << depth << " m under "
                << stress << " Pa";
        }
    };

    // The gust, at DeGray's three-hour met spacing: rise, peak, fall, calm.
    sweep(0.02, 10800.0);
    sweep(kGust, 10800.0);
    sweep(0.01, 10800.0);
    sweep(0.0, 21600.0);

    // And it must settle rather than ring indefinitely.
    EXPECT_LT(maximumSpeed(solver), 0.20)
        << "the basin was still ringing long after the gust passed";
}

TEST(Barotropic, SteadyWindSetupIsHeldAgainstTheAnalyticSlope)
{
    // Steady-state wind setup balances wind stress against the pressure
    // gradient:  d(eta)/dx = tau / (rho g H).
    constexpr int    kColumns = 100;
    constexpr double kLength  = 1000.0;
    constexpr double kDepth   = 10.0;
    constexpr double kStress  = 1.0; // Pa

    BarotropicOptions options;
    options.gravity         = kGravity;
    options.windStressX     = kStress;
    options.frictionEnabled = true;

    const double slope = kStress / (options.waterDensity * kGravity * kDepth);

    const HorizontalMesh mesh =
        channel(kColumns, kLength, 20.0, [](double) { return -kDepth; });

    // Start from the analytic tilt; a correct wind stress holds it.
    std::vector<double> tilted(static_cast<size_t>(kColumns));
    for (int64_t c = 0; c < mesh.faceCount(); ++c)
    {
        const double x = mesh.faceCentroidX()[static_cast<size_t>(c)];
        tilted[static_cast<size_t>(c)] = (x - 0.5 * kLength) * slope;
    }

    const auto measuredSlope = [&](const BarotropicSolver &solver) {
        const size_t last = static_cast<size_t>(solver.columnCount() - 1);
        const double dx   = mesh.faceCentroidX()[last] - mesh.faceCentroidX()[0];
        return (solver.surfaceElevation()[last] - solver.surfaceElevation()[0]) / dx;
    };

    BarotropicSolver withWind;
    std::string      message;

    ASSERT_TRUE(withWind.initialize(mesh, options, message)) << message;
    ASSERT_TRUE(withWind.setSurface(tilted, message)) << message;
    ASSERT_TRUE(withWind.advance(300.0, message)) << message;

    EXPECT_NEAR(measuredSlope(withWind), slope, 0.10 * slope)
        << "wind-supported slope drifted";
    EXPECT_LT(maximumSpeed(withWind), 0.05) << "the balanced state should stay near rest";

    // Control: the same tilt without wind is not an equilibrium — it must
    // collapse. Without this, the test above would pass on a solver that
    // simply never moves.
    BarotropicOptions noWind = options;
    noWind.windStressX       = 0.0;

    BarotropicSolver without;
    ASSERT_TRUE(without.initialize(mesh, noWind, message)) << message;
    ASSERT_TRUE(without.setSurface(tilted, message)) << message;
    ASSERT_TRUE(without.advance(300.0, message)) << message;

    EXPECT_LT(measuredSlope(without), 0.5 * slope)
        << "an unforced tilt should slosh away, not persist";
}

// ── Shock capturing ───────────────────────────────────────────────────────

TEST(Barotropic, RitterDamBreakConvergesToTheAnalyticSolution)
{
    // Ritter's dry-bed dam break. Rather than pin a single tolerance — which
    // says nothing about whether an error is first-order accuracy or a bug —
    // this runs a refinement study and requires the error to fall like the
    // scheme's formal order.
    constexpr double kLength   = 1000.0;
    constexpr double kDamX     = 500.0;
    constexpr double kUpstream = 1.0;
    constexpr double kEndTime  = 20.0;

    const double celerity = std::sqrt(kGravity * kUpstream);

    struct Run
    {
        double relativeError = 0.0;
        double frontX        = 0.0;
        double volumeDrift   = 0.0;
        double upstreamDepth = 0.0;
        double minimumDepth  = 0.0;
    };

    const auto solve = [&](int columns) {
        const HorizontalMesh mesh =
            channel(columns, kLength, 10.0, [](double) { return 0.0; });

        BarotropicOptions options = frictionlessOptions();
        options.dryDepth          = 1.0e-6;

        BarotropicSolver solver;
        std::string      message;

        EXPECT_TRUE(solver.initialize(mesh, options, message)) << message;

        std::vector<double> surface(static_cast<size_t>(columns));
        for (int64_t c = 0; c < mesh.faceCount(); ++c)
        {
            surface[static_cast<size_t>(c)] =
                mesh.faceCentroidX()[static_cast<size_t>(c)] < kDamX ? kUpstream : 0.0;
        }

        EXPECT_TRUE(solver.setSurface(surface, message)) << message;

        const double initialVolume = solver.totalVolume();
        EXPECT_TRUE(solver.advance(kEndTime, message)) << message;

        Run run;
        run.volumeDrift =
            std::abs(solver.totalVolume() - initialVolume) / initialVolume;
        run.upstreamDepth = solver.depth()[0];
        run.minimumDepth  = solver.depth()[0];
        run.frontX        = kDamX;

        int64_t nearest = 0;
        double  gap     = 1.0e30;

        for (int64_t c = 0; c < mesh.faceCount(); ++c)
        {
            const size_t i = static_cast<size_t>(c);
            const double x = mesh.faceCentroidX()[i];

            run.minimumDepth = std::min(run.minimumDepth, solver.depth()[i]);

            if (solver.depth()[i] > 1.0e-3)
                run.frontX = std::max(run.frontX, x);

            const double candidate = std::abs(x - kDamX);
            if (candidate < gap)
            {
                gap     = candidate;
                nearest = c;
            }
        }

        // Evaluate the analytic profile at the cell's own centre, not at the
        // nominal dam position — half a cell of offset is itself a percent
        // of the answer at these resolutions.
        const double xi = (mesh.faceCentroidX()[static_cast<size_t>(nearest)] - kDamX) /
                          kEndTime;
        const double analytic = std::pow(2.0 * celerity - xi, 2.0) / (9.0 * kGravity);

        run.relativeError =
            std::abs(solver.depth()[static_cast<size_t>(nearest)] - analytic) / analytic;

        return run;
    };

    const Run coarse = solve(200);
    const Run medium = solve(400);
    const Run fine   = solve(800);

    // First-order convergence: halving the cell size must roughly halve the
    // error. A bug in the flux or the source term shows up here as an error
    // that stalls instead of falling.
    EXPECT_LT(medium.relativeError, coarse.relativeError);
    EXPECT_LT(fine.relativeError, medium.relativeError);

    const double ratio = coarse.relativeError / fine.relativeError;
    EXPECT_GT(ratio, 2.0) << "error fell only " << ratio << "x over a 4x refinement";

    EXPECT_LT(fine.relativeError, 0.05)
        << "depth over the dam is off by " << 100.0 * fine.relativeError << "%";

    // Positivity at the wetting front, and mass closure across it.
    EXPECT_GE(fine.minimumDepth, 0.0);
    EXPECT_LT(fine.volumeDrift, 1.0e-9);

    // The undisturbed upstream reservoir must be untouched at this time.
    EXPECT_NEAR(fine.upstreamDepth, kUpstream, 1.0e-6);

    // Front position. The analytic tip has zero depth, so a finite threshold
    // necessarily measures a point behind it; the reference is the analytic
    // location of the same contour,
    //     x(h) = x0 + t (2 c0 - 3 sqrt(g h)).
    // The front must never run *ahead* of theory — that would signal an
    // unstable wave-speed estimate — but it does trail, and unlike the depth
    // field the lag does not shrink with refinement: the dry-front speed of a
    // first-order scheme converges at a reduced rate because the tip is a
    // degenerate point where the depth and the gradient both vanish. MUSCL
    // reconstruction is the fix, and lands with the transport engines in
    // Phase 5; until then this bound is documentation of a known limitation,
    // not a target.
    constexpr double kFrontThreshold = 1.0e-3;

    const double contourSpeed =
        2.0 * celerity - 3.0 * std::sqrt(kGravity * kFrontThreshold);
    const double expectedContour = kDamX + contourSpeed * kEndTime;

    EXPECT_LE(fine.frontX, expectedContour + 2.5) << "front ran ahead of theory";
    EXPECT_GE(fine.frontX, expectedContour - 25.0) << "front lagged excessively";
}

// ── Time stepping ─────────────────────────────────────────────────────────

TEST(Barotropic, StableTimeStepFollowsTheCflCondition)
{
    constexpr double kDepth = 25.0;
    constexpr double kDx    = 10.0;

    const HorizontalMesh mesh =
        channel(50, 50.0 * kDx, 20.0, [](double) { return -kDepth; });

    BarotropicOptions options = frictionlessOptions();
    options.maxTimeStep       = 1000.0;

    BarotropicSolver solver;
    std::string      message;

    ASSERT_TRUE(solver.initialize(mesh, options, message)) << message;
    solver.setUniformSurface(0.0);

    // Still water in a uniform channel: every wave speed is sqrt(gH), so the
    // cell limit is area / (perimeter * c). Wall faces are included in the
    // sum — they transmit no wave, so counting them is conservative, and a
    // conservative external mode is the right default when the internal mode
    // will subcycle it thousands of times per run.
    constexpr double kWidth = 20.0;

    const double celerity  = std::sqrt(kGravity * kDepth);
    const double area      = kDx * kWidth;
    const double perimeter = 2.0 * (kDx + kWidth);
    const double expected  = options.cflNumber * area / (perimeter * celerity);

    EXPECT_NEAR(solver.stableTimeStep(), expected, 1.0e-9 * expected);
}

TEST(Barotropic, AdvanceLandsExactlyOnTheTargetTime)
{
    // The internal mode will subcycle this solver, so the two clocks have to
    // agree exactly — a drifting external clock desynchronizes mode
    // splitting.
    const HorizontalMesh mesh =
        channel(40, 400.0, 20.0, [](double) { return -8.0; });

    BarotropicSolver solver;
    std::string      message;

    ASSERT_TRUE(solver.initialize(mesh, frictionlessOptions(), message)) << message;
    solver.setUniformSurface(0.1);

    ASSERT_TRUE(solver.advance(7.5, message)) << message;
    EXPECT_DOUBLE_EQ(solver.time(), 7.5);

    ASSERT_TRUE(solver.advance(2.5, message)) << message;
    EXPECT_DOUBLE_EQ(solver.time(), 10.0);

    ASSERT_TRUE(solver.advance(0.0, message)) << message;
    EXPECT_DOUBLE_EQ(solver.time(), 10.0);
}

TEST(Barotropic, InvalidConfigurationIsRejected)
{
    const HorizontalMesh mesh =
        channel(10, 100.0, 20.0, [](double) { return -5.0; });

    BarotropicSolver solver;
    std::string      message;

    BarotropicOptions badCfl = frictionlessOptions();
    badCfl.cflNumber         = 1.5;
    EXPECT_FALSE(solver.initialize(mesh, badCfl, message));
    EXPECT_FALSE(message.empty());

    BarotropicOptions badSteps = frictionlessOptions();
    badSteps.minTimeStep       = 10.0;
    badSteps.maxTimeStep       = 1.0;
    EXPECT_FALSE(solver.initialize(mesh, badSteps, message));

    ASSERT_TRUE(solver.initialize(mesh, frictionlessOptions(), message)) << message;
    EXPECT_FALSE(solver.step(0.0, message));
    EXPECT_FALSE(solver.step(-1.0, message));
}

// ── Friction ──────────────────────────────────────────────────────────────

TEST(Barotropic, FrictionDissipatesFlowAndStaysStable)
{
    // Note on what this does and does not test. The semi-implicit drag
    // factor 1/(1 + dt * drag) lies in (0, 1] by construction, so the
    // operator itself can never reverse or amplify a velocity — that is a
    // property of the formula, not something a system-level test can isolate,
    // because in a closed basin the walls legitimately reverse the flow all
    // on their own. What is worth testing is the system behaviour: strong
    // friction must dissipate the flow relative to a frictionless control,
    // and must stay stable while doing it.
    const HorizontalMesh mesh =
        channel(20, 200.0, 20.0, [](double) { return -0.5; });

    BarotropicOptions options;
    options.gravity         = kGravity;
    options.frictionEnabled = true;
    options.manningN        = 0.20; // deliberately extreme

    const std::vector<double> u(static_cast<size_t>(mesh.faceCount()), 2.0);
    const std::vector<double> v(static_cast<size_t>(mesh.faceCount()), 0.0);

    const auto run = [&](const BarotropicOptions &opts, BarotropicSolver &solver) {
        std::string message;
        ASSERT_TRUE(solver.initialize(mesh, opts, message)) << message;
        solver.setUniformSurface(0.0);
        ASSERT_TRUE(solver.setVelocity(u, v, message)) << message;
        ASSERT_TRUE(solver.advance(50.0, message)) << message;
    };

    BarotropicSolver damped;
    run(options, damped);

    BarotropicOptions frictionless = options;
    frictionless.frictionEnabled   = false;

    BarotropicSolver undamped;
    run(frictionless, undamped);

    const double dampedPeak   = maximumSpeed(damped);
    const double undampedPeak = maximumSpeed(undamped);

    EXPECT_LT(dampedPeak, 0.30 * undampedPeak)
        << "friction " << dampedPeak << " vs frictionless " << undampedPeak;

    // Stability: finite, non-negative depths and conserved volume throughout.
    for (int64_t c = 0; c < damped.columnCount(); ++c)
    {
        const size_t i = static_cast<size_t>(c);
        EXPECT_TRUE(std::isfinite(damped.depth()[i])) << "column " << c;
        EXPECT_GE(damped.depth()[i], 0.0) << "column " << c;
        EXPECT_TRUE(std::isfinite(damped.velocityX()[i])) << "column " << c;
    }

    const double expectedVolume = 0.5 * 200.0 * 20.0;
    EXPECT_NEAR(damped.totalVolume(), expectedVolume, 1.0e-9 * expectedVolume);
}
