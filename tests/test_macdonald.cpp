// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_macdonald.cpp
 * \brief MacDonald-type manufactured steady-flow benchmarks (F2.4).
 *
 * The Ritter dam break in Phase 2 tests an unsteady shock against a closed
 * form. It says nothing about the balance that dominates real riverine
 * hydraulics: a steady discharge, a sloping bed, and bed friction holding each
 * other in equilibrium. A model can capture a dam break perfectly and still get
 * a backwater profile wrong.
 *
 * MacDonald's construction inverts the usual problem. Rather than solving for
 * the depth on a given bed — which has no closed form once friction is present
 * — it *chooses* a smooth depth profile `h(x)` and integrates the steady
 * momentum equation to find the bed that makes it an exact solution:
 *
 *     dz_b/dx = ( q^2/(g h^3) - 1 ) h'(x) - S_f(h) .
 *
 * Any smooth `h` will do, so the profile here is manufactured by that same
 * construction rather than copied from a particular paper's constants. What
 * matters, and what is easy to get wrong, is that `S_f` must be *exactly* the
 * friction law the model applies. FVQual uses the wide-channel Manning form
 * `S_f = n^2 u|u| / h^(4/3)`, so with `u = q/h` the manufactured bed uses
 * `S_f = n^2 q^2 / h^(10/3)`. Using the hydraulic-radius form that MacDonald's
 * paper writes would manufacture a bed for a *different* model, and the
 * benchmark would measure the mismatch between two friction laws while
 * appearing to measure the solver.
 */

#include "fvqual/hydro/barotropicsolver.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

using FVQual::Hydro::BarotropicOptions;
using FVQual::Hydro::BarotropicSolver;
using FVQual::Mesh::HorizontalMesh;
namespace Testing = FVQual::Testing;

namespace {

constexpr double kGravity = 9.80665;

/*!
 * \brief The manufactured solution: a smooth bump in depth along a channel.
 */
struct Exact
{
    static constexpr double kLength  = 1000.0; //!< Channel length (m).
    static constexpr double kWidth   = 10.0;   //!< Channel width (m).
    static constexpr double kUnitQ   = 2.0;    //!< Discharge per unit width (m²/s).
    static constexpr double kManning = 0.03;
    static constexpr double kBase    = 2.0; //!< Depth away from the bump (m).
    static constexpr double kBump    = 0.4; //!< Relative bump height.
    static constexpr double kSharp   = 16.0;

    [[nodiscard]] static double depth(double x)
    {
        const double s = x / kLength - 0.5;
        return kBase * (1.0 + kBump * std::exp(-kSharp * s * s));
    }

    [[nodiscard]] static double depthSlope(double x)
    {
        const double s = x / kLength - 0.5;
        return kBase * kBump * std::exp(-kSharp * s * s) * (-2.0 * kSharp * s) / kLength;
    }

    /*!
     * \brief Friction slope, in the form the solver itself applies.
     */
    [[nodiscard]] static double frictionSlope(double h)
    {
        return kManning * kManning * kUnitQ * kUnitQ / std::pow(h, 10.0 / 3.0);
    }

    [[nodiscard]] static double bedSlope(double x)
    {
        const double h = depth(x);

        return (kUnitQ * kUnitQ / (kGravity * h * h * h) - 1.0) * depthSlope(x) -
               frictionSlope(h);
    }

    /*!
     * \brief Bed elevation, with the datum at the downstream end.
     *
     * Integrated by Simpson's rule rather than differenced, because the bed is
     * the *input* to the model and any error in it is not a solver error — it
     * silently becomes a different problem with a different exact answer.
     */
    [[nodiscard]] static double bed(double x)
    {
        constexpr int intervals = 2000;

        const double step = (kLength - x) / intervals;

        if (step <= 0.0)
            return 0.0;

        double total = bedSlope(x) + bedSlope(kLength);

        for (int i = 1; i < intervals; ++i)
        {
            const double xi = x + i * step;
            total += (i % 2 == 0 ? 2.0 : 4.0) * bedSlope(xi);
        }

        // z_b(L) = 0, so z_b(x) = -integral from x to L.
        return -total * step / 3.0;
    }

    [[nodiscard]] static double surface(double x) { return bed(x) + depth(x); }

    [[nodiscard]] static double velocity(double x) { return kUnitQ / depth(x); }
};

/*!
 * \brief The channel, its manufactured bed, and its boundaries.
 */
struct Benchmark
{
    HorizontalMesh   mesh;
    BarotropicSolver solver;

    void build(int columns)
    {
        std::string message;

        ASSERT_TRUE(mesh.build(Testing::channelStrip(columns, Exact::kLength,
                                                     Exact::kWidth),
                               message))
            << message;

        for (int64_t f = 0; f < mesh.faceCount(); ++f)
        {
            const double x = mesh.faceCentroidX()[static_cast<size_t>(f)];
            mesh.faceBed()[static_cast<size_t>(f)] = Exact::bed(x);
        }

        BarotropicOptions options;
        options.gravity         = kGravity;
        options.manningN        = Exact::kManning;
        options.cflNumber       = 0.4;
        options.frictionEnabled = true;

        ASSERT_TRUE(solver.initialize(mesh, options, message)) << message;
    }

    void openEnds()
    {
        std::string message;

        std::vector<int64_t> inlet;
        std::vector<int64_t> outlet;

        for (int64_t e = 0; e < mesh.edgeCount(); ++e)
        {
            if (mesh.edgeRightFace()[static_cast<size_t>(e)] != FVQual::Mesh::kNoFace)
                continue;

            const double x = mesh.edgeMidX()[static_cast<size_t>(e)];

            if (std::abs(x) < 1.0e-6)
                inlet.push_back(e);
            else if (std::abs(x - Exact::kLength) < 1.0e-6)
                outlet.push_back(e);
        }

        ASSERT_FALSE(inlet.empty());
        ASSERT_FALSE(outlet.empty());

        BarotropicSolver::BoundaryCondition upstream;
        upstream.id                = "inflow";
        upstream.type              = BarotropicSolver::BoundaryType::Inflow;
        upstream.specificDischarge = Exact::kUnitQ;

        ASSERT_GE(solver.addBoundary(upstream, inlet, message), 0) << message;

        BarotropicSolver::BoundaryCondition downstream;
        downstream.id    = "stage";
        downstream.type  = BarotropicSolver::BoundaryType::Stage;
        downstream.stage = Exact::surface(Exact::kLength);

        ASSERT_GE(solver.addBoundary(downstream, outlet, message), 0) << message;
    }

    void startFromExact()
    {
        std::string message;

        const size_t        faces = static_cast<size_t>(mesh.faceCount());
        std::vector<double> eta(faces, 0.0);
        std::vector<double> u(faces, 0.0);
        std::vector<double> v(faces, 0.0);

        for (int64_t f = 0; f < mesh.faceCount(); ++f)
        {
            const double x = mesh.faceCentroidX()[static_cast<size_t>(f)];

            eta[static_cast<size_t>(f)] = Exact::surface(x);
            u[static_cast<size_t>(f)]   = Exact::velocity(x);
        }

        ASSERT_TRUE(solver.setSurface(eta, message)) << message;
        ASSERT_TRUE(solver.setVelocity(u, v, message)) << message;
    }

    void startFromUniform()
    {
        std::string message;

        const size_t        faces = static_cast<size_t>(mesh.faceCount());
        std::vector<double> eta(faces, Exact::surface(Exact::kLength));

        ASSERT_TRUE(solver.setSurface(eta, message)) << message;
    }

    /*!
     * \brief Mean absolute depth error against the exact profile (m).
     *
     * The end cells are excluded: a finite-volume cell adjacent to a boundary
     * is where the ghost state is felt most directly, and including them would
     * report a boundary-treatment error as a solver error.
     */
    [[nodiscard]] double depthError() const
    {
        double total = 0.0;
        int    count = 0;

        for (int64_t f = 0; f < mesh.faceCount(); ++f)
        {
            const double x = mesh.faceCentroidX()[static_cast<size_t>(f)];

            if (x < 0.03 * Exact::kLength || x > 0.97 * Exact::kLength)
                continue;

            const double modelled = solver.depth()[static_cast<size_t>(f)];

            total += std::abs(modelled - Exact::depth(x));
            ++count;
        }

        return count > 0 ? total / count : 0.0;
    }
};

} // namespace

TEST(MacDonald, TheManufacturedBedIsConsistentWithItsOwnSolution)
{
    // Before testing the solver, test the fixture. The bed is produced by
    // integrating a slope; if that integration is wrong then everything
    // downstream measures the wrong problem, and it would do so quietly.
    // Differencing the integrated bed must return the slope it came from.
    for (double x : {100.0, 300.0, 500.0, 700.0, 900.0})
    {
        const double delta      = 0.05;
        const double numerical  = (Exact::bed(x + delta) - Exact::bed(x - delta)) /
                                 (2.0 * delta);
        const double analytical = Exact::bedSlope(x);

        EXPECT_NEAR(numerical, analytical, 1.0e-6 + 1.0e-3 * std::abs(analytical))
            << "at x = " << x;
    }

    // And the flow it describes must be subcritical throughout, or the
    // downstream stage boundary would not be the right condition to impose.
    for (double x = 0.0; x <= Exact::kLength; x += 25.0)
    {
        const double h  = Exact::depth(x);
        const double fr = Exact::velocity(x) / std::sqrt(kGravity * h);

        EXPECT_LT(fr, 0.9) << "Froude number " << fr << " at x = " << x;
    }
}

TEST(MacDonald, TheExactProfileIsAFixedPointOfTheSolver)
{
    // Started at the exact solution, a correct steady solver should leave it
    // alone. This is the sharpest of the three: it isolates the balance
    // between the bed-slope source, the friction source, and the flux
    // divergence, with no transient to hide behind.
    Benchmark benchmark;
    benchmark.build(200);
    benchmark.openEnds();
    benchmark.startFromExact();

    const double initial = benchmark.depthError();

    std::string message;
    ASSERT_TRUE(benchmark.solver.advance(2000.0, message)) << message;

    const double settled = benchmark.depthError();

    // The bound is set from what the scheme actually achieves — 0.57 mm on a
    // 2 m depth at this resolution — not from what feels safe. A percent-level
    // tolerance here would be satisfied by a solver that had lost the balance
    // entirely and was merely in the right neighbourhood.
    EXPECT_LT(settled, 1.0e-3)
        << "the exact profile drifted to a mean depth error of " << settled
        << " m (started at " << initial << " m)";
}

TEST(MacDonald, TheSolverFindsTheProfileFromAFlatStart)
{
    // The harder direction: given only the bed, the inflow, and the downstream
    // stage, recover the depth profile. A model that merely preserved the exact
    // state could still fail this.
    Benchmark benchmark;
    benchmark.build(200);
    benchmark.openEnds();
    benchmark.startFromUniform();

    std::string message;
    ASSERT_TRUE(benchmark.solver.advance(20000.0, message)) << message;

    const double error = benchmark.depthError();

    // It must reach the *same* discrete steady state the fixed-point test
    // holds, not merely somewhere close to the exact profile.
    EXPECT_LT(error, 1.0e-3)
        << "converged to a mean depth error of " << error << " m";

    // At steady state the discharge must be uniform along the channel, which
    // is the statement of continuity the profile depends on.
    double minimum = 1.0e30;
    double maximum = 0.0;

    for (int64_t f = 0; f < benchmark.mesh.faceCount(); ++f)
    {
        const double x = benchmark.mesh.faceCentroidX()[static_cast<size_t>(f)];

        if (x < 0.05 * Exact::kLength || x > 0.95 * Exact::kLength)
            continue;

        const double unit = benchmark.solver.depth()[static_cast<size_t>(f)] *
                            benchmark.solver.velocityX()[static_cast<size_t>(f)];

        minimum = std::min(minimum, unit);
        maximum = std::max(maximum, unit);
    }

    EXPECT_NEAR(minimum, Exact::kUnitQ, 0.05 * Exact::kUnitQ);
    EXPECT_NEAR(maximum, Exact::kUnitQ, 0.05 * Exact::kUnitQ);
}

TEST(MacDonald, RefinementConvergesTowardTheExactProfile)
{
    // A single tolerance cannot distinguish a first-order scheme from a bug.
    // Refinement can: the error must fall, and fall at a rate consistent with
    // the discretisation rather than stalling at some fixed floor.
    std::vector<double> errors;

    for (int columns : {50, 100, 200})
    {
        Benchmark benchmark;
        benchmark.build(columns);
        benchmark.openEnds();
        benchmark.startFromExact();

        std::string message;
        ASSERT_TRUE(benchmark.solver.advance(2000.0, message)) << message;

        errors.push_back(benchmark.depthError());
    }

    // Each doubling of resolution must roughly halve the error. That is first
    // order, which is what the barotropic scheme is: MUSCL reconstruction was
    // added to scalar transport, not to the shallow-water depth and momentum.
    // Asserting the *rate* rather than a floor is what separates a converging
    // scheme from one that has stalled at a fixed error and would pass any
    // single tolerance.
    for (size_t i = 1; i < errors.size(); ++i)
    {
        const double ratio = errors[i - 1] / errors[i];

        EXPECT_GT(ratio, 1.7) << "halving dx cut the error by only " << ratio
                              << "x (from " << errors[i - 1] << " m to " << errors[i]
                              << " m), which is short of first order";
    }

    EXPECT_LT(errors.back(), 0.35 * errors.front())
        << "quadrupling the resolution moved the error only from " << errors.front()
        << " m to " << errors.back() << " m";
}
