// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_transport.cpp
 * \brief Verification of the Eulerian MUSCL + FCT transport engine (F5.1).
 *
 * Three properties are non-negotiable for a water-quality engine, and each is
 * pinned here:
 *
 * - **Conservation.** A budget model that loses mass cannot answer the
 *   question it exists to answer.
 * - **Boundedness.** An overshoot in a concentration is a negative dissolved
 *   oxygen, which nonlinear kinetics will turn into nonsense downstream.
 * - **Accuracy.** The scheme must actually beat the first-order baseline it
 *   replaces, measured rather than asserted.
 */

#include "fvqual/transport/euleriantransport.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <numeric>
#include <vector>

using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;
using FVQual::Transport::AdvectionScheme;
using FVQual::Transport::EulerianTransport;
namespace Testing = FVQual::Testing;

namespace {

/*!
 * \brief Advection fixtures on a channel.
 *
 * Two flow fields, because the properties under test need different things.
 *
 * A *closed circulation* — rightward along the surface layer, leftward along
 * the bed, joined by vertical flow at the two end columns — is exactly
 * divergence-free. That matters: uniform-field preservation, boundedness, and
 * conservation are only meaningful when the discrete volume identity holds,
 * and with static cell volumes that requires zero net flux in every cell. A
 * wall-bounded channel with a uniform interior flux does *not* satisfy this —
 * the first cell drains and the last accumulates — so a uniform field is
 * legitimately destroyed there, and a test asserting otherwise would be
 * testing the fixture rather than the code.
 *
 * A *translating* field, uniform through the interior, is the right fixture
 * for the accuracy comparison: away from the two end cells the flow is
 * divergence-free, so a pulse kept clear of the walls translates cleanly and
 * the surviving peak measures numerical diffusion and nothing else.
 */
struct Advector
{
    HorizontalMesh horizontal;
    LayeredMesh    layered;

    std::vector<double> edgeFlux;     //!< m³/s per (edge, layer), left to right.
    std::vector<double> verticalFlux; //!< m/s upward, per interface slot.
    std::vector<double> volume;

    double cellSize  = 0.0;
    double velocity  = 0.0;
    double thickness = 0.0;

    void makeMesh(int columns, int layers, double length, double width, double depth)
    {
        std::string message;

        ASSERT_TRUE(horizontal.build(Testing::channelStrip(columns, length, width),
                                     message))
            << message;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = -depth;

        SigmaSpec spec;
        spec.layerCount        = layers;
        spec.minLayerThickness = 0.0;

        ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

        cellSize  = length / static_cast<double>(columns);
        thickness = depth / static_cast<double>(layers);

        edgeFlux.assign(static_cast<size_t>(horizontal.edgeCount()) *
                            static_cast<size_t>(layers),
                        0.0);
        verticalFlux.assign(static_cast<size_t>(layered.interfaceSlotCount()), 0.0);
        volume = layered.cellVolume();
    }

    /*!
     * \brief Is this edge one of the streamwise faces?
     */
    [[nodiscard]] bool isStreamwise(int64_t edge) const
    {
        const size_t ei = static_cast<size_t>(edge);

        return horizontal.edgeRightFace()[ei] != FVQual::Mesh::kNoFace &&
               std::abs(horizontal.edgeNormalX()[ei]) > 0.5;
    }

    /*!
     * \brief Uniform interior flow at \p speed, for accuracy measurements.
     */
    void buildTranslation(int columns, int layers, double length, double width,
                          double depth, double speed)
    {
        makeMesh(columns, layers, length, width, depth);
        velocity = speed;

        for (int64_t e = 0; e < horizontal.edgeCount(); ++e)
        {
            if (!isStreamwise(e))
                continue;

            const size_t ei   = static_cast<size_t>(e);
            const double sign = horizontal.edgeNormalX()[ei] > 0.0 ? 1.0 : -1.0;

            for (int k = 0; k < layers; ++k)
            {
                edgeFlux[ei * static_cast<size_t>(layers) + static_cast<size_t>(k)] =
                    sign * speed * thickness * horizontal.edgeLength()[ei];
            }
        }
    }

    /*!
     * \brief A closed, exactly divergence-free circulation.
     */
    void buildCirculation(int columns, int layers, double length, double width,
                          double depth, double speed)
    {
        makeMesh(columns, layers, length, width, depth);
        velocity = speed;

        const int last = layers - 1;

        double discharge = 0.0;

        for (int64_t e = 0; e < horizontal.edgeCount(); ++e)
        {
            if (!isStreamwise(e))
                continue;

            const size_t ei   = static_cast<size_t>(e);
            const double sign = horizontal.edgeNormalX()[ei] > 0.0 ? 1.0 : -1.0;

            discharge = speed * thickness * horizontal.edgeLength()[ei];

            // Rightward along the surface, leftward along the bed.
            edgeFlux[ei * static_cast<size_t>(layers)] = sign * discharge;
            edgeFlux[ei * static_cast<size_t>(layers) + static_cast<size_t>(last)] =
                -sign * discharge;
        }

        // The two end columns carry the circuit vertically: up at the left,
        // down at the right, through every interior interface.
        const int64_t leftColumn  = 0;
        const int64_t rightColumn = layered.columnCount() - 1;

        for (int k = 1; k < layers; ++k)
        {
            verticalFlux[static_cast<size_t>(layered.interfaceSlot(leftColumn, k))] =
                discharge / horizontal.faceArea()[static_cast<size_t>(leftColumn)];

            verticalFlux[static_cast<size_t>(layered.interfaceSlot(rightColumn, k))] =
                -discharge / horizontal.faceArea()[static_cast<size_t>(rightColumn)];
        }
    }

    [[nodiscard]] double totalMass(const std::vector<double> &values) const
    {
        double total = 0.0;

        for (int64_t i = 0; i < layered.cellCount(); ++i)
            total += values[static_cast<size_t>(i)] * volume[static_cast<size_t>(i)];

        return total;
    }

    /*!
     * \brief Largest net volume imbalance over the cells (m³/s) — zero for a
     *        divergence-free field, and worth asserting so the fixture cannot
     *        quietly stop being one.
     */
    [[nodiscard]] double peakDivergence() const
    {
        const int layers = layered.layerCount();

        double peak = 0.0;

        for (int64_t c = 0; c < layered.columnCount(); ++c)
        {
            const int64_t begin = horizontal.faceEdgeOffsets()[static_cast<size_t>(c)];
            const int64_t end   = horizontal.faceEdgeOffsets()[static_cast<size_t>(c) + 1];

            for (int k = 0; k < layers; ++k)
            {
                double net = 0.0;

                for (int64_t index = begin; index < end; ++index)
                {
                    const size_t ei = static_cast<size_t>(
                        horizontal.faceEdges()[static_cast<size_t>(index)]);

                    if (horizontal.edgeRightFace()[ei] == FVQual::Mesh::kNoFace)
                        continue;

                    net += (horizontal.edgeLeftFace()[ei] == c ? 1.0 : -1.0) *
                           edgeFlux[ei * static_cast<size_t>(layers) +
                                    static_cast<size_t>(k)];
                }

                const double area = horizontal.faceArea()[static_cast<size_t>(c)];

                net += verticalFlux[static_cast<size_t>(layered.interfaceSlot(c, k))] * area;
                net -= verticalFlux[static_cast<size_t>(layered.interfaceSlot(c, k + 1))] *
                       area;

                peak = std::max(peak, std::abs(net));
            }
        }

        return peak;
    }
};

/*!
 * \brief Advances a field for a fixed duration under a scheme.
 */
void translate(Advector &setup, EulerianTransport &engine, AdvectionScheme scheme,
               std::vector<double> &values, double duration, double courant = 0.4)
{
    const double dt = courant * setup.cellSize / std::abs(setup.velocity);

    for (double t = 0.0; t < duration - 1.0e-12; t += dt)
    {
        const double step = std::min(dt, duration - t);

        engine.advect(setup.layered, setup.volume, setup.edgeFlux, setup.verticalFlux,
                      step, scheme, values);
    }
}

} // namespace

TEST(EulerianTransport, UniformFieldIsPreservedExactly)
{
    // The reconstruction must produce zero gradient for a uniform field, so
    // the high-order and low-order fluxes coincide and nothing moves.
    Advector setup;
    setup.buildCirculation(24, 6, 240.0, 20.0, 12.0, 0.5);

    // The fixture must actually be divergence-free, or this proves nothing.
    ASSERT_LT(setup.peakDivergence(), 1.0e-9);

    EulerianTransport engine;
    engine.initialize(setup.layered);

    std::vector<double> values(static_cast<size_t>(setup.layered.cellCount()), 7.25);

    translate(setup, engine, AdvectionScheme::MusclFct, values, 200.0);

    for (int64_t i = 0; i < setup.layered.cellCount(); ++i)
        EXPECT_NEAR(values[static_cast<size_t>(i)], 7.25, 1.0e-12) << "cell " << i;
}

TEST(EulerianTransport, MassIsConservedByBothSchemes)
{
    for (AdvectionScheme scheme : {AdvectionScheme::Upwind, AdvectionScheme::MusclFct})
    {
        Advector setup;
        setup.buildCirculation(40, 4, 400.0, 20.0, 8.0, 0.5);

        EulerianTransport engine;
        engine.initialize(setup.layered);

        std::vector<double> values(static_cast<size_t>(setup.layered.cellCount()), 0.0);

        // A blob well away from the walls.
        for (int64_t c = 8; c < 16; ++c)
        {
            for (int k = 0; k < setup.layered.layerCount(); ++k)
                values[static_cast<size_t>(setup.layered.cell(c, k))] = 4.0;
        }

        const double before = setup.totalMass(values);
        ASSERT_GT(before, 0.0);

        translate(setup, engine, scheme, values, 120.0);

        EXPECT_LT(std::abs(setup.totalMass(values) - before) / before, 1.0e-12)
            << "scheme " << static_cast<int>(scheme);
    }
}

TEST(EulerianTransport, SharpFrontStaysBounded)
{
    // The property that makes flux correction worth its cost: a second-order
    // reconstruction alone overshoots at a discontinuity, and an overshoot in
    // a concentration is a negative oxygen or a negative biomass.
    Advector setup;
    setup.buildCirculation(60, 6, 600.0, 20.0, 12.0, 0.5);

    ASSERT_LT(setup.peakDivergence(), 1.0e-9);

    EulerianTransport engine;
    engine.initialize(setup.layered);

    std::vector<double> values(static_cast<size_t>(setup.layered.cellCount()), 0.0);

    for (int64_t c = 0; c < 20; ++c)
    {
        for (int k = 0; k < setup.layered.layerCount(); ++k)
            values[static_cast<size_t>(setup.layered.cell(c, k))] = 1.0;
    }

    translate(setup, engine, AdvectionScheme::MusclFct, values, 300.0);

    for (int64_t i = 0; i < setup.layered.cellCount(); ++i)
    {
        const double value = values[static_cast<size_t>(i)];

        EXPECT_GE(value, -1.0e-12) << "undershoot at cell " << i;
        EXPECT_LE(value, 1.0 + 1.0e-12) << "overshoot at cell " << i;
    }
}

TEST(EulerianTransport, MusclIsMoreAccurateThanUpwind)
{
    // The accuracy claim, measured against the exact solution rather than
    // asserted. A square pulse is translated 20 cells and compared with where
    // it should be.
    //
    // Note the metric: the *peak* is useless here, because a 16-cell plateau
    // keeps its middle value under any scheme — first-order upwind erodes the
    // edges while the centre sits at 1.0. Integrated error is what actually
    // measures the smearing.
    constexpr int    kColumns  = 80;
    constexpr double kLength   = 800.0;
    constexpr double kSpeed    = 0.5;
    constexpr double kDuration = 400.0;
    constexpr int    kStart    = 20;
    constexpr int    kEnd      = 36;

    const auto errorFor = [&](AdvectionScheme scheme) {
        Advector setup;
        setup.buildTranslation(kColumns, 4, kLength, 20.0, 8.0, kSpeed);

        EulerianTransport engine;
        engine.initialize(setup.layered);

        std::vector<double> values(static_cast<size_t>(setup.layered.cellCount()), 0.0);

        for (int64_t c = kStart; c < kEnd; ++c)
        {
            for (int k = 0; k < setup.layered.layerCount(); ++k)
                values[static_cast<size_t>(setup.layered.cell(c, k))] = 1.0;
        }

        translate(setup, engine, scheme, values, kDuration);

        // Exact solution: the same pulse, shifted by speed * duration.
        const int shift = static_cast<int>(std::lround(kSpeed * kDuration /
                                                       setup.cellSize));

        double error = 0.0;
        int64_t counted = 0;

        for (int64_t c = 0; c < setup.layered.columnCount(); ++c)
        {
            const double exact =
                (c >= kStart + shift && c < kEnd + shift) ? 1.0 : 0.0;

            for (int k = 0; k < setup.layered.layerCount(); ++k)
            {
                error += std::abs(
                    values[static_cast<size_t>(setup.layered.cell(c, k))] - exact);
                ++counted;
            }
        }

        return error / static_cast<double>(counted);
    };

    const double upwindError = errorFor(AdvectionScheme::Upwind);
    const double musclError  = errorFor(AdvectionScheme::MusclFct);

    EXPECT_LT(musclError, 0.5 * upwindError)
        << "MUSCL L1 error " << musclError << " vs upwind " << upwindError;
}

TEST(EulerianTransport, VerticalTransportIsAlsoBoundedAndConservative)
{
    // The vertical direction has its own reconstruction and its own half of
    // the flux correction, so a purely horizontal test would hide a bug
    // there. The circulation carries the tracer up one end column and down
    // the other, which exercises both.
    Advector setup;
    setup.buildCirculation(12, 20, 240.0, 20.0, 20.0, 0.4);

    ASSERT_LT(setup.peakDivergence(), 1.0e-9);

    EulerianTransport engine;
    engine.initialize(setup.layered);

    std::vector<double> values(static_cast<size_t>(setup.layered.cellCount()), 0.0);

    // A band in the middle layers, which the circuit must carry vertically.
    for (int64_t c = 0; c < setup.layered.columnCount(); ++c)
    {
        for (int k = 8; k < 12; ++k)
            values[static_cast<size_t>(setup.layered.cell(c, k))] = 2.0;
    }

    const double before = setup.totalMass(values);
    ASSERT_GT(before, 0.0);

    translate(setup, engine, AdvectionScheme::MusclFct, values, 600.0);

    EXPECT_LT(std::abs(setup.totalMass(values) - before) / before, 1.0e-12);

    for (int64_t i = 0; i < setup.layered.cellCount(); ++i)
    {
        EXPECT_GE(values[static_cast<size_t>(i)], -1.0e-12) << "cell " << i;
        EXPECT_LE(values[static_cast<size_t>(i)], 2.0 + 1.0e-12) << "cell " << i;
    }

    // The tracer must have actually moved out of its initial band.
    double outsideBand = 0.0;
    for (int64_t c = 0; c < setup.layered.columnCount(); ++c)
    {
        for (int k = 0; k < 6; ++k)
            outsideBand = std::max(outsideBand,
                                   values[static_cast<size_t>(setup.layered.cell(c, k))]);
    }

    EXPECT_GT(outsideBand, 1.0e-3) << "nothing was transported vertically";
}

TEST(EulerianTransport, FluxCorrectionStaysOutOfTheWayInSmoothFlow)
{
    // Second-order accuracy survives only if the correction is admitted where
    // the solution is smooth; a limiter that clips everywhere would silently
    // degrade the scheme back to first order.
    //
    // Worth recording what this diagnostic does *not* show. It averages over
    // faces that carry a nonzero antidiffusive flux, and at a discontinuity
    // the Barth-Jespersen limiter has usually already driven the gradient to
    // zero — so those faces drop out of the average entirely and Zalesak has
    // little left to do. The admitted fraction is therefore not a monotone
    // indicator of front sharpness, and boundedness at fronts is pinned
    // separately by SharpFrontStaysBounded.
    Advector setup;
    setup.buildTranslation(60, 4, 600.0, 20.0, 8.0, 0.5);

    EulerianTransport engine;
    engine.initialize(setup.layered);

    std::vector<double> values(static_cast<size_t>(setup.layered.cellCount()), 0.0);

    for (int64_t c = 0; c < setup.layered.columnCount(); ++c)
    {
        const double x = setup.horizontal.faceCentroidX()[static_cast<size_t>(c)];
        const double value =
            1.0 + 0.5 * std::sin(2.0 * 3.14159265358979 * x / 600.0);

        for (int k = 0; k < setup.layered.layerCount(); ++k)
            values[static_cast<size_t>(setup.layered.cell(c, k))] = value;
    }

    const double dt = 0.4 * setup.cellSize / setup.velocity;

    engine.advect(setup.layered, setup.volume, setup.edgeFlux, setup.verticalFlux, dt,
                  AdvectionScheme::MusclFct, values);

    EXPECT_GT(engine.lastAdmittedFraction(), 0.9)
        << "the limiter is clipping a smooth field";
}
