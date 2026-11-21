// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_boundarytransport.cpp
 * \brief Three-dimensional scalar transport across open boundaries (F2.3b).
 *
 * Until now the internal mode refused to start with an open boundary, because
 * the external mode could pass water through one but nothing said what that
 * water was made of. This is that story, and the gates are about the two claims
 * it rests on.
 *
 * **The modes still agree.** Mode splitting is only trustworthy while the
 * vertical velocity diagnosed from continuity vanishes at the free surface —
 * that residual is what proves the internal layer fluxes moved the same water
 * the external mode did. Opening a boundary is precisely the change most likely
 * to break it, so it is checked with the boundary open and flowing.
 *
 * **An unconfigured boundary is transparent, not clean.** Inflow carries the
 * receiving column's own values unless told otherwise, so a uniform field
 * crosses unchanged. The alternative default — zero — would import distilled
 * water through every boundary nobody had got round to configuring.
 *
 * **Three defects these gates caught during F2.3b**, each of which would have
 * produced a model that ran and looked plausible:
 *
 * - Boundary water was reported to the external mode as a column source *and*
 *   crossed the boundary edge, so every open boundary moved its water twice.
 *   Ablating the `externallyAccounted` flag reproduces it.
 * - The density placement was applied to outflow as well as inflow. A
 *   density-matched search has no meaning for water leaving, and it placed the
 *   withdrawal somewhere the volume change was not: the mode-consistency
 *   residual sat at 1.2e-2 against 1e-15 for uniform.
 * - The ambient composition was read before `distribute()` had computed this
 *   step's placement weights, so it evaluated to zero — an open boundary
 *   quietly importing distilled water, which is precisely the failure the
 *   ambient default exists to prevent. It showed up as the bed layer of the
 *   inflow column diluting by exactly the ratio of inflow volume to layer
 *   volume.
 */

#include "fvqual/hydro/internalmode.h"
#include "fvqual/mesh/horizontalmesh.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/mesh/sigmaspec.h"
#include "fvqual/state/speciesregistry.h"
#include "fvqual/state/statestore.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

using FVQual::Forcing::Placement;
using FVQual::Hydro::BarotropicSolver;
using FVQual::Hydro::InternalModeOptions;
using FVQual::Hydro::InternalModeSolver;
using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::kNoFace;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;
using FVQual::State::SpeciesDefinition;
using FVQual::State::SpeciesRegistry;
using FVQual::State::StateStore;
namespace Testing = FVQual::Testing;

namespace {

/*!
 * \brief A channel with an inflow boundary at one end and a stage at the other.
 */
struct Channel
{
    HorizontalMesh     horizontal;
    LayeredMesh        layered;
    SpeciesRegistry    registry;
    StateStore         state;
    InternalModeSolver solver;

    int temperature = -1;
    int dye         = -1;

    int inlet  = -1;
    int outlet = -1;

    double depth     = 0.0;
    double thickness = 0.0;
    int    layers    = 0;

    //! Edges on the left (x = 0) and right (x = length) walls.
    std::vector<int64_t> leftEdges;
    std::vector<int64_t> rightEdges;

    void build(int columns, int layerCount, double length, double width, double totalDepth,
               double waterTemperature = 12.0)
    {
        std::string message;

        ASSERT_TRUE(horizontal.build(Testing::channelStrip(columns, length, width),
                                     message))
            << message;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = -totalDepth;

        SigmaSpec spec;
        spec.layerCount        = layerCount;
        spec.minLayerThickness = 0.0;

        ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

        ASSERT_GE(registry.addTemperature(waterTemperature, message), 0) << message;
        temperature = registry.temperatureIndex();

        SpeciesDefinition tracer;
        tracer.id    = "DYE";
        tracer.units = "g/m3";
        ASSERT_GE(registry.add(tracer, message), 0) << message;
        dye = registry.indexOf("DYE");

        state.allocate(layered, registry);

        for (int64_t cell = 0; cell < layered.cellCount(); ++cell)
            state.speciesAt(temperature, cell) = waterTemperature;

        InternalModeOptions options;
        options.gravity         = 9.80665;
        options.manningN        = 0.0;
        options.cflNumber       = 0.5;
        options.surfaceHeatFlux = 0.0;

        ASSERT_TRUE(solver.initialize(layered, state, registry, options, message))
            << message;

        depth     = totalDepth;
        layers    = layerCount;
        thickness = totalDepth / static_cast<double>(layerCount);

        // Boundary edges are the vertical walls at each end of the strip: the
        // ones with no neighbor whose centroid sits at x = 0 or x = length.
        for (int64_t e = 0; e < horizontal.edgeCount(); ++e)
        {
            if (horizontal.edgeRightFace()[static_cast<size_t>(e)] != kNoFace)
                continue;

            const double x = horizontal.edgeMidX()[static_cast<size_t>(e)];

            if (x < 1.0e-6)
                leftEdges.push_back(e);
            else if (x > length - 1.0e-6)
                rightEdges.push_back(e);
        }

        ASSERT_FALSE(leftEdges.empty());
        ASSERT_FALSE(rightEdges.empty());
    }

    //! Adds a prescribed-discharge inlet and a fixed-stage outlet.
    void openEnds(double specificDischarge, double stage,
                  Placement placement = Placement::DensityMatched)
    {
        std::string message;

        BarotropicSolver::BoundaryCondition inflow;
        inflow.type              = BarotropicSolver::BoundaryType::Inflow;
        inflow.specificDischarge = specificDischarge;

        inlet = solver.external().addBoundary(inflow, leftEdges, message);
        ASSERT_GE(inlet, 0) << message;

        BarotropicSolver::BoundaryCondition level;
        level.type  = BarotropicSolver::BoundaryType::Stage;
        level.stage = stage;

        outlet = solver.external().addBoundary(level, rightEdges, message);
        ASSERT_GE(outlet, 0) << message;

        ASSERT_TRUE(solver.bindOpenBoundaries(placement, message)) << message;
    }

    void setUniform(int species, double value)
    {
        for (int64_t cell = 0; cell < layered.cellCount(); ++cell)
            state.speciesAt(species, cell) = value;
    }

    [[nodiscard]] double at(int species, int64_t column, int k) const
    {
        return state.speciesAt(species, layered.cell(column, k));
    }

    [[nodiscard]] double totalMass(int species) const
    {
        double total = 0.0;

        for (int64_t cell = 0; cell < layered.cellCount(); ++cell)
        {
            total += state.speciesAt(species, cell) *
                     layered.cellVolume()[static_cast<size_t>(cell)];
        }

        return total;
    }

    void run(double dt, int steps)
    {
        std::string message;

        for (int i = 0; i < steps; ++i)
            ASSERT_TRUE(solver.step(dt, message)) << message;
    }
};

} // namespace

// ── The refusal is gone, and the modes still agree ────────────────────────

TEST(BoundaryTransport, TheInternalModeNowAcceptsAnOpenBoundary)
{
    Channel channel;
    channel.build(10, 5, 500.0, 40.0, 10.0);

    std::string message;

    // Before F2.3b this configuration was refused outright.
    channel.openEnds(0.2, 0.0);

    EXPECT_TRUE(channel.solver.external().hasOpenBoundaries());

    ASSERT_TRUE(channel.solver.step(10.0, message)) << message;
}

TEST(BoundaryTransport, ModeConsistencySurvivesAnOpenBoundary)
{
    // The gate that matters most. The vertical velocity is integrated up from a
    // no-flux bed, so its value at the free surface must vanish if the internal
    // fluxes and sources together account for exactly what the external mode
    // moved. Opening a boundary adds a whole new way for the two modes to
    // disagree, so the residual is checked while water is actually crossing.
    Channel channel;
    channel.build(12, 6, 600.0, 40.0, 12.0);
    channel.openEnds(0.3, 0.0);

    std::string message;

    double worst = 0.0;

    for (int i = 0; i < 30; ++i)
    {
        ASSERT_TRUE(channel.solver.step(5.0, message)) << message;

        worst = std::max(worst, std::abs(channel.solver.diagnostics().surfaceFluxResidual));
    }

    // The same tolerance the closed-basin gates hold to: this is a round-off
    // identity, not an approximation.
    EXPECT_LT(worst, 1.0e-10) << "surface flux residual reached " << worst;

    // And water really was crossing, or the residual is trivially zero.
    EXPECT_GT(channel.solver.boundaryVolume(channel.inlet), 0.0);
    EXPECT_LT(channel.solver.boundaryVolume(channel.outlet), 0.0);
}

// ── An unset boundary is transparent ──────────────────────────────────────

TEST(BoundaryTransport, AnAmbientBoundaryLeavesAUniformFieldAlone)
{
    // The default. Water crosses, volume changes, and the constituent field does
    // not move — which is what "ambient" has to mean for it to be a neutral
    // choice rather than a hidden one.
    constexpr double c0 = 7.0;

    Channel channel;
    channel.build(10, 5, 500.0, 40.0, 10.0);
    channel.openEnds(0.25, 0.0);
    channel.setUniform(channel.dye, c0);

    channel.run(5.0, 40);

    for (int64_t c = 0; c < channel.layered.columnCount(); ++c)
    {
        for (int k = 0; k < channel.layers; ++k)
        {
            EXPECT_NEAR(channel.at(channel.dye, c, k), c0, 1.0e-8 * c0)
                << "column " << c << " layer " << k;
        }
    }

    // Volume did change, so the field held still through real throughflow
    // rather than through nothing happening.
    EXPECT_GT(channel.solver.boundaryVolume(channel.inlet), 0.0);
}

// ── A prescribed inflow ───────────────────────────────────────────────────

TEST(BoundaryTransport, APrescribedInflowDeliversTheMassItsFluxImplies)
{
    // The constituent budget, closed against the boundary rather than assumed.
    // With a clean domain and an inlet carrying a stated concentration, the mass
    // that accumulates must equal the concentration times the volume that
    // entered — minus whatever left through the outlet at ambient.
    constexpr double inflowValue = 5.0;

    Channel channel;
    channel.build(12, 4, 600.0, 40.0, 10.0);
    channel.openEnds(0.3, 0.0);
    channel.setUniform(channel.dye, 0.0);

    std::string message;
    ASSERT_TRUE(channel.solver.setBoundaryValue(channel.inlet, "DYE", inflowValue,
                                                message))
        << message;

    double delivered = 0.0;
    double exported  = 0.0;

    const double before = channel.totalMass(channel.dye);

    for (int i = 0; i < 40; ++i)
    {
        // The outlet exports at ambient, so what leaves has to be measured as
        // it happens — the concentration there is changing.
        const double outletBefore = channel.at(channel.dye,
                                               channel.layered.columnCount() - 1, 0);

        ASSERT_TRUE(channel.solver.step(5.0, message)) << message;

        delivered += channel.solver.boundaryVolume(channel.inlet) * inflowValue;

        // A first-order estimate of the export is enough for the bound below;
        // the exact accounting is the gate two tests down.
        exported += -channel.solver.boundaryVolume(channel.outlet) * outletBefore;
    }

    const double gained = channel.totalMass(channel.dye) - before;

    ASSERT_GT(delivered, 0.0);

    // Everything that arrived either stayed or left.
    EXPECT_NEAR(gained, delivered - exported, 0.05 * delivered)
        << "gained " << gained << ", delivered " << delivered << ", exported "
        << exported;

    // The dye reached the domain and is strongest near the inlet.
    EXPECT_GT(channel.at(channel.dye, 0, 0), 0.1 * inflowValue);
    EXPECT_GT(channel.at(channel.dye, 0, 0),
              channel.at(channel.dye, channel.layered.columnCount() - 1, 0));
}

TEST(BoundaryTransport, WithNoOutletTheBudgetClosesExactly)
{
    // A closed-ended channel with one inflow boundary: nothing leaves, so the
    // mass gained must equal concentration times volume entered, to round-off.
    // This is the gate the previous one approximates.
    constexpr double inflowValue = 4.0;

    Channel channel;
    channel.build(8, 4, 400.0, 40.0, 10.0);

    std::string message;

    BarotropicSolver::BoundaryCondition inflow;
    inflow.type              = BarotropicSolver::BoundaryType::Inflow;
    inflow.specificDischarge = 0.15;

    channel.inlet = channel.solver.external().addBoundary(inflow, channel.leftEdges,
                                                          message);
    ASSERT_GE(channel.inlet, 0) << message;
    ASSERT_TRUE(channel.solver.bindOpenBoundaries(message)) << message;

    channel.setUniform(channel.dye, 0.0);

    ASSERT_TRUE(channel.solver.setBoundaryValue(channel.inlet, "DYE", inflowValue,
                                                message))
        << message;

    double entered = 0.0;

    for (int i = 0; i < 30; ++i)
    {
        ASSERT_TRUE(channel.solver.step(4.0, message)) << message;
        entered += channel.solver.boundaryVolume(channel.inlet);
    }

    ASSERT_GT(entered, 0.0);

    EXPECT_NEAR(channel.totalMass(channel.dye), entered * inflowValue,
                1.0e-6 * entered * inflowValue);
}

// ── Placement ─────────────────────────────────────────────────────────────

TEST(BoundaryTransport, ADenseInflowPlungesAndAUniformOneDoesNot)
{
    // The reason placement is selectable, and the reason density is the
    // default. A cold tributary entering a stratified reservoir sinks to its own
    // density rather than mixing across the surface, and the two placements have
    // to give visibly different answers or the choice is decorative.
    const auto runWith = [](Placement placement) {
        Channel channel;

        // Warm at the surface, cold at depth — a real thermocline for the
        // inflow to find its level in.
        channel.build(10, 8, 500.0, 40.0, 16.0, /*waterTemperature=*/20.0);

        for (int64_t c = 0; c < channel.layered.columnCount(); ++c)
        {
            for (int k = 0; k < channel.layers; ++k)
            {
                const double z =
                    -0.5 * (channel.layered.z(c, k) + channel.layered.z(c, k + 1));

                channel.state.speciesAt(channel.temperature, channel.layered.cell(c, k)) =
                    20.0 - 10.0 * z / 16.0;
            }
        }

        channel.openEnds(0.25, 0.0, placement);
        channel.setUniform(channel.dye, 0.0);

        std::string message;

        // Cold, dense inflow, marked with dye so its path is visible.
        EXPECT_TRUE(channel.solver.setBoundaryValue(channel.inlet, "DYE", 1.0, message))
            << message;
        EXPECT_TRUE(channel.solver.setBoundaryValue(
            channel.inlet, SpeciesRegistry::kTemperature, 8.0, message))
            << message;

        channel.run(4.0, 25);

        // Where did the dye go? Report the depth-weighted centroid in the
        // receiving column.
        double moment = 0.0;
        double mass   = 0.0;

        for (int k = 0; k < channel.layers; ++k)
        {
            const double z =
                -0.5 * (channel.layered.z(0, k) + channel.layered.z(0, k + 1));
            const double value = channel.at(channel.dye, 0, k);

            moment += value * z;
            mass += value;
        }

        EXPECT_GT(mass, 0.0);

        return mass > 0.0 ? moment / mass : 0.0;
    };

    const double densityMatched = runWith(Placement::DensityMatched);
    const double uniform        = runWith(Placement::Uniform);

    // A uniform placement puts the inflow at the column's mid-depth; a dense
    // one puts it deeper. If these agreed, the placement would not be reaching
    // the boundary source at all.
    EXPECT_GT(densityMatched, uniform + 1.0)
        << "density-matched centroid " << densityMatched << " m, uniform " << uniform
        << " m";
}

// ── Refusals ──────────────────────────────────────────────────────────────

TEST(BoundaryTransport, MisuseIsRefused)
{
    Channel channel;
    channel.build(8, 4, 400.0, 40.0, 10.0);

    std::string message;

    // Nothing to bind yet.
    EXPECT_FALSE(channel.solver.bindOpenBoundaries(message));
    EXPECT_NE(message.find("no open boundaries"), std::string::npos) << message;

    // Values cannot be set before binding, because there is nothing to set them
    // on and silently remembering them would be worse.
    EXPECT_FALSE(channel.solver.setBoundaryValue(0, "DYE", 1.0, message));

    channel.openEnds(0.2, 0.0);

    // Binding twice would create a second source per boundary and double the
    // inflow.
    EXPECT_FALSE(channel.solver.bindOpenBoundaries(message));
    EXPECT_NE(message.find("already bound"), std::string::npos) << message;

    EXPECT_FALSE(channel.solver.setBoundaryValue(99, "DYE", 1.0, message));
    EXPECT_FALSE(channel.solver.setBoundaryValue(channel.inlet, "GHOST", 1.0, message));

    EXPECT_TRUE(channel.solver.setBoundaryValue(channel.inlet, "DYE", 1.0, message))
        << message;
    EXPECT_TRUE(channel.solver.clearBoundaryValue(channel.inlet, "DYE", message))
        << message;
}


