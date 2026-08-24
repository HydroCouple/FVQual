// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_settling.cpp
 * \brief Verification of particulate settling and the benthic handoff (F6.2c).
 *
 * Three claims, in descending order of how much trouble getting them wrong
 * would cause.
 *
 * **Nothing is created or destroyed.** What leaves the water column arrives on
 * the sediment, exactly. This is an engine identity — unlike reaction
 * conservation, which belongs to the modeller — so it is gated to round-off, and
 * the budget includes the benthic store rather than stopping at the waterline.
 *
 * **Nothing goes negative, at any time step.** This is the whole reason the
 * solve is implicit. An explicit upwind settling term needs `w dt < h`, and the
 * layer that binds is the thin one at the surface, not the fast-settling
 * species. The gate below runs at a hundred times the explicit stability limit
 * and requires the answer to stay positive and conservative — which an explicit
 * scheme cannot do, and which is checked by ablation.
 *
 * **The velocity means what it says.** The profile translates downward at the
 * settling velocity, measured as L1 convergence against the exact translated
 * solution with the rate measured rather than asserted.
 *
 * That gate began as a centroid displacement and was wrong twice, both mine:
 * the blob travelled 40 m in a 20 m column, so it was measuring the remnant
 * left behind, and the refinement halved the layer thickness while holding the
 * step fixed, which leaves implicit upwind's O(w dt) error as a constant floor
 * and drives the convergence ratio to 1. The step is now refined with the grid
 * at a fixed Courant number.
 */

#include "fvqual/hydro/internalmode.h"
#include "fvqual/mesh/horizontalmesh.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/mesh/sigmaspec.h"
#include "fvqual/state/speciesregistry.h"
#include "fvqual/state/statestore.h"
#include "fvqual/transport/settling.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

using FVQual::Hydro::InternalModeOptions;
using FVQual::Hydro::InternalModeSolver;
using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;
using FVQual::State::SpeciesDefinition;
using FVQual::State::SpeciesKind;
using FVQual::State::SpeciesRegistry;
using FVQual::State::StateStore;
using FVQual::Transport::SettlingSet;
namespace Testing = FVQual::Testing;

namespace {

//! A settling velocity of one metre per day, in SI.
constexpr double kMetrePerDay = 1.0 / 86400.0;

/*!
 * \brief A column with a settling particulate and a bed store for it.
 */
struct Column
{
    HorizontalMesh  horizontal;
    LayeredMesh     layered;
    SpeciesRegistry registry;
    StateStore      state;
    SettlingSet     settling;

    int particulate = -1;
    int sediment    = -1;
    int dissolved   = -1;

    double depth     = 0.0;
    double thickness = 0.0;

    void build(int layers, double totalDepth, bool withSediment = true,
               bool withDissolved = false)
    {
        std::string message;

        ASSERT_TRUE(horizontal.build(Testing::channelStrip(3, 300.0, 40.0), message))
            << message;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = -totalDepth;

        SigmaSpec spec;
        spec.layerCount        = layers;
        spec.minLayerThickness = 0.0;

        ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

        SpeciesDefinition solids;
        solids.id    = "TSS";
        solids.units = "g/m3";
        ASSERT_GE(registry.add(solids, message), 0) << message;
        particulate = registry.indexOf("TSS");

        if (withDissolved)
        {
            SpeciesDefinition tracer;
            tracer.id    = "CL";
            tracer.units = "g/m3";
            ASSERT_GE(registry.add(tracer, message), 0) << message;
            dissolved = registry.indexOf("CL");
        }

        if (withSediment)
        {
            SpeciesDefinition bed;
            bed.id    = "TSS_BED";
            bed.kind  = SpeciesKind::MsxBenthic;
            bed.units = "g/m2";
            ASSERT_GE(registry.add(bed, message), 0) << message;
            sediment = registry.indexOf("TSS_BED");
        }

        state.allocate(layered, registry);
        settling.initialize(registry);

        depth     = totalDepth;
        thickness = totalDepth / static_cast<double>(layers);
    }

    void setUniform(int species, double value)
    {
        for (int64_t i = 0; i < layered.cellCount(); ++i)
            state.speciesAt(species, i) = value;
    }

    //! Volume-weighted mass of \p species in the water column (g).
    [[nodiscard]] double waterMass(int species) const
    {
        double total = 0.0;

        for (int64_t cell = 0; cell < layered.cellCount(); ++cell)
        {
            total += state.speciesAt(species, cell) *
                     layered.cellVolume()[static_cast<size_t>(cell)];
        }

        return total;
    }

    //! Area-weighted mass in the benthic store (g).
    [[nodiscard]] double bedMass() const
    {
        if (sediment < 0)
            return 0.0;

        double total = 0.0;

        for (int64_t c = 0; c < layered.columnCount(); ++c)
        {
            total += state.speciesAt(sediment, layered.cell(c, layered.layerCount() - 1)) *
                     horizontal.faceArea()[static_cast<size_t>(c)];
        }

        return total;
    }

    [[nodiscard]] double at(int species, int k) const
    {
        return state.speciesAt(species, layered.cell(0, k));
    }
};

} // namespace

// ── Configuration ─────────────────────────────────────────────────────────

TEST(Settling, MisconfigurationIsRefused)
{
    Column column;
    column.build(5, 10.0);

    std::string message;

    // Upward motion is buoyancy, and it would make the triangular solve below
    // wrong rather than merely inaccurate — the substitution assumes each cell's
    // inflow comes from the one above.
    EXPECT_FALSE(column.settling.setVelocity(column.particulate, -1.0e-5,
                                             column.registry, message));
    EXPECT_FALSE(message.empty());

    EXPECT_FALSE(column.settling.setVelocity(99, 1.0e-5, column.registry, message));

    // A benthic species is already on the bed.
    EXPECT_FALSE(column.settling.setVelocity(column.sediment, 1.0e-5, column.registry,
                                             message));

    // A transported receiver would be advected away, and the sediment would flow
    // downstream.
    ASSERT_TRUE(column.settling.setVelocity(column.particulate, kMetrePerDay,
                                            column.registry, message))
        << message;

    EXPECT_FALSE(column.settling.setBenthicReceiver(column.particulate,
                                                    column.particulate,
                                                    column.registry, message));
    EXPECT_TRUE(column.settling.setBenthicReceiver(column.particulate, column.sediment,
                                                   column.registry, message))
        << message;
}

TEST(Settling, ASpeciesWithNoVelocityDoesNotMove)
{
    // Dissolved constituents must be untouched, and it must cost nothing to
    // leave them so — the alternative is every model needing an explicit zero
    // for every solute.
    Column column;
    column.build(5, 10.0, /*withSediment=*/true, /*withDissolved=*/true);
    column.setUniform(column.particulate, 4.0);
    column.setUniform(column.dissolved, 7.0);

    std::string message;
    ASSERT_TRUE(column.settling.setVelocity(column.particulate, 10.0 * kMetrePerDay,
                                            column.registry, message))
        << message;

    EXPECT_TRUE(column.settling.any());

    column.settling.apply(column.layered, column.state, 3600.0);

    for (int k = 0; k < 5; ++k)
        EXPECT_EQ(column.at(column.dissolved, k), 7.0) << "layer " << k;

    EXPECT_EQ(column.settling.deliveredToBed(column.dissolved), 0.0);
}

// ── Conservation ──────────────────────────────────────────────────────────

TEST(Settling, WaterPlusSedimentIsConservedExactly)
{
    // The engine identity. Not "approximately conserved": the telescoping sum of
    // the substitution leaves exactly `w dt C_bed` per unit area, and that is
    // exactly what is added to the store.
    constexpr double c0    = 12.0;
    constexpr double dt    = 1800.0;
    constexpr int    steps = 40;

    Column column;
    column.build(8, 16.0);
    column.setUniform(column.particulate, c0);

    std::string message;
    ASSERT_TRUE(column.settling.setVelocity(column.particulate, 2.0 * kMetrePerDay,
                                            column.registry, message))
        << message;
    ASSERT_TRUE(column.settling.setBenthicReceiver(column.particulate, column.sediment,
                                                   column.registry, message))
        << message;

    const double initial = column.waterMass(column.particulate);
    ASSERT_GT(initial, 0.0);

    for (int i = 0; i < steps; ++i)
    {
        column.settling.apply(column.layered, column.state, dt);

        // Checked every step. A drift that cancelled at the final time would
        // otherwise pass.
        const double total = column.waterMass(column.particulate) + column.bedMass();

        EXPECT_NEAR(total, initial, 1.0e-12 * initial) << "after step " << i;
    }

    // And something actually moved, or conservation is trivial.
    EXPECT_GT(column.bedMass(), 0.05 * initial);
    EXPECT_LT(column.waterMass(column.particulate), 0.95 * initial);

    // The tally agrees with the store it filled.
    EXPECT_NEAR(column.settling.deliveredToBed(column.particulate), column.bedMass(),
                1.0e-12 * initial);
    EXPECT_EQ(column.settling.lostToBed(column.particulate), 0.0);
}

TEST(Settling, WithNoReceiverTheLossIsReportedRatherThanHidden)
{
    // A model may legitimately not track sediment. What it may not do is lose
    // mass without saying so — an unreported loss is indistinguishable from a
    // conservation bug, and this budget is the engine's to keep.
    constexpr double c0 = 5.0;
    constexpr double dt = 3600.0;

    Column column;
    column.build(6, 12.0, /*withSediment=*/false);
    column.setUniform(column.particulate, c0);

    std::string message;
    ASSERT_TRUE(column.settling.setVelocity(column.particulate, 5.0 * kMetrePerDay,
                                            column.registry, message))
        << message;

    const double initial = column.waterMass(column.particulate);

    for (int i = 0; i < 20; ++i)
        column.settling.apply(column.layered, column.state, dt);

    const double remaining = column.waterMass(column.particulate);
    const double lost      = column.settling.lostToBed(column.particulate);

    EXPECT_GT(lost, 0.0);
    EXPECT_NEAR(initial - remaining, lost, 1.0e-12 * initial)
        << "the reported loss does not account for the missing mass";
}

// ── Stability: the reason the solve is implicit ────────────────────────────

TEST(Settling, PositiveAndConservativeFarBeyondTheExplicitStabilityLimit)
{
    // The gate that earns the implicit solve. The step here is a hundred times
    // the explicit limit `h/w`, which is the situation a thin surface sigma layer
    // creates without anyone choosing it. An explicit update would produce large
    // negative concentrations; this must stay positive and closed.
    constexpr double c0         = 8.0;
    constexpr double totalDepth = 10.0;
    constexpr int    layers     = 10;
    constexpr double w          = 20.0 * kMetrePerDay;

    Column column;
    column.build(layers, totalDepth);
    column.setUniform(column.particulate, c0);

    const double h             = totalDepth / static_cast<double>(layers);
    const double explicitLimit = h / w;
    const double dt            = 100.0 * explicitLimit;

    std::string message;
    ASSERT_TRUE(column.settling.setVelocity(column.particulate, w, column.registry,
                                            message))
        << message;
    ASSERT_TRUE(column.settling.setBenthicReceiver(column.particulate, column.sediment,
                                                   column.registry, message))
        << message;

    const double initial = column.waterMass(column.particulate);

    for (int i = 0; i < 5; ++i)
    {
        column.settling.apply(column.layered, column.state, dt);

        for (int64_t cell = 0; cell < column.layered.cellCount(); ++cell)
        {
            ASSERT_GE(column.state.speciesAt(column.particulate, cell), 0.0)
                << "negative concentration in cell " << cell << " at step " << i
                << ", Courant number " << w * dt / h;
        }

        EXPECT_NEAR(column.waterMass(column.particulate) + column.bedMass(), initial,
                    1.0e-12 * initial);
    }

    // At this Courant number the column should be nearly swept out, so the gate
    // is exercising the regime it claims to.
    EXPECT_GT(column.bedMass(), 0.9 * initial);
}

TEST(Settling, ASingleLayerFollowsTheSchemesClosedForm)
{
    // One layer reduces the substitution to a scalar recurrence with an exact
    // solution, `C0 / (1 + w dt/h)^N`. This tests the arithmetic of the update
    // rather than the physics, which is worth separating: an error in the
    // coefficient would otherwise hide inside a plausible-looking profile.
    constexpr double c0         = 20.0;
    constexpr double totalDepth = 4.0;
    constexpr double w          = 6.0 * kMetrePerDay;
    constexpr double dt         = 7200.0;
    constexpr int    steps      = 12;

    Column column;
    column.build(1, totalDepth);
    column.setUniform(column.particulate, c0);

    std::string message;
    ASSERT_TRUE(column.settling.setVelocity(column.particulate, w, column.registry,
                                            message))
        << message;
    ASSERT_TRUE(column.settling.setBenthicReceiver(column.particulate, column.sediment,
                                                   column.registry, message))
        << message;

    for (int i = 0; i < steps; ++i)
        column.settling.apply(column.layered, column.state, dt);

    const double courant = w * dt / totalDepth;
    const double exact   = c0 / std::pow(1.0 + courant, steps);

    EXPECT_NEAR(column.at(column.particulate, 0), exact, 1.0e-12 * c0);
}

// ── The velocity means what it says ───────────────────────────────────────

TEST(Settling, TheProfileTranslatesAtTheSettlingVelocity)
{
    // The physical claim. A blob's centroid should move down by `w t`, and
    // first-order upwinding approaches that as the discretization refines. The
    // rate is measured, not asserted: if it came out zeroth order the velocity
    // would be wrong by a constant factor and a fixed-tolerance test could
    // still pass.
    //
    // Two things this fixture got wrong at first, both mine. The blob travelled
    // 40 m in a 20 m column, so what was being measured was the centroid of the
    // numerical remnant after the material had left. And the refinement halved
    // `h` while holding `dt` fixed — implicit upwind carries an O(w dt) error
    // alongside its O(w h) one, so the constant term dominated and the ratio
    // came out 1.005. The step is now refined with the grid at a fixed Courant
    // number, which is the only way a first-order rate is observable.
    constexpr double totalDepth = 20.0;
    constexpr double w          = 10.0 * kMetrePerDay;
    constexpr double elapsed    = 0.5 * 86400.0; // 5 m of travel
    constexpr double centre     = 6.0;           // initial centroid depth
    constexpr double spread     = 2.0;
    constexpr double courant    = 0.25;

    double previousError = 0.0;
    double ratio         = 0.0;

    for (int layers : {40, 80, 160})
    {
        Column column;
        column.build(layers, totalDepth, /*withSediment=*/true);

        const double h = totalDepth / static_cast<double>(layers);

        // A Gaussian, kept well clear of both the surface and the bed so the
        // boundaries do not clip it over the interval.
        for (int k = 0; k < layers; ++k)
        {
            const double z = (static_cast<double>(k) + 0.5) * h;
            const double x = (z - centre) / spread;

            column.state.speciesAt(column.particulate, column.layered.cell(0, k)) =
                std::exp(-0.5 * x * x);
        }

        std::string message;
        ASSERT_TRUE(column.settling.setVelocity(column.particulate, w, column.registry,
                                                message))
            << message;
        ASSERT_TRUE(column.settling.setBenthicReceiver(column.particulate,
                                                       column.sediment,
                                                       column.registry, message))
            << message;

        const double dt    = courant * h / w;
        const int    steps = static_cast<int>(std::lround(elapsed / dt));

        // The refinement is only clean if the interval divides exactly.
        ASSERT_NEAR(static_cast<double>(steps) * dt, elapsed, 1.0e-9 * elapsed);

        for (int i = 0; i < steps; ++i)
            column.settling.apply(column.layered, column.state, dt);

        // L1 error against the exact solution of pure downward advection, which
        // is the initial profile translated by w t. Measured this way rather
        // than as a centroid displacement: the centroid of what *remains* is
        // biased by whatever has left through the bed, so it converges to the
        // wrong thing near the boundary. The blob is kept several standard
        // deviations clear of the bed so the translated profile is exact here.
        double error = 0.0;
        double norm  = 0.0;

        for (int k = 0; k < layers; ++k)
        {
            const double z = (static_cast<double>(k) + 0.5) * h;
            const double x = (z - (centre + w * elapsed)) / spread;

            const double exact = std::exp(-0.5 * x * x);

            error += std::abs(column.at(column.particulate, k) - exact) * h;
            norm += exact * h;
        }

        ASSERT_GT(norm, 0.0);

        const double relative = error / norm;

        if (previousError > 0.0)
            ratio = previousError / relative;

        previousError = relative;
    }

    // First order: halving the discretization halves the error.
    EXPECT_GT(ratio, 1.6) << "convergence ratio " << ratio;
    EXPECT_LT(ratio, 2.6) << "convergence ratio " << ratio;
}

TEST(Settling, MaterialMovesDownwardAndNotUpward)
{
    // Direction, pinned separately because a sign error in the sweep would
    // produce a perfectly conservative, perfectly stable, upside-down model.
    Column column;
    column.build(6, 12.0);

    // All the material in the middle, nothing above or below it.
    column.state.speciesAt(column.particulate, column.layered.cell(0, 2)) = 10.0;

    std::string message;
    ASSERT_TRUE(column.settling.setVelocity(column.particulate, 20.0 * kMetrePerDay,
                                            column.registry, message))
        << message;
    ASSERT_TRUE(column.settling.setBenthicReceiver(column.particulate, column.sediment,
                                                   column.registry, message))
        << message;

    column.settling.apply(column.layered, column.state, 7200.0);

    // Layers above the blob must still be empty: upwinding from above means the
    // surface layer has no source at all.
    EXPECT_EQ(column.at(column.particulate, 0), 0.0);
    EXPECT_EQ(column.at(column.particulate, 1), 0.0);

    // And the layer below must have received some.
    EXPECT_GT(column.at(column.particulate, 3), 0.0);
    EXPECT_LT(column.at(column.particulate, 2), 10.0);
}

// ── Through the full solver ───────────────────────────────────────────────

TEST(SettlingSolver, SettlingRunsInTheStepAndKeepsTheBudgetClosed)
{
    // The wiring gate. Everything above exercises SettlingSet directly; this
    // checks the solver actually calls it, times it, and that the water-plus-
    // sediment budget still closes when real transport is moving the field
    // around underneath.
    constexpr double c0         = 6.0;
    constexpr double totalDepth = 12.0;
    constexpr int    layers     = 8;
    constexpr double dt         = 60.0;
    constexpr int    steps      = 25;

    HorizontalMesh     horizontal;
    LayeredMesh        layered;
    SpeciesRegistry    registry;
    StateStore         state;
    InternalModeSolver solver;

    std::string message;

    ASSERT_TRUE(horizontal.build(Testing::channelStrip(20, 400.0, 40.0), message))
        << message;

    for (int64_t f = 0; f < horizontal.faceCount(); ++f)
        horizontal.faceBed()[static_cast<size_t>(f)] = -totalDepth;

    SigmaSpec spec;
    spec.layerCount        = layers;
    spec.minLayerThickness = 0.0;

    ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

    ASSERT_GE(registry.addTemperature(10.0, message), 0) << message;
    const int temperature = registry.temperatureIndex();

    SpeciesDefinition solids;
    solids.id    = "TSS";
    solids.units = "g/m3";
    ASSERT_GE(registry.add(solids, message), 0) << message;
    const int particulate = registry.indexOf("TSS");

    SpeciesDefinition bed;
    bed.id    = "TSS_BED";
    bed.kind  = SpeciesKind::MsxBenthic;
    bed.units = "g/m2";
    ASSERT_GE(registry.add(bed, message), 0) << message;
    const int sediment = registry.indexOf("TSS_BED");

    state.allocate(layered, registry);

    // A temperature front, so the basin is stirring and settling has to stay
    // conservative while advection moves material between columns.
    for (int64_t c = 0; c < layered.columnCount(); ++c)
    {
        const double x = horizontal.faceCentroidX()[static_cast<size_t>(c)];

        for (int k = 0; k < layers; ++k)
        {
            const int64_t cell = layered.cell(c, k);

            state.speciesAt(temperature, cell)  = x < 200.0 ? 8.0 : 18.0;
            state.speciesAt(particulate, cell)  = c0;
        }
    }

    InternalModeOptions options;
    options.gravity         = 9.80665;
    options.manningN        = 0.0;
    options.cflNumber       = 0.5;
    options.surfaceHeatFlux = 0.0;

    ASSERT_TRUE(solver.initialize(layered, state, registry, options, message)) << message;

    ASSERT_TRUE(solver.settling().setVelocity(particulate, 30.0 * kMetrePerDay, registry,
                                              message))
        << message;
    ASSERT_TRUE(solver.settling().setBenthicReceiver(particulate, sediment, registry,
                                                     message))
        << message;

    const auto totalMass = [&]() {
        double water = 0.0;

        for (int64_t cell = 0; cell < layered.cellCount(); ++cell)
        {
            water += state.speciesAt(particulate, cell) *
                     layered.cellVolume()[static_cast<size_t>(cell)];
        }

        double onBed = 0.0;

        for (int64_t c = 0; c < layered.columnCount(); ++c)
        {
            onBed += state.speciesAt(sediment, layered.cell(c, layers - 1)) *
                     horizontal.faceArea()[static_cast<size_t>(c)];
        }

        return water + onBed;
    };

    const double initial = totalMass();
    ASSERT_GT(initial, 0.0);

    for (int i = 0; i < steps; ++i)
        ASSERT_TRUE(solver.step(dt, message)) << message;

    // Transport conserves and settling conserves, so the sum must, through a
    // step that is doing both at once.
    EXPECT_NEAR(totalMass(), initial, 1.0e-9 * initial);

    // Settling was actually reached — without this the gate would pass on a
    // step() that never called it.
    EXPECT_GT(solver.profile().settling, 0)
        << "no time was spent in the settling operator";

    EXPECT_GT(solver.settling().deliveredToBed(particulate), 0.01 * initial)
        << "nothing reached the bed";
}
