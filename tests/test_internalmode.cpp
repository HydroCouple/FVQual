// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_internalmode.cpp
 * \brief The Phase 3 verification gate for the 3-D hydrostatic internal mode.
 *
 * These pin the properties whose loss is silent: mode-splitting consistency,
 * exactness of the rest state, and the accuracy of the baroclinic pressure
 * gradient — the quantity risk R1 is about.
 *
 * **What is deliberately not here.** Baroclinic *dynamics* — the internal
 * seiche and lock-exchange gates — are not yet passing and are recorded as an
 * open item rather than as skipped tests. The three
 * rest-state tests below are what localize that failure: the solver holds a
 * homogeneous basin and a stratified basin motionless to round-off, so the
 * coupling machinery is exact and the defect lies in the response to a
 * baroclinic forcing, not in the plumbing.
 */

#include "fvqual/hydro/equationofstate.h"
#include "fvqual/hydro/internalmode.h"
#include "fvqual/hydro/pressureclosure.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <utility>
#include <vector>

using FVQual::Hydro::EquationOfState;
using FVQual::Hydro::HydrostaticPressureClosure;
using FVQual::Hydro::InternalModeOptions;
using FVQual::Hydro::InternalModeSolver;
using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;
using FVQual::State::SpeciesDefinition;
using FVQual::State::SpeciesRegistry;
using FVQual::State::StateStore;
namespace Testing = FVQual::Testing;

namespace {

constexpr double kGravity          = 9.80665;
constexpr double kReferenceDensity = 1000.0;

/*!
 * \brief Builds a channel mesh with a prescribed bed.
 */
template <typename BedFunction>
HorizontalMesh channel(int columns, double length, double width, BedFunction bed)
{
    HorizontalMesh mesh;
    std::string    message;

    EXPECT_TRUE(mesh.build(Testing::channelStrip(columns, length, width), message))
        << message;

    for (int64_t f = 0; f < mesh.faceCount(); ++f)
    {
        mesh.faceBed()[static_cast<size_t>(f)] =
            bed(mesh.faceCentroidX()[static_cast<size_t>(f)]);
    }

    return mesh;
}

/*!
 * \brief A complete stratified channel model.
 */
struct Model
{
    HorizontalMesh     horizontal;
    LayeredMesh        layered;
    SpeciesRegistry    registry;
    StateStore         state;
    InternalModeSolver solver;

    int temperature = -1;
    int tracer      = -1;

    template <typename BedFunction, typename TemperatureFunction>
    void build(int columns, int layers, double length, double width, BedFunction bed,
               TemperatureFunction temp, const InternalModeOptions &options,
               bool withTracer = false)
    {
        std::string message;

        horizontal = channel(columns, length, width, bed);

        SigmaSpec spec;
        spec.layerCount        = layers;
        spec.minLayerThickness = 0.0;

        ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

        ASSERT_GE(registry.addTemperature(10.0, message), 0) << message;
        temperature = registry.temperatureIndex();

        if (withTracer)
        {
            SpeciesDefinition dye;
            dye.id           = "DYE";
            dye.units        = "g/m3";
            dye.initialValue = 3.0;
            ASSERT_GE(registry.add(dye, message), 0) << message;
            tracer = registry.indexOf("DYE");
        }

        state.allocate(layered, registry);

        for (int64_t c = 0; c < layered.columnCount(); ++c)
        {
            const double x = horizontal.faceCentroidX()[static_cast<size_t>(c)];

            for (int k = 0; k < layered.layerCount(); ++k)
            {
                const double z = 0.5 * (layered.z(c, k) + layered.z(c, k + 1));
                state.speciesAt(temperature, layered.cell(c, k)) = temp(x, z);
            }
        }

        ASSERT_TRUE(solver.initialize(layered, state, registry, options, message))
            << message;
    }

    [[nodiscard]] double maximumSpeed() const
    {
        double peak = 0.0;

        for (int64_t i = 0; i < state.cellCount(); ++i)
        {
            peak = std::max(peak, std::hypot(state.u[static_cast<size_t>(i)],
                                             state.v[static_cast<size_t>(i)]));
        }

        return peak;
    }
};

InternalModeOptions quietOptions()
{
    InternalModeOptions options;
    options.gravity   = kGravity;
    options.manningN  = 0.0;
    options.cflNumber = 0.5;
    return options;
}

/*!
 * \brief A horizontally uniform exponential thermocline.
 */
double stratifiedProfile(double, double z)
{
    return 5.0 + 20.0 * std::exp(z / 8.0);
}

/*!
 * \brief Largest baroclinic acceleration a density field produces, evaluated
 *        directly from the closure with no time stepping involved.
 */
double peakBaroclinicAcceleration(const LayeredMesh &mesh,
                                  const std::vector<double> &density,
                                  bool linearReconstruction = true)
{
    HydrostaticPressureClosure closure(
        HydrostaticPressureClosure::Options{linearReconstruction});

    std::vector<double> ax;
    std::vector<double> ay;

    closure.baroclinicAcceleration(mesh, density, kReferenceDensity, kGravity, ax, ay);

    double peak = 0.0;
    for (double value : ax)
        peak = std::max(peak, std::abs(value));

    return peak;
}

/*!
 * \brief Fills a density field from a horizontally uniform T(z) profile.
 */
std::vector<double> uniformStratification(const LayeredMesh &mesh)
{
    std::vector<double> density(static_cast<size_t>(mesh.cellCount()), 0.0);

    for (int64_t c = 0; c < mesh.columnCount(); ++c)
    {
        for (int k = 0; k < mesh.layerCount(); ++k)
        {
            const double z = 0.5 * (mesh.z(c, k) + mesh.z(c, k + 1));
            density[static_cast<size_t>(mesh.cell(c, k))] =
                EquationOfState::density(stratifiedProfile(0.0, z));
        }
    }

    return density;
}

} // namespace

// ── The pressure gradient itself (risk R1) ────────────────────────────────

TEST(PressureClosure, FlatBedUniformStratificationIsExactlyZero)
{
    // The design target: a horizontally uniform ocean must produce no
    // baroclinic force at all. With every face evaluated at the cell's own
    // elevation, the identical face pressures multiply a closed sum of
    // outward normals and cancel exactly — not approximately.
    const HorizontalMesh horizontal = channel(20, 400.0, 30.0, [](double) { return -25.0; });

    SigmaSpec spec;
    spec.layerCount        = 16;
    spec.minLayerThickness = 0.0;

    LayeredMesh layered;
    std::string message;
    ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

    EXPECT_EQ(peakBaroclinicAcceleration(layered, uniformStratification(layered)), 0.0);
}

TEST(PressureClosure, SlopingBedErrorConvergesUnderVerticalRefinement)
{
    // Over a sloping bed the cancellation is no longer exact: two columns
    // represent the same profile on layers at different elevations. What
    // matters is that the residual converges, and quickly.
    const auto bed = [](double x) { return -30.0 + 10.0 * std::sin(x / 200.0); };

    double previous = 0.0;

    for (int layers : {10, 20, 40})
    {
        const HorizontalMesh horizontal = channel(40, 600.0, 50.0, bed);

        SigmaSpec spec;
        spec.layerCount        = layers;
        spec.minLayerThickness = 0.0;

        LayeredMesh layered;
        std::string message;
        ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

        const double error =
            peakBaroclinicAcceleration(layered, uniformStratification(layered));

        EXPECT_LT(error, 1.0e-4) << "at " << layers << " layers";

        if (previous > 0.0)
        {
            EXPECT_LT(error, 0.5 * previous)
                << "error did not at least halve when refining to " << layers
                << " layers";
        }

        previous = error;
    }
}

TEST(PressureClosure, TiltedSurfaceOverHomogeneousWaterIsExactlyZero)
{
    // Regression test for a bug that cost most of a development session.
    //
    // Homogeneous water under a tilted free surface is a purely *barotropic*
    // situation: the entire pressure gradient is the surface slope, which
    // belongs to the external mode. The baroclinic closure must return
    // nothing at all.
    //
    // Integrating the full density rather than the anomaly returns
    // rho0 g d(eta)/dx here instead of zero — the external mode then receives
    // a duplicate of its own restoring force, lagged across a subcycling
    // interval, and the free surface grows without bound. Every rest-state
    // test still passed, because a flat surface hides it completely. Only a
    // tilt exposes it.
    const HorizontalMesh horizontal = channel(20, 400.0, 30.0, [](double) { return -25.0; });

    SigmaSpec spec;
    spec.layerCount        = 12;
    spec.minLayerThickness = 0.0;

    LayeredMesh layered;
    std::string message;
    ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

    // A pronounced surface tilt.
    std::vector<double> surface(static_cast<size_t>(horizontal.faceCount()));
    for (int64_t c = 0; c < horizontal.faceCount(); ++c)
    {
        surface[static_cast<size_t>(c)] =
            0.5 * (horizontal.faceCentroidX()[static_cast<size_t>(c)] / 400.0 - 0.5);
    }

    ASSERT_TRUE(layered.updateSurface(surface, message)) << message;

    const double slope = 0.5 / 400.0;

    // Water exactly at the reference density: the anomaly is identically zero,
    // so the answer must be identically zero. With the bug this returns
    // g * d(eta)/dx — about 0.012 m/s², nine orders of magnitude adrift.
    const std::vector<double> atReference(static_cast<size_t>(layered.cellCount()),
                                          kReferenceDensity);

    EXPECT_EQ(peakBaroclinicAcceleration(layered, atReference), 0.0);

    // Water at a realistic density leaves a residual, and it is worth knowing
    // exactly what that residual is: the Boussinesq split approximates
    // 1/rho by 1/rho0, so a uniform anomaly rho' contributes
    // g (rho'/rho0) d(eta)/dx. That is the approximation's own error, not a
    // discretization artefact, and it is three orders of magnitude below the
    // barotropic term it sits beside.
    const double density = EquationOfState::density(12.0);
    const std::vector<double> realistic(static_cast<size_t>(layered.cellCount()), density);

    const double residual = peakBaroclinicAcceleration(layered, realistic);
    const double expected =
        kGravity * std::abs(density - kReferenceDensity) / kReferenceDensity * slope;

    EXPECT_NEAR(residual, expected, 0.02 * expected);
    EXPECT_LT(residual, 1.0e-3 * kGravity * slope);
}

TEST(PressureClosure, RespondsToAHorizontalDensityContrast)
{
    // The complement of the tests above: having established that the closure
    // is silent when it should be, confirm it is not simply always silent.
    const HorizontalMesh horizontal = channel(20, 400.0, 30.0, [](double) { return -20.0; });

    SigmaSpec spec;
    spec.layerCount        = 8;
    spec.minLayerThickness = 0.0;

    LayeredMesh layered;
    std::string message;
    ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

    std::vector<double> density(static_cast<size_t>(layered.cellCount()), 0.0);

    for (int64_t c = 0; c < layered.columnCount(); ++c)
    {
        const double x = horizontal.faceCentroidX()[static_cast<size_t>(c)];

        for (int k = 0; k < layered.layerCount(); ++k)
        {
            density[static_cast<size_t>(layered.cell(c, k))] =
                EquationOfState::density(x < 200.0 ? 5.0 : 25.0);
        }
    }

    // A 20 °C front over one cell is a strong forcing; the response should be
    // of order g/rho0 * d(rho)/dx * depth.
    EXPECT_GT(peakBaroclinicAcceleration(layered, density), 1.0e-3);
}

// ── Rest states: the mode-splitting machinery ─────────────────────────────

TEST(InternalMode, HomogeneousBasinOverSteepBedStaysAtRest)
{
    // Homogeneous water has no baroclinic forcing whatsoever, so every part
    // of the splitting — the body force handed to the external mode, the mode
    // adjustment, the flux adjustment, the continuity integration — must
    // conspire to produce exactly nothing. This is the test that localizes
    // any failure of baroclinic dynamics to the forcing response rather than
    // to the plumbing.
    Model model;
    model.build(
        20, 10, 600.0, 50.0,
        [](double x) { return -10.0 - 40.0 / (1.0 + std::exp(-(x - 300.0) / 25.0)); },
        [](double, double) { return 10.0; }, quietOptions());

    std::string message;
    ASSERT_TRUE(model.solver.advance(200.0, message)) << message;

    EXPECT_LT(model.maximumSpeed(), 1.0e-9);
    EXPECT_LT(model.solver.diagnostics().surfaceFluxResidual, 1.0e-10);
}

TEST(InternalMode, StratifiedBasinOnAFlatBedStaysAtRest)
{
    Model model;
    model.build(20, 16, 600.0, 50.0, [](double) { return -20.0; }, stratifiedProfile,
                quietOptions());

    const auto range = [&model]() {
        double coldest = 1.0e30;
        double warmest = -1.0e30;

        for (int64_t i = 0; i < model.state.cellCount(); ++i)
        {
            const double value = model.state.speciesAt(model.temperature, i);
            coldest            = std::min(coldest, value);
            warmest            = std::max(warmest, value);
        }

        return std::pair{coldest, warmest};
    };

    const auto before = range();

    std::string message;
    ASSERT_TRUE(model.solver.advance(200.0, message)) << message;

    EXPECT_LT(model.maximumSpeed(), 1.0e-9);

    // The stratification must survive: with no motion, the only thing acting
    // on it is the background vertical diffusivity, which over this interval
    // is worth a few hundredths of a millikelvin. Anything larger would be
    // spurious mixing.
    const auto after = range();

    EXPECT_NEAR(after.second, before.second, 1.0e-3);
    EXPECT_NEAR(after.first, before.first, 1.0e-3);
    EXPECT_GT(after.first, 4.99);
}

// ── Mode-splitting consistency ────────────────────────────────────────────

TEST(InternalMode, SurfaceContinuityResidualVanishes)
{
    // The vertical velocity is integrated upward from a no-flux bed, so its
    // value at the free surface is a pure diagnostic of whether the layer
    // fluxes agree with the transport the external mode delivered. A
    // non-negligible residual means the two modes have drifted apart and
    // scalar transport can no longer be trusted, however plausible the fields
    // still look.
    InternalModeOptions options = quietOptions();
    options.windStressX         = 0.5; // stir it, so the test is not trivial

    Model model;
    model.build(20, 10, 600.0, 50.0, [](double) { return -15.0; }, stratifiedProfile,
                options);

    std::string message;

    for (int i = 0; i < 20; ++i)
    {
        ASSERT_TRUE(model.solver.step(2.0, message)) << message;
        EXPECT_LT(model.solver.diagnostics().surfaceFluxResidual, 1.0e-10)
            << "at step " << i;
    }

    // The forcing must actually have done something.
    EXPECT_GT(model.maximumSpeed(), 1.0e-5);
}

TEST(InternalMode, UniformTracerStaysUniformUnderForcing)
{
    // Conservative transport divided by the new cell volume reproduces a
    // uniform field only if the volume change equals the net flux exactly.
    // Any inconsistency in the splitting shows up here first.
    InternalModeOptions options = quietOptions();
    options.windStressX         = 0.5;

    Model model;
    model.build(16, 8, 400.0, 40.0, [](double) { return -12.0; }, stratifiedProfile,
                options, /*withTracer=*/true);

    ASSERT_GE(model.tracer, 0);

    std::string message;
    ASSERT_TRUE(model.solver.advance(40.0, message)) << message;

    for (int64_t i = 0; i < model.state.cellCount(); ++i)
    {
        EXPECT_NEAR(model.state.speciesAt(model.tracer, i), 3.0, 1.0e-9)
            << "cell " << i;
    }
}

TEST(InternalMode, TracerMassIsConserved)
{
    InternalModeOptions options = quietOptions();
    options.windStressX         = 0.5;

    Model model;
    model.build(16, 8, 400.0, 40.0, [](double) { return -12.0; }, stratifiedProfile,
                options, /*withTracer=*/true);

    // A non-uniform tracer, so conservation is not trivially satisfied.
    for (int64_t c = 0; c < model.layered.columnCount(); ++c)
    {
        for (int k = 0; k < model.layered.layerCount(); ++k)
        {
            model.state.speciesAt(model.tracer, model.layered.cell(c, k)) =
                c < model.layered.columnCount() / 2 ? 1.0 : 5.0;
        }
    }

    const double before = model.solver.totalSpeciesMass(model.tracer);
    ASSERT_GT(before, 0.0);

    std::string message;
    ASSERT_TRUE(model.solver.advance(40.0, message)) << message;

    const double after = model.solver.totalSpeciesMass(model.tracer);

    EXPECT_LT(std::abs(after - before) / before, 1.0e-10)
        << "tracer mass drifted by " << (after - before);
}

TEST(InternalMode, WaterVolumeIsConserved)
{
    Model model;
    model.build(20, 10, 600.0, 50.0, [](double) { return -15.0; }, stratifiedProfile,
                quietOptions());

    const double before = model.layered.totalVolume();

    std::string message;
    ASSERT_TRUE(model.solver.advance(100.0, message)) << message;

    EXPECT_LT(std::abs(model.layered.totalVolume() - before) / before, 1.0e-10);
}

// ── Baroclinic dynamics ───────────────────────────────────────────────────

TEST(InternalMode, QuiescentStratifiedBasinOverSteepBathymetryStaysNearlyAtRest)
{
    // The risk R1 gate. Horizontally uniform stratification over steeply
    // sloping bathymetry has an exact solution — no motion — and a sigma
    // coordinate cannot represent it exactly, so truncation error appears as
    // spurious currents. What must not happen is growth: the coupled system
    // has to stay bounded.
    Model model;
    model.build(
        30, 16, 600.0, 50.0,
        [](double x) { return -10.0 - 40.0 / (1.0 + std::exp(-(x - 300.0) / 25.0)); },
        stratifiedProfile, quietOptions());

    // The bathymetry passes the conventional Haney screening, which is worth
    // recording: the gate is necessary, not sufficient.
    EXPECT_TRUE(model.layered.haneyReport(0.2).passed());

    std::string message;
    ASSERT_TRUE(model.solver.advance(600.0, message)) << message;

    // Measured at ~1.1e-2 m/s for this configuration. The bound is what the
    // gate defends: a runaway would exceed it by orders of magnitude.
    EXPECT_LT(model.maximumSpeed(), 5.0e-2)
        << "spurious velocity " << model.maximumSpeed() << " m/s";

    EXPECT_LT(model.solver.diagnostics().surfaceFluxResidual, 1.0e-10);
}

TEST(InternalMode, TwoLayerSeichePeriodMatchesLinearTheory)
{
    // A tilted interface between two layers oscillates at
    //     T = 2 L / sqrt(g' h1 h2 / H),   g' = g (rho2 - rho1) / rho0.
    // The internal wave is a thousand times slower than the surface wave the
    // external mode carries, so getting this period right is the clearest
    // evidence that mode splitting separates the two correctly.
    constexpr double kLength = 1000.0;
    constexpr double kDepth  = 20.0;
    constexpr double kUpper  = 5.0;
    constexpr double kTilt   = 1.5;

    const double light   = EquationOfState::density(25.0);
    const double heavy   = EquationOfState::density(5.0);
    const double reduced = kGravity * (heavy - light) / kReferenceDensity;
    const double speed   = std::sqrt(reduced * kUpper * (kDepth - kUpper) / kDepth);
    const double period  = 2.0 * kLength / speed;

    // Resolution matters to this measurement: the tilted interface must span
    // enough layers to be represented as a slope rather than a staircase.
    // Twenty layers over 20 m puts about three layers across the tilt.
    Model model;
    model.build(40, 20, kLength, 50.0, [](double) { return -kDepth; },
                [&](double x, double z) {
                    const double interface =
                        -kUpper + kTilt * (x / kLength - 0.5) * 2.0;
                    return z > interface ? 25.0 : 5.0;
                },
                quietOptions());

    // Track the depth of the 15 °C isotherm at the closed end.
    const auto interfaceDepth = [&]() {
        for (int k = 0; k < model.layered.layerCount(); ++k)
        {
            if (model.state.speciesAt(model.temperature, model.layered.cell(0, k)) < 15.0)
                return 0.5 * (model.layered.z(0, k) + model.layered.z(0, k + 1));
        }

        return model.layered.z(0, model.layered.layerCount());
    };

    std::vector<double> crossings;
    double              previous     = interfaceDepth() + kUpper;
    double              previousTime = 0.0;

    constexpr double kSample = 20.0;
    std::string      message;

    while (model.solver.time() < 2.2 * period)
    {
        ASSERT_TRUE(model.solver.advance(kSample, message)) << message;

        const double current = interfaceDepth() + kUpper;

        if ((previous > 0.0 && current <= 0.0) || (previous < 0.0 && current >= 0.0))
        {
            const double fraction = previous / (previous - current);
            crossings.push_back(previousTime + fraction * kSample);
        }

        previous     = current;
        previousTime = model.solver.time();
    }

    ASSERT_GE(crossings.size(), 2u) << "the interface never oscillated";

    double sum = 0.0;
    for (size_t i = 1; i < crossings.size(); ++i)
        sum += crossings[i] - crossings[i - 1];

    const double measured = 2.0 * sum / static_cast<double>(crossings.size() - 1);

    // Numerical diffusion thickens the interface, which lowers the effective
    // internal wave speed and lengthens the period. With MUSCL + flux
    // correction the bias is about +8%; it was +19% under first-order upwind,
    // which is the clearest single measurement of what the transport scheme
    // buys.
    EXPECT_NEAR(measured, period, 0.15 * period)
        << "measured " << measured << " s vs analytic " << period << " s";

    // The velocity scale must be physical, not merely oscillatory.
    EXPECT_LT(model.maximumSpeed(), 0.5);
}

TEST(InternalMode, LockExchangeFrontSpeedMatchesGravityCurrentTheory)
{
    // Release dense water against light water over a flat bed: the dense
    // front runs along the bed at close to 0.5 sqrt(g' H).
    constexpr double kLength = 200.0;
    constexpr double kDepth  = 10.0;
    constexpr double kEnd    = 60.0;

    const double light   = EquationOfState::density(25.0);
    const double heavy   = EquationOfState::density(5.0);
    const double reduced = kGravity * (heavy - light) / kReferenceDensity;
    const double analytic = 0.5 * std::sqrt(reduced * kDepth);

    Model model;
    model.build(80, 10, kLength, 20.0, [](double) { return -kDepth; },
                [](double x, double) { return x < 100.0 ? 5.0 : 25.0; }, quietOptions());

    std::string message;
    ASSERT_TRUE(model.solver.advance(kEnd, message)) << message;

    // Furthest point the dense water has reached along the bed.
    double front = 0.5 * kLength;

    for (int64_t c = 0; c < model.layered.columnCount(); ++c)
    {
        const double bottom = model.state.speciesAt(
            model.temperature, model.layered.cell(c, model.layered.layerCount() - 1));

        if (bottom < 15.0)
            front = std::max(front, model.horizontal.faceCentroidX()[static_cast<size_t>(c)]);
    }

    const double measured = (front - 0.5 * kLength) / kEnd;

    // Measured at about -19% of theory: a first-order scheme smears the
    // density interface, and a smeared front runs slow. Same character as the
    // Ritter dry-front lag in the external mode.
    EXPECT_NEAR(measured, analytic, 0.35 * analytic)
        << "measured " << measured << " m/s vs theory " << analytic << " m/s";

    // Velocities must stay on the gravity-current scale.
    EXPECT_LT(model.maximumSpeed(), 3.0 * analytic);
    EXPECT_LT(model.solver.diagnostics().surfaceFluxResidual, 1.0e-10);
}

// ── Vertical mixing ───────────────────────────────────────────────────────

TEST(InternalMode, StratificationSuppressesVerticalMixing)
{
    // Pacanowski-Philander damping is what lets a thermocline survive: mixing
    // must collapse where the Richardson number is large.
    Model stratified;
    stratified.build(10, 12, 300.0, 30.0, [](double) { return -20.0; },
                     stratifiedProfile, quietOptions());

    Model homogeneous;
    homogeneous.build(10, 12, 300.0, 30.0, [](double) { return -20.0; },
                      [](double, double) { return 10.0; }, quietOptions());

    std::string message;
    ASSERT_TRUE(stratified.solver.step(1.0, message)) << message;
    ASSERT_TRUE(homogeneous.solver.step(1.0, message)) << message;

    const auto interiorPeak = [](const Model &model) {
        double peak = 0.0;

        for (int64_t c = 0; c < model.layered.columnCount(); ++c)
        {
            for (int k = 1; k < model.layered.layerCount(); ++k)
            {
                peak = std::max(peak,
                                model.solver.eddyDiffusivity()[static_cast<size_t>(
                                    model.layered.interfaceSlot(c, k))]);
            }
        }

        return peak;
    };

    EXPECT_LT(interiorPeak(stratified), interiorPeak(homogeneous));
    EXPECT_GT(stratified.solver.diagnostics().maxRichardson, 1.0);
}

// ── Contracts ─────────────────────────────────────────────────────────────

TEST(InternalMode, DefaultClosureIsHydrostaticAndSwappable)
{
    Model model;
    model.build(8, 6, 200.0, 20.0, [](double) { return -10.0; }, stratifiedProfile,
                quietOptions());

    EXPECT_FALSE(model.solver.pressureClosure().isNonHydrostatic());
    EXPECT_NE(model.solver.pressureClosure().name().find("hydrostatic"),
              std::string::npos);

    // The seam: swapping the closure is how a
    // non-hydrostatic correction arrives without touching momentum assembly.
    model.solver.setPressureClosure(std::make_unique<HydrostaticPressureClosure>(
        HydrostaticPressureClosure::Options{false}));

    EXPECT_NE(model.solver.pressureClosure().name().find("piecewise"),
              std::string::npos);
}

TEST(InternalMode, TimeStepAccountsForHorizontalDensityContrast)
{
    // A lock exchange begins with every column internally uniform, so a
    // column-local density range sees no baroclinic wave and hands back the
    // maximum step — which is unstable from the very first step. The estimate
    // must therefore use the domain-wide contrast.
    Model lock;
    lock.build(20, 8, 200.0, 20.0, [](double) { return -10.0; },
               [](double x, double) { return x < 100.0 ? 5.0 : 25.0; }, quietOptions());

    const double step = lock.solver.stableTimeStep();

    EXPECT_GT(step, 0.0);
    EXPECT_LT(step, 30.0) << "step " << step
                          << " s ignores the horizontal density contrast";
}

TEST(InternalMode, AdvanceLandsExactlyOnTheTargetTime)
{
    Model model;
    model.build(10, 8, 300.0, 20.0, [](double) { return -12.0; }, stratifiedProfile,
                quietOptions());

    std::string message;

    ASSERT_TRUE(model.solver.advance(60.0, message)) << message;
    EXPECT_DOUBLE_EQ(model.solver.time(), 60.0);

    ASSERT_TRUE(model.solver.advance(0.0, message)) << message;
    EXPECT_DOUBLE_EQ(model.solver.time(), 60.0);

    EXPECT_FALSE(model.solver.step(0.0, message));
    EXPECT_FALSE(model.solver.step(-1.0, message));
}

TEST(InternalMode, RejectsAMismatchedStateStore)
{
    const HorizontalMesh horizontal =
        channel(6, 60.0, 10.0, [](double) { return -5.0; });

    SigmaSpec spec;
    spec.layerCount        = 4;
    spec.minLayerThickness = 0.0;

    LayeredMesh layered;
    std::string message;
    ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

    SigmaSpec other = spec;
    other.layerCount = 7;

    LayeredMesh mismatched;
    ASSERT_TRUE(mismatched.build(horizontal, other, 0.0, message)) << message;

    SpeciesRegistry registry;
    ASSERT_GE(registry.addTemperature(10.0, message), 0) << message;

    StateStore state;
    state.allocate(mismatched, registry);

    InternalModeSolver solver;
    EXPECT_FALSE(solver.initialize(layered, state, registry, quietOptions(), message));
    EXPECT_FALSE(message.empty());
}

// ── Wind on a stratified, tapering reservoir (the DeGray failure) ──────────

/*!
 * \brief The gust that breaks DeGray, in a controlled coupled model.
 *
 * DeGray fails at day 96.75 and nowhere else. The met record there reads 3.4,
 * **12.1**, 2.3, 0.0 m/s at three-hour spacing — one isolated gust, a 12.7-fold
 * jump in stress and back within nine hours — and a shallow column drains to its
 * bed. Which resolutions survive reshuffles when the formulation changes, so it
 * is a knife edge rather than under-resolution.
 *
 * The barotropic solver is not the culprit: `Barotropic.WindSetupIsReachedFrom
 * RestOverAShoalingBed` and `Barotropic.AnImpulsiveGustOverATaperingBasinDoes
 * NotEmptyIt` put exactly this forcing on exactly this geometry — shoaling bed,
 * storage curve tapering to a tenth of plan area, the same gust sequence through
 * the same setWindStress path — and both pass.
 *
 * What those lack is stratification and the mode coupling. This adds them.
 */
TEST(InternalMode, AGustOverAStratifiedTaperingReservoirDoesNotDrainAColumn)
{
    constexpr int    kColumns = 24;
    constexpr int    kLayers  = 20;
    constexpr double kLength  = 6000.0;
    constexpr double kDeep    = 30.0;
    constexpr double kShallow = 5.0;

    InternalModeOptions options = quietOptions();
    options.manningN    = 0.03;
    options.maxTimeStep = 300.0;

    Model model;

    ASSERT_NO_FATAL_FAILURE(model.build(
        kColumns, kLayers, kLength, 300.0,
        [](double x) {
            const double f = x / kLength;
            return -(kDeep + (kShallow - kDeep) * f);
        },
        stratifiedProfile, options));

    // Storage tapering to a tenth of plan area at the bed — DeGray's shape, and
    // the ingredient the constant-area gates cannot have.
    {
        std::string message;

        for (int64_t c = 0; c < model.horizontal.faceCount(); ++c)
        {
            const size_t i   = static_cast<size_t>(c);
            const double bed = model.horizontal.faceBed()[i];

            std::vector<double> elevation;
            std::vector<double> planArea;

            constexpr int kSamples = 12;

            for (int s = 0; s < kSamples; ++s)
            {
                const double f = static_cast<double>(s) / (kSamples - 1);

                elevation.push_back(bed + f * (1.0 - bed));
                planArea.push_back(model.horizontal.faceArea()[i] *
                                   (0.1 + 0.9 * f));
            }

            ASSERT_TRUE(model.layered.setHypsography(c, elevation, planArea, message))
                << message;
        }
    }

    model.solver.enableMeteorology(/*attenuation=*/0.45);

    // No heat exchange: this is about momentum. Only the wind varies.
    FVQual::Hydro::Meteorology weather;
    weather.airTemperature = 15.0;
    weather.dewPoint       = 15.0;
    weather.cloudCover     = 1.0;
    weather.shortwave      = 0.0;
    weather.windDirection  = 0.0;

    const auto blow = [&](double speed, double seconds) {
        weather.windSpeed          = speed;
        model.solver.meteorology() = weather;

        std::string  message;
        const double until = seconds;

        for (double elapsed = 0.0; elapsed < until;)
        {
            const double dt = std::min(model.solver.stableTimeStep(), 300.0);

            ASSERT_TRUE(model.solver.step(dt, message))
                << "at " << speed << " m/s: " << message;

            elapsed += dt;

            for (int64_t c = 0; c < model.layered.columnCount(); ++c)
            {
                const double depth = model.layered.surfaceElevation()
                                         [static_cast<size_t>(c)] -
                                     model.horizontal.faceBed()[static_cast<size_t>(c)];

                ASSERT_GT(depth, 0.2)
                    << "column " << c << " drained to " << depth << " m under "
                    << speed << " m/s";
            }
        }
    };

    // DeGray's day 96.75 sequence, at its three-hour met spacing.
    blow(3.4, 10800.0);
    blow(12.1, 10800.0);
    blow(2.3, 10800.0);
    blow(0.0, 21600.0);
}
