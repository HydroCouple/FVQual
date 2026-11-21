// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_lightfield.cpp
 * \brief Verification of the Beer-Lambert irradiance march (F6.2a).
 *
 * The unusual thing about this kernel is that the discrete answer is not an
 * approximation. Attenuation is piecewise constant within a layer, and the
 * exponential composes exactly, so the interface irradiances equal
 * `I0 exp(-gamma z)` to round-off no matter how coarse the layering is. The
 * layer means have a closed form too. So every gate below is written against an
 * exact value at 1e-12 rather than a discretization tolerance — which means any
 * real defect has nowhere to hide, and there is no loose threshold absorbing it.
 *
 * The one quantity that *is* an approximation is the relationship between the
 * layer mean and the point value at the cell center; that converges at second
 * order, and the refinement gate measures the rate rather than asserting it.
 */

#include "fvqual/kinetics/lightfield.h"
#include "fvqual/mesh/horizontalmesh.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/mesh/sigmaspec.h"
#include "fvqual/state/speciesregistry.h"
#include "fvqual/state/statestore.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

using FVQual::Kinetics::LightField;
using FVQual::Kinetics::LightOptions;
using FVQual::Kinetics::meanAttenuationFactor;
using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;
using FVQual::State::SpeciesDefinition;
using FVQual::State::SpeciesRegistry;
using FVQual::State::StateStore;
namespace Testing = FVQual::Testing;

namespace {

/*!
 * \brief A flat column of uniform layers, with an optional shading species.
 */
struct Column
{
    HorizontalMesh  horizontal;
    LayeredMesh     layered;
    SpeciesRegistry registry;
    StateStore      state;
    LightField      light;

    double depth     = 0.0;
    double thickness = 0.0;
    int    shader    = -1;

    void build(int layers, double totalDepth, bool withShader = false)
    {
        std::string message;

        ASSERT_TRUE(horizontal.build(Testing::channelStrip(2, 100.0, 40.0), message))
            << message;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = -totalDepth;

        SigmaSpec spec;
        spec.layerCount        = layers;
        spec.minLayerThickness = 0.0;

        ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

        if (withShader)
        {
            SpeciesDefinition definition;
            definition.id           = "ALGAE";
            definition.units        = "g/m3";
            definition.initialValue = 0.0;

            ASSERT_GE(registry.add(definition, message), 0) << message;
            shader = registry.indexOf("ALGAE");
        }

        state.allocate(layered, registry);
        light.initialize(layered, registry);

        depth     = totalDepth;
        thickness = totalDepth / static_cast<double>(layers);
    }

    //! Irradiance at the top of layer \p k in column 0.
    [[nodiscard]] double atInterface(int k) const
    {
        return light.interfaceIrradiance()[static_cast<size_t>(
            layered.interfaceSlot(0, k))];
    }

    //! Layer-mean irradiance in layer \p k of column 0.
    [[nodiscard]] double inLayer(int k) const
    {
        return light.meanIrradiance()[static_cast<size_t>(layered.cell(0, k))];
    }
};

} // namespace

// ── The layer-mean factor, alone ──────────────────────────────────────────

TEST(LightField, TheMeanFactorIsAccurateAcrossTheWholeRange)
{
    // (1 - e^-tau)/tau, which is where the naive expression fails: at tau = 0 it
    // divides by zero, and just above it subtracts two nearly equal numbers.
    // An unlit or vanishingly thin layer is ordinary, not exotic — the top layer
    // of a sigma grid gets thin whenever the surface drops.
    EXPECT_NEAR(meanAttenuationFactor(0.0), 1.0, 1.0e-15);

    // Against the series, where the closed form cannot be trusted.
    for (double tau : {1.0e-12, 1.0e-9, 1.0e-6})
    {
        const double series = 1.0 - tau / 2.0 + tau * tau / 6.0;

        EXPECT_NEAR(meanAttenuationFactor(tau), series, 1.0e-15) << "tau = " << tau;
    }

    // Against the closed form, where it can.
    for (double tau : {0.01, 0.5, 1.0, 5.0, 50.0})
    {
        const double exact = (1.0 - std::exp(-tau)) / tau;

        EXPECT_NEAR(meanAttenuationFactor(tau), exact, 1.0e-14) << "tau = " << tau;
    }

    // And it is a mean of values in (0, 1], so it must be too.
    for (double tau : {0.0, 0.1, 1.0, 10.0, 1000.0})
    {
        EXPECT_GT(meanAttenuationFactor(tau), 0.0) << "tau = " << tau;
        EXPECT_LE(meanAttenuationFactor(tau), 1.0) << "tau = " << tau;
    }
}

// ── The march ─────────────────────────────────────────────────────────────

TEST(LightField, InterfaceIrradianceIsExactlyBeerLambert)
{
    constexpr double incident   = 400.0;
    constexpr double beta       = 0.45;
    constexpr double kd         = 0.6;
    constexpr double totalDepth = 12.0;
    constexpr int    layers     = 8;

    Column column;
    column.build(layers, totalDepth);

    LightOptions options;
    options.incidentShortwave    = incident;
    options.surfaceAbsorption    = beta;
    options.backgroundAttenuation = kd;

    column.light.update(column.layered, column.state, options);

    const double penetrating = incident * (1.0 - beta);

    // Exact, not approximate: exp(-kd h) composed k times is exp(-kd k h).
    for (int k = 0; k <= layers; ++k)
    {
        const double z     = static_cast<double>(k) * column.thickness;
        const double exact = penetrating * std::exp(-kd * z);

        EXPECT_NEAR(column.atInterface(k), exact, 1.0e-12 * penetrating)
            << "interface " << k << " at z = " << z;
    }

    // The surface value is what got past the surface, not the incident value —
    // conflating the two overstates available light by nearly a factor of two.
    EXPECT_NEAR(column.atInterface(0), penetrating, 1.0e-12 * incident);
    EXPECT_NEAR(column.light.surfaceAbsorbed()[0], incident * beta, 1.0e-12 * incident);
}

TEST(LightField, CellValuesAreLayerMeansAndNotTopOfLayerValues)
{
    // The distinction this test exists for. In a 2 m layer at kd = 0.6 the mean
    // is 30% below the top value, and using the top value would bias production
    // high everywhere — worst in the surface layer, where most of it happens.
    constexpr double incident   = 300.0;
    constexpr double kd         = 0.6;
    constexpr double totalDepth = 12.0;
    constexpr int    layers     = 6;

    Column column;
    column.build(layers, totalDepth);

    LightOptions options;
    options.incidentShortwave     = incident;
    options.surfaceAbsorption     = 0.0; // simpler arithmetic to check by hand
    options.backgroundAttenuation = kd;

    column.light.update(column.layered, column.state, options);

    const double h   = column.thickness;
    const double tau = kd * h;

    for (int k = 0; k < layers; ++k)
    {
        const double top = incident * std::exp(-kd * static_cast<double>(k) * h);

        // The analytic mean over the layer, computed independently here.
        const double exact = top * (1.0 - std::exp(-tau)) / tau;

        EXPECT_NEAR(column.inLayer(k), exact, 1.0e-12 * incident) << "layer " << k;

        // And it must be strictly between the interface values that bound it.
        EXPECT_LT(column.inLayer(k), column.atInterface(k)) << "layer " << k;
        EXPECT_GT(column.inLayer(k), column.atInterface(k + 1)) << "layer " << k;
    }

    // The gate with teeth: a top-of-layer implementation would pass every
    // "irradiance decreases with depth" check ever written. It differs from the
    // mean by this much, which is far outside the tolerance above.
    const double topOfLayer = incident;

    EXPECT_LT(column.inLayer(0), 0.9 * topOfLayer)
        << "the surface cell is holding a top-of-layer value";
}

TEST(LightField, AbsorbedEnergyAccountsForEverythingIncident)
{
    // An energy budget. Whatever is not absorbed in the water column or at the
    // surface reaches the bed, and the three must sum to the incident flux. This
    // is the identity a solar heat source will depend on, so it is pinned before
    // anything depends on it.
    constexpr double incident   = 500.0;
    constexpr double beta       = 0.35;
    constexpr double kd         = 0.4;
    constexpr double totalDepth = 15.0;
    constexpr int    layers     = 10;

    Column column;
    column.build(layers, totalDepth);

    LightOptions options;
    options.incidentShortwave     = incident;
    options.surfaceAbsorption     = beta;
    options.backgroundAttenuation = kd;

    column.light.update(column.layered, column.state, options);

    double absorbedInWater = 0.0;

    for (int k = 0; k < layers; ++k)
        absorbedInWater += column.atInterface(k) - column.atInterface(k + 1);

    const double atSurface = column.light.surfaceAbsorbed()[0];
    const double toTheBed  = column.atInterface(layers);

    EXPECT_NEAR(absorbedInWater + atSurface + toTheBed, incident, 1.0e-12 * incident);
}

// ── Self-shading ──────────────────────────────────────────────────────────

TEST(LightField, ConstituentsShadeWhatIsBelowThem)
{
    // The coupling that makes light a 3-D problem rather than a boundary
    // condition: algae attenuate the light they are competing for.
    constexpr double incident   = 400.0;
    constexpr double kd         = 0.3;
    constexpr double alpha      = 0.02; // m²/g
    constexpr double algae      = 25.0; // g/m³
    constexpr double totalDepth = 10.0;
    constexpr int    layers     = 5;

    Column column;
    column.build(layers, totalDepth, /*withShader=*/true);

    std::string message;
    ASSERT_TRUE(column.light.setSpecificAttenuation(column.shader, alpha, message))
        << message;

    // A bloom confined to the top two layers, so the layers beneath are shaded
    // by constituents they do not contain — which a per-cell formula that
    // forgot to carry the march downward would get wrong.
    for (int k = 0; k < 2; ++k)
        column.state.speciesAt(column.shader, column.layered.cell(0, k)) = algae;

    LightOptions options;
    options.incidentShortwave     = incident;
    options.surfaceAbsorption     = 0.0;
    options.backgroundAttenuation = kd;

    column.light.update(column.layered, column.state, options);

    const double h      = column.thickness;
    const double shaded = kd + alpha * algae;

    // Computed here from the two optical depths, independently of the march.
    const double belowBloom = incident * std::exp(-shaded * 2.0 * h);

    EXPECT_NEAR(column.atInterface(2), belowBloom, 1.0e-12 * incident);

    for (int k = 2; k <= layers; ++k)
    {
        const double clearBelow = static_cast<double>(k - 2) * h;
        const double exact      = belowBloom * std::exp(-kd * clearBelow);

        EXPECT_NEAR(column.atInterface(k), exact, 1.0e-12 * incident) << "interface " << k;
    }

    // And the shading matters: without the bloom, the bed sees this instead.
    const double unshaded = incident * std::exp(-kd * totalDepth);

    EXPECT_LT(column.atInterface(layers), 0.5 * unshaded)
        << "the bloom had almost no effect; is the attenuation table being read?";
}

TEST(LightField, ASpeciesOnlyShadesIfItWasGivenACoefficient)
{
    // Otherwise every registered constituent would silently darken the water,
    // and a model would get dimmer as it got more detailed.
    constexpr double incident   = 200.0;
    constexpr double kd         = 0.25;
    constexpr double totalDepth = 8.0;

    Column column;
    column.build(4, totalDepth, /*withShader=*/true);

    for (int64_t i = 0; i < column.state.cellCount(); ++i)
        column.state.speciesAt(column.shader, i) = 100.0;

    LightOptions options;
    options.incidentShortwave     = incident;
    options.surfaceAbsorption     = 0.0;
    options.backgroundAttenuation = kd;

    column.light.update(column.layered, column.state, options);

    EXPECT_NEAR(column.atInterface(4), incident * std::exp(-kd * totalDepth),
                1.0e-12 * incident)
        << "a species with no coefficient attenuated the light";

    EXPECT_EQ(column.light.specificAttenuation(column.shader), 0.0);
}

TEST(LightField, ANegativeCoefficientIsRefused)
{
    // It would be a constituent that creates light on the way down, and it would
    // let the total attenuation go negative and the march grow without bound.
    Column column;
    column.build(4, 8.0, /*withShader=*/true);

    std::string message;

    EXPECT_FALSE(column.light.setSpecificAttenuation(column.shader, -0.01, message));
    EXPECT_FALSE(message.empty());

    EXPECT_FALSE(column.light.setSpecificAttenuation(99, 0.01, message));
    EXPECT_FALSE(column.light.setSpecificAttenuation(-1, 0.01, message));
}

// ── Refinement ────────────────────────────────────────────────────────────

TEST(LightField, LayerMeansConvergeToThePointValueAtSecondOrder)
{
    // The one genuinely approximate relationship here: a cell holds the mean
    // over its layer, and that approaches the value at the cell center as the
    // layers thin. Second order, and measured rather than asserted — if it came
    // out first order, the mean would be being computed by something cruder
    // than the closed form.
    constexpr double incident   = 400.0;
    constexpr double kd         = 0.5;
    constexpr double totalDepth = 10.0;
    constexpr double probe      = 5.0; // depth at which to compare

    double previousError = 0.0;
    double previousRatio = 0.0;

    for (int layers : {10, 20, 40, 80})
    {
        Column column;
        column.build(layers, totalDepth);

        LightOptions options;
        options.incidentShortwave     = incident;
        options.surfaceAbsorption     = 0.0;
        options.backgroundAttenuation = kd;

        column.light.update(column.layered, column.state, options);

        // The layer whose center sits at the probe depth.
        const double h     = totalDepth / static_cast<double>(layers);
        const int    index = static_cast<int>(probe / h) - 1;
        const double centre = (static_cast<double>(index) + 0.5) * h;

        const double point = incident * std::exp(-kd * centre);
        const double error = std::abs(column.inLayer(index) - point);

        if (previousError > 0.0)
            previousRatio = previousError / error;

        previousError = error;
    }

    // Halving the layer thickness must quarter the error.
    EXPECT_GT(previousRatio, 3.5) << "convergence is slower than second order: ratio "
                                  << previousRatio;
    EXPECT_LT(previousRatio, 4.5) << "ratio " << previousRatio;
}

TEST(LightField, ADarkSurfaceLeavesTheColumnDark)
{
    // Night. Worth pinning because a multiplicative march is exactly the kind of
    // code that grows a stray additive term, and because zero incident light
    // drives tau-independent paths through the mean factor.
    Column column;
    column.build(6, 12.0);

    LightOptions options;
    options.incidentShortwave     = 0.0;
    options.surfaceAbsorption     = 0.45;
    options.backgroundAttenuation = 0.5;

    column.light.update(column.layered, column.state, options);

    for (int k = 0; k < 6; ++k)
        EXPECT_EQ(column.inLayer(k), 0.0) << "layer " << k;

    EXPECT_EQ(column.light.surfaceAbsorbed()[0], 0.0);
}
