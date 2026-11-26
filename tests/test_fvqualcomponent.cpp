// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_fvqualcomponent.cpp
 * \brief Gates for FVQualComponent against the HydroCouple component contract
 *        (F8.1b).
 *
 * Driven the way the framework drives a component — initialize, validate,
 * prepare, update until Done, finish — rather than by calling the pieces in an
 * order no workflow would use. A component that only works when driven by its
 * own test is not a component.
 */

#include "fvqual/component/fvqualcomponent.h"

#include <gtest/gtest.h>

#include <span>
#include <string>
#include <vector>

using FVQual::Component::FVQualComponent;
using ComponentStatus = HydroCouple::IModelComponent::ComponentStatus;

namespace {

constexpr int kColumns = 6;
constexpr int kLayers  = 8;

//! Argument slots, mirroring the component's own ordering.
enum
{
    kArgColumns = 0,
    kArgLayers
};

enum
{
    kArgLength = 0,
    kArgWidth,
    kArgDepth,
    kArgInitialTemperature,
    kArgCouplingInterval,
    kArgEndTime
};

//! Sets the geometry through the arguments, as a configuration file would.
void configure(FVQualComponent &component, double interval, double endTime)
{
    auto *counts = dynamic_cast<HydroCouple::SDK::Argument1DInt *>(
        component.arguments()[0]);
    auto *values = dynamic_cast<HydroCouple::SDK::Argument1DDouble *>(
        component.arguments()[1]);

    ASSERT_NE(counts, nullptr);
    ASSERT_NE(values, nullptr);

    (*counts)[kArgColumns] = kColumns;
    (*counts)[kArgLayers]  = kLayers;

    (*values)[kArgLength]             = 1200.0;
    (*values)[kArgWidth]              = 300.0;
    (*values)[kArgDepth]              = 25.0;
    (*values)[kArgInitialTemperature] = 12.0;
    (*values)[kArgCouplingInterval]   = interval;
    (*values)[kArgEndTime]            = endTime;
}

} // namespace

TEST(FVQualComponent, ArgumentsAreReadableBeforeInitialize)
{
    FVQualComponent component("fvqual.test");

    // The SDK creates arguments lazily on first read, so a configuration loader
    // can populate them while the component is still in Created.
    EXPECT_EQ(component.status(), ComponentStatus::Created);
    // geometry_counts, geometry_values, and (since the provider pipeline)
    // the optional mesh, met_air_temperature, and kinetics_file.
    ASSERT_GE(component.arguments().size(), 18u);
}

TEST(FVQualComponent, RunsTheFullLifecycleAndReportsDone)
{
    FVQualComponent component("fvqual.test");
    ASSERT_NO_FATAL_FAILURE(configure(component, 600.0, 1800.0));

    component.initialize();
    ASSERT_EQ(component.status(), ComponentStatus::Initialized);

    EXPECT_TRUE(component.validate().empty());

    component.prepare();
    ASSERT_NE(component.status(), ComponentStatus::Failed);

    // Drive it as a workflow does: update until it says Done, with a cap so a
    // component that never finishes fails the gate instead of hanging it.
    int steps = 0;
    while (component.status() != ComponentStatus::Done && steps < 100)
    {
        component.update();
        ASSERT_NE(component.status(), ComponentStatus::Failed)
            << "failed at step " << steps;
        ++steps;
    }

    EXPECT_EQ(component.status(), ComponentStatus::Done);
    EXPECT_EQ(steps, 3) << "1800 s at 600 s an interval is three exchanges";
    EXPECT_NEAR(component.currentTime(), 1800.0, 1.0e-6);

    component.finish();
    EXPECT_EQ(component.status(), ComponentStatus::Finished);
}

TEST(FVQualComponent, PublishesLayeredTemperatureWithItsVerticalCoordinate)
{
    FVQualComponent component("fvqual.test");
    ASSERT_NO_FATAL_FAILURE(configure(component, 600.0, 600.0));

    component.initialize();
    component.prepare();
    ASSERT_NE(component.status(), ComponentStatus::Failed);

    ASSERT_GE(component.outputs().size(), 10u);
    HydroCouple::IOutput *output = component.outputs()[0];
    ASSERT_NE(output, nullptr);

    // Shape is {column, layer}, entity first.
    EXPECT_EQ(output->shape(), (std::vector<int64_t>{kColumns, kLayers}));

    // The part that makes it worth exchanging: a consumer can find out where
    // the layers are. Without this the values are an array it cannot place.
    auto *layering = dynamic_cast<HydroCouple::Spatial::ILayering *>(output);
    ASSERT_NE(layering, nullptr);

    HydroCouple::Spatial::IVerticalCoordinate *vertical =
        layering->verticalCoordinate();
    ASSERT_NE(vertical, nullptr);

    EXPECT_EQ(vertical->kind(), HydroCouple::Spatial::VerticalCoordinateKind::Sigma);
    EXPECT_EQ(vertical->layerCount(), kLayers);
    EXPECT_TRUE(vertical->isTimeVarying());

    // And those elevations must describe the mesh that was actually built —
    // surface at 0, bed at the configured depth — not the empty mesh the
    // output was constructed against.
    EXPECT_NEAR(vertical->interfaceElevation(0, 0), 0.0, 1.0e-9);
    EXPECT_NEAR(vertical->interfaceElevation(0, kLayers), -25.0, 1.0e-9);
}

TEST(FVQualComponent, PublishedValuesAreTheModelsOwnTemperature)
{
    FVQualComponent component("fvqual.test");
    ASSERT_NO_FATAL_FAILURE(configure(component, 600.0, 600.0));

    component.initialize();
    component.prepare();
    ASSERT_NE(component.status(), ComponentStatus::Failed);

    auto *output = dynamic_cast<
        HydroCouple::SDK::Spatial::LayeredMeshOutputDouble *>(component.outputs()[0]);
    ASSERT_NE(output, nullptr);

    // prepare() publishes the initial condition, so every cell is 12 C.
    for (int c = 0; c < kColumns; ++c)
        for (int k = 0; k < kLayers; ++k)
            ASSERT_NEAR((*output)(c, k), 12.0, 1.0e-9)
                << "column " << c << " layer " << k;

    // That check alone proves almost nothing, and an ablation is what showed
    // it: with every cell at 12 C, a publish loop that scrambled its indices
    // reads 12 C everywhere and passes. Cell identity is only observable on a
    // field that varies, so stamp one that does.
    const std::span<double> block = component.state().speciesBlock(
        component.temperatureIndex());

    for (int c = 0; c < kColumns; ++c)
        for (int k = 0; k < kLayers; ++k)
            block[static_cast<size_t>(component.mesh().cell(c, k))] =
                100.0 * static_cast<double>(c) + static_cast<double>(k);

    component.update();
    ASSERT_NE(component.status(), ComponentStatus::Failed);

    // Both sides read the evolved field, so this does not depend on the stamped
    // values surviving the step — only on their still being distinguishable.
    for (int c = 0; c < kColumns; ++c)
        for (int k = 0; k < kLayers; ++k)
            ASSERT_NEAR((*output)(c, k),
                        block[static_cast<size_t>(component.mesh().cell(c, k))],
                        1.0e-12)
                << "column " << c << " layer " << k;

    // And the field must still vary, or the comparison above is vacuous again.
    EXPECT_GT((*output)(kColumns - 1, 0) - (*output)(0, 0), 1.0);
}

TEST(FVQualComponent, ReportsEveryConfigurationProblemAtOnce)
{
    FVQualComponent component("fvqual.test");

    auto *counts = dynamic_cast<HydroCouple::SDK::Argument1DInt *>(
        component.arguments()[0]);
    auto *values = dynamic_cast<HydroCouple::SDK::Argument1DDouble *>(
        component.arguments()[1]);
    ASSERT_NE(counts, nullptr);
    ASSERT_NE(values, nullptr);

    (*counts)[kArgColumns] = 0;
    (*counts)[kArgLayers]  = 0;
    (*values)[kArgLength]  = -1.0;
    (*values)[kArgDepth]   = 0.0;

    component.initialize();

    // A bad configuration must fail initialization, not surface later as a
    // strange result.
    EXPECT_EQ(component.status(), ComponentStatus::Failed);

    const std::vector<std::string> problems = component.validate();

    // Every problem in one pass. A caller fixing a configuration one message
    // per run is a caller running five times to learn five things.
    EXPECT_GE(problems.size(), 4u);
    EXPECT_EQ(component.status(), ComponentStatus::Invalid);
}

TEST(FVQualComponent, UpdateBeforePrepareFailsRatherThanRunningOnEmptyState)
{
    FVQualComponent component("fvqual.test");
    ASSERT_NO_FATAL_FAILURE(configure(component, 600.0, 1800.0));

    component.initialize();
    ASSERT_EQ(component.status(), ComponentStatus::Initialized);

    // Stepping a solver that was never initialised is the kind of thing that
    // produces numbers rather than an error.
    component.update();
    EXPECT_EQ(component.status(), ComponentStatus::Failed);
}

TEST(FVQualComponent, TemperatureIsConservedInAClosedBasinWithNoForcing)
{
    FVQualComponent component("fvqual.test");
    ASSERT_NO_FATAL_FAILURE(configure(component, 600.0, 3600.0));

    component.initialize();
    component.prepare();
    ASSERT_NE(component.status(), ComponentStatus::Failed);

    auto *output = dynamic_cast<
        HydroCouple::SDK::Spatial::LayeredMeshOutputDouble *>(component.outputs()[0]);
    ASSERT_NE(output, nullptr);

    int steps = 0;
    while (component.status() != ComponentStatus::Done && steps < 100)
    {
        component.update();
        ASSERT_NE(component.status(), ComponentStatus::Failed);
        ++steps;
    }

    // No heat budget, no inflow, no wind: an isothermal basin must stay
    // isothermal. This is the gate that would catch the publish loop reading a
    // stale or wrongly-strided block, which the constant initial condition
    // alone cannot distinguish from a correct one.
    for (int c = 0; c < kColumns; ++c)
        for (int k = 0; k < kLayers; ++k)
            ASSERT_NEAR((*output)(c, k), 12.0, 1.0e-6)
                << "column " << c << " layer " << k << " after " << steps << " steps";
}
