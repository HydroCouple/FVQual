// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/state/statestore.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;
using FVQual::State::SpeciesDefinition;
using FVQual::State::SpeciesKind;
using FVQual::State::SpeciesRegistry;
using FVQual::State::StateStore;
namespace Testing = FVQual::Testing;

namespace {

struct Fixture
{
    HorizontalMesh horizontal;
    LayeredMesh    layered;

    explicit Fixture(int columns = 4, int layers = 5, double bed = -10.0)
    {
        auto definition = Testing::quadStrip(columns);
        definition.nodeZ.assign(definition.nodeX.size(), bed);

        std::string message;
        EXPECT_TRUE(horizontal.build(definition, message)) << message;

        SigmaSpec spec;
        spec.layerCount        = layers;
        spec.minLayerThickness = 0.0;

        EXPECT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;
    }
};

} // namespace

TEST(SpeciesRegistry, ReservedSpeciesAreTracked)
{
    SpeciesRegistry registry;
    std::string     message;

    const int temperature = registry.addTemperature(12.5, message);
    ASSERT_GE(temperature, 0) << message;

    const int age = registry.addAge(message);
    ASSERT_GE(age, 0) << message;

    EXPECT_EQ(registry.temperatureIndex(), temperature);
    EXPECT_EQ(registry.ageIndex(), age);
    EXPECT_EQ(registry.indexOf(SpeciesRegistry::kTemperature), temperature);
    EXPECT_EQ(registry.count(), 2);

    EXPECT_EQ(registry.definition(temperature).units, "degC");
    EXPECT_DOUBLE_EQ(registry.definition(temperature).initialValue, 12.5);
    EXPECT_TRUE(registry.definition(temperature).isTransported());
}

TEST(SpeciesRegistry, DuplicatesAreRejected)
{
    SpeciesRegistry registry;
    std::string     message;

    ASSERT_GE(registry.addTemperature(10.0, message), 0) << message;
    EXPECT_LT(registry.addTemperature(11.0, message), 0);
    EXPECT_FALSE(message.empty());

    SpeciesDefinition oxygen;
    oxygen.id    = "DO";
    oxygen.units = "g/m3";

    ASSERT_GE(registry.add(oxygen, message), 0) << message;
    EXPECT_LT(registry.add(oxygen, message), 0);
}

TEST(SpeciesRegistry, BenthicSpeciesAreNotTransported)
{
    SpeciesDefinition sediment;
    sediment.id   = "SOD_POOL";
    sediment.kind = SpeciesKind::MsxBenthic;

    EXPECT_FALSE(sediment.isTransported());
}

TEST(StateStore, BlocksAreSizedForTheMesh)
{
    Fixture         fixture(4, 5, -10.0);
    SpeciesRegistry registry;
    std::string     message;

    ASSERT_GE(registry.addTemperature(15.0, message), 0) << message;
    ASSERT_GE(registry.addAge(message), 0) << message;

    StateStore state;
    state.allocate(fixture.layered, registry);

    EXPECT_EQ(state.columnCount(), 4);
    EXPECT_EQ(state.cellCount(), 20);
    EXPECT_EQ(state.speciesCount(), 2);

    EXPECT_EQ(state.eta.size(), 4u);
    EXPECT_EQ(state.u.size(), 20u);
    EXPECT_EQ(state.w.size(), static_cast<size_t>(fixture.layered.interfaceSlotCount()));
    EXPECT_EQ(state.species.size(), 40u);
}

TEST(StateStore, InitialValuesAreAppliedPerSpecies)
{
    Fixture         fixture;
    SpeciesRegistry registry;
    std::string     message;

    const int temperature = registry.addTemperature(18.0, message);
    const int age         = registry.addAge(message);
    ASSERT_GE(temperature, 0);
    ASSERT_GE(age, 0);

    StateStore state;
    state.allocate(fixture.layered, registry);

    for (int64_t cell = 0; cell < state.cellCount(); ++cell)
    {
        EXPECT_DOUBLE_EQ(state.speciesAt(temperature, cell), 18.0);
        EXPECT_DOUBLE_EQ(state.speciesAt(age, cell), 0.0);
    }
}

TEST(StateStore, SpeciesBlocksAreContiguousAndDisjoint)
{
    Fixture         fixture;
    SpeciesRegistry registry;
    std::string     message;

    const int temperature = registry.addTemperature(4.0, message);
    const int age         = registry.addAge(message);

    StateStore state;
    state.allocate(fixture.layered, registry);

    std::span<double> block = state.speciesBlock(age);
    ASSERT_EQ(block.size(), static_cast<size_t>(state.cellCount()));

    // Writing through one species' block must not disturb the other's.
    for (double &value : block)
        value = 7.0;

    for (int64_t cell = 0; cell < state.cellCount(); ++cell)
    {
        EXPECT_DOUBLE_EQ(state.speciesAt(age, cell), 7.0);
        EXPECT_DOUBLE_EQ(state.speciesAt(temperature, cell), 4.0);
    }

    EXPECT_EQ(&block[0], &state.speciesAt(age, 0));
}

TEST(StateStore, TotalMassIsVolumeWeighted)
{
    // Four unit-area columns, 10 m deep: 40 m³ of water.
    Fixture         fixture(4, 5, -10.0);
    SpeciesRegistry registry;
    std::string     message;

    SpeciesDefinition tracer;
    tracer.id           = "TRACER";
    tracer.units        = "g/m3";
    tracer.initialValue = 2.0;

    const int index = registry.add(tracer, message);
    ASSERT_GE(index, 0) << message;

    StateStore state;
    state.allocate(fixture.layered, registry);

    EXPECT_NEAR(fixture.layered.totalVolume(), 40.0, 1e-9);
    EXPECT_NEAR(state.totalMass(index, fixture.layered), 80.0, 1e-9);
    EXPECT_NEAR(state.volumeWeightedMean(index, fixture.layered), 2.0, 1e-12);
}

TEST(StateStore, SnapshotSlicesCarryTheLayeredShape)
{
    Fixture         fixture(3, 6, -8.0);
    SpeciesRegistry registry;
    std::string     message;

    const int temperature = registry.addTemperature(9.0, message);
    ASSERT_GE(temperature, 0) << message;

    StateStore state;
    state.allocate(fixture.layered, registry);

    const auto layered = state.layeredSlice("temperature", state.speciesBlock(temperature));

    ASSERT_EQ(layered.shape.size(), 2u);
    EXPECT_EQ(layered.shape[0], 3);
    EXPECT_EQ(layered.shape[1], 6);
    EXPECT_EQ(layered.elementCount(), 18);
    EXPECT_EQ(layered.kind, HydroCouple::DataKind::Float64);
    EXPECT_DOUBLE_EQ(layered.values<double>()[0], 9.0);

    const auto column = state.columnSlice("eta", state.eta);

    ASSERT_EQ(column.shape.size(), 1u);
    EXPECT_EQ(column.shape[0], 3);
    EXPECT_EQ(column.elementCount(), 3);
}
