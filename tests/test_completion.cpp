// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
#include "fvqual/component/fvqualcomponent.h"
#include "fvqual/coupling/intervaltransfer.h"
#include "hydrocouplesdk/io/meshpackage.h"
#include "testmeshes.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <limits>
#include <fstream>
#include <stdexcept>

namespace {
struct Basin {
    FVQual::Mesh::HorizontalMesh horizontal;
    FVQual::Mesh::LayeredMesh layered;
    FVQual::State::SpeciesRegistry registry;
    FVQual::State::StateStore state;
    FVQual::Hydro::InternalModeSolver solver;
    int temperature = -1, tracer = -1;
    void build(HydroCouple::SDK::IO::MeshDefinition definition =
                   FVQual::Testing::twoTriangleSquare()) {
        definition.nodeZ.assign(definition.nodeX.size(), -5.0);
        std::string message;
        ASSERT_TRUE(horizontal.build(definition, message)) << message;
        FVQual::Mesh::SigmaSpec spec;
        spec.layerCount = 2;
        spec.minLayerThickness = 0.0;
        ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;
        temperature = registry.addTemperature(10.0, message);
        ASSERT_GE(temperature, 0) << message;
        FVQual::State::SpeciesDefinition species;
        species.id = "tracer";
        species.units = "kg/m3";
        tracer = registry.add(species, message);
        ASSERT_GE(tracer, 0) << message;
        state.allocate(layered, registry);
        FVQual::Hydro::InternalModeOptions options;
        options.maxTimeStep = 0.05;
        ASSERT_TRUE(solver.initialize(layered, state, registry, options, message)) << message;
    }
};
}

TEST(CompletionMesh, ImportedRegionalMeshRetainsCellsAndConservesAClosedReservoir) {
    HydroCouple::SDK::IO::MeshPackage package;
    HydroCouple::SDK::IO::MeshImportOptions options;
    options.horizontalCrs = "LOCAL";
    options.verticalDatum = "synthetic-terrain";
    options.elevationKind = "terrain-elevation";
    std::string message;
    ASSERT_TRUE(HydroCouple::SDK::IO::importSurfaceMesh(
        std::string(FVQUAL_FIXTURES_DIR) + "/mesh/swmm_regional_triquad.2dm",
        options, package, message)) << message;
    ASSERT_EQ(package.mesh.faceCount(), 1058);
    ASSERT_EQ(package.mesh.nodeCount(), 565);
    int quads = 0;
    for (size_t i = 1; i < package.mesh.faceNodeOffsets.size(); ++i)
        quads += package.mesh.faceNodeOffsets[i] - package.mesh.faceNodeOffsets[i-1] == 4;
    EXPECT_EQ(quads, 24);
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(basin.build(package.mesh));
    std::fill(basin.state.speciesBlock(basin.tracer).begin(),
              basin.state.speciesBlock(basin.tracer).end(), 0.25);
    const double volume = basin.layered.totalVolume();
    const double mass = basin.solver.totalSpeciesMass(basin.tracer);
    ASSERT_TRUE(basin.solver.advance(0.1, message)) << message;
    EXPECT_NEAR(basin.layered.totalVolume(), volume, volume * 1e-12);
    EXPECT_NEAR(basin.solver.totalSpeciesMass(basin.tracer), mass, mass * 1e-12);
    for (double c : basin.state.speciesBlock(basin.tracer)) EXPECT_NEAR(c, 0.25, 1e-12);
}

TEST(CompletionMesh, RejectsNonfiniteAndConcaveQuadsAndPreservesProjectedGeometry) {
    std::string message;
    FVQual::Mesh::HorizontalMesh mesh;
    auto definition = FVQual::Testing::twoTriangleSquare();
    definition.nodeX[0] = std::numeric_limits<double>::quiet_NaN();
    EXPECT_FALSE(mesh.build(definition, message));
    definition = FVQual::Testing::quadStrip(1);
    definition.nodeX = {0,1,0.2,0};
    definition.nodeY = {0,0,0.2,1};
    definition.faceNodes = {0,1,2,3};
    EXPECT_FALSE(mesh.build(definition, message));
    definition.nodeX = {1e9,1e9+1,1e9+1,1e9};
    definition.nodeY = {1e9,1e9,1e9+1,1e9+1};
    ASSERT_TRUE(mesh.build(definition, message)) << message;
    EXPECT_DOUBLE_EQ(mesh.faceArea()[0], 1.0);
    EXPECT_DOUBLE_EQ(mesh.faceCentroidX()[0], 1e9+0.5);
    EXPECT_DOUBLE_EQ(mesh.faceCentroidY()[0], 1e9+0.5);
}

TEST(CompletionReplay, RestoresHistorySourcesAndStableStateViewsExactly) {
    Basin basin, foreign;
    ASSERT_NO_FATAL_FAILURE(basin.build());
    ASSERT_NO_FATAL_FAILURE(foreign.build());
    std::string message;
    auto checkpoint = basin.solver.checkpoint(message);
    ASSERT_NE(checkpoint, nullptr) << message;
    const auto *view = basin.state.species.data();
    const auto epoch = basin.layered.geometryEpoch();
    ASSERT_TRUE(basin.solver.advance(0.1, message)) << message;
    const auto species = basin.state.species;
    const auto eta = basin.state.eta;
    const auto u = basin.state.u;
    ASSERT_TRUE(basin.solver.restoreCheckpoint(*checkpoint, message)) << message;
    EXPECT_EQ(basin.state.species.data(), view);
    EXPECT_GT(basin.layered.geometryEpoch(), epoch);
    EXPECT_DOUBLE_EQ(basin.solver.time(), 0.0);
    ASSERT_TRUE(basin.solver.advance(0.1, message)) << message;
    EXPECT_EQ(basin.state.species, species);
    EXPECT_EQ(basin.state.eta, eta);
    EXPECT_EQ(basin.state.u, u);
    EXPECT_FALSE(foreign.solver.restoreCheckpoint(*checkpoint, message));
}

TEST(CompletionCoupling, ConservesHeatAndSpeciesAndNeverReplaysADonorDebit) {
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(basin.build());
    FVQual::Coupling::IntervalTransferReceiver receiver(basin.solver, basin.registry, basin.layered);
    FVQual::Coupling::TransferFrame frame;
    frame.id = "SWMM/OF1/1";
    frame.endSeconds = 0.1;
    FVQual::Coupling::InflowTransfer inflow;
    inflow.column = 0;
    inflow.volumeM3 = 0.001;
    inflow.heatJ = 1000.0 * 4186.0 * inflow.volumeM3 * 12.0;
    inflow.speciesMassKg["tracer"] = 0.0002;
    frame.inflows.push_back(inflow);
    const double volume = basin.layered.totalVolume();
    const double thermalInventory = basin.solver.totalSpeciesMass(basin.temperature);
    std::string message;
    FVQual::Coupling::TransferReceipt receipt;
    ASSERT_TRUE(receiver.apply(frame, receipt, message)) << message;
    EXPECT_NEAR(basin.layered.totalVolume() - volume, inflow.volumeM3, 1e-12);
    EXPECT_NEAR(basin.solver.totalSpeciesMass(basin.tracer), 0.0002, 1e-12);
    EXPECT_NEAR(basin.solver.totalSpeciesMass(basin.temperature) - thermalInventory, 0.012, 1e-10);
    EXPECT_FALSE(receipt.replayed);
    EXPECT_DOUBLE_EQ(receipt.acceptedVolumeM3, 0.001);
    EXPECT_EQ(basin.solver.sources().count(), 0u);
    const auto species = basin.state.species;
    ASSERT_TRUE(receiver.apply(frame, receipt, message)) << message;
    EXPECT_TRUE(receipt.replayed);
    EXPECT_EQ(basin.state.species, species);
    EXPECT_DOUBLE_EQ(basin.solver.time(), 0.1);
    frame.inflows[0].volumeM3 *= 2.0;
    EXPECT_FALSE(receiver.apply(frame, receipt, message));
    EXPECT_EQ(basin.state.species, species);
}

TEST(CompletionCoupling, InvalidFramesLeaveStateAndLedgerUntouched) {
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(basin.build());
    FVQual::Coupling::IntervalTransferReceiver receiver(basin.solver, basin.registry, basin.layered);
    FVQual::Coupling::TransferFrame frame;
    frame.id = "bad";
    frame.endSeconds = 0.1;
    FVQual::Coupling::InflowTransfer inflow;
    inflow.volumeM3 = -1.0;
    frame.inflows.push_back(inflow);
    FVQual::Coupling::TransferReceipt receipt;
    std::string message;
    const auto state = basin.state.species;
    EXPECT_FALSE(receiver.apply(frame, receipt, message));
    EXPECT_EQ(state, basin.state.species);
    EXPECT_TRUE(receiver.receipts().empty());
    EXPECT_DOUBLE_EQ(basin.solver.time(), 0.0);
    frame.inflows[0].volumeM3 = 1e-300;
    frame.inflows[0].heatJ = 1e300;
    frame.inflows[0].speciesMassKg["tracer"] = 0.0;
    EXPECT_FALSE(receiver.apply(frame, receipt, message));
    EXPECT_NE(message.find("overflows"), std::string::npos) << message;
    EXPECT_EQ(state, basin.state.species);
    EXPECT_TRUE(receiver.receipts().empty());
}

TEST(CompletionReplay, ComponentCheckpointLifecycleAndTemperatureUnits) {
    FVQual::Component::FVQualComponent component("completion.test");
    component.initialize();
    component.prepare();
    ASSERT_EQ(component.status(), HydroCouple::IModelComponent::ComponentStatus::Updated);
    auto *quantity = dynamic_cast<HydroCouple::IQuantity *>(component.outputs()[0]->valueDefinition());
    ASSERT_NE(quantity, nullptr);
    ASSERT_NE(quantity->unit(), nullptr);
    EXPECT_DOUBLE_EQ(quantity->unit()->offsetToSI(), 273.15);
    std::string token, message;
    ASSERT_TRUE(component.saveState(token, message)) << message;
    component.update();
    const auto species = component.state().species;
    ASSERT_TRUE(component.restoreState(token, message)) << message;
    EXPECT_DOUBLE_EQ(component.currentTime(), 0.0);
    component.update();
    EXPECT_EQ(component.state().species, species);
    EXPECT_TRUE(component.releaseState(token, message));
    EXPECT_FALSE(component.releaseState(token, message));
}

TEST(CompletionTime, PublishesUtcEpochAndRecordedNextExchangeAfterRestore) {
    FVQual::Component::FVQualComponent component("time");
    for (auto *argument : component.arguments())
        if (argument->id() == "epoch_julian_day") {
            std::string message;
            ASSERT_TRUE(argument->initialize(R"({"values":[2461323.5]})",
                HydroCouple::IArgument::ArgumentInputType::JSON, message)) << message;
        }
    component.initialize();
    EXPECT_DOUBLE_EQ(component.currentDateTime()->julianDay(), 2461323.5);
    EXPECT_DOUBLE_EQ(component.nextDateTimeJulianDay(), 2461323.5);
    EXPECT_DOUBLE_EQ(component.simulationPeriod()->duration(), 1);
    component.prepare();
    EXPECT_DOUBLE_EQ(component.nextDateTimeJulianDay(), 2461323.5 + 3600.0/86400.0);
    std::string token, message;
    ASSERT_TRUE(component.saveState(token, message)) << message;
    component.update();
    EXPECT_DOUBLE_EQ(component.currentDateTime()->julianDay(), 2461323.5 + 3600.0/86400.0);
    ASSERT_TRUE(component.restoreState(token, message)) << message;
    EXPECT_DOUBLE_EQ(component.currentDateTime()->julianDay(), 2461323.5);
    EXPECT_DOUBLE_EQ(component.nextDateTimeJulianDay(), 2461323.5 + 3600.0/86400.0);
}

TEST(CompletionCoupling, NumericalFailureRollsBackPartialAdvanceAndClosureHistory) {
    struct FailingClosure : FVQual::Hydro::IPressureClosure {
        int calls = 0;
        std::string name() const override { return "failure-injection"; }
        bool isNonHydrostatic() const noexcept override { return false; }
        std::unique_ptr<FVQual::Hydro::IPressureClosure> clone() const override {
            return std::make_unique<FailingClosure>(*this);
        }
        void baroclinicAcceleration(const FVQual::Mesh::LayeredMesh &mesh,
            std::span<const double>, double, double,
            std::vector<double> &x, std::vector<double> &y) override {
            if (++calls == 2) throw std::runtime_error("injected pressure failure");
            x.assign(static_cast<size_t>(mesh.cellCount()), 0.0);
            y=x;
        }
    };
    Basin basin;
    ASSERT_NO_FATAL_FAILURE(basin.build());
    basin.solver.setPressureClosure(std::make_unique<FailingClosure>());
    const auto state = basin.state.species;
    const auto eta = basin.state.eta;
    const auto *view = basin.state.species.data();
    FVQual::Coupling::IntervalTransferReceiver receiver(basin.solver,basin.registry,basin.layered);
    FVQual::Coupling::TransferFrame frame{"fails-after-one-substep",0,.2,{}};
    frame.inflows.push_back({0,.001,50232,FVQual::Forcing::Placement::Uniform,0,{{"tracer",.0002}}});
    std::string message;
    FVQual::Coupling::TransferReceipt receipt;
    EXPECT_FALSE(receiver.apply(frame,receipt,message));
    EXPECT_NE(message.find("rolled back"),std::string::npos) << message;
    EXPECT_DOUBLE_EQ(basin.solver.time(),0);
    EXPECT_EQ(basin.state.species,state);
    EXPECT_EQ(basin.state.eta,eta);
    EXPECT_EQ(basin.state.species.data(),view);
    EXPECT_TRUE(receiver.receipts().empty());
    EXPECT_EQ(dynamic_cast<const FailingClosure &>(basin.solver.pressureClosure()).calls,0);
}

TEST(CompletionReplay, RecordedComponentRestoresTransferLedgerAndCanReinitialize) {
    FVQual::Component::FVQualComponent component("recorded");
    auto *counts=dynamic_cast<HydroCouple::SDK::Argument1DInt *>(component.arguments()[0]);
    auto *values=dynamic_cast<HydroCouple::SDK::Argument1DDouble *>(component.arguments()[1]);
    (*counts)[0]=3; (*counts)[1]=4;
    (*values)[0]=20; (*values)[1]=10; (*values)[2]=5;
    (*values)[4]=1; (*values)[5]=4;
    FVQual::Component::RecordedInflowArgument *inflows=nullptr;
    for (auto *arg : component.arguments()) if (arg->id() == "recorded_inflows") inflows=dynamic_cast<FVQual::Component::RecordedInflowArgument *>(arg);
    ASSERT_NE(inflows,nullptr);
    const auto path=std::filesystem::path(FVQUAL_FIXTURES_DIR).parent_path().parent_path() /
        "examples/swmm_reservoir/inflow_series.json";
    std::string message, token;
    ASSERT_TRUE(inflows->initialize(path.string(),HydroCouple::IArgument::ArgumentInputType::File,message)) << message;
    component.initialize(); component.prepare();
    ASSERT_EQ(component.status(),HydroCouple::IModelComponent::ComponentStatus::Updated);
    auto *spatial=dynamic_cast<HydroCouple::Spatial::ILayeredMeshComponentDataItem *>(component.outputs()[0]);
    ASSERT_NE(spatial,nullptr);
    ASSERT_NE(spatial->polyhedralSurface(),nullptr);
    EXPECT_EQ(spatial->polyhedralSurface()->meshView()->faceCount(),3);
    ASSERT_TRUE(component.saveState(token,message));
    component.update();
    std::string errors;
    for (const auto &error : component.errors()) errors += error.message + "\n";
    ASSERT_NE(component.status(),HydroCouple::IModelComponent::ComponentStatus::Failed) << errors;
    const auto expected=component.state().species;
    EXPECT_DOUBLE_EQ(component.currentTime(),1);
    ASSERT_TRUE(component.restoreState(token,message)) << message;
    component.update();
    EXPECT_EQ(component.state().species,expected);
    EXPECT_DOUBLE_EQ(component.currentTime(),1);
    ASSERT_TRUE(component.releaseState(token,message));
    component.finish(); component.initialize(); component.prepare();
    EXPECT_EQ(component.status(),HydroCouple::IModelComponent::ComponentStatus::Updated);
    EXPECT_DOUBLE_EQ(component.currentTime(),0);
    EXPECT_EQ(component.registry().count(),2);
}
