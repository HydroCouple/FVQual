// SPDX-License-Identifier: MIT
#include "fvqual/hydro/internalmode.h"
#include "testmeshes.h"
#include <gtest/gtest.h>
#include <cmath>
#include <fstream>
#include <limits>

namespace {
struct Tank {
    FVQual::Mesh::HorizontalMesh horizontal;
    FVQual::Mesh::LayeredMesh mesh;
    FVQual::State::SpeciesRegistry registry;
    FVQual::State::StateStore state;
    FVQual::Hydro::InternalModeSolver solver;
    void build(double dt = 300, int layers = 1) {
        auto definition = FVQual::Testing::channelStrip(1, 10, 10);
        definition.nodeZ.assign(definition.nodeX.size(), -5);
        std::string message;
        ASSERT_TRUE(horizontal.build(definition, message)) << message;
        FVQual::Mesh::SigmaSpec sigma; sigma.layerCount = layers;
        ASSERT_TRUE(mesh.build(horizontal, sigma, 0, message)) << message;
        FVQual::State::SpeciesDefinition tracer; tracer.id = "tracer"; tracer.units = "kg/m3"; tracer.initialValue = 10;
        ASSERT_GE(registry.add(tracer, message), 0) << message;
        state.allocate(mesh, registry);
        FVQual::Hydro::InternalModeOptions options; options.maxTimeStep = dt;
        options.backgroundDiffusivity = 0; options.neutralViscosity = 0;
        ASSERT_TRUE(solver.initialize(mesh, state, registry, options, message)) << message;
    }
    void exchange(double q) {
        std::string message;
        FVQual::Forcing::Source inlet; inlet.discharge = q; inlet.concentration = {0};
        ASSERT_GE(solver.sources().add(inlet, message), 0);
        inlet.discharge = -q;
        ASSERT_GE(solver.sources().add(inlet, message), 0);
    }
};
}
TEST(GrossTurnover, NetZeroFirstStepStaysPositiveAndRefinesToTankOracle) {
    std::ofstream evidence(std::string(FVQual_TEST_OUTPUT_DIR) + "/gross_turnover.csv");
    evidence << "dt_s,tracer_kg_m3,exact,error,water_residual,mass_residual\n";
    double previous = 1e30;
    for (double dt : {.25, .125, .0625, .03125}) {
        Tank tank; ASSERT_NO_FATAL_FAILURE(tank.build(dt));
        ASSERT_NO_FATAL_FAILURE(tank.exchange(1000));
        std::string message;
        EXPECT_NEAR(tank.solver.stableTimeStep(), dt, 1e-12);
        ASSERT_TRUE(tank.solver.step(1, message)) << message; // direct callers are substepped too
        const double value = tank.state.speciesAt(0,0), exact = 10 * std::exp(-1000. / 500.);
        EXPECT_GE(value, 0); EXPECT_LE(value, 10);
        const double error = std::abs(value - exact);
        EXPECT_LT(error, previous * .6); previous = error;
        EXPECT_NEAR(tank.mesh.totalVolume(), 500, 1e-10);
        EXPECT_NEAR(tank.solver.waterBudgetResidual(), 0, 1e-10);
        EXPECT_NEAR(tank.solver.scalarBudgetResidual(0), 0, 1e-9);
        EXPECT_NEAR(tank.solver.budget().flows[0].waterOut, 1000, 1e-10);
        evidence << dt << ',' << value << ',' << exact << ',' << error << ','
                 << tank.solver.waterBudgetResidual() << ',' << tank.solver.scalarBudgetResidual(0) << '\n';
    }
}
TEST(GrossTurnover, SummedSinksAndLayerPlacementBindBelowHydrodynamicFloor) {
    Tank tank; ASSERT_NO_FATAL_FAILURE(tank.build(300, 3));
    std::string message;
    FVQual::Forcing::Source source; source.discharge = 1e7; source.placement = FVQual::Forcing::Placement::Surface;
    ASSERT_GE(tank.solver.sources().add(source, message), 0);
    source.discharge = -5e6;
    ASSERT_GE(tank.solver.sources().add(source, message), 0);
    ASSERT_GE(tank.solver.sources().add(source, message), 0);
    EXPECT_NEAR(tank.solver.stableTimeStep(), .5 * (500./3.) / 1e7, 1e-15);
    ASSERT_TRUE(tank.solver.advance(.0001, message)) << message;
    for (double c : tank.state.speciesBlock(0)) { EXPECT_GE(c, 0); EXPECT_LE(c, 10); }
    EXPECT_NEAR(tank.solver.scalarBudgetResidual(0), 0, 1e-9);
}
TEST(GrossTurnover, NewlyActivatedPumpAndHeldScheduleRespectGrossExchange) {
    Tank tank; ASSERT_NO_FATAL_FAILURE(tank.build());
    std::string message;
    FVQual::Forcing::Source inlet; inlet.discharge = 0; inlet.concentration = {0};
    const int index = tank.solver.sources().add(inlet, message); ASSERT_GE(index, 0);
    FVQual::Forcing::Structure pump; pump.type = FVQual::Forcing::StructureType::Pump;
    pump.column = 0; pump.capacity = -1000; pump.placement = FVQual::Forcing::Placement::Uniform; pump.scheduled = true;
    const int port = tank.solver.structures().add(pump, tank.solver.sources(), message); ASSERT_GE(port, 0) << message;
    FVQual::Forcing::OperationControl in, release;
    in.kind = FVQual::Forcing::OperationControl::Kind::SourceDischarge; in.index = index;
    in.times = {0,.1}; in.values = {0,1000};
    release.kind = FVQual::Forcing::OperationControl::Kind::StructureDemand; release.index = port;
    release.times = {0,.1}; release.values = {0,1000};
    tank.solver.setOperationControls({in, release});
    ASSERT_TRUE(tank.solver.advance(1.1, message)) << message;
    EXPECT_GE(tank.state.speciesAt(0,0), 0); EXPECT_LT(tank.state.speciesAt(0,0), 10);
    EXPECT_NEAR(tank.mesh.totalVolume(), 500, 1e-9);
    EXPECT_NEAR(tank.solver.budget().flows[0].waterIn, 1000, 1e-9);
    EXPECT_NEAR(tank.solver.budget().flows[1].waterOut, 1000, 1e-9);
    EXPECT_NEAR(tank.solver.scalarBudgetResidual(0), 0, 1e-8);
}

TEST(GrossTurnover, NonfiniteDurationsRefuseBeforeMutation) {
    Tank tank;
    ASSERT_NO_FATAL_FAILURE(tank.build());
    std::string message;
    const auto values = tank.state.species;
    EXPECT_FALSE(tank.solver.step(std::numeric_limits<double>::infinity(), message));
    EXPECT_FALSE(tank.solver.advance(std::numeric_limits<double>::infinity(), message));
    EXPECT_FALSE(tank.solver.advance(std::numeric_limits<double>::quiet_NaN(), message));
    EXPECT_EQ(tank.state.species, values);
    EXPECT_DOUBLE_EQ(tank.solver.time(), 0);
}

TEST(GrossTurnover, ActualBoundaryWithdrawalRetriesWithoutCountingRejectedWater) {
    Tank tank;
    ASSERT_NO_FATAL_FAILURE(tank.build());
    std::string message;
    FVQual::Forcing::Source inlet; inlet.discharge = 1000; inlet.concentration = {0};
    ASSERT_GE(tank.solver.sources().add(inlet, message), 0);
    std::vector<int64_t> outletEdges;
    for (int64_t e=0; e<tank.horizontal.edgeCount(); ++e)
        if (tank.horizontal.edgeRightFace()[e] == FVQual::Mesh::kNoFace && tank.horizontal.edgeMidX()[e] == 10)
            outletEdges.push_back(e);
    ASSERT_EQ(outletEdges.size(), 1);
    FVQual::Hydro::BarotropicSolver::BoundaryCondition outlet;
    outlet.type = FVQual::Hydro::BarotropicSolver::BoundaryType::Stage;
    outlet.stage = 0; // rising inlet water drives a dynamically rated outlet
    ASSERT_GE(tank.solver.external().addBoundary(outlet, outletEdges, message), 0) << message;
    ASSERT_TRUE(tank.solver.bindOpenBoundaries(FVQual::Forcing::Placement::Uniform, message)) << message;
    auto saved = tank.solver.checkpoint(message); ASSERT_NE(saved, nullptr) << message;
    ASSERT_TRUE(tank.solver.step(10, message)) << message;
    const auto values = tank.state.species;
    EXPECT_GE(values[0], 0); EXPECT_LE(values[0], 10);
    EXPECT_NEAR(tank.solver.budget().flows[0].waterIn, 10000, 1e-8);
    EXPECT_GT(tank.solver.budget().flows[2].waterOut, 500);
    EXPECT_NEAR(tank.solver.waterBudgetResidual(), 0, 1e-8);
    EXPECT_NEAR(tank.solver.scalarBudgetResidual(0), 0, 1e-8);
    ASSERT_TRUE(tank.solver.restoreCheckpoint(*saved, message)) << message;
    ASSERT_TRUE(tank.solver.step(10, message)) << message;
    EXPECT_EQ(tank.state.species, values);
}
