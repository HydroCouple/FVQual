// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
#include "fvqual/hydro/internalmode.h"
#include "testmeshes.h"
#include <gtest/gtest.h>
#include <cmath>

namespace {
struct Column {
    FVQual::Mesh::HorizontalMesh horizontal;
    FVQual::Mesh::LayeredMesh mesh;
    FVQual::State::SpeciesRegistry registry;
    FVQual::State::StateStore state;
    FVQual::Hydro::InternalModeSolver solver;
    void build(FVQual::Hydro::InternalModeOptions options, bool tapered = true) {
        std::string message;
        auto definition = FVQual::Testing::channelStrip(1, 10, 10);
        definition.nodeZ.assign(4, -4);
        ASSERT_TRUE(horizontal.build(definition, message)) << message;
        FVQual::Mesh::SigmaSpec spec;
        spec.layerCount = 2; spec.minLayerThickness = 0;
        ASSERT_TRUE(mesh.build(horizontal, spec, 0, message)) << message;
        if (tapered) {
            ASSERT_TRUE(mesh.setHypsography(0, std::vector<double>{-4,0},
                                               std::vector<double>{20,100}, message)) << message;
            ASSERT_TRUE(mesh.updateSurface(std::vector<double>{0}, message)) << message;
        }
        ASSERT_GE(registry.addTemperature(10, message), 0) << message;
        state.allocate(mesh, registry);
        options.neutralViscosity = 0;
        options.backgroundDiffusivity = 0;
        ASSERT_TRUE(solver.initialize(mesh,state,registry,options,message)) << message;
    }
};
}

TEST(MomentumGeometry, TaperedViscosityMatchesIndependentTwoBoxImpulse) {
    // A(z)=20+20*(z+4): Vtop=160, Vbottom=80, Acontact=60, spacing=2.
    // G=nu*Acontact/spacing=3. Backward Euler damps the velocity difference
    // by 1/[1+dt*G*(1/160+1/80)]. A uniform mode correction cannot change it.
    Column column;
    FVQual::Hydro::InternalModeOptions options;
    options.manningN=0; options.backgroundViscosity=.1; options.auditMomentum=true;
    ASSERT_NO_FATAL_FAILURE(column.build(options));
    column.state.u={1,-1};
    std::string message;
    ASSERT_TRUE(column.solver.step(10,message)) << message;
    EXPECT_NEAR(column.state.u[0]-column.state.u[1],1.28,1e-12);
    const auto &audit=column.solver.momentumAudit().at(0);
    EXPECT_NEAR(audit.verticalResidualX,0,1e-10);
    EXPECT_NEAR(audit.viscosityDissipation,49152,1e-8);
    EXPECT_NEAR(audit.implicitDissipation,13824,1e-8);
    EXPECT_NEAR(audit.verticalEnergyResidual,0,1e-8);
    // The rectangular conveyance mean, reconciled to zero, differs from the
    // volume mean in a taper. Measure this impulse rather than hiding it.
    EXPECT_NEAR(audit.modeImpulseX,-28800,1e-8);
    auto checkpoint=column.solver.checkpoint(message);
    ASSERT_TRUE(checkpoint) << message;
    ASSERT_TRUE(column.solver.step(10,message)) << message;
    const auto energy=column.solver.momentumAudit()[0].viscosityDissipation;
    ASSERT_TRUE(column.solver.restoreCheckpoint(*checkpoint,message)) << message;
    EXPECT_NEAR(column.solver.momentumAudit()[0].viscosityDissipation,49152,1e-8);
    ASSERT_TRUE(column.solver.step(10,message)) << message;
    EXPECT_DOUBLE_EQ(column.solver.momentumAudit()[0].viscosityDissipation,energy);
}

TEST(MomentumGeometry, WindShearUsesSurfaceContactAreaAndActualTopStorage) {
    Column column;
    FVQual::Hydro::InternalModeOptions options;
    options.manningN=0; options.backgroundViscosity=0;
    options.windStressX=1; options.windStressY=-2;
    ASSERT_NO_FATAL_FAILURE(column.build(options));
    std::string message;
    ASSERT_TRUE(column.solver.step(10,message)) << message;
    // 10 s * tau * 100 m2 / (1000 kg/m3 * 160 m3).
    EXPECT_NEAR(column.state.u[0]-column.state.u[1],.00625,1e-12);
    EXPECT_NEAR(column.state.v[0]-column.state.v[1],-.0125,1e-12);
}

TEST(MomentumGeometry, BottomDragUsesBedContactAreaAndPreservesVectorDirection) {
    Column column;
    FVQual::Hydro::InternalModeOptions options;
    options.manningN=.3; options.backgroundViscosity=0;
    ASSERT_NO_FATAL_FAILURE(column.build(options));
    column.state.u={0,.6}; column.state.v={0,.8};
    const double drag=9.80665*.3*.3/std::cbrt(4.);
    const double decay=1./(1.+10*drag*20/80); // speed=1, Abed=20, Vbottom=80.
    std::string message;
    ASSERT_TRUE(column.solver.step(10,message)) << message;
    EXPECT_NEAR(column.state.u[1]-column.state.u[0],.6*decay,1e-12);
    EXPECT_NEAR(column.state.v[1]-column.state.v[0],.8*decay,1e-12);
}

TEST(MomentumGeometry, PrismaticViscosityRetainsTheAnalyticLimit) {
    Column column;
    FVQual::Hydro::InternalModeOptions options;
    options.manningN=0; options.backgroundViscosity=.1;
    ASSERT_NO_FATAL_FAILURE(column.build(options,false));
    column.state.u={1,-1};
    std::string message;
    ASSERT_TRUE(column.solver.step(10,message)) << message;
    EXPECT_NEAR(column.state.u[0]-column.state.u[1],4./3.,1e-12);
}

TEST(MomentumGeometry, StableScalarAndViscosityScalesAreIndependent) {
    Column baseline, noScalar, noViscosity;
    FVQual::Hydro::InternalModeOptions options;
    options.manningN=0; options.backgroundViscosity=.1;
    ASSERT_NO_FATAL_FAILURE(baseline.build(options));
    options.stableDiffusivityScale=0;
    ASSERT_NO_FATAL_FAILURE(noScalar.build(options));
    options.stableDiffusivityScale=1; options.stableViscosityScale=0;
    ASSERT_NO_FATAL_FAILURE(noViscosity.build(options));
    for (auto *column : {&baseline,&noScalar,&noViscosity}) {
        column->state.u={1,-1};
        std::string message;
        ASSERT_TRUE(column->solver.step(10,message)) << message;
    }
    EXPECT_EQ(baseline.state.u,noScalar.state.u);
    EXPECT_EQ(baseline.solver.eddyViscosity(),noScalar.solver.eddyViscosity());
    EXPECT_EQ(baseline.solver.eddyDiffusivity(),noViscosity.solver.eddyDiffusivity());
    EXPECT_NEAR(noViscosity.state.u[0]-noViscosity.state.u[1],2,1e-12);
    EXPECT_DOUBLE_EQ(noScalar.solver.eddyDiffusivity()[1],0);
    EXPECT_DOUBLE_EQ(noViscosity.solver.eddyViscosity()[1],0);
    // Convective adjustment survives both stable-mixing interventions.
    noScalar.state.speciesAt(0,0)=5; noScalar.state.speciesAt(0,1)=20;
    std::string message;
    ASSERT_TRUE(noScalar.solver.step(.1,message)) << message;
    EXPECT_DOUBLE_EQ(noScalar.solver.eddyDiffusivity()[1],options.convectiveMixing);
}

TEST(MomentumGeometry, InvalidAndUnsupportedDiagnosticScalesAreRefused) {
    FVQual::Hydro::InternalModeOptions options;
    Column column;
    ASSERT_NO_FATAL_FAILURE(column.build(options));
    std::string message;
    options.stableDiffusivityScale=-1;
    EXPECT_FALSE(column.solver.initialize(column.mesh,column.state,column.registry,options,message));
    options.stableDiffusivityScale=0; options.useTurbulenceClosure=true;
    EXPECT_FALSE(column.solver.initialize(column.mesh,column.state,column.registry,options,message));
    options.useTurbulenceClosure=false; options.stableViscosityScale=INFINITY;
    EXPECT_FALSE(column.solver.initialize(column.mesh,column.state,column.registry,options,message));
}

TEST(PressureGeometry, LinearDensityAtRestOverASteepBedStepIsBalanced) {
    for (int layers : {4,8,16,32}) {
        FVQual::Mesh::HorizontalMesh horizontal;
        FVQual::Mesh::LayeredMesh mesh;
        auto definition=FVQual::Testing::channelStrip(2,100,10);
        definition.nodeZ.assign(6,-20);
        std::string message;
        ASSERT_TRUE(horizontal.build(definition,message)) << message;
        horizontal.faceBed()[0]=-4;
        FVQual::Mesh::SigmaSpec spec;
        spec.layerCount=layers; spec.minLayerThickness=0;
        ASSERT_TRUE(mesh.build(horizontal,spec,0,message)) << message;
        std::vector<double> density(static_cast<size_t>(mesh.cellCount()));
        for (int64_t c=0;c<mesh.columnCount();++c)
            for (int k=0;k<layers;++k)
                density[static_cast<size_t>(mesh.cell(c,k))]=1000-.05*.5*(mesh.z(c,k)+mesh.z(c,k+1));
        FVQual::Hydro::HydrostaticPressureClosure closure;
        std::vector<double> ax,ay;
        closure.baroclinicAcceleration(mesh,density,1000,9.80665,ax,ay);
        for (double acceleration:ax) EXPECT_NEAR(acceleration,0,1e-12) << "layers=" << layers;
        for (double acceleration:ay) EXPECT_NEAR(acceleration,0,1e-12);
    }
}
