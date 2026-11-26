// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
#include "fvqual/hydro/internalmode.h"
#include "fvqual/coupling/intervaltransfer.h"
#include "testmeshes.h"
#include <gtest/gtest.h>
#include <cmath>
#include <numeric>

namespace {
using namespace FVQual;
struct Case {
    Mesh::HorizontalMesh horizontal;
    Mesh::LayeredMesh mesh;
    State::SpeciesRegistry registry;
    State::StateStore state;
    Hydro::InternalModeSolver solver;
    int tracer=-1;
    void build(double scale=1, bool mixed=true, bool flow=true) {
        std::string error;
        auto d=mixed?Testing::quadPlusTriangle():Testing::channelStrip(1,10,10);
        d.nodeZ.assign(d.nodeX.size(),-4);
        if (mixed) {for (auto &x:d.nodeX) x*=10;for (auto &y:d.nodeY) y*=10;}
        ASSERT_TRUE(horizontal.build(d,error)) << error;
        Mesh::SigmaSpec sigma;sigma.layerCount=2;sigma.minLayerThickness=0;
        ASSERT_TRUE(mesh.build(horizontal,sigma,0,error)) << error;
        for (int64_t c=0;c<mesh.columnCount();++c) {
            const double a=horizontal.faceArea()[c];
            ASSERT_TRUE(mesh.setHypsography(c,std::vector<double>{-4,0},std::vector<double>{.2*a,a},error)) << error;
        }
        ASSERT_TRUE(mesh.updateSurface(std::vector<double>(mesh.columnCount(),0),error)) << error;
        ASSERT_GE(registry.addTemperature(10,error),0);
        tracer=registry.add({"tracer",State::SpeciesKind::Pollutant,"kg/m3","",0},error);ASSERT_GE(tracer,0);
        state.allocate(mesh,registry);
        for (size_t i=0;i<static_cast<size_t>(mesh.cellCount());++i) state.speciesAt(tracer,i)=1.2-.19*i;
        Hydro::InternalModeOptions options;options.advectionScheme=Transport::AdvectionScheme::Upwind;
        options.neutralViscosity=0;options.backgroundViscosity=0;options.backgroundDiffusivity=.1;
        options.stableDiffusivityScale=scale;options.maxTimeStep=10;options.minTimeStep=10;
        options.manningN=0;options.windStressX=flow?.02:0;
        ASSERT_TRUE(solver.initialize(mesh,state,registry,options,error)) << error;
        if (flow) {
            Forcing::Source in;in.id="incoming";in.column=0;in.discharge=.01;in.temperature=10;
            in.placement=Forcing::Placement::Surface;in.concentration={0,.3};
            ASSERT_GE(solver.sources().add(in,error),0);
            auto out=in;out.id="withdrawal";out.discharge=-.007;out.placement=Forcing::Placement::Bed;
            ASSERT_GE(solver.sources().add(out,error),0);
        }
    }
    std::vector<double> run(int steps=8) {
        std::string error;std::vector<double> result;
        const auto record=[&] {const auto x=state.speciesBlock(tracer);result.insert(result.end(),x.begin(),x.end());};
        record();for (int k=0;k<steps;++k) {EXPECT_TRUE(solver.step(10,error)) << error;record();}
        return result;
    }
};
double dot(std::span<const double> a,std::span<const double> b) {
    return std::inner_product(a.begin(),a.end(),b.begin(),0.0);
}
}

TEST(ScalarTrace, ActualMixedMeshForwardAndWholeRunDirectionalSweep) {
    Case baseline;ASSERT_NO_FATAL_FAILURE(baseline.build());std::string error;
    ASSERT_TRUE(baseline.solver.beginScalarTrace(baseline.tracer,100,error)) << error;
    const auto forward=baseline.run();const auto trace=baseline.solver.scalarTrace();
    ASSERT_EQ(trace.count,8u);
    const auto sources=trace.sourceValues();const auto replay=trace.evaluate(trace.initial,sources);
    for (size_t i=0;i<forward.size();++i) EXPECT_NEAR(replay[i],forward[i],2e-15);
    std::vector<double> v={.2,-.1,.4,.05},sourceDirection(sources.size(),0);
    for (size_t i=0;i<sourceDirection.size();i+=2) sourceDirection[i]=.13;
    const auto tangent=trace.jvp(v,sourceDirection,.37);
    for (double epsilon:{1e-2,1e-3,1e-4,1e-5}) {
        std::vector<std::vector<double>> outputs;
        for (double sign:{1.,-1.}) {
            Case perturbed;ASSERT_NO_FATAL_FAILURE(perturbed.build(1+sign*epsilon*.37));
            for (size_t i=0;i<v.size();++i) perturbed.state.speciesAt(perturbed.tracer,i)+=sign*epsilon*v[i];
            ASSERT_TRUE(perturbed.solver.sources().setInflowValue(0,perturbed.tracer,0,.3+sign*epsilon*.13,error));
            outputs.push_back(perturbed.run());
        }
        double maxError=0;
        for (size_t i=0;i<tangent.size();++i) maxError=std::max(maxError,std::abs((outputs[0][i]-outputs[1][i])/(2*epsilon)-tangent[i]));
        EXPECT_LT(maxError,epsilon*epsilon*.03+2e-9) << epsilon;
    }
    std::vector<double> w(forward.size());for(size_t i=0;i<w.size();++i) w[i]=std::sin(.7*i);
    const auto bar=trace.vjp(w);
    EXPECT_NEAR(dot(w,tangent),dot(v,bar.initial)+dot(sourceDirection,bar.sources)+.37*bar.diffusionMultiplier,2e-14);
    // Deliberate bad sign in recorded source mapping must fail the independent
    // whole-run oracle, even though a consistently wrong JVP/VJP can pass a dot test.
    auto bad=trace;auto node=std::make_shared<Transport::ScalarTraceNode>(*trace.head);
    ASSERT_FALSE(node->step.inflow.empty());node->step.inflow[0].value*=-1;bad.head=node;
    EXPECT_GT(std::abs(dot(w,bad.jvp(v,sourceDirection,.37))-dot(w,tangent)),1e-7);
}

TEST(ScalarTrace, IndependentTwoBoxImplicitDerivativeAndTangentMass) {
    Case c;ASSERT_NO_FATAL_FAILURE(c.build(1,false,false));std::string error;
    c.state.speciesAt(c.tracer,0)=1;c.state.speciesAt(c.tracer,1)=0;
    ASSERT_TRUE(c.solver.beginScalarTrace(c.tracer,10,error));c.run(1);
    const auto &t=c.solver.scalarTrace();
    const auto tangent=t.jvp(std::vector<double>{0,0},{},1);
    const double damping=10*3*(1./160+1./80);
    const double difference=-damping/std::pow(1+damping,2);
    EXPECT_NEAR(tangent[2],difference/3,2e-15);
    EXPECT_NEAR(tangent[3],-2*difference/3,2e-15);
    EXPECT_NEAR(160*tangent[2]+80*tangent[3],0,2e-13);
    const auto v=t.jvp(std::vector<double>{.2,-.1},{},0);
    EXPECT_NEAR(160*v[2]+80*v[3],160*.2-80*.1,2e-13);
    // Inventory observation is volume-weighted, requiring a true transpose
    // on the unequal-volume (nonsymmetric concentration) diffusion matrix.
    const auto b=t.vjp(std::vector<double>{0,0,160,80});
    EXPECT_NEAR(b.initial[0],160,2e-13);EXPECT_NEAR(b.initial[1],80,2e-13);
    EXPECT_NEAR(b.diffusionMultiplier,0,2e-13);
}

TEST(ScalarTrace, CheckpointHistoryReplayAndSnapshotLifetime) {
    Transport::ScalarTrajectory snapshot;
    {
        Case c;ASSERT_NO_FATAL_FAILURE(c.build());std::string error;
        ASSERT_TRUE(c.solver.beginScalarTrace(c.tracer,100,error));c.run(3);
        auto saved=c.solver.checkpoint(error);ASSERT_TRUE(saved);
        const auto head=c.solver.scalarTrace().head;c.run(2);snapshot=c.solver.scalarTrace();
        const auto result=snapshot.jvp(std::vector<double>{1,0,0,0},std::vector<double>(10,0),.1);
        ASSERT_TRUE(c.solver.restoreCheckpoint(*saved,error));EXPECT_EQ(c.solver.scalarTrace().head,head);
        EXPECT_EQ(c.solver.scalarTrace().count,3u);c.run(2);
        EXPECT_EQ(result,c.solver.scalarTrace().jvp(std::vector<double>{1,0,0,0},std::vector<double>(10,0),.1));
        EXPECT_EQ(snapshot.head->step.after,c.solver.scalarTrace().head->step.after);
        const auto state=c.state.species;
        EXPECT_THROW((void)snapshot.jvp(std::vector<double>{1},{},0),std::invalid_argument);
        EXPECT_EQ(c.state.species,state);c.solver.endScalarTrace();
    }
    EXPECT_EQ(snapshot.count,5u);
    EXPECT_FALSE(snapshot.evaluate(snapshot.initial,snapshot.sourceValues()).empty());
}

TEST(ScalarTrace, UnsupportedMutationAndStepLimitRefusedBeforeAdvance) {
    Case c;ASSERT_NO_FATAL_FAILURE(c.build());std::string error;
    EXPECT_FALSE(c.solver.beginScalarTrace(0,10,error));
    ASSERT_TRUE(c.solver.beginScalarTrace(c.tracer,1,error));c.run(1);
    const auto state=c.state.species;const double time=c.solver.time();
    EXPECT_FALSE(c.solver.step(10,error));EXPECT_EQ(c.state.species,state);EXPECT_EQ(c.solver.time(),time);
    ASSERT_TRUE(c.solver.beginScalarTrace(c.tracer,10,error));c.state.speciesAt(c.tracer,0)+=.01;
    EXPECT_FALSE(c.solver.step(10,error));EXPECT_NE(error.find("changed outside"),std::string::npos);
    ASSERT_TRUE(c.solver.beginScalarTrace(c.tracer,10,error));
    ASSERT_TRUE(c.solver.settling().setVelocity(c.tracer,.01,c.registry,error));
    const auto changed=c.state.species;EXPECT_FALSE(c.solver.step(10,error));EXPECT_EQ(c.state.species,changed);
}

TEST(ScalarTrace, CoupledIntervalLoadGradientAndRollbackHistory) {
    Case c;ASSERT_NO_FATAL_FAILURE(c.build(1,true,false));std::string error;
    ASSERT_TRUE(c.solver.beginScalarTrace(c.tracer,100,error));
    Coupling::IntervalTransferReceiver receiver(c.solver,c.registry,c.mesh);
    Coupling::TransferFrame frame;frame.id="event";frame.startSeconds=0;frame.endSeconds=20;
    Coupling::InflowTransfer inflow;inflow.column=0;inflow.volumeM3=.1;
    inflow.heatJ=.1*1000*4186*10;inflow.speciesMassKg["tracer"]=.02;frame.inflows={inflow};
    Coupling::TransferReceipt receipt;ASSERT_TRUE(receiver.apply(frame,receipt,error)) << error;
    const auto trace=c.solver.scalarTrace();ASSERT_EQ(trace.count,2u);
    EXPECT_EQ(trace.steps()[0]->sourceIds[0],"event:0");
    std::vector<double> zeros(trace.initial.size(),0);
    const auto tangent=trace.jvp(zeros,trace.sourceValues());
    const size_t n=zeros.size();double mass=0;
    for(size_t i=0;i<n;++i) mass+=tangent[2*n+i]*c.mesh.cellVolume()[i];
    EXPECT_NEAR(mass,.02,1e-13);
    ASSERT_TRUE(receiver.apply(frame,receipt,error));EXPECT_TRUE(receipt.replayed);
    EXPECT_EQ(c.solver.scalarTrace().head,trace.head);
}
