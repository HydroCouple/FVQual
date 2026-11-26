// SPDX-License-Identifier: MIT
#include "fvqual/hydro/internalmode.h"
#include "fvqual/kinetics/reactiontrace.h"
#include "testmeshes.h"
#include <gtest/gtest.h>
#include <cmath>
#include <numeric>
namespace {
using namespace FVQual;
double dot(std::span<const double> a,std::span<const double> b){return std::inner_product(a.begin(),a.end(),b.begin(),0.);}
struct ReactiveCase {
    Mesh::HorizontalMesh horizontal;Mesh::LayeredMesh mesh;State::SpeciesRegistry registry;
    State::StateStore state;Hydro::InternalModeSolver solver;int a=-1,b=-1;
    void build(double rate=.04) {
        std::string error;auto d=Testing::quadPlusTriangle();for(auto &x:d.nodeX)x*=10;for(auto &y:d.nodeY)y*=10;d.nodeZ.assign(d.nodeX.size(),-4);
        ASSERT_TRUE(horizontal.build(d,error))<<error;Mesh::SigmaSpec sigma;sigma.layerCount=2;sigma.minLayerThickness=0;
        ASSERT_TRUE(mesh.build(horizontal,sigma,0,error))<<error;
        for(int64_t c=0;c<mesh.columnCount();++c){const double area=horizontal.faceArea()[c];ASSERT_TRUE(mesh.setHypsography(c,std::vector<double>{-4,0},std::vector<double>{.3*area,area},error));}
        ASSERT_TRUE(mesh.updateSurface(std::vector<double>(mesh.columnCount(),0),error));
        ASSERT_GE(registry.addTemperature(12,error),0);
        a=registry.add({"A",State::SpeciesKind::Pollutant,"kg/m3","",1},error);
        b=registry.add({"B",State::SpeciesKind::MsxBulk,"kg/m3","",.2},error);state.allocate(mesh,registry);
        for(size_t i=0;i<state.cellCount();++i){state.speciesAt(a,i)=1.7-.22*i;state.speciesAt(b,i)=.25+.07*i;}
        Hydro::InternalModeOptions o;o.advectionScheme=Transport::AdvectionScheme::Upwind;o.backgroundDiffusivity=.02;
        o.neutralViscosity=0;o.backgroundViscosity=0;o.manningN=0;o.windStressX=.02;o.minTimeStep=10;o.maxTimeStep=10;
        ASSERT_TRUE(solver.initialize(mesh,state,registry,o,error))<<error;
        ASSERT_TRUE(solver.reactions().addConstant("k",rate,error));
        ASSERT_TRUE(solver.reactions().addRate("A","-k*A*(1 + B)",error));
        ASSERT_TRUE(solver.reactions().addRate("B","k*A*(1 + B)",error));ASSERT_TRUE(solver.compileReactions(error))<<error;
        Forcing::Source in;in.id="in";in.column=0;in.discharge=.01;in.temperature=12;in.concentration={0,.5,.2};ASSERT_GE(solver.sources().add(in,error),0);
        in.id="out";in.discharge=-.007;in.placement=Forcing::Placement::Bed;ASSERT_GE(solver.sources().add(in,error),0);
    }
    std::vector<double> run(int n=7) {
        std::vector<double> result;const auto record=[&]{for(int s:{a,b}){auto x=state.speciesBlock(s);result.insert(result.end(),x.begin(),x.end());}};
        std::string error;record();for(int i=0;i<n;++i){EXPECT_TRUE(solver.step(10,error))<<error;record();}return result;
    }
};
}
TEST(ReactiveTrace, NativeMixedMeshMultispeciesEpsilonSweepAndAdjoint) {
    ReactiveCase c;ASSERT_NO_FATAL_FAILURE(c.build());std::string error;ASSERT_TRUE(c.solver.beginReactiveTrace(100,error))<<error;
    const auto native=c.run();const auto trace=c.solver.reactiveTrace();const auto replay=trace.evaluate(trace.initial,trace.constants);
    ReactiveCase ordinary;ASSERT_NO_FATAL_FAILURE(ordinary.build());EXPECT_EQ(native,ordinary.run());
    ASSERT_EQ(native.size(),replay.size());for(size_t i=0;i<native.size();++i)EXPECT_NEAR(native[i],replay[i],2e-13);
    std::vector<double> v(trace.initial.size());for(size_t i=0;i<v.size();++i)v[i]=.1*std::sin(i+.2);
    const std::vector<double> p{.007};const auto dy=trace.jvp(v,p);
    for(double eps:{.01,.001,.0001,.00001}) {
        std::vector<std::vector<double>> runs;
        for(double sign:{1.,-1.}){ReactiveCase perturbed;ASSERT_NO_FATAL_FAILURE(perturbed.build(.04+sign*eps*p[0]));
            for(size_t i=0;i<v.size();++i)perturbed.state.speciesAt(i<4?perturbed.a:perturbed.b,i%4)+=sign*eps*v[i];runs.push_back(perturbed.run());}
        double worst=0;for(size_t i=0;i<dy.size();++i)worst=std::max(worst,std::abs(dy[i]-(runs[0][i]-runs[1][i])/(2*eps)));
        EXPECT_LT(worst,3e-3*eps*eps+2e-8)<<eps;
    }
    std::vector<double> w(dy.size());for(size_t i=0;i<w.size();++i)w[i]=std::cos(.7*i);
    const auto bar=trace.vjp(w);EXPECT_NEAR(dot(w,dy),dot(v,bar.initial)+dot(p,bar.constants),2e-12);
    EXPECT_GT(bar.peakTapeNodes,1000u); // adaptive stages are present in the actual reverse operator.
}
TEST(ReactionTrace, IndependentDecayAndScopedGeometryDerivatives) {
    ReactiveCase c;ASSERT_NO_FATAL_FAILURE(c.build());std::string error;
    Kinetics::ReactionNetwork network;ASSERT_TRUE(network.addConstant("k",.1,error));
    ASSERT_TRUE(network.addRate("A","-k*A",error));ASSERT_TRUE(network.addRate("B","k*A",error));ASSERT_TRUE(network.compile(c.registry,error));
    const auto trace=Kinetics::ReactionStep::capture(network,c.registry,c.state,50,{},{});
    std::vector<double> w(c.state.species.size(),0);w[c.a*4]=1;const auto bar=trace.vjp(w);
    EXPECT_NEAR(bar.species[c.a*4],std::exp(-5),3e-9);
    EXPECT_NEAR(bar.constants[0],-50*c.state.speciesAt(c.a,0)*std::exp(-5),2e-7);
    Kinetics::ReactionNetwork flux;ASSERT_TRUE(flux.addConstant("j",.02,error));
    ASSERT_TRUE(flux.addRate("A","j",Kinetics::RateScope::Bed,error));ASSERT_TRUE(flux.compile(c.registry,error));
    std::vector<double> thickness(4,2),ratio(4,.3);std::vector<unsigned char> flags(4,Kinetics::CellAtBed);
    Kinetics::DerivedFields fields;fields.thickness=thickness;fields.bedAreaPerVolume=ratio;
    auto f=Kinetics::ReactionStep::capture(flux,c.registry,c.state,10,fields,flags);const auto fb=f.vjp(w);
    EXPECT_NEAR(fb.geometry[4*4],.2,2e-14);EXPECT_NEAR(fb.constants[0],3,2e-14);EXPECT_DOUBLE_EQ(fb.geometry[2*4],0);
}
TEST(ReactiveTrace, CheckpointRestoreRefusalsAndUnchangedForwardState) {
    ReactiveCase c;ASSERT_NO_FATAL_FAILURE(c.build());std::string error;ASSERT_TRUE(c.solver.beginReactiveTrace(10,error));c.run(2);
    const auto saved=c.solver.checkpoint(error);ASSERT_TRUE(saved);const auto head=c.solver.reactiveTrace().head;c.run(2);
    const auto trace=c.solver.reactiveTrace();const auto forward=c.state.species;std::vector<double> w(trace.initial.size()*(trace.count+1),1);
    ASSERT_NO_THROW((void)trace.vjp(w));EXPECT_EQ(forward,c.state.species);
    ASSERT_TRUE(c.solver.restoreCheckpoint(*saved,error));EXPECT_EQ(c.solver.reactiveTrace().head,head);c.run(2);EXPECT_EQ(c.state.species,forward);
    c.solver.endReactiveTrace();ASSERT_TRUE(c.solver.reactions().addRate(State::SpeciesRegistry::kTemperature,"A",error));ASSERT_TRUE(c.solver.compileReactions(error));
    EXPECT_FALSE(c.solver.beginReactiveTrace(10,error));EXPECT_NE(error.find("temperature"),std::string::npos);
    Kinetics::ReactionNetwork nonsmooth;ASSERT_TRUE(nonsmooth.addRate("A","abs(B - .25)",error));ASSERT_TRUE(nonsmooth.compile(c.registry,error));
    State::StateStore branch=c.state;branch.speciesAt(c.b,0)=.25;
    auto bad=Kinetics::ReactionStep::capture(nonsmooth,c.registry,branch,1,{},{});
    EXPECT_THROW((void)bad.jvp(std::vector<double>(branch.species.size(),0),{}),std::invalid_argument);
    EXPECT_NO_THROW((void)bad.evaluate(bad.before(),{}));
}

TEST(ReactionTrace, AdaptiveControllerTieRefusesLocalGradientAndPreservesNativeEvaluation) {
    ReactiveCase c;ASSERT_NO_FATAL_FAILURE(c.build());std::string error;
    Kinetics::ReactionNetwork network;ASSERT_TRUE(network.addConstant("k",.1,error));
    ASSERT_TRUE(network.addRate("A","-k*A",error));ASSERT_TRUE(network.addRate("B","-k*B",error));ASSERT_TRUE(network.compile(c.registry,error));
    for(size_t i=0;i<c.state.cellCount();++i)c.state.speciesAt(c.b,i)=c.state.speciesAt(c.a,i);
    const auto trace=Kinetics::ReactionStep::capture(network,c.registry,c.state,50,{},{});
    EXPECT_NO_THROW((void)trace.evaluate(trace.before(),trace.constants()));
    EXPECT_THROW((void)trace.jvp(std::vector<double>(trace.before().size(),0),std::vector<double>{1}),std::invalid_argument);
}
