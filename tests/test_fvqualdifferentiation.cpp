// SPDX-License-Identifier: MIT
#include "fvqual/component/fvqualcomponent.h"
#include "hydrocouplesdk/differentiation/differentiableworkflow.h"
#include "hydrocouplesdk/differentiation/differentials.h"
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <numeric>
namespace {
using namespace HydroCouple;using namespace HydroCouple::SDK;using namespace HydroCouple::SDK::Differentiation;
using FVQual::Component::FVQualComponent;using R=DifferentialRole;
IArgument *arg(FVQualComponent &c,const std::string &id){for(auto *a:c.arguments())if(a->id()==id)return a;return nullptr;}
IOutput *out(FVQualComponent &c,const std::string &id){for(auto *a:c.outputs())if(a->id()==id)return a;return nullptr;}
void assign(FVQualComponent &c,const std::string &id,const nlohmann::json &j){std::string m;ASSERT_NE(arg(c,id),nullptr);ASSERT_TRUE(arg(c,id)->initialize(j.dump(),IArgument::ArgumentInputType::JSON,m))<<m;}
std::vector<double> values(IComponentDataItem *item){const auto shape=item->shape();std::vector<double> v(elementCount(shape));std::string m;EXPECT_TRUE(item->getValuesInto(hostBuffer(v,shape),std::vector<int64_t>(shape.size(),0),shape,&m))<<m;return v;}
void configure(FVQualComponent &c,double k=.015,double initial=1.4,bool enabled=true) {
    const auto *info=::testing::UnitTest::GetInstance()->current_test_info();
    const auto dir=std::filesystem::path(FVQual_TEST_OUTPUT_DIR)/"sdk_reactive"/(info?info->name():"review");std::filesystem::create_directories(dir);
    const auto file=dir/(c.id()+".rxn");std::ofstream deck(file);deck<<std::setprecision(17)<<"[SPECIES]\nA kg/m3 "<<initial<<"\nB kg/m3 .25\n[CONSTANTS]\nk "<<k<<"\n[RATES]\nA WATER -k*A*(1+B)\nB WATER k*A*(1+B)\n";deck.close();
    assign(c,"geometry_counts",{2,2});assign(c,"geometry_values",{20.,10.,4.,12.,5.,80.});
    assign(c,"maximum_time_step_s",{1.});assign(c,"kinetics_file",{file.string()});assign(c,"conditional_species_derivatives",{enabled?1:0});
}
struct Rig {
    FVQualComponent c;DifferentiableWorkflow workflow;
    Rig(std::string id,double k=.015,double initial=1.4,CheckpointPolicy policy=CheckpointPolicy::EveryStep):c(id),workflow("workflow-"+id) {
        configure(c,k,initial);std::string m;EXPECT_TRUE(workflow.addModelComponent(&c,nullptr,&m))<<m;workflow.setCheckpointPolicy(policy);
        workflow.initialize();EXPECT_NE(workflow.status(),IWorkflowComponent::WorkflowStatus::Failed);
        EXPECT_TRUE(workflow.validate().empty());workflow.prepare();EXPECT_NE(c.status(),IModelComponent::ComponentStatus::Failed);
    }
    double run() {
        double loss=0;for(int step=1;step<=16;++step){workflow.update();EXPECT_NE(c.status(),IModelComponent::ComponentStatus::Failed);
            const auto a=values(out(c,"A")),b=values(out(c,"B"));for(size_t i=0;i<a.size();++i)loss+=(i+1)*a[i]*a[i]+.3*b[i];}
        return loss;
    }
    SeedFunction seeds(){return [this](int,Seeds &s){auto a=values(out(c,"A"));std::vector<double>b(a.size(),.3);for(size_t i=0;i<a.size();++i)a[i]*=2*(i+1);EXPECT_TRUE(s.add(out(c,"A"),a));EXPECT_TRUE(s.add(out(c,"B"),b));};}
};
double loss(double k,double initial){Rig r("finite-difference",k,initial);return r.run();}
}
TEST(FVQualDifferentiation, NativeSdkEveryStepAndLogarithmicReplayMatchIndependentRuns) {
    std::vector<double> reference;double initialGradient=0;
    for(auto policy:{CheckpointPolicy::EveryStep,CheckpointPolicy::Logarithmic}) {
        Rig r(policy==CheckpointPolicy::EveryStep?"every":"log",.015,1.4,policy);r.run();
        ASSERT_TRUE(r.workflow.differentiabilityProblems().empty());ASSERT_EQ(r.workflow.recordedSteps(),16);
        EXPECT_TRUE(r.c.capabilities().count(Capability::Differentiable));const auto before=r.c.state().species;
        const auto checkpoints=r.workflow.storedCheckpoints();if(policy==CheckpointPolicy::Logarithmic)EXPECT_LT(checkpoints,10);
        Gradients g;std::string m;ASSERT_TRUE(r.workflow.backward(r.seeds(),g,m))<<m;EXPECT_EQ(before,r.c.state().species);EXPECT_NEAR(r.c.currentTime(),80,1e-12);
        const auto *parameter=g.find(&r.c,arg(r.c,"reaction_constants"));ASSERT_NE(parameter,nullptr);ASSERT_EQ(parameter->values.size(),1u);
        const auto *initial=g.findState(&r.c,out(r.c,"A"));ASSERT_NE(initial,nullptr);const double sum=std::accumulate(initial->values.begin(),initial->values.end(),0.);
        if(reference.empty()){reference=parameter->values;initialGradient=sum;}else{EXPECT_NEAR(reference[0],parameter->values[0],1e-10);EXPECT_NEAR(initialGradient,sum,1e-11);}
        EXPECT_GT(r.workflow.replays(),0);EXPECT_EQ(r.workflow.refusedReleases(),0);
        Gradients again;ASSERT_TRUE(r.workflow.backward(r.seeds(),again,m))<<m;EXPECT_NEAR(again.find(&r.c,arg(r.c,"reaction_constants"))->values[0],reference[0],1e-10);
    }
    for(double eps:{1e-3,1e-4,1e-5}) {
        const double fd=(loss(.015+eps*.005,1.4+eps*.2)-loss(.015-eps*.005,1.4-eps*.2))/(2*eps);
        EXPECT_NEAR(fd,.005*reference[0]+.2*initialGradient,1e-4*eps+3e-7)<<eps;
    }
}
TEST(FVQualDifferentiation, StepProductsRespectRolesStridesAndAtomicFailure) {
    FVQualComponent c("products");configure(c);c.initialize();c.prepare();c.update();
    ASSERT_EQ(c.differentiableInputs().size(),0u);ASSERT_EQ(c.differentiableOutputs().size(),2u);
    auto *a=out(c,"A");auto *k=arg(c,"reaction_constants");const auto shape=a->shape(),pshape=k->shape();
    std::vector<double> v(4,.2),dp(1,.005),padded(8,91),w(4,1),bar(4,0),pk(1,0);
    std::vector<int64_t> stride(shape.size());int64_t bytes=2*sizeof(double);for(size_t i=shape.size();i>0;--i){stride[i-1]=bytes;bytes*=shape[i-1];}
    auto outputBuffer=hostBuffer(padded,shape);outputBuffer.stridesBytes=stride.data();
    std::vector<DifferentialEntry> seeds{{a,R::StateBefore,hostBuffer(v,shape)},{k,R::Argument,hostBuffer(dp,pshape)}};
    std::vector<DifferentialEntry> results{{a,R::Output,outputBuffer}};std::string message;ASSERT_TRUE(c.jvp(seeds,results,&message))<<message;
    for(size_t i=1;i<padded.size();i+=2)EXPECT_EQ(padded[i],91);
    std::vector<DifferentialEntry> reverseSeeds{{a,R::Output,hostBuffer(w,shape)}};
    std::vector<DifferentialEntry> reverseResults{{a,R::StateBefore,hostBuffer(bar,shape)},{k,R::Argument,hostBuffer(pk,pshape)}};
    ASSERT_TRUE(c.vjp(reverseSeeds,reverseResults,&message))<<message;double tangent=0;for(size_t i=0;i<4;++i)tangent+=padded[2*i];
    EXPECT_NEAR(tangent,.2*std::accumulate(bar.begin(),bar.end(),0.)+.005*pk[0],2e-13);
    const auto saved=padded;results.push_back({out(c,"temperature"),R::Output,outputBuffer});EXPECT_FALSE(c.jvp(seeds,results,&message));EXPECT_EQ(padded,saved);
    EXPECT_EQ(c.currentTime(),5);
}
TEST(FVQualDifferentiation, DefaultOffUnsupportedConfigurationsAndReinitialization) {
    FVQualComponent off("off");configure(off,.015,1.4,false);off.initialize();off.prepare();EXPECT_FALSE(off.capabilities().count(Capability::Differentiable));
    std::string message;EXPECT_FALSE(off.jvp({}, {},&message));EXPECT_TRUE(off.differentiableOutputs().empty());
    FVQualComponent bad("diffusion");configure(bad);assign(bad,"horizontal_diffusivity_m2_s",{.1});bad.initialize();bad.prepare();EXPECT_EQ(bad.status(),IModelComponent::ComponentStatus::Failed);
    FVQualComponent unknown("unknown");configure(unknown);assign(unknown,"reaction_constants",{{"typo",1}});unknown.initialize();unknown.prepare();EXPECT_EQ(unknown.status(),IModelComponent::ComponentStatus::Failed);
    FVQualComponent c("reset");configure(c);c.initialize();c.prepare();c.update();std::string token;ASSERT_TRUE(c.saveState(token,message));c.update();ASSERT_TRUE(c.restoreState(token,message));
    EXPECT_FALSE(c.jvp({}, {},&message));c.update();EXPECT_TRUE(c.jvp({}, {},&message))<<message;ASSERT_TRUE(c.releaseState(token,message));
}
