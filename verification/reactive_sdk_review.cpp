// SPDX-License-Identifier: MIT
#include "fvqual/component/fvqualcomponent.h"
#include "hydrocouplesdk/differentiation/differentiableworkflow.h"
#include "hydrocouplesdk/differentiation/differentials.h"
#include <nlohmann/json.hpp>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <cmath>
using namespace HydroCouple;using namespace HydroCouple::SDK;using namespace HydroCouple::SDK::Differentiation;
using FVQual::Component::FVQualComponent;using nlohmann::json;
namespace {
void require(bool ok,const std::string &m){if(!ok)throw std::runtime_error(m);}
IArgument *argument(FVQualComponent &c,const std::string &id){for(auto *a:c.arguments())if(a->id()==id)return a;throw std::runtime_error("missing argument "+id);}
IOutput *output(FVQualComponent &c,const std::string &id){for(auto *o:c.outputs())if(o->id()==id)return o;throw std::runtime_error("missing output "+id);}
void assign(FVQualComponent &c,const std::string &id,const json &j){std::string m;if(!argument(c,id)->initialize(j.dump(),IArgument::ArgumentInputType::JSON,m))throw std::runtime_error(id+": "+m+" ("+j.dump()+")");}
std::vector<double> values(IComponentDataItem *item){const auto shape=item->shape();std::vector<double> v(elementCount(shape));std::string m;require(item->getValuesInto(hostBuffer(v,shape),std::vector<int64_t>(shape.size(),0),shape,&m),m);return v;}
json run(const json &spec,const std::filesystem::path &dir,CheckpointPolicy policy,double k,bool reverse) {
    FVQualComponent c("reservoir");DifferentiableWorkflow workflow("reactive-review");std::string m;
    const auto &g=spec.at("geometry");HydroCouple::SDK::IO::MeshDefinition mesh;mesh.meshName="synthetic-bellinge-receiver";
    mesh.nodeX=g.at("x").get<std::vector<double>>();mesh.nodeY=g.at("y").get<std::vector<double>>();mesh.nodeZ=g.at("z").get<std::vector<double>>();
    mesh.faceNodeOffsets=g.at("offsets").get<std::vector<int64_t>>();mesh.faceNodes=g.at("nodes").get<std::vector<int64_t>>();
    require(dynamic_cast<HydroCouple::SDK::Spatial::PolyhedralSurfaceArgument *>(argument(c,"mesh"))->setMesh(std::move(mesh),m),m);
    assign(c,"geometry_counts",{3,3});assign(c,"geometry_values",{100.,100.,5.,10.,1080.,spec.at("duration_s")});
    assign(c,"epoch_julian_day",json::array({spec.at("epoch_julian_day")}));assign(c,"conditional_species_derivatives",{1});
    assign(c,"maximum_time_step_s",{300.});assign(c,"kinetics_file",{(dir/"reaction.rxn").string()});
    assign(c,"recorded_inflows",{{"@uri",(dir/"native_inflows.json").string()}});assign(c,"reaction_constants",{{"k",k}});
    require(workflow.addModelComponent(&c,nullptr,&m),m);workflow.setCheckpointPolicy(policy);workflow.initialize();
    require(c.status()==IModelComponent::ComponentStatus::Initialized,"component initialization failed");
    const auto problems=workflow.validate();require(problems.empty(),problems.empty()?"":problems.front());workflow.prepare();
    require(c.status()!=IModelComponent::ComponentStatus::Failed,"component prepare failed");
    double loss=0;json history=json::array();
    for(size_t i=0;i<spec.at("frames").size();++i) {
        workflow.update();require(c.status()!=IModelComponent::ComponentStatus::Failed,"native Bellinge update failed");
        auto a=values(output(c,"A")),b=values(output(c,"B"));for(size_t cell=0;cell<a.size();++cell)loss+=a[cell]+.3*b[cell];
        history.push_back({{"time_s",c.currentTime()},{"A",a},{"B",b}});
    }
    json result={{"loss",loss},{"history",history},{"stored_checkpoints",workflow.storedCheckpoints()},{"steps",workflow.recordedSteps()}};
    if(reverse) {
        const auto before=c.state().species;Gradients gradients;
        require(workflow.backward([&](int,Seeds &seeds){require(seeds.add(output(c,"A"),std::vector<double>(9,1)),"A seeds");require(seeds.add(output(c,"B"),std::vector<double>(9,.3)),"B seeds");},gradients,m),m);
        require(before==c.state().species,"SDK reverse mutated final native state");
        const auto *p=gradients.find(&c,argument(c,"reaction_constants"));require(p,"missing SDK kinetic gradient");
        const auto names=dynamic_cast<IIdBasedComponentDataItem *>(argument(c,"reaction_constants"))->identifiers();
        for(size_t i=0;i<names.size();++i)result["constant_gradient"][names[i]]=p->values[i];
        for(const std::string id:{"A","B"}){const auto *s=gradients.findState(&c,output(c,id));require(s,"missing SDK initial gradient");result["initial_gradient"][id]=s->values;}
        result["replays"]=workflow.replays();result["peak_checkpoints"]=workflow.peakCheckpoints();result["refused_releases"]=workflow.refusedReleases();
        require(workflow.refusedReleases()==0,"SDK checkpoint release failure");
    }
    workflow.finish();return result;
}
}
int main(int argc,char **argv) {
    try {
        if(argc!=2)throw std::invalid_argument("Usage: fvqual-reactive-sdk-review review-directory");
        const auto dir=std::filesystem::absolute(argv[1]);json spec;std::ifstream(dir/"sdk_case.json")>>spec;
        const double k=spec.at("rate_s_1");const auto every=run(spec,dir,CheckpointPolicy::EveryStep,k,true),log=run(spec,dir,CheckpointPolicy::Logarithmic,k,true);
        const double a=every.at("constant_gradient").at("k"),b=log.at("constant_gradient").at("k");require(std::abs(a-b)<1e-10,"SDK replay gradients disagree");
        require(every.at("history")==log.at("history"),"SDK forward results differ between policies");json fd=json::array();
        for(double eps:{.01,.001,.0001}){const auto plus=run(spec,dir,CheckpointPolicy::EveryStep,k+eps*3e-6,false),minus=run(spec,dir,CheckpointPolicy::EveryStep,k-eps*3e-6,false);
            const double error=std::abs((plus.at("loss").get<double>()-minus.at("loss").get<double>())/(2*eps)-3e-6*a);
            require(error<1e-8,"independent SDK whole-run directional derivative failed");fd.push_back({{"epsilon",eps},{"absolute_loss_direction_error",error}});}
        const json report={{"scope","native FVQual conditional species differentiation through HydroCoupleSDK; prescribed Bellinge water/heat"},{"every_step",every},{"logarithmic",log},{"finite_difference_sweep",fd},{"maximum_policy_gradient_difference",std::abs(a-b)}};
        std::ofstream(dir/"sdk_qualification.json")<<report.dump(2)<<'\n';std::cout<<"Native SDK review passed; "<<every.at("stored_checkpoints")<<" versus "<<log.at("stored_checkpoints")<<" stored composition checkpoints\n";return 0;
    }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
