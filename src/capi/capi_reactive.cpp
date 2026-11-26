// SPDX-License-Identifier: MIT
#include "capi_internal.hpp"
#include "fvqual/capi/fvqual_reactive.h"
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <algorithm>
struct FVQ_ReactiveTrace_ {FVQual::Transport::ReactiveTrajectory trajectory;std::string metadata;};
namespace {
template<class F> FVQ_ErrorCode checked(F f){return FVQual::CApi::guard([&]()->FVQ_ErrorCode{try{f();return FVQual::CApi::succeed();}catch(const std::invalid_argument &e){return FVQual::CApi::fail(FVQ_ERR_INVALID_ARGUMENT,e.what());}});}
void array(const void *p,int64_t count,size_t n){if(count<0 || static_cast<size_t>(count)!=n || (!p&&n))throw std::invalid_argument("reactive trace array is null or has the wrong length");}
const auto &trace(const FVQ_ReactiveTrace *p){if(!p)throw std::invalid_argument("reactive trace is null");return p->trajectory;}
void dimensions(const FVQ_ReactiveTrace *p,const void *initial,int64_t ni,const void *parameters,int64_t np,const void *states,int64_t ns){const auto &t=trace(p);array(initial,ni,t.initial.size());array(parameters,np,t.constants.size());array(states,ns,t.initial.size()*(t.count+1));}
}
extern "C" {
FVQ_ErrorCode fvq_sim_begin_reactive_trace(FVQ_Simulation *s,int64_t limit){return checked([&]{
    if(!s||!s->initialized||limit<1||limit>10000)throw std::invalid_argument("reactive recording requires initialized simulation and limit 1..10000");
    std::string m;if(!s->solver.beginReactiveTrace(static_cast<size_t>(limit),m))throw std::invalid_argument(m);
});}
FVQ_ErrorCode fvq_sim_end_reactive_trace(FVQ_Simulation *s){return checked([&]{if(!s||!s->initialized)throw std::invalid_argument("initialized simulation required");s->solver.endReactiveTrace();});}
FVQ_ErrorCode fvq_sim_capture_reactive_trace(const FVQ_Simulation *s,FVQ_ReactiveTrace **out){return checked([&]{
    if(!s||!s->initialized||!out||s->solver.reactiveTrace().initial.empty())throw std::invalid_argument("begin reactive recording before capture");
    auto p=std::make_unique<FVQ_ReactiveTrace>();p->trajectory=s->solver.reactiveTrace();std::vector<std::string> ids,units;
    for(int i:p->trajectory.species){ids.push_back(s->registry.definition(i).id);units.push_back(s->registry.definition(i).units);}
    p->metadata=nlohmann::json{{"schema","fvqual.reactive-trace/1.0"},{"species",ids},{"species_units",units},{"parameter_names",p->trajectory.constantNames},
        {"scope","conditional on prescribed hydraulics, heat, geometry, mixing and source composition"},{"kinetics","native adaptive Cash-Karp including smooth controller dependence"}}.dump();*out=p.release();
});}
void fvq_reactive_trace_destroy(FVQ_ReactiveTrace *p){delete p;}
FVQ_ErrorCode fvq_reactive_trace_size(const FVQ_ReactiveTrace *p,FVQ_ReactiveSize *out){return checked([&]{const auto &t=trace(p);if(!out)throw std::invalid_argument("reactive size output is null");*out={static_cast<int64_t>(t.cells),static_cast<int64_t>(t.species.size()),static_cast<int64_t>(t.count),static_cast<int64_t>(t.constants.size())};});}
const char *fvq_reactive_trace_metadata(const FVQ_ReactiveTrace *p){if(!p){FVQual::CApi::fail(FVQ_ERR_INVALID_ARGUMENT,"reactive trace is null");return nullptr;}FVQual::CApi::succeed();return p->metadata.c_str();}
FVQ_ErrorCode fvq_reactive_trace_baseline(const FVQ_ReactiveTrace *p,double *states,int64_t ns,double *parameters,int64_t np,double *times,int64_t nt){return checked([&]{
    const auto &t=trace(p);array(states,ns,t.initial.size()*(t.count+1));array(parameters,np,t.constants.size());array(times,nt,t.count+1);
    const auto history=t.steps();
    std::copy(t.initial.begin(),t.initial.end(),states);std::copy(t.constants.begin(),t.constants.end(),parameters);times[0]=t.startTime;size_t offset=t.initial.size(),i=1;
    for(const auto *step:history){std::copy(step->after.begin(),step->after.end(),states+offset);offset+=t.initial.size();times[i++]=step->transport[0].time+step->transport[0].dt;}
});}
FVQ_ErrorCode fvq_reactive_trace_evaluate(const FVQ_ReactiveTrace *p,const double *initial,int64_t ni,const double *parameters,int64_t np,double *states,int64_t ns){return checked([&]{dimensions(p,initial,ni,parameters,np,states,ns);const auto y=p->trajectory.evaluate({initial,static_cast<size_t>(ni)},{parameters,static_cast<size_t>(np)});std::copy(y.begin(),y.end(),states);});}
FVQ_ErrorCode fvq_reactive_trace_jvp(const FVQ_ReactiveTrace *p,const double *initial,int64_t ni,const double *parameters,int64_t np,double *states,int64_t ns){return checked([&]{dimensions(p,initial,ni,parameters,np,states,ns);const auto y=p->trajectory.jvp({initial,static_cast<size_t>(ni)},{parameters,static_cast<size_t>(np)});std::copy(y.begin(),y.end(),states);});}
FVQ_ErrorCode fvq_reactive_trace_vjp(const FVQ_ReactiveTrace *p,const double *states,int64_t ns,double *initial,int64_t ni,double *parameters,int64_t np,int64_t *peak){return checked([&]{dimensions(p,initial,ni,parameters,np,states,ns);if(!peak)throw std::invalid_argument("peak tape node output is null");const auto bar=p->trajectory.vjp({states,static_cast<size_t>(ns)});std::copy(bar.initial.begin(),bar.initial.end(),initial);std::copy(bar.constants.begin(),bar.constants.end(),parameters);*peak=static_cast<int64_t>(bar.peakTapeNodes);});}
FVQ_ErrorCode fvq_sim_set_reaction_constant(FVQ_Simulation *s,const char *name,double value){return checked([&]{if(!s||!s->initialized||!name)throw std::invalid_argument("initialized simulation and constant name required");std::string m;if(!s->solver.reactions().setConstant(name,value,m))throw std::invalid_argument(m);});}
const char *fvq_sim_species_id(const FVQ_Simulation *s,int species){if(!s||species<0||species>=s->registry.count()){FVQual::CApi::fail(FVQ_ERR_INVALID_ARGUMENT,"invalid species index");return nullptr;}FVQual::CApi::succeed();return s->registry.definition(species).id.c_str();}
}
