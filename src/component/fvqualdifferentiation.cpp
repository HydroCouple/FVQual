// SPDX-License-Identifier: MIT
#include "fvqual/component/fvqualcomponent.h"
#include "hydrocouplesdk/differentiation/differentials.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
namespace FVQual::Component {
using namespace HydroCouple;
namespace D=HydroCouple::SDK::Differentiation;
bool FVQualComponent::configureDifferentiation(std::string &message) {
    auto &network=m_solver.reactions();const auto names=network.constantNames();
    const auto keys=m_reactionParameters->identifiers();std::map<std::string,double> overrides;
    for(size_t i=0;i<keys.size();++i) {
        if(std::find(names.begin(),names.end(),keys[i])==names.end() || !std::isfinite((*m_reactionParameters)[i])) {
            message="reaction_constants names must exist in the deck and values must be finite";return false;
        }
        overrides[keys[i]]=(*m_reactionParameters)[i];
    }
    m_reactionParameters->setIdentifiers(names);
    for(size_t i=0;i<names.size();++i)(*m_reactionParameters)[i]=overrides.count(names[i])?overrides.at(names[i]):network.constantValues()[i];
    if(!applyReactionParameters(message))return false;
    if(!(*m_differentiationArgument)[0])return true;
    if(!m_solver.beginReactiveTrace(10000,message))return false;
    for(int s:m_solver.reactiveTrace().species) {
        const auto id=m_registry.definition(s).id;
        const auto it=std::find_if(m_fields.begin(),m_fields.end(),[&](const auto &f){return f->id()==id;});
        if(it==m_fields.end()){m_solver.endReactiveTrace();message="missing differentiable species output";return false;}
        m_differentiableFields.push_back(it->get());
    }
    m_solver.endReactiveTrace();m_differentiationReady=true;message.clear();return true;
}
bool FVQualComponent::applyReactionParameters(std::string &message) {
    auto &network=m_solver.reactions();const auto names=network.constantNames();
    if(m_reactionParameters->identifiers()!=names){message="reaction_constants identifiers changed; reinitialize the component";return false;}
    for(size_t i=0;i<names.size();++i)if(!std::isfinite((*m_reactionParameters)[i])){message="reaction constant must be finite";return false;}
    for(size_t i=0;i<names.size();++i)if(!network.setConstant(names[i],(*m_reactionParameters)[i],message))return false;
    message.clear();return true;
}
std::vector<IArgument *> FVQualComponent::differentiableArguments() const {
    return m_differentiationReady && !m_reactionParameters->identifiers().empty()?std::vector<IArgument *>{m_reactionParameters.get()}:std::vector<IArgument *>{};
}
std::vector<IOutput *> FVQualComponent::differentiableOutputs() const {return m_differentiationReady?m_differentiableFields:std::vector<IOutput *>{};}
std::vector<IComponentDataItem *> FVQualComponent::differentiableStates() const {
    std::vector<IComponentDataItem *> result;for(auto *f:differentiableOutputs())result.push_back(f);return result;
}
bool FVQualComponent::jvp(DifferentialSet seeds,DifferentialSet results,std::string *message){return differential(false,seeds,results,message);}
bool FVQualComponent::vjp(DifferentialSet seeds,DifferentialSet results,std::string *message){return differential(true,seeds,results,message);}
bool FVQualComponent::differential(bool reverse,DifferentialSet seeds,DifferentialSet results,std::string *message) {
    try {
        if(!m_prepared || !m_differentiationReady || !m_linearized)throw std::invalid_argument("no supported most-recent FVQual step to differentiate");
        const auto &trace=m_solver.reactiveTrace();const size_t n=trace.initial.size(),cells=trace.cells;
        size_t j=0;for(int s:trace.species)for(double v:m_state.speciesBlock(s))
            if(v!=trace.head->step.after[j++])throw std::invalid_argument("state changed since linearization; run update again");
        using R=DifferentialRole;
        const auto index=[&](const DifferentialEntry &e,bool seed) -> int {
            if(e.item==m_reactionParameters.get() && !trace.constants.empty() && e.role==R::Argument && (reverse?!seed:seed))return -1;
            const bool forwardSide=(reverse==seed);
            if(forwardSide && e.role!=R::Output && e.role!=R::StateAfter)throw std::invalid_argument("invalid derivative result/seed role");
            if(!forwardSide && e.role!=R::StateBefore)throw std::invalid_argument("only initial species and kinetic constants are differentiated");
            const auto it=std::find_if(m_differentiableFields.begin(),m_differentiableFields.end(),[&](auto *f){return e.item==static_cast<const IComponentDataItem *>(f);});
            if(it==m_differentiableFields.end())throw std::invalid_argument("item is not a declared differentiable FVQual species");
            return static_cast<int>(it-m_differentiableFields.begin());
        };
        // Validate every role and descriptor before overwriting any result.
        for(bool seed:{true,false}) {
            std::set<std::pair<const IComponentDataItem *,R>> seen;
            for(const auto &e:seed?seeds:results) {
                (void)index(e,seed);
                if(!seen.emplace(e.item,e.role).second)throw std::invalid_argument("duplicate differential entry");
                const auto shape=e.item->shape();const std::vector<int64_t> start(shape.size(),0);
                const auto error=HydroCouple::SDK::detail::validateSelection(shape,e.value,DataKind::Float64,start,shape);
                if(!error.empty())throw std::invalid_argument(error);
            }
        }
        std::vector<double> state(n,0),parameters(trace.constants.size(),0);
        for(const auto &e:seeds) {
            std::vector<double> values;std::string error;
            if(!D::readDifferential(e,e.item->shape(),values,&error))throw std::invalid_argument(error);
            if(!std::all_of(values.begin(),values.end(),[](double v){return std::isfinite(v);}))throw std::invalid_argument("nonfinite derivative seed");
            const int s=index(e,true);auto &target=s<0?parameters:state;const size_t offset=s<0?0:static_cast<size_t>(s)*cells;
            for(size_t i=0;i<values.size();++i)target[offset+i]+=values[i];
        }
        if(reverse) {
            std::vector<double> w(n*(trace.count+1),0);std::copy(state.begin(),state.end(),w.end()-n);
            const auto bar=trace.vjp(w);state=bar.initial;parameters=bar.constants;
        } else {
            const auto tangent=trace.jvp(state,parameters);state.assign(tangent.end()-n,tangent.end());
        }
        for(const auto &e:results) {
            const int s=index(e,false);std::vector<double> values;
            if(s<0)values=parameters;else values.assign(state.begin()+static_cast<size_t>(s)*cells,state.begin()+(static_cast<size_t>(s)+1)*cells);
            std::string error;if(!D::writeDifferential(e,e.item->shape(),values,&error))throw std::invalid_argument(error);
        }
        if(message)message->clear();return true;
    } catch(const std::exception &e){if(message)*message=e.what();return false;}
}
}
