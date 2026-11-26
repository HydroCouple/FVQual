// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
#include "fvqual/transport/scalartrace.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace FVQual::Transport {
namespace {
void finite(std::span<const double> x, size_t n) {
    if (x.size()!=n || !std::all_of(x.begin(),x.end(),[](double v){return std::isfinite(v);}))
        throw std::invalid_argument("scalar trace dimension mismatch or nonfinite value");
}
std::vector<double> solve(const ScalarStep &s, std::span<const double> rhs,
                          bool transpose=false, double scale=1) {
    std::vector<double> y(rhs.begin(),rhs.end()), sweep(y.size());
    const size_t layers=static_cast<size_t>(s.layers);
    for (size_t base=0;base<y.size();base+=layers) {
        double beta=1+scale*(s.diagonal[base]-1);
        y[base]/=beta;
        for (size_t k=1;k<layers;++k) {
            const size_t i=base+k;
            const double upper=scale*(transpose?s.lower[i]:s.upper[i-1]);
            const double lower=scale*(transpose?s.upper[i-1]:s.lower[i]);
            sweep[i]=upper/beta;
            beta=1+scale*(s.diagonal[i]-1)-lower*sweep[i];
            y[i]=(y[i]-lower*y[i-1])/beta;
        }
        for (size_t k=layers-1;k>0;--k) y[base+k-1]-=sweep[base+k]*y[base+k];
    }
    finite(y,rhs.size()); return y;
}
std::vector<double> multiply(const std::vector<ScalarEntry> &entries,
                            std::span<const double> x, size_t n) {
    std::vector<double> y(n,0);
    for (const auto &a:entries) y[a.row]+=a.value*x[a.column];
    return y;
}
std::vector<double> diffusionAction(const ScalarStep &s) {
    std::vector<double> y(s.after.size());
    for (size_t i=0;i<y.size();++i) {
        y[i]=(s.diagonal[i]-1)*s.after[i];
        if (i%static_cast<size_t>(s.layers)) y[i]+=s.lower[i]*s.after[i-1];
        if ((i+1)%static_cast<size_t>(s.layers)) y[i]+=s.upper[i]*s.after[i+1];
    }
    return y;
}
}

ScalarStep ScalarStep::capture(const Mesh::LayeredMesh &mesh,
    std::span<const double> oldVolume, std::span<const double> edgeFlux,
    std::span<const double> verticalFlux, const Forcing::SourceSet &sources,
    int species, std::span<const double> values, double time, double dt) {
    ScalarStep s; s.time=time; s.dt=dt; s.layers=mesh.layerCount();
    const size_t n=static_cast<size_t>(mesh.cellCount());
    s.before.assign(values.begin(),values.end()); s.after.resize(n);
    s.lower.resize(n);s.diagonal.resize(n);s.upper.resize(n);
    finite(values,n); finite(oldVolume,n);
    const auto &vol=mesh.cellVolume();
    for (size_t i=0;i<n;++i) {
        if (!(vol[i]>0) || !std::isfinite(vol[i]))
            throw std::invalid_argument("scalar differentiation requires positive finite cell volumes");
        s.advection.push_back({i,i,oldVolume[i]/vol[i]});
    }
    // Add the two conservative contributions of each donor flux. Boundary
    // edges have zero advective scalar flux in EulerianTransport.
    const auto connection=[&](size_t left,size_t right,double flow) {
        const size_t donor=flow>=0?left:right;
        s.advection.push_back({left,donor,-dt*flow/vol[left]});
        s.advection.push_back({right,donor,dt*flow/vol[right]});
    };
    const auto &h=mesh.horizontal();
    for (size_t e=0;e<static_cast<size_t>(h.edgeCount());++e) {
        if (h.edgeRightFace()[e]==Mesh::kNoFace) continue;
        for (int k=0;k<s.layers;++k)
            connection(mesh.cell(h.edgeLeftFace()[e],k),mesh.cell(h.edgeRightFace()[e],k),
                       edgeFlux[e*s.layers+k]);
    }
    for (int64_t c=0;c<mesh.columnCount();++c)
        for (int k=1;k<s.layers;++k)
            connection(mesh.cell(c,k),mesh.cell(c,k-1),
                       verticalFlux[mesh.interfaceSlot(c,k)]*h.faceArea()[c]);
    std::vector<double> withdrawal;
    std::vector<Forcing::ScalarSourceEntry> incoming;
    sources.scalarLinearization(mesh,withdrawal,incoming);
    for (size_t i=0;i<n;++i) s.advection.push_back({i,i,dt*withdrawal[i]/vol[i]});
    for (int i=0;i<sources.count();++i) {
        const auto &source=sources.source(i);
        s.sourceIds.push_back(source.id.empty()?"source:"+std::to_string(i):source.id);
        s.sourceValues.push_back(static_cast<size_t>(species)<source.concentration.size()?source.concentration[species]:0);
    }
    for (const auto &a:incoming) s.inflow.push_back({a.cell,a.source,dt*a.rate/vol[a.cell]});
    finite(s.sourceValues,s.sourceValues.size());
    return s;
}

std::vector<const ScalarStep *> ScalarTrajectory::steps() const {
    std::vector<const ScalarStep *> result;
    for (auto node=head;node;node=node->previous) result.push_back(&node->step);
    if (result.size()!=count) throw std::invalid_argument("inconsistent scalar history");
    std::reverse(result.begin(),result.end());return result;
}
std::vector<double> ScalarStep::apply(std::span<const double> values,std::span<const double> sources) const {
    finite(values,before.size());finite(sources,sourceValues.size());
    auto rhs=multiply(advection,values,values.size());const auto load=multiply(inflow,sources,values.size());
    for(size_t i=0;i<rhs.size();++i)rhs[i]+=load[i];return solve(*this,rhs);
}
std::pair<std::vector<double>,std::vector<double>> ScalarStep::transpose(std::span<const double> cotangents) const {
    finite(cotangents,before.size());const auto lambda=solve(*this,cotangents,true);
    std::vector<double> state(before.size(),0),sources(sourceValues.size(),0);
    for(const auto &a:advection)state[a.column]+=a.value*lambda[a.row];
    for(const auto &a:inflow)sources[a.column]+=a.value*lambda[a.row];return {state,sources};
}
std::vector<double> ScalarTrajectory::sourceValues() const {
    std::vector<double> result;
    for (const auto *s:steps()) result.insert(result.end(),s->sourceValues.begin(),s->sourceValues.end());
    return result;
}
std::vector<double> ScalarTrajectory::evaluate(std::span<const double> initialValues,
    std::span<const double> sources,double scale) const {
    finite(initialValues,initial.size());finite(sources,sourceValues().size());
    if (!(scale>0) || !std::isfinite(scale)) throw std::invalid_argument("diffusion multiplier must be positive and finite");
    std::vector<double> state(initialValues.begin(),initialValues.end()), result=state;
    size_t offset=0;
    for (const auto *s:steps()) {
        auto rhs=multiply(s->advection,state,state.size());
        auto load=multiply(s->inflow,sources.subspan(offset,s->sourceValues.size()),state.size());
        offset+=s->sourceValues.size();
        for (size_t i=0;i<rhs.size();++i) rhs[i]+=load[i];
        state=solve(*s,rhs,false,scale);
        result.insert(result.end(),state.begin(),state.end());
    }
    return result;
}
std::vector<double> ScalarTrajectory::jvp(std::span<const double> initialDirection,
    std::span<const double> sourceDirection,double diffusionDirection) const {
    finite(initialDirection,initial.size());finite(sourceDirection,sourceValues().size());
    if (!std::isfinite(diffusionDirection)) throw std::invalid_argument("nonfinite diffusion direction");
    std::vector<double> state(initialDirection.begin(),initialDirection.end()), result=state;
    size_t offset=0;
    for (const auto *s:steps()) {
        auto rhs=multiply(s->advection,state,state.size());
        auto load=multiply(s->inflow,sourceDirection.subspan(offset,s->sourceValues.size()),state.size());
        const auto change=diffusionAction(*s);
        offset+=s->sourceValues.size();
        for (size_t i=0;i<rhs.size();++i) rhs[i]+=load[i]-diffusionDirection*change[i];
        state=solve(*s,rhs);
        result.insert(result.end(),state.begin(),state.end());
    }
    return result;
}
ScalarPullback ScalarTrajectory::vjp(std::span<const double> cotangents) const {
    const size_t n=initial.size();finite(cotangents,n*(count+1));
    ScalarPullback result;result.sources.assign(sourceValues().size(),0);
    std::vector<double> carry(n,0);
    const auto history=steps();size_t offset=result.sources.size();
    for (size_t k=count;k>0;--k) {
        const auto &s=*history[k-1];offset-=s.sourceValues.size();
        for (size_t i=0;i<n;++i) carry[i]+=cotangents[k*n+i];
        const auto lambda=solve(s,carry,true);
        const auto change=diffusionAction(s);
        for (size_t i=0;i<n;++i) result.diffusionMultiplier-=lambda[i]*change[i];
        for (const auto &a:s.inflow) result.sources[offset+a.column]+=a.value*lambda[a.row];
        std::fill(carry.begin(),carry.end(),0);
        for (const auto &a:s.advection) carry[a.column]+=a.value*lambda[a.row];
    }
    for (size_t i=0;i<n;++i) carry[i]+=cotangents[i];
    finite(carry,n);finite(result.sources,result.sources.size());
    if (!std::isfinite(result.diffusionMultiplier)) throw std::invalid_argument("scalar pullback overflow");
    result.initial=std::move(carry);return result;
}
size_t ScalarTrajectory::bytes() const {
    size_t result=initial.size()*sizeof(double);
    for (const auto *s:steps()) {
        result+=sizeof(ScalarTraceNode)+(s->advection.size()+s->inflow.size())*sizeof(ScalarEntry);
        result+=(s->sourceValues.size()+s->lower.size()+s->diagonal.size()+s->upper.size()+
                 s->before.size()+s->after.size())*sizeof(double);
        for (const auto &id:s->sourceIds) result+=sizeof(std::string)+id.size();
    }
    return result;
}
}
