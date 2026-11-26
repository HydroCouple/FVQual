// SPDX-License-Identifier: MIT
#include "fvqual/transport/reactivetrace.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace FVQual::Transport {
namespace {
void finite(std::span<const double> a,size_t n) {
    if(a.size()!=n || !std::all_of(a.begin(),a.end(),[](double v){return std::isfinite(v);}))
        throw std::invalid_argument("reactive trace dimensions differ or values are nonfinite");
}
std::vector<double> selected(std::span<const double> full,const ReactiveTrajectory &t) {
    std::vector<double> result;for(int s:t.species)result.insert(result.end(),full.begin()+s*t.cells,full.begin()+(s+1)*t.cells);return result;
}
void scatter(std::span<const double> x,std::vector<double> &full,const ReactiveTrajectory &t) {
    for(size_t i=0;i<t.species.size();++i)std::copy_n(x.begin()+i*t.cells,t.cells,full.begin()+t.species[i]*t.cells);
}
}
std::vector<const ReactiveStep *> ReactiveTrajectory::steps() const {
    std::vector<const ReactiveStep *> result;for(auto n=head;n;n=n->previous)result.push_back(&n->step);
    if(result.size()!=count)throw std::invalid_argument("inconsistent reactive history");std::reverse(result.begin(),result.end());return result;
}
std::vector<double> ReactiveTrajectory::evaluate(std::span<const double> x,std::span<const double> p) const {
    finite(x,initial.size());finite(p,constants.size());std::vector<double> current(x.begin(),x.end()),result=current;
    for(const auto *step:steps()) {
        std::vector<double> transported;
        for(size_t s=0;s<species.size();++s) {
            auto v=step->transport[s].apply(std::span<const double>(current.data()+s*cells,cells),step->transport[s].sourceValues);
            transported.insert(transported.end(),v.begin(),v.end());
        }
        auto full=step->reaction.before();scatter(transported,full,*this);
        current=selected(step->reaction.evaluate(full,p),*this);result.insert(result.end(),current.begin(),current.end());
    }
    finite(result,initial.size()*(count+1));return result;
}
std::vector<double> ReactiveTrajectory::jvp(std::span<const double> x,std::span<const double> p) const {
    finite(x,initial.size());finite(p,constants.size());std::vector<double> tangent(x.begin(),x.end()),result=tangent;
    for(const auto *step:steps()) {
        std::vector<double> transported;
        for(size_t s=0;s<species.size();++s) {
            const std::vector<double> zero(step->transport[s].sourceValues.size(),0);
            auto v=step->transport[s].apply(std::span<const double>(tangent.data()+s*cells,cells),zero);
            transported.insert(transported.end(),v.begin(),v.end());
        }
        std::vector<double> full(step->reaction.before().size(),0);scatter(transported,full,*this);
        tangent=selected(step->reaction.jvp(full,p),*this);result.insert(result.end(),tangent.begin(),tangent.end());
    }
    finite(result,initial.size()*(count+1));return result;
}
ReactivePullback ReactiveTrajectory::vjp(std::span<const double> w) const {
    const size_t n=initial.size();finite(w,n*(count+1));ReactivePullback result;result.constants.assign(constants.size(),0);
    std::vector<double> carry(n,0);const auto history=steps();
    for(size_t k=count;k>0;--k) {
        const auto &step=*history[k-1];for(size_t i=0;i<n;++i)carry[i]+=w[k*n+i];
        std::vector<double> full(step.reaction.before().size(),0);scatter(carry,full,*this);
        const auto bar=step.reaction.vjp(full);carry=selected(bar.species,*this);
        result.peakTapeNodes=std::max(result.peakTapeNodes,bar.peakTapeNodes);
        for(size_t i=0;i<result.constants.size();++i)result.constants[i]+=bar.constants[i];
        for(size_t s=0;s<species.size();++s) {
            const auto b=step.transport[s].transpose(std::span<const double>(carry.data()+s*cells,cells));
            std::copy(b.first.begin(),b.first.end(),carry.begin()+s*cells);
        }
    }
    for(size_t i=0;i<n;++i)carry[i]+=w[i];result.initial=std::move(carry);
    finite(result.initial,n);finite(result.constants,constants.size());return result;
}
}
