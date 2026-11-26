// SPDX-License-Identifier: MIT
#include "fvqual/kinetics/reactiontrace.h"
#include "fvqual/kinetics/odesolver.h"
#include "../hydro/differential.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
namespace FVQual::Kinetics {
namespace {
using Differentiation::Scalar;
using Differentiation::Tape;
void finite(std::span<const double> a,size_t n) {
    if(a.size()!=n || !std::all_of(a.begin(),a.end(),[](double v){return std::isfinite(v);}))
        throw std::invalid_argument("reaction derivative dimensions differ or values are nonfinite");
}
void boundary(Scalar x,double at,const char *name) {
    if(x.index>=0 && std::abs(x.value-at)<=1e-10*std::max(1.,std::abs(at)))
        throw std::invalid_argument(std::string("reaction derivative refused near ")+name+" branch");
}
Scalar expression(const Expression &e,std::span<const Scalar> vars) {
    using namespace Differentiation;
    std::vector<Scalar> stack;
    for(const auto &t:e.postfix) {
        if(t.type==TokenType::NUMBER){stack.emplace_back(t.value);continue;}
        if(t.type==TokenType::VARIABLE){stack.push_back(vars[t.var_idx]);continue;}
        if(stack.empty())throw std::invalid_argument("invalid derivative expression stack");
        Scalar b=stack.back(),a;
        const bool binary=t.type==TokenType::ADD || t.type==TokenType::SUB || t.type==TokenType::MUL ||
            t.type==TokenType::DIV || t.type==TokenType::POW || t.type==TokenType::FUNC_MIN || t.type==TokenType::FUNC_MAX;
        if(binary){stack.pop_back();if(stack.empty())throw std::invalid_argument("invalid binary expression");a=stack.back();}
        Scalar y;
        switch(t.type) {
        case TokenType::ADD:y=a+b;break;
        case TokenType::SUB:y=a-b;break;
        case TokenType::MUL:y=a*b;break;
        case TokenType::DIV:boundary(b,0,"division");if(b.value==0)throw std::invalid_argument("zero reaction denominator");y=a/b;break;
        case TokenType::POW:y=pow(a,b);break;
        case TokenType::NEG:y=-b;break;
        case TokenType::FUNC_ABS:boundary(b,0,"absolute value");y=abs(b);break;
        case TokenType::FUNC_SQRT:boundary(b,0,"square root");if(b.value<0)throw std::invalid_argument("negative reaction square root");y=pow(b,.5);break;
        case TokenType::FUNC_LOG:boundary(b,0,"logarithm");if(b.value<=0)throw std::invalid_argument("nonpositive reaction logarithm");y=log(b);break;
        case TokenType::FUNC_EXP:y=exp(b);break;
        case TokenType::FUNC_SIN:y=sin(b);break;
        case TokenType::FUNC_COS:y=cos(b);break;
        case TokenType::FUNC_TAN:y=Tape::op(b,{},std::tan(b.value),1/std::pow(std::cos(b.value),2),0);break;
        case TokenType::FUNC_ASIN:y=Tape::op(b,{},std::asin(b.value),1/std::sqrt(1-b.value*b.value),0);break;
        case TokenType::FUNC_ACOS:y=Tape::op(b,{},std::acos(b.value),-1/std::sqrt(1-b.value*b.value),0);break;
        case TokenType::FUNC_ATAN:y=Tape::op(b,{},std::atan(b.value),1/(1+b.value*b.value),0);break;
        case TokenType::FUNC_SGN:boundary(b,0,"sign");y=Scalar(b.value>0?1:b.value<0?-1:0);break;
        case TokenType::FUNC_STEP:boundary(b,0,"step");y=Scalar(b.value>=0?1:0);break;
        case TokenType::FUNC_MIN:case TokenType::FUNC_MAX:
            if(a.index!=b.index)boundary(a-b,0,"min/max");
            y=t.type==TokenType::FUNC_MIN?(a.value<b.value?a:b):(a.value>b.value?a:b);break;
        default:throw std::invalid_argument("reaction expression function is not qualified for differentiation");
        }
        stack.back()=y;
    }
    if(stack.size()!=1)throw std::invalid_argument("invalid derivative expression result");return stack.front();
}
// These operations follow the native adaptive controller, including changing
// accepted step lengths. Only discrete accept/reject decisions are held local.
template<class F> std::vector<Scalar> integrateLocal(std::vector<Scalar> y,double duration,F rhs) {
    Scalar x=0,h=duration;
    const size_t n=y.size();if(!n)return y;
    std::vector<Scalar> d(n),scale(n),trial(n),err(n),stage(n);
    std::array<std::vector<Scalar>,5> k;for(auto &v:k)v.resize(n);
    for(int count=0;count<MAXSTP;++count) {
        d=rhs(y);
        bool scaleBoundary=false;
        for(size_t i=0;i<n;++i) {
            const auto rateStep=d[i]*h;
            if((y[i].index>=0 && y[i].value==0) || (rateStep.index>=0 && rateStep.value==0))scaleBoundary=true;
            scale[i]=Differentiation::abs(y[i])+Differentiation::abs(rateStep)+TINY;
        }
        if(h.index>=0 || x.index>=0)boundary(x+h-duration,0,"adaptive final-step clipping");
        if((x.value+h.value-duration)*(x.value+h.value)>0)h=duration-x;
        Scalar error=0;
        for(int attempt=0;;++attempt) {
            if(attempt>=MAXSTP)throw std::invalid_argument("reaction derivative trial-step limit");
            for(size_t i=0;i<n;++i)stage[i]=y[i]+.2*h*d[i];k[0]=rhs(stage);
            for(size_t i=0;i<n;++i)stage[i]=y[i]+h*(3./40*d[i]+9./40*k[0][i]);k[1]=rhs(stage);
            for(size_t i=0;i<n;++i)stage[i]=y[i]+h*(.3*d[i]-.9*k[0][i]+1.2*k[1][i]);k[2]=rhs(stage);
            for(size_t i=0;i<n;++i)stage[i]=y[i]+h*(-11./54*d[i]+2.5*k[0][i]-70./27*k[1][i]+35./27*k[2][i]);k[3]=rhs(stage);
            for(size_t i=0;i<n;++i)stage[i]=y[i]+h*(1631./55296*d[i]+175./512*k[0][i]+575./13824*k[1][i]+44275./110592*k[2][i]+253./4096*k[3][i]);k[4]=rhs(stage);
            error=0;
            bool tiedMaximum=false;
            for(size_t i=0;i<n;++i) {
                trial[i]=y[i]+h*(37./378*d[i]+250./621*k[1][i]+125./594*k[2][i]+512./1771*k[4][i]);
                err[i]=h*((37./378-2825./27648)*d[i]+(250./621-18575./48384)*k[1][i]+(125./594-13525./55296)*k[2][i]-277./14336*k[3][i]+(512./1771-.25)*k[4][i]);
                const auto candidate=Differentiation::abs(err[i]/scale[i]);
                const double margin=1e-10*std::max(std::abs(error.value),std::abs(candidate.value));
                if(candidate.value>error.value+margin){error=candidate;tiedMaximum=false;}
                else if(candidate.value>0 && std::abs(candidate.value-error.value)<=margin && candidate.index!=error.index) {
                    tiedMaximum=true;if(candidate.value>error.value)error=candidate;
                }
            }
            error=error/1e-8;boundary(error,1,"adaptive acceptance");
            if((error.value>1 || x.value+h.value<duration) && error.value>ERRCON && (scaleBoundary || tiedMaximum))
                throw std::invalid_argument("reaction derivative refused near an adaptive scaling/max-error branch");
            if(error.value<=1)break;
            const auto next=SAFETY*h*Differentiation::pow(error,PSHRNK);
            boundary(next-.1*h,0,"adaptive shrink limit");h=next.value>.1*h.value?next:.1*h;
            if(x.value+h.value==x.value)throw std::invalid_argument("reaction derivative step underflow");
        }
        x=x+h;y=trial;
        if((x.value-duration)*duration>=0)return y;
        boundary(error,ERRCON,"adaptive growth");
        h=error.value>ERRCON?SAFETY*h*Differentiation::pow(error,PGROW):5.*h;
        if(!(h.value>0))throw std::invalid_argument("reaction derivative invalid step length");
    }
    throw std::invalid_argument("reaction derivative integration limit");
}
}
struct ReactionStep::Data {
    struct Law {Expression expression;int slot;RateScope scope;bool areal;};
    ReactionNetwork network;
    State::SpeciesRegistry registry;
    State::StateStore state;
    std::vector<int> reacting;
    std::vector<Law> laws;
    std::vector<double> constants;
    std::array<std::vector<double>,5> fields;
    std::vector<unsigned char> positions;
    double dt=0;
    std::vector<Scalar> map(size_t cell,std::span<const Scalar> input) const {
        const size_t species=registry.count(),nc=constants.size();
        std::vector<Scalar> vars(input.begin(),input.begin()+species);
        // Native binding order is species, three derived fields, constants.
        for(size_t d=0;d<3;++d)vars.push_back(input[species+nc+d]);
        for(size_t c=0;c<nc;++c)vars.push_back(input[species+c]);
        std::vector<Scalar> y;for(int s:reacting)y.push_back(vars[s]);
        const unsigned char flags=positions[cell];
        const auto derivatives=[&](const std::vector<Scalar> &current) {
            for(size_t r=0;r<reacting.size();++r)vars[reacting[r]]=current[r];
            std::vector<Scalar> rates(reacting.size(),Scalar(0));
            for(const auto &law:laws) {
                if((law.scope==RateScope::Bed || law.areal) && !(flags&CellAtBed))continue;
                if(law.scope==RateScope::Surface && !(flags&CellAtSurface))continue;
                auto rate=expression(law.expression,vars);
                if(law.scope!=RateScope::WaterColumn && !law.areal) {
                    const size_t d=law.scope==RateScope::Surface?3:4;
                    rate=rate*(fields[d].empty()?1./input[species+nc+2]:input[species+nc+d]);
                }
                rates[law.slot]=rates[law.slot]+rate;
            }
            return rates;
        };
        y=integrateLocal(std::move(y),dt,derivatives);
        std::vector<Scalar> result(input.begin(),input.begin()+species);
        for(size_t r=0;r<reacting.size();++r)result[reacting[r]]=y[r];return result;
    }
    std::vector<Scalar> inputs(Tape &tape,size_t cell) const {
        std::vector<Scalar> x;for(int s=0;s<registry.count();++s)x.push_back(tape.input(state.speciesAt(s,cell)));
        for(double c:constants)x.push_back(tape.input(c));
        for(const auto &field:fields)x.push_back(tape.input(field.empty()?0:field[cell]));return x;
    }
};
ReactionStep ReactionStep::capture(const ReactionNetwork &network,const State::SpeciesRegistry &registry,
    const State::StateStore &state,double dt,const DerivedFields &fields,std::span<const unsigned char> position) {
    if(!(dt>0) || !std::isfinite(dt) || (!network.empty()&&!network.compiled()))
        throw std::invalid_argument("reaction trace requires compiled laws and a finite positive interval");
    auto d=std::make_shared<Data>();d->network=network;d->registry=registry;d->state=state;d->dt=dt;
    d->reacting=network.m_reacting;d->constants=network.m_constantValues;
    for(const auto &law:network.m_rates)d->laws.push_back({law.expression,law.slot,law.scope,law.arealSpecies});
    const std::span<const double> source[]={fields.light,fields.depth,fields.thickness,fields.surfaceAreaPerVolume,fields.bedAreaPerVolume};
    const size_t cells=state.cellCount();
    for(size_t i=0;i<5;++i) {
        if(!source[i].empty()){finite(source[i],cells);d->fields[i].assign(source[i].begin(),source[i].end());}
        if(i<3 && network.m_derivedUsed[i] && source[i].empty())throw std::invalid_argument("reaction trace missing a required derived field");
        if(i>=3 && std::any_of(source[i].begin(),source[i].end(),[](double v){return v<0;}))throw std::invalid_argument("reaction area/volume factors must be nonnegative");
    }
    if(!position.empty() && position.size()!=cells)throw std::invalid_argument("reaction trace position size mismatch");
    d->positions.assign(cells,CellInterior);if(!position.empty())d->positions.assign(position.begin(),position.end());
    if(network.usesCellPosition() && position.empty())throw std::invalid_argument("reaction trace missing position flags");
    if(network.usesCellPosition() && (fields.thickness.empty() || std::any_of(fields.thickness.begin(),fields.thickness.end(),[](double v){return v<=0;})))
        throw std::invalid_argument("scoped reaction trace requires positive thickness");
    ReactionStep result;result.m_data=std::move(d);return result;
}
const std::vector<double> &ReactionStep::before() const {if(!m_data)throw std::invalid_argument("empty reaction trace");return m_data->state.species;}
const std::vector<double> &ReactionStep::constants() const {if(!m_data)throw std::invalid_argument("empty reaction trace");return m_data->constants;}
std::vector<double> ReactionStep::evaluate(std::span<const double> values,std::span<const double> constants) const {
    finite(values,before().size());finite(constants,m_data->constants.size());
    auto network=m_data->network;network.m_constantValues.assign(constants.begin(),constants.end());
    auto state=m_data->state;state.species.assign(values.begin(),values.end());
    DerivedFields f{m_data->fields[0],m_data->fields[1],m_data->fields[2],m_data->fields[3],m_data->fields[4]};
    std::string error;
    if(!network.integrate(m_data->registry,state,state.cellCount(),m_data->dt,error,f,m_data->positions) || network.diagnostics().unconvergedCells)
        throw std::invalid_argument("native reaction replay failed: "+error);
    finite(state.species,values.size());return state.species;
}
std::vector<double> ReactionStep::jvp(std::span<const double> direction,std::span<const double> constants,
    std::span<const double> geometry) const {
    finite(direction,before().size());finite(constants,m_data->constants.size());const size_t cells=m_data->state.cellCount();
    if(!geometry.empty())finite(geometry,5*cells);
    const auto native=evaluate(before(),m_data->constants);std::vector<double> result(direction.size());
    for(size_t cell=0;cell<cells;++cell) {
        Tape tape;const auto x=m_data->inputs(tape,cell),y=m_data->map(cell,x);std::vector<double> v;
        for(int s=0;s<m_data->registry.count();++s)v.push_back(direction[s*cells+cell]);
        v.insert(v.end(),constants.begin(),constants.end());
        for(size_t d=0;d<5;++d)v.push_back(geometry.empty()?0:geometry[d*cells+cell]);
        const auto dy=tape.jvp(v,y);
        for(size_t s=0;s<dy.size();++s) {
            if(std::abs(y[s].value-native[s*cells+cell])>2e-11*std::max(1.,std::abs(native[s*cells+cell])))
                throw std::invalid_argument("reaction derivative operator differs from native solve");
            result[s*cells+cell]=dy[s];
        }
    }
    finite(result,direction.size());return result;
}
ReactionPullback ReactionStep::vjp(std::span<const double> cotangents) const {
    finite(cotangents,before().size());const size_t cells=m_data->state.cellCount(),ns=m_data->registry.count(),nc=m_data->constants.size();
    const auto native=evaluate(before(),m_data->constants);ReactionPullback result;
    result.species.resize(cotangents.size());result.constants.assign(nc,0);result.geometry.resize(5*cells);
    for(size_t cell=0;cell<cells;++cell) {
        Tape tape;const auto x=m_data->inputs(tape,cell),y=m_data->map(cell,x);std::vector<double> w;
        for(size_t s=0;s<ns;++s) {
            if(std::abs(y[s].value-native[s*cells+cell])>2e-11*std::max(1.,std::abs(native[s*cells+cell])))
                throw std::invalid_argument("reaction derivative operator differs from native solve");
            w.push_back(cotangents[s*cells+cell]);
        }
        const auto bar=tape.vjp(y,w);result.peakTapeNodes=std::max(result.peakTapeNodes,tape.nodes.size());
        for(size_t s=0;s<ns;++s)result.species[s*cells+cell]=bar[s];
        for(size_t c=0;c<nc;++c)result.constants[c]+=bar[ns+c];
        for(size_t d=0;d<5;++d)result.geometry[d*cells+cell]=bar[ns+nc+d];
    }
    finite(result.species,cotangents.size());finite(result.constants,nc);finite(result.geometry,5*cells);return result;
}
}
