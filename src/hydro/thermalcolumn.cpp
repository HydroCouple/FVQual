// SPDX-License-Identifier: MIT
#include "fvqual/hydro/thermalcolumn.h"
#include "fvqual/hydro/columnphysics.h"
#include "differential.hpp"
#include <numeric>
#include <cmath>
#include <limits>
#include <stdexcept>
namespace FVQual::Hydro {
namespace {
using Differentiation::Scalar;
using Differentiation::Tape;
enum Parameter {Wind,Direction,Air,Dew,Shortwave,Cloud,Neutral,BackgroundViscosity,
    BackgroundDiffusivity,ViscosityScale,DiffusivityScale,DepthScale,Attenuation,
    HeatCapacity,EvaporationA,EvaporationB,Drag,AreaScale,ParameterCount};
void finite(std::span<const double> a,size_t n) {
    if(a.size()!=n || !std::all_of(a.begin(),a.end(),[](double v){return std::isfinite(v);}))
        throw std::invalid_argument("thermal column dimensions differ or values are nonfinite");
}
void validateParameters(std::span<const double> p) {
    finite(p,ParameterCount);
    for(size_t i:{Wind,Shortwave,Neutral,BackgroundViscosity,BackgroundDiffusivity,
                  ViscosityScale,DiffusivityScale,Attenuation,EvaporationA,EvaporationB,Drag})
        if(p[i]<0)throw std::invalid_argument("thermal column coefficients must be nonnegative");
    if(p[DepthScale]<=0 || p[AreaScale]<=0 || p[HeatCapacity]<=0 || p[Cloud]<0 || p[Cloud]>1 ||
       p[Air]<=-100 || p[Air]>=100 || p[Dew]<=-100 || p[Dew]>=100)
        throw std::invalid_argument("thermal column parameter outside its supported physical range");
}
template<class T> std::vector<T> project(std::span<const T> state,int layers) {
    std::vector<T> result(state.begin(),state.end());
    for(int block=1;block<3;++block) {
        T mean=0;for(int k=0;k<layers;++k)mean=mean+state[block*layers+k];
        mean=mean/layers;for(int k=0;k<layers;++k)result[block*layers+k]=result[block*layers+k]-mean;
    }
    return result;
}
template<class T> std::vector<T> solve(std::vector<T> lower,std::vector<T> diagonal,
    std::vector<T> upper,std::vector<T> rhs) {
    std::vector<T> sweep(rhs.size());T beta=diagonal[0];rhs[0]=rhs[0]/beta;
    for(size_t i=1;i<rhs.size();++i) {
        sweep[i]=upper[i-1]/beta;beta=diagonal[i]-lower[i]*sweep[i];
        rhs[i]=(rhs[i]-lower[i]*rhs[i-1])/beta;
    }
    for(size_t i=rhs.size()-1;i>0;--i)rhs[i-1]=rhs[i-1]-sweep[i]*rhs[i];return rhs;
}
struct Regime {double buoyancyMargin=1;bool smooth=true;};
template<class T> std::vector<T> step(const ThermalColumnConfiguration &c,
    std::span<const T> x,std::span<const T> p,double dt,Regime *regime=nullptr) {
    using std::exp;using ColumnPhysics::number;
    if(regime && (number(p[Cloud])>=1-1e-10 ||
       (number(p[Wind])<=1e-10 && c.weather.windExponent==1 && number(p[EvaporationB])>0)))regime->smooth=false;
    const int n=c.layers;const T depth=c.depth*p[DepthScale],area=c.area*p[AreaScale],dz=depth/n;
    std::vector<T> volume(n),contact(n+1),rho(n),km(n+1,T(0)),kz(n+1,T(0));
    for(int k=0;k<=n;++k)contact[k]=area*(1.0-(1.0-c.bedAreaFraction)*k/n);
    for(int k=0;k<n;++k) {
        volume[k]=.5*(contact[k]+contact[k+1])*dz;rho[k]=ColumnPhysics::freshwaterDensity(x[k]);
        if(regime && (number(x[k])<=-2+1e-8 || number(x[k])>=100-1e-8))regime->smooth=false;
    }
    for(int k=1;k<n;++k) {
        const T buoyancy=c.options.gravity/c.options.referenceDensity*(rho[k]-rho[k-1])/dz;
        const T su=(x[n+k-1]-x[n+k])/dz,sv=(x[2*n+k-1]-x[2*n+k])/dz;
        const T shear=su*su+sv*sv;
        if(regime) {
            regime->buoyancyMargin=std::min(regime->buoyancyMargin,std::abs(number(buoyancy)));
            if(std::abs(number(buoyancy))<1e-11 || std::abs(number(shear)-1e-12)<1e-14)regime->smooth=false;
        }
        const auto mixing=ColumnPhysics::richardson(buoyancy,shear,p[Neutral],p[BackgroundViscosity],
            p[BackgroundDiffusivity],p[ViscosityScale],p[DiffusivityScale],c.options.convectiveMixing);
        km[k]=mixing.viscosity;kz[k]=mixing.diffusivity;
    }
    const auto stress=ColumnPhysics::wind(c.weather,p[Wind],p[Direction],p[Drag]);
    std::vector<T> result(x.begin(),x.end());
    for(int field=0;field<3;++field) {
        std::vector<T> lo(n,T(0)),diagonal(n,T(1)),up(n,T(0)),rhs(n);
        const auto &coefficient=field==0?kz:km;
        for(int k=0;k<n;++k) {
            if(k>0)lo[k]=ColumnPhysics::exchange(dt,coefficient[k],contact[k],volume[k],dz);
            if(k+1<n)up[k]=ColumnPhysics::exchange(dt,coefficient[k+1],contact[k+1],volume[k],dz);
            diagonal[k]=1.0-lo[k]-up[k];rhs[k]=x[field*n+k];
            if(field==0)rhs[k]=rhs[k]*volume[k]/volume[k]; // native conservative no-flow transport
        }
        if(field>0)rhs[0]=rhs[0]+dt*(field==1?stress.first:stress.second)/c.options.referenceDensity*contact[0]/volume[0];
        auto solution=solve(std::move(lo),std::move(diagonal),std::move(up),std::move(rhs));
        for(int k=0;k<n;++k)result[field*n+k]=solution[k];
    }
    result=project<T>(result,n);
    const auto heat=ColumnPhysics::surface(c.weather,result[0],p[Wind],p[Air],p[Dew],p[Shortwave],p[Cloud],p[EvaporationA],p[EvaporationB]);
    T irradiance=heat.netShortwave*(1.0-c.weather.surfaceAbsorption);
    for(int k=0;k<n;++k) {
        const T incoming=irradiance;irradiance=irradiance*exp(-p[Attenuation]*dz);
        T watts=incoming*contact[k]-irradiance*contact[k+1];
        if(k==0)watts=watts+heat.surfaceNet()*contact[0];
        result[k]=result[k]+dt*watts/(c.options.referenceDensity*p[HeatCapacity]*volume[k]);
    }
    result[n-1]=result[n-1]+dt*irradiance/(c.options.referenceDensity*p[HeatCapacity]*dz)*(contact[n]*dz/volume[n-1]);
    return result;
}
std::pair<double,double> check(std::span<const double> state) {
    double a=0,b=0;for(size_t i=0;i<state.size();++i){a+=(i+1)*state[i];b+=state[i]*state[i];}return {a,b};
}
std::vector<Scalar> inputs(Tape &tape,std::span<const double> x,std::span<const double> p) {
    std::vector<Scalar> result;for(double v:x)result.push_back(tape.input(v));for(double v:p)result.push_back(tape.input(v));return result;
}
std::vector<double> values(const std::vector<Scalar> &x) {
    std::vector<double> result;for(const auto &v:x)result.push_back(v.value);return result;
}
}

ThermalColumnConfiguration::ThermalColumnConfiguration() {
    options.manningN=0;options.coriolis=0;options.advectionScheme=Transport::AdvectionScheme::Upwind;
    options.backgroundDiffusivity=1e-5;weather.shortwave=250;weather.windSpeed=3;
}
const std::vector<std::string> &ThermalColumnTrace::parameterNames() {
    static const std::vector<std::string> names={"wind_speed_m_s","wind_direction_rad","air_temperature_c","dew_point_c",
        "shortwave_w_m2","cloud_fraction","neutral_viscosity_m2_s","background_viscosity_m2_s","background_diffusivity_m2_s",
        "stable_viscosity_scale","stable_diffusivity_scale","depth_scale","solar_attenuation_m_1","heat_capacity_j_kg_k",
        "evaporation_a","evaporation_b","drag_coefficient","area_scale"};return names;
}
std::vector<double> ThermalColumnTrace::defaultParameters(const ThermalColumnConfiguration &c) {
    const auto &w=c.weather;const auto &o=c.options;
    return {w.windSpeed,w.windDirection,w.airTemperature,w.dewPoint,w.shortwave,w.cloudCover,
        o.neutralViscosity,o.backgroundViscosity,o.backgroundDiffusivity,o.stableViscosityScale,o.stableDiffusivityScale,
        1,.45,o.heatCapacity,w.evaporationA,w.evaporationB,w.dragCoefficient,1};
}
std::vector<double> ThermalColumnTrace::defaultInitial(const ThermalColumnConfiguration &c) {
    std::vector<double> x(c.layers*3,0);for(int k=0;k<c.layers;++k) {
        x[k]=22.-10.*(k+.5)/c.layers;x[c.layers+k]=.03*std::cos(3.141592653589793*(k+.5)/c.layers);
    }return x;
}
std::vector<double> ThermalColumnTrace::evaluate(std::span<const double> initialValues,std::span<const double> p) const {
    finite(initialValues,m_config.layers*3);validateParameters(p);
    const auto initial=project<double>(initialValues,m_config.layers);
    const double depth=m_config.depth*p[DepthScale],area=m_config.area*p[AreaScale],side=std::sqrt(area);
    Mesh::HorizontalMesh horizontal;Mesh::LayeredMesh mesh;State::SpeciesRegistry registry;
    State::StateStore state;InternalModeSolver solver;std::string error;
    const auto require=[&](bool ok){if(!ok)throw std::invalid_argument(error);};
    HydroCouple::SDK::IO::MeshDefinition d;d.meshName="qualified-thermal-column";
    d.nodeX={0,side,side,0};d.nodeY={0,0,side,side};d.nodeZ={-depth,-depth,-depth,-depth};d.faceNodeOffsets={0,4};d.faceNodes={0,1,2,3};
    require(horizontal.build(d,error));Mesh::SigmaSpec sigma;sigma.layerCount=m_config.layers;sigma.minLayerThickness=0;
    require(mesh.build(horizontal,sigma,0,error));
    require(mesh.setHypsography(0,std::vector<double>{-depth,0},std::vector<double>{area*m_config.bedAreaFraction,area},error));
    require(mesh.updateSurface(std::vector<double>{0},error));
    if(registry.addTemperature(10,error)<0)throw std::invalid_argument(error);state.allocate(mesh,registry);
    auto options=m_config.options;options.neutralViscosity=p[Neutral];options.backgroundViscosity=p[BackgroundViscosity];
    options.backgroundDiffusivity=p[BackgroundDiffusivity];options.stableViscosityScale=p[ViscosityScale];options.stableDiffusivityScale=p[DiffusivityScale];
    options.heatCapacity=p[HeatCapacity];require(solver.initialize(mesh,state,registry,options,error));
    Meteorology weather=m_config.weather;weather.windSpeed=p[Wind];weather.windDirection=p[Direction];
    weather.airTemperature=p[Air];weather.dewPoint=p[Dew];weather.shortwave=p[Shortwave];weather.cloudCover=p[Cloud];
    weather.evaporationA=p[EvaporationA];weather.evaporationB=p[EvaporationB];weather.dragCoefficient=p[Drag];
    solver.meteorology()=weather;solver.enableMeteorology(p[Attenuation]);
    const int n=m_config.layers;
    for(int k=0;k<n;++k){state.speciesAt(0,k)=initial[k];state.u[k]=initial[n+k];state.v[k]=initial[2*n+k];}
    std::vector<double> result=initial;
    for(size_t t=1;t<m_times.size();++t) {
        require(solver.step(m_times[t]-m_times[t-1],error));
        std::vector<double> x;for(double v:state.speciesBlock(0))x.push_back(v);
        x.insert(x.end(),state.u.begin(),state.u.end());x.insert(x.end(),state.v.begin(),state.v.end());x=project<double>(x,n);
        result.insert(result.end(),x.begin(),x.end());
    }
    finite(result,m_times.size()*initial.size());return result;
}
ThermalColumnTrace::ThermalColumnTrace(ThermalColumnConfiguration c,std::vector<double> x,std::vector<double> p):m_config(std::move(c)) {
    if(m_config.layers<1 || m_config.layers>64 || !std::isfinite(m_config.dt) || m_config.dt<=0 ||
       !std::isfinite(m_config.options.gravity) || m_config.options.gravity<=0 ||
       !std::isfinite(m_config.options.referenceDensity) || m_config.options.referenceDensity<=0 ||
       !std::isfinite(m_config.options.convectiveMixing) || m_config.options.convectiveMixing<0 ||
       !std::isfinite(m_config.duration) || m_config.duration<=0 || m_config.duration/m_config.dt>10000 ||
       !std::isfinite(m_config.depth) || m_config.depth<=0 || !std::isfinite(m_config.area) || m_config.area<=0 ||
       !(m_config.bedAreaFraction>0 && m_config.bedAreaFraction<=1) || m_config.checkpointStride<0 ||
       m_config.options.manningN!=0 || m_config.options.coriolis!=0 || m_config.options.useTurbulenceClosure ||
       m_config.options.horizontalDiffusivity!=0 || m_config.options.horizontalViscosity!=0 ||
       m_config.options.advectionScheme!=Transport::AdvectionScheme::Upwind)
        throw std::invalid_argument("thermal column requires a finite closed fixed-step column, upwind, Richardson mixing, zero bed friction/Coriolis/horizontal mixing");
    m_initial=x.empty()?defaultInitial(m_config):std::move(x);m_parameters=p.empty()?defaultParameters(m_config):std::move(p);
    m_times={0};const size_t steps=static_cast<size_t>(std::ceil(m_config.duration/m_config.dt));
    for(size_t i=1;i<=steps;++i)m_times.push_back(std::min(m_config.duration,i*m_config.dt));
    const auto actual=evaluate(m_initial,m_parameters);const size_t n=m_initial.size();
    const size_t stride=m_config.checkpointStride?m_config.checkpointStride:static_cast<size_t>(std::ceil(std::sqrt(m_times.size()-1)));
    for(size_t t=0;t<m_times.size();++t) {
        const std::span<const double> row(actual.data()+t*n,n);m_checks.push_back(check(row));
        if(t%stride==0 || t+1==m_times.size())m_checkpoints.emplace(t,std::vector<double>(row.begin(),row.end()));
        if(t+1<m_times.size()) {
            Regime regime;const auto predicted=step<double>(m_config,row,m_parameters,m_times[t+1]-m_times[t],&regime);
            for(size_t i=0;i<n;++i)m_forwardDifference=std::max(m_forwardDifference,std::abs(predicted[i]-actual[(t+1)*n+i]));
            m_supported=m_supported&&regime.smooth;m_buoyancyMargin=std::min(m_buoyancyMargin,regime.buoyancyMargin);
        }
    }
    if(m_forwardDifference>2e-10)throw std::invalid_argument("native thermal-column forward qualification failed");
}
void ThermalColumnTrace::validateReplay(size_t t,std::span<const double> state) const {
    const auto c=check(state),expected=m_checks.at(t);
    if(std::abs(c.first-expected.first)>2e-9*std::max(1.,std::abs(expected.first)) ||
       std::abs(c.second-expected.second)>2e-9*std::max(1.,std::abs(expected.second)))
        throw std::invalid_argument("thermal column replay does not match its native baseline");
}
std::vector<double> ThermalColumnTrace::jvp(std::span<const double> initialDirection,std::span<const double> parameterDirection) const {
    finite(initialDirection,m_initial.size());finite(parameterDirection,m_parameters.size());
    if(!m_supported)throw std::invalid_argument("thermal derivative refused near a density/mixing/clamp branch; use full-run ensembles");
    const size_t n=m_initial.size();Tape first;auto input=inputs(first,m_initial,m_parameters);
    auto state=project<Scalar>(std::span<const Scalar>(input.data(),n),m_config.layers);
    std::vector<double> direction(initialDirection.begin(),initialDirection.end());direction.insert(direction.end(),parameterDirection.begin(),parameterDirection.end());
    auto tangent=first.jvp(direction,state);auto current=values(state);std::vector<double> result=tangent;
    for(size_t t=1;t<m_times.size();++t) {
        Tape tape;input=inputs(tape,current,m_parameters);
        state=step<Scalar>(m_config,{input.data(),n},{input.data()+n,m_parameters.size()},m_times[t]-m_times[t-1]);
        direction=tangent;direction.insert(direction.end(),parameterDirection.begin(),parameterDirection.end());
        tangent=tape.jvp(direction,state);current=values(state);validateReplay(t,current);
        result.insert(result.end(),tangent.begin(),tangent.end());
    }
    finite(result,n*m_times.size());return result;
}
ColumnPullback ThermalColumnTrace::vjp(std::span<const double> cotangents) const {
    const size_t n=m_initial.size();finite(cotangents,n*m_times.size());
    if(!m_supported)throw std::invalid_argument("thermal derivative refused near a density/mixing/clamp branch; use full-run ensembles");
    ColumnPullback result;result.parameters.assign(m_parameters.size(),0);std::vector<double> carry(n,0);
    auto end=m_checkpoints.rbegin();
    while(std::next(end)!=m_checkpoints.rend()) {
        const auto start=std::next(end);std::vector<std::vector<double>> block{start->second};
        for(size_t t=start->first;t<end->first;++t) {
            block.push_back(step<double>(m_config,block.back(),m_parameters,m_times[t+1]-m_times[t]));
            validateReplay(t+1,block.back());++result.recomputedSteps;
        }
        result.peakBlockStates=std::max(result.peakBlockStates,block.size());
        for(size_t t=end->first;t>start->first;--t) {
            for(size_t i=0;i<n;++i)carry[i]+=cotangents[t*n+i];
            Tape tape;auto input=inputs(tape,block[t-start->first-1],m_parameters);
            auto output=step<Scalar>(m_config,{input.data(),n},{input.data()+n,m_parameters.size()},m_times[t]-m_times[t-1]);
            auto bar=tape.vjp(output,carry);result.peakTapeNodes=std::max(result.peakTapeNodes,tape.nodes.size());
            carry.assign(bar.begin(),bar.begin()+n);for(size_t i=0;i<m_parameters.size();++i)result.parameters[i]+=bar[n+i];
        }
        end=start;
    }
    for(size_t i=0;i<n;++i)carry[i]+=cotangents[i];
    Tape tape;auto input=inputs(tape,m_initial,m_parameters);
    const auto normalized=project<Scalar>({input.data(),n},m_config.layers);auto bar=tape.vjp(normalized,carry);
    result.initial.assign(bar.begin(),bar.begin()+n);for(size_t i=0;i<m_parameters.size();++i)result.parameters[i]+=bar[n+i];
    finite(result.initial,n);finite(result.parameters,m_parameters.size());return result;
}
}
