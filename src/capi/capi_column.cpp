// SPDX-License-Identifier: MIT
#include "capi_internal.hpp"
#include "fvqual/capi/fvqual_column.h"
#include "fvqual/hydro/thermalcolumn.h"
#include <set>
using FVQual::Hydro::ThermalColumnTrace;
using FVQual::Hydro::ThermalColumnConfiguration;
struct FVQ_ColumnTrace_ {std::unique_ptr<ThermalColumnTrace> column;std::string metadata;};
namespace {
template<class F> FVQ_ErrorCode checked(F f) {
    return FVQual::CApi::guard([&]() -> FVQ_ErrorCode {
        try {f();return FVQual::CApi::succeed();}
        catch(const std::invalid_argument &e){return FVQual::CApi::fail(FVQ_ERR_INVALID_ARGUMENT,e.what());}
        catch(const nlohmann::json::exception &e){return FVQual::CApi::fail(FVQ_ERR_INVALID_ARGUMENT,e.what());}
    });
}
void array(const void *p,int64_t n,size_t expected) {
    if(!p || n<0 || static_cast<size_t>(n)!=expected)throw std::invalid_argument("thermal column array or dimensions invalid");
}
const ThermalColumnTrace &column(const FVQ_ColumnTrace *t) {
    if(!t)throw std::invalid_argument("thermal column handle is null");return *t->column;
}
void dimensions(const FVQ_ColumnTrace *t,const void *x,int64_t nx,const void *p,int64_t np,const void *y,int64_t ny) {
    const auto &c=column(t);array(x,nx,c.initial().size());array(p,np,c.parameters().size());array(y,ny,c.initial().size()*c.times().size());
}
void known(const nlohmann::json &j,const std::set<std::string> &names) {
    if(!j.is_object())throw std::invalid_argument("thermal configuration must be an object");
    for(const auto &[key,value]:j.items())if(!names.contains(key))throw std::invalid_argument("unsupported thermal setting: "+key);
}
}
extern "C" {
FVQ_ErrorCode fvq_column_create(const char *configuration,FVQ_ColumnTrace **out) {
    return checked([&] {
        if(!configuration || !out)throw std::invalid_argument("thermal configuration and output handle are required");
        auto j=nlohmann::json::parse(configuration);
        known(j,{"layers","depth_m","surface_area_m2","bed_area_fraction","dt_s","duration_s","checkpoint_stride","initial","parameters","weather_constants"});
        ThermalColumnConfiguration c;
        for(const char *key:{"layers","checkpoint_stride"})if(j.contains(key)&&!j[key].is_number_integer())throw std::invalid_argument("column counts must be integers");
        for(const char *key:{"depth_m","surface_area_m2","bed_area_fraction","dt_s","duration_s"})
            if(j.contains(key)&&!j[key].is_number())throw std::invalid_argument("column geometry and time settings must be numeric");
        c.layers=j.value("layers",c.layers);c.depth=j.value("depth_m",c.depth);c.area=j.value("surface_area_m2",c.area);
        c.bedAreaFraction=j.value("bed_area_fraction",c.bedAreaFraction);c.dt=j.value("dt_s",c.dt);
        c.duration=j.value("duration_s",c.duration);c.checkpointStride=j.value("checkpoint_stride",c.checkpointStride);
        if(j.contains("weather_constants")) {
            const auto &w=j["weather_constants"];
            known(w,{"albedo","surface_absorption","water_emissivity","longwave_reflection","pressure_pa","wind_exponent",
                     "wind_height_m","evaporation_height_m","wind_roughness_m","psychrometric","air_density_kg_m3"});
            c.weather.albedo=w.value("albedo",c.weather.albedo);c.weather.surfaceAbsorption=w.value("surface_absorption",c.weather.surfaceAbsorption);
            c.weather.waterEmissivity=w.value("water_emissivity",c.weather.waterEmissivity);c.weather.longwaveReflection=w.value("longwave_reflection",c.weather.longwaveReflection);
            c.weather.pressure=w.value("pressure_pa",c.weather.pressure);c.weather.windExponent=w.value("wind_exponent",c.weather.windExponent);
            c.weather.windHeight=w.value("wind_height_m",c.weather.windHeight);c.weather.evaporationHeight=w.value("evaporation_height_m",c.weather.evaporationHeight);
            c.weather.windRoughness=w.value("wind_roughness_m",c.weather.windRoughness);c.weather.psychrometric=w.value("psychrometric",c.weather.psychrometric);
            c.weather.airDensity=w.value("air_density_kg_m3",c.weather.airDensity);
            for(const auto &[key,value]:w.items())if(!value.is_number() || !std::isfinite(value.get<double>()))throw std::invalid_argument("nonfinite weather constant");
            if(c.weather.albedo<0 || c.weather.albedo>1 || c.weather.surfaceAbsorption<0 || c.weather.surfaceAbsorption>1 ||
               c.weather.waterEmissivity<0 || c.weather.waterEmissivity>1 || c.weather.longwaveReflection<0 || c.weather.longwaveReflection>1 ||
               c.weather.pressure<=0 || c.weather.airDensity<=0 || c.weather.windExponent<1 || c.weather.psychrometric<0 ||
               c.weather.windHeight<=0 || c.weather.evaporationHeight<=0 || c.weather.windRoughness<0)
                throw std::invalid_argument("weather constant outside physical range");
        }
        std::vector<double> x,p=ThermalColumnTrace::defaultParameters(c);
        if(j.contains("initial"))x=j["initial"].get<std::vector<double>>();
        if(j.contains("parameters")) {
            const auto &names=ThermalColumnTrace::parameterNames();known(j["parameters"],std::set<std::string>(names.begin(),names.end()));
            for(size_t i=0;i<names.size();++i)if(j["parameters"].contains(names[i])) {
                if(!j["parameters"][names[i]].is_number())throw std::invalid_argument("column parameters must be numeric");
                p[i]=j["parameters"][names[i]].get<double>();
            }
        }
        auto t=std::make_unique<FVQ_ColumnTrace>();t->column=std::make_unique<ThermalColumnTrace>(c,std::move(x),std::move(p));
        t->metadata=nlohmann::json{{"configuration",j},{"parameters",ThermalColumnTrace::parameterNames()},
            {"state_order",{"temperature_c","shear_u_m_s","shear_v_m_s"}},
            {"coverage","closed thermal column with density/Richardson/heat/wind feedback; depth/area scaling; blocked local-adjoint replay"},
            {"limits","zero bed friction and Coriolis; no sources, horizontal exchange, salinity, kinetics or SWMM gradients; branch-local derivatives"}}.dump();
        *out=t.release();
    });
}
void fvq_column_destroy(FVQ_ColumnTrace *t){delete t;}
FVQ_ErrorCode fvq_column_size(const FVQ_ColumnTrace *t,FVQ_ColumnSize *out) {
    return checked([&]{const auto &c=column(t);if(!out)throw std::invalid_argument("column size output is null");
        *out={static_cast<int64_t>(c.initial().size()),static_cast<int64_t>(c.parameters().size()),static_cast<int64_t>(c.times().size()),
              static_cast<int64_t>(c.checkpointCount()),c.maximumForwardDifference(),c.minimumBuoyancyMargin(),c.derivativeSupported()?1:0};});
}
const char *fvq_column_metadata(const FVQ_ColumnTrace *t) {
    if(!t){FVQual::CApi::fail(FVQ_ERR_NULL_ARGUMENT,"column is null");return nullptr;}FVQual::CApi::succeed();return t->metadata.c_str();
}
FVQ_ErrorCode fvq_column_baseline(const FVQ_ColumnTrace *t,double *x,int64_t nx,double *p,int64_t np,double *times,int64_t nt) {
    return checked([&]{const auto &c=column(t);array(x,nx,c.initial().size());array(p,np,c.parameters().size());array(times,nt,c.times().size());
        std::copy(c.initial().begin(),c.initial().end(),x);std::copy(c.parameters().begin(),c.parameters().end(),p);std::copy(c.times().begin(),c.times().end(),times);});
}
FVQ_ErrorCode fvq_column_evaluate(const FVQ_ColumnTrace *t,const double *x,int64_t nx,const double *p,int64_t np,double *y,int64_t ny) {
    return checked([&]{dimensions(t,x,nx,p,np,y,ny);const auto result=t->column->evaluate({x,static_cast<size_t>(nx)},{p,static_cast<size_t>(np)});std::copy(result.begin(),result.end(),y);});
}
FVQ_ErrorCode fvq_column_jvp(const FVQ_ColumnTrace *t,const double *x,int64_t nx,const double *p,int64_t np,double *y,int64_t ny) {
    return checked([&]{dimensions(t,x,nx,p,np,y,ny);const auto result=t->column->jvp({x,static_cast<size_t>(nx)},{p,static_cast<size_t>(np)});std::copy(result.begin(),result.end(),y);});
}
FVQ_ErrorCode fvq_column_vjp(const FVQ_ColumnTrace *t,const double *y,int64_t ny,double *x,int64_t nx,double *p,int64_t np,FVQ_ColumnAdjointStats *stats) {
    return checked([&]{dimensions(t,x,nx,p,np,y,ny);const auto result=t->column->vjp({y,static_cast<size_t>(ny)});
        std::copy(result.initial.begin(),result.initial.end(),x);std::copy(result.parameters.begin(),result.parameters.end(),p);
        if(stats)*stats={static_cast<int64_t>(result.recomputedSteps),static_cast<int64_t>(result.peakTapeNodes),static_cast<int64_t>(result.peakBlockStates)};});
}
}
