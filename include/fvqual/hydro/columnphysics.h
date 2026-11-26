// SPDX-License-Identifier: MIT
// Shared production physics for double evaluation and discrete differentiation.
#ifndef FVQUAL_HYDRO_COLUMNPHYSICS_H
#define FVQUAL_HYDRO_COLUMNPHYSICS_H
#include "fvqual/hydro/heatbudget.h"
#include <cmath>
#include <algorithm>
namespace FVQual::Hydro::ColumnPhysics {
template<class T> double number(const T &x) { return static_cast<double>(x); }
template<class T> T freshwaterDensity(T t) {
    if(number(t)<-2) t=T(-2);if(number(t)>100)t=T(100);
    return 999.842594 + 6.793952e-2*t - 9.095290e-3*t*t +
        1.001685e-4*t*t*t - 1.120083e-6*t*t*t*t + 6.536332e-9*t*t*t*t*t;
}
template<class T> T vapour(const T &t) {
    using std::exp;return 610.78*exp(17.27*t/(t+237.3));
}
template<class T> struct Mixing { T viscosity,diffusivity,richardson; };
template<class T> Mixing<T> richardson(const T &buoyancy,const T &shear,
    const T &neutral,const T &backgroundViscosity,const T &backgroundDiffusivity,
    const T &viscosityScale,const T &diffusivityScale,double convective) {
    if(number(buoyancy)<0)return {T(convective),T(convective),T(0)};
    const T ri=buoyancy/(number(shear)>1e-12?shear:T(1e-12));
    const T damping=1.0+5.0*ri;
    const T viscosity=neutral/(damping*damping)+backgroundViscosity;
    return {viscosity*viscosityScale,(viscosity/damping+backgroundDiffusivity)*diffusivityScale,ri};
}
template<class T> T exchange(double dt,const T &coefficient,const T &area,
    const T &volume,const T &spacing) {
    return -dt*coefficient*area/(volume*(number(spacing)>1e-12?spacing:T(1e-12)));
}
template<class T> struct HeatTerms {
    T netShortwave,surfaceSolar,penetratingSolar,longwaveDown,longwaveUp,latent,sensible;
    T surfaceNet() const {return surfaceSolar+longwaveDown+longwaveUp+latent+sensible;}
};
template<class T> HeatTerms<T> surface(const Meteorology &w,const T &temperature,
    const T &windSpeed,const T &airTemperature,const T &dewPoint,const T &shortwave,
    T cloud,const T &evaporationA,const T &evaporationB) {
    using std::log;using std::pow;
    HeatTerms<T> t;
    t.netShortwave=(1.0-w.albedo)*shortwave;
    t.surfaceSolar=w.surfaceAbsorption*t.netShortwave;
    t.penetratingSolar=t.netShortwave-t.surfaceSolar;
    const T airK=airTemperature+273.15;
    const T emissivity=.937e-5*airK*airK;
    if(number(cloud)<0)cloud=T(0);if(number(cloud)>1)cloud=T(1);
    const T cover=1.0+.17*cloud*cloud;
    t.longwaveDown=(1.0-w.longwaveReflection)*emissivity*5.670374419e-8*airK*airK*airK*airK*cover;
    const T waterK=temperature+273.15;
    t.longwaveUp=-w.waterEmissivity*5.670374419e-8*waterK*waterK*waterK*waterK;
    T speed=number(windSpeed)>0?windSpeed:T(0);
    if(w.evaporationHeight!=w.windHeight && w.windRoughness>0 &&
       w.windHeight>w.windRoughness && w.evaporationHeight>w.windRoughness)
        speed=speed*log(w.evaporationHeight/w.windRoughness)/log(w.windHeight/w.windRoughness);
    const T wind=evaporationA+evaporationB*(w.windExponent==1?speed:pow(speed,w.windExponent));
    t.latent=-wind*(vapour(temperature)-vapour(dewPoint));
    const double gamma=w.psychrometric*w.pressure/101325.0;
    t.sensible=-wind*gamma*(temperature-airTemperature);return t;
}
template<class T> std::pair<T,T> wind(const Meteorology &w,const T &speed,
    const T &direction,const T &drag) {
    using std::cos;using std::sin;
    const T positive=number(speed)>0?speed:T(0);
    const T magnitude=w.airDensity*drag*positive*positive;
    return {magnitude*cos(direction),magnitude*sin(direction)};
}
}
#endif
