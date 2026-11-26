// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
#include "fvqual/hydro/internalmode.h"
#include <nlohmann/json.hpp>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
using namespace FVQual;
using nlohmann::json;
namespace {
void require(bool ok,const std::string &message) {if(!ok)throw std::runtime_error(message);}
bool transportedTemperature=false;
struct Result {std::vector<double> values;Transport::ScalarTrajectory trace;};
Result run(int layers,double dt,double scale=1,bool record=true) {
    Mesh::HorizontalMesh horizontal;Mesh::LayeredMesh mesh;State::SpeciesRegistry registry;
    State::StateStore state;Hydro::InternalModeSolver solver;std::string error;
    HydroCouple::SDK::IO::MeshDefinition d;d.meshName="diffusion-column";
    d.nodeX={0,10,10,0};d.nodeY={0,0,10,10};d.nodeZ={-4,-4,-4,-4};d.faceNodeOffsets={0,4};d.faceNodes={0,1,2,3};
    require(horizontal.build(d,error),error);Mesh::SigmaSpec sigma;sigma.layerCount=layers;sigma.minLayerThickness=0;
    // Homogeneous density without a transported temperature fixes Kz exactly.
    // A nominally uniform temperature may acquire roundoff stratification and
    // trigger convective mixing; that would be a different executed regime.
    require(mesh.build(horizontal,sigma,0,error),error);
    if(transportedTemperature)require(registry.addTemperature(10,error)>=0,error);
    const int s=registry.add({"tracer",State::SpeciesKind::Pollutant,"kg/m3","",1},error);
    state.allocate(mesh,registry);
    const double pi=std::acos(-1.0),dz=4./layers;
    // Independent exact layer averages of C(z,0)=1+0.2*cos(pi*depth/4).
    for(int k=0;k<layers;++k)state.speciesAt(s,k)=1+.2*4/(pi*dz)*(std::sin(pi*(k+1)/layers)-std::sin(pi*k/layers));
    Hydro::InternalModeOptions o;o.advectionScheme=Transport::AdvectionScheme::Upwind;o.manningN=0;
    o.neutralViscosity=0;o.backgroundViscosity=0;o.backgroundDiffusivity=.01;o.stableDiffusivityScale=scale;
    o.maxTimeStep=dt;o.minTimeStep=dt;require(solver.initialize(mesh,state,registry,o,error),error);
    if(record)require(solver.beginScalarTrace(s,1000,error),error);
    require(solver.advance(120,error),error);
    Result r;r.values.assign(state.speciesBlock(s).begin(),state.speciesBlock(s).end());
    if(record)r.trace=solver.scalarTrace();return r;
}
double sample(const std::vector<double> &values,double depth) {
    const double cell=depth*values.size()/4-.5;
    const size_t i=static_cast<size_t>(std::floor(cell));const double fraction=cell-i;
    return values.at(i)*(1-fraction)+values.at(i+1)*fraction;
}
}
int main(int argc,char **argv) {
    try {
        if(argc!=2 && argc!=3)throw std::invalid_argument("Usage: fvqual-scalar-derivative-review output.json [--temperature-feedback-diagnostic]");
        if(argc==3) {
            require(std::string(argv[2])=="--temperature-feedback-diagnostic","Unknown diagnostic option");
            transportedTemperature=true;
            double difference=0,kzDeparture=0;json cases=json::array();
            for(const auto &[layers,dt]:std::vector<std::pair<int,double>>{{16,10},{16,5},{16,2.5},{4,2.5},{8,2.5},{32,2.5}}) {
                auto r=run(layers,dt);const auto tangent=r.trace.jvp(std::vector<double>(layers,0),{},1);
                auto plus=run(layers,dt,1.0001),minus=run(layers,dt,.9999);
                double fd=0,departure=0;const double dz=4./layers;
                for(int k=0;k<layers;++k)fd=std::max(fd,std::abs((plus.values[k]-minus.values[k])/.0002-tangent[tangent.size()-layers+k]));
                for(const auto &[trajectory,scale]:std::vector<std::pair<const Transport::ScalarTrajectory *,double>>{{&r.trace,1},{&plus.trace,1.0001},{&minus.trace,.9999}})
                    for(const auto *s:trajectory->steps())for(int k=1;k<layers;++k)
                        departure=std::max(departure,std::abs(-s->lower[k]*dz*dz/s->dt-.01*scale));
                difference=std::max(difference,fd);kzDeparture=std::max(kzDeparture,departure);
                cases.push_back({{"layers",layers},{"dt_s",dt},{"max_gradient_difference",fd},{"max_kz_departure",departure}});
            }
            json diagnostic={{"status","diagnostic only; thermal mixing feedback is outside the qualified multiplier derivative"},
                {"max_whole_run_gradient_difference",difference},{"max_kz_departure_from_prescribed_0_01",kzDeparture},{"cases",cases},
                {"interpretation","Transported nominally uniform temperature can acquire roundoff stratification and change Richardson/convective branches. Scaling a recorded Kz trajectory does not differentiate those branches or the thermal feedback."}};
            std::ofstream stream(argv[1]);require(bool(stream),"Cannot write diagnostic");stream<<diagnostic.dump(2)<<'\n';
            std::cout<<diagnostic.dump(2)<<'\n';return 0;
        }
        json cases=json::array();const double pi=std::acos(-1.),a=.01*pi*pi*120/16;
        for(const auto &[layers,dt]:std::vector<std::pair<int,double>>{{16,10},{16,5},{16,2.5},{4,2.5},{8,2.5},{32,2.5}}) {
            auto r=run(layers,dt);const auto full=r.trace.jvp(std::vector<double>(layers,0),{},1);
            std::vector<double> gradient(full.end()-layers,full.end());
            const auto plus=run(layers,dt,1.0001,false),minus=run(layers,dt,.9999,false);
            double fd=0,mass=0,rms=0;
            for(int k=0;k<layers;++k) {
                fd=std::max(fd,std::abs((plus.values[k]-minus.values[k])/.0002-gradient[k]));mass+=gradient[k]*400/layers;
            }
            json samples=json::array();
            // Common physical depths remain within all centre supports.
            for(double depth:{.75,1.25,2.75,3.25}) {
                const double analytic=-.2*a*std::exp(-a)*std::cos(pi*depth/4);
                const double value=sample(gradient,depth);rms+=(value-analytic)*(value-analytic);
                samples.push_back({{"depth_m",depth},{"discrete_derivative",value},{"analytic_derivative",analytic}});
            }
            require(fd<1e-8,"diffusion derivative finite difference gate failed: "+std::to_string(fd));require(std::abs(mass)<1e-10,"tangent mass gate failed");
            cases.push_back({{"layers",layers},{"dt_s",dt},{"steps",r.trace.count},{"gradient_rms_error",std::sqrt(rms/4)},
                {"max_fd_error",fd},{"tangent_mass",mass},{"trace_payload_bytes",r.trace.bytes()},{"samples",samples}});
        }
        auto baseline=run(16,2.5);auto jvp=baseline.trace.jvp(std::vector<double>(16,0),{},1);
        const double sensitivity=sample(std::vector<double>(jvp.end()-16,jvp.end()),.75);
        const double centre=sample(baseline.values,.75);
        std::mt19937 rng(60213);std::uniform_real_distribution<double> uniform(-1,1);
        std::vector<double> draws(512);for(auto &v:draws)v=uniform(rng);
        json ensembles=json::array();
        for(double radius:{.05,.5,.9}) {
            double mean=0,square=0,linearSquare=0,maxNonlinear=0;int failures=0;
            for(double z:draws) {
                try {
                    const double value=sample(run(16,2.5,1+radius*z,false).values,.75);
                    mean+=value;square+=value*value;
                    const double linear=centre+sensitivity*radius*z;linearSquare+=(linear-centre)*(linear-centre);
                    maxNonlinear=std::max(maxNonlinear,std::abs(value-linear));
                }catch(const std::exception &){++failures;}
            }
            require(failures==0,"diffusion ensemble failure");mean/=draws.size();
            ensembles.push_back({{"uniform_radius",radius},{"samples",draws.size()},{"failed",failures},
                {"actual_mean",mean},{"baseline",centre},{"actual_sd",std::sqrt(std::max(0.,square/draws.size()-mean*mean))},
                {"local_sd",std::abs(sensitivity)*radius/std::sqrt(3.)},
                {"max_nonlinear_difference",maxNonlinear},{"linear_rms_about_baseline",std::sqrt(linearSquare/draws.size())}});
        }
        json report={{"scope","Actual no-flow vertical diffusion with uniform prescribed Kz; initial layer averages and observation depths independent of implementation"},
            {"seed",60213},{"refinement",cases},{"ensembles",ensembles},
            {"interpretation","Uniform positive diffusivity multiplier uncertainty; local SD is a covariance approximation, not a 95% probability interval. Wider uncertainty exposes nonlinear bias."}};
        std::ofstream stream(argv[1]);if(!stream)throw std::runtime_error("Cannot write report");stream<<report.dump(2)<<'\n';
        std::cout<<report.dump(2)<<'\n';return 0;
    }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
