// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
// Reproducible closed-basin diagnostics. Storage and conveying width are
// specified separately; changing one never infers a change to the other.
#include "fvqual/hydro/internalmode.h"
#include <nlohmann/json.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>

using nlohmann::json;
using namespace FVQual;
namespace {
void require(bool ok,const std::string &message) {
    if (!ok) throw std::runtime_error(message);
}
HydroCouple::SDK::IO::MeshDefinition channel(int columns,double width) {
    HydroCouple::SDK::IO::MeshDefinition result;
    result.meshName="closed-channel";
    for (int i=0;i<=columns;++i) {
        const double x=1000.*i/columns;
        for (double y:{0.,width}) {
            result.nodeX.push_back(x); result.nodeY.push_back(y); result.nodeZ.push_back(-20);
        }
    }
    result.faceNodeOffsets.push_back(0);
    for (int i=0;i<columns;++i) {
        for (int n:{2*i,2*i+2,2*i+3,2*i+1}) result.faceNodes.push_back(n);
        result.faceNodeOffsets.push_back(result.faceNodes.size());
    }
    return result;
}
struct Case {
    std::string id;
    int columns=20,layers=20;
    double dt=15,width=100,storageWidth=100,wind=.02;
    bool taper=true;
    double diffusivityScale=1,viscosityScale=1;
    Transport::AdvectionScheme scheme=Transport::AdvectionScheme::MusclFct;
};
json momentum(const std::vector<Hydro::MomentumColumnAudit> &columns) {
    Hydro::MomentumColumnAudit total;
    for (const auto &a:columns) {
        total.windWork+=a.windWork; total.viscosityDissipation+=a.viscosityDissipation;
        total.bedDissipation+=a.bedDissipation; total.implicitDissipation+=a.implicitDissipation;
        total.modeEnergy+=a.modeEnergy; total.modeImpulseX+=a.modeImpulseX;
        total.verticalResidualX+=a.verticalResidualX; total.verticalEnergyResidual+=a.verticalEnergyResidual;
        total.maxRelativeImpulseResidual=std::max(total.maxRelativeImpulseResidual,a.maxRelativeImpulseResidual);
        total.maxRelativeEnergyResidual=std::max(total.maxRelativeEnergyResidual,a.maxRelativeEnergyResidual);
        total.maxAbsoluteImpulseResidual=std::max(total.maxAbsoluteImpulseResidual,a.maxAbsoluteImpulseResidual);
        total.maxAbsoluteEnergyResidual=std::max(total.maxAbsoluteEnergyResidual,a.maxAbsoluteEnergyResidual);
    }
    return { {"wind_work_J",total.windWork}, {"viscosity_dissipation_J",total.viscosityDissipation},
        {"bed_dissipation_J",total.bedDissipation}, {"implicit_dissipation_J",total.implicitDissipation},
        {"mode_energy_J",total.modeEnergy}, {"mode_impulse_x_kg_m_s",total.modeImpulseX},
        {"vertical_residual_x_kg_m_s",total.verticalResidualX}, {"vertical_energy_residual_J",total.verticalEnergyResidual},
        {"max_relative_impulse_residual",total.maxRelativeImpulseResidual},
        {"max_relative_energy_residual",total.maxRelativeEnergyResidual},
        {"max_absolute_impulse_residual_kg_m_s",total.maxAbsoluteImpulseResidual},
        {"max_absolute_energy_residual_J",total.maxAbsoluteEnergyResidual} };
}
json run(const Case &test) {
    const auto definition=channel(test.columns,test.width);
    Mesh::HorizontalMesh horizontal;
    Mesh::LayeredMesh mesh;
    State::SpeciesRegistry registry;
    State::StateStore state;
    Hydro::InternalModeSolver solver;
    std::string message;
    require(horizontal.build(definition,message),message);
    Mesh::SigmaSpec spec;
    spec.layerCount=test.layers; spec.minLayerThickness=0;
    require(mesh.build(horizontal,spec,0,message),message);
    const double storageArea=1000.*test.storageWidth/test.columns;
    // A(z)=Asurface*(1+0.04*z), with 20% contact area at z=-20.
    if (test.taper) {
        for (int c=0;c<test.columns;++c)
            require(mesh.setHypsography(c,std::vector<double>{-20,0},
                std::vector<double>{.2*storageArea,storageArea},message),message);
        require(mesh.updateSurface(std::vector<double>(test.columns,0),message),message);
    }
    const int temperature=registry.addTemperature(10,message);
    require(temperature>=0,message);
    State::SpeciesDefinition dye; dye.id="tracer"; dye.units="kg/m3";
    const int tracer=registry.add(dye,message); require(tracer>=0,message);
    state.allocate(mesh,registry);
    for (int c=0;c<test.columns;++c)
        for (int k=0;k<test.layers;++k) {
            const double bottom=mesh.z(c,k+1),top=mesh.z(c,k);
            // Independent exact first moment of the linear storage curve.
            const double moment=test.taper ? storageArea*((top*top-bottom*bottom)/2+
                .04*(top*top*top-bottom*bottom*bottom)/3) :
                horizontal.faceArea()[c]*(top*top-bottom*bottom)/2;
            const auto cell=mesh.cell(c,k);
            const double z=moment/mesh.cellVolume()[cell];
            state.speciesAt(temperature,cell)=5+.75*(z+20);
            const double x=horizontal.faceCentroidX()[c];
            state.speciesAt(tracer,cell)=.25+.25*std::cos(2*3.141592653589793*x/1000);
        }
    Hydro::InternalModeOptions options;
    options.manningN=.03; options.maxTimeStep=test.dt;
    options.windStressX=test.wind; options.surfaceHeatFlux=100;
    options.stableViscosityScale=test.viscosityScale;
    options.stableDiffusivityScale=test.diffusivityScale;
    options.advectionScheme=test.scheme; options.auditHeatTransport=true; options.auditMomentum=true;
    require(solver.initialize(mesh,state,registry,options,message),message);
    const double volume=mesh.totalVolume(),heat=solver.totalSpeciesMass(temperature);
    const double mass=solver.totalSpeciesMass(tracer);
    double maxWaterResidual=0,maxHeatResidual=0,maxTracerResidual=0,maxContinuity=0;
    double maxModeMeanResidual=0;
    json frames=json::array();
    constexpr double duration=7200;
    for (double target=0;target<=duration;target+=600) {
        while (solver.time()<target-1e-10) {
            const double dt=std::min({test.dt,solver.stableTimeStep(),target-solver.time()});
            require(solver.step(dt,message),message);
            maxWaterResidual=std::max(maxWaterResidual,std::abs(mesh.totalVolume()-volume)/volume);
            maxHeatResidual=std::max(maxHeatResidual,std::abs(solver.scalarBudgetResidual(temperature))/heat);
            maxTracerResidual=std::max(maxTracerResidual,std::abs(solver.totalSpeciesMass(tracer)-mass)/mass);
            maxContinuity=std::max(maxContinuity,solver.diagnostics().surfaceFluxResidual);
            for (int c=0;c<test.columns;++c) {
                double sum=0;
                for (int k=0;k<test.layers;++k) sum+=mesh.thickness()[mesh.cell(c,k)]*state.u[mesh.cell(c,k)];
                maxModeMeanResidual=std::max(maxModeMeanResidual,
                    std::abs(sum/mesh.columnDepth()[c]-solver.external().velocityX()[c]));
            }
        }
        double speed=0;
        for (size_t i=0;i<state.u.size();++i) speed=std::max(speed,std::hypot(state.u[i],state.v[i]));
        json frame={{"time_s",solver.time()},{"temperature",state.speciesBlock(temperature)},
            {"tracer",state.speciesBlock(tracer)},{"velocity_x",state.u},{"velocity_y",state.v},
            {"interfaces",mesh.interfaceZ()},{"volumes",mesh.cellVolume()},
            {"viscosity",solver.eddyViscosity()},{"diffusivity",solver.eddyDiffusivity()},
            {"pressure_x",solver.baroclinicAccelerationX()},
            {"scalars",{{"max_speed",speed},{"heat_inventory_C_m3",solver.totalSpeciesMass(temperature)},
                {"tracer_inventory_kg",solver.totalSpeciesMass(tracer)},
                {"water_relative_residual",maxWaterResidual},{"heat_relative_residual",maxHeatResidual},
                {"tracer_relative_residual",maxTracerResidual},{"continuity_residual_m_s",maxContinuity},
                {"mode_mean_residual_m_s",maxModeMeanResidual}}}};
        frame["scalars"].update(momentum(solver.momentumAudit()));
        frames.push_back(std::move(frame));
    }
    require(maxWaterResidual<1e-8 && maxHeatResidual<1e-8 && maxTracerResidual<1e-8,
            test.id+": closed-basin conservation gate failed");
    const auto finalMomentum=momentum(solver.momentumAudit());
    require(finalMomentum["max_relative_impulse_residual"].get<double>()<1e-8 &&
            finalMomentum["max_relative_energy_residual"].get<double>()<1e-8 &&
            maxModeMeanResidual<1e-12 && maxContinuity<1e-10,
            test.id+": vertical momentum or split-mode consistency gate failed");
    return {{"id",test.id},{"columns",test.columns},{"layers",test.layers},
        {"dt_max_s",test.dt},{"conveyance_width_m",test.width},{"storage_width_m",test.storageWidth},
        {"tapered_storage",test.taper},{"wind_stress_Pa",test.wind},
        {"scalar_scale",test.diffusivityScale},{"viscosity_scale",test.viscosityScale},
        {"advection",test.scheme==Transport::AdvectionScheme::MusclFct ? "MUSCL-FCT" : "upwind"},
        {"geometry",{{"node_x",definition.nodeX},{"node_y",definition.nodeY},{"node_z",definition.nodeZ},
                     {"face_nodes",definition.faceNodes},{"face_offsets",definition.faceNodeOffsets}}},
        {"frames",frames}};
}
json pressureReview() {
    json results=json::array();
    for (int layers:{10,20,40,80}) {
        Mesh::HorizontalMesh horizontal; Mesh::LayeredMesh mesh;
        std::string message;
        require(horizontal.build(channel(20,100),message),message);
        for (int c=0;c<20;++c) horizontal.faceBed()[c]=c<10 ? -4 : -20;
        Mesh::SigmaSpec spec;spec.layerCount=layers;spec.minLayerThickness=0;
        require(mesh.build(horizontal,spec,0,message),message);
        for (bool linear:{true,false}) {
            std::vector<double> density(mesh.cellCount());
            for (int c=0;c<20;++c) for (int k=0;k<layers;++k) {
                const double z=.5*(mesh.z(c,k)+mesh.z(c,k+1));
                density[mesh.cell(c,k)]=linear ? 1000-.05*z : 1001-std::exp(z/8);
            }
            Hydro::HydrostaticPressureClosure closure;
            std::vector<double> ax,ay;
            closure.baroclinicAcceleration(mesh,density,1000,9.80665,ax,ay);
            double peak=0;
            for (size_t i=0;i<ax.size();++i) peak=std::max(peak,std::hypot(ax[i],ay[i]));
            results.push_back({{"layers",layers},{"profile",linear ? "linear" : "exponential"},
                {"max_acceleration_m_s2",peak}});
        }
    }
    return results;
}
}
int main(int argc,char **argv) {
    if (argc!=2) {std::cerr<<"Usage: fvqual-circulation-review <output.json>\n";return 2;}
    try {
        std::vector<Case> cases;
        auto add=[&](std::string id) -> Case& { cases.push_back(Case{});cases.back().id=std::move(id);return cases.back(); };
        add("baseline"); add("dt7_5").dt=7.5; add("dt3_75").dt=3.75;
        auto &l10=add("layers10"); l10.layers=10;l10.dt=3.75;
        auto &l40=add("layers40"); l40.layers=40;l40.dt=3.75;
        auto &l80=add("layers80"); l80.layers=80;l80.dt=3.75;
        auto &n10=add("columns10");n10.columns=10;n10.dt=3.75;
        auto &n40=add("columns40");n40.columns=40;n40.dt=3.75;
        add("scalar_mixing_off").diffusivityScale=0;
        add("viscosity_off").viscosityScale=0;
        add("no_wind").wind=0;
        add("prismatic").taper=false;
        add("narrow_conveyance").width=50;
        add("upwind").scheme=Transport::AdvectionScheme::Upwind;
        json result={{"schema","fvqual.circulation-review/1.0"},
            {"scope","Synthetic closed tapered basin; storage and rectangular conveying width prescribed independently; not field validation or full subgrid momentum closure"},
            {"momentum_scope","Implicit vertical momentum and energy balances; mode correction recorded separately"},
            {"cases",json::array()},{"pressure",pressureReview()}};
        for (const auto &test:cases) {
            std::cout<<"Running "<<test.id<<std::endl;
            result["cases"].push_back(run(test));
        }
        const auto path=std::filesystem::path(argv[1]);
        std::filesystem::create_directories(path.parent_path());
        std::ofstream stream(path); stream<<result.dump(2)<<'\n';
        require(static_cast<bool>(stream),"Could not write review JSON");
        return 0;
    } catch (const std::exception &error) {std::cerr<<error.what()<<'\n';return 1;}
}
