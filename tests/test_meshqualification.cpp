// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
#include "fvqual/hydro/internalmode.h"
#include "fvqual/hydro/equationofstate.h"
#include "fvqual/component/fvqualcomponent.h"
#include "fvqual/transport/euleriantransport.h"
#include "hydrocouplesdk/io/meshpackage.h"
#include "testmeshes.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <numeric>
#include <sstream>

namespace {
using HydroCouple::SDK::IO::MeshDefinition;
constexpr double pi = 3.14159265358979323846;

// Shared nodes and fixed seams. Alternating cell types exercise every edge
// orientation; the hole and L-shaped corridor are intentional domain features.
MeshDefinition corpus(int n, const std::string &family) {
    auto mesh = FVQual::Testing::basinGrid(n,n,1,1);
    mesh.nodeZ.assign(mesh.nodeCount(), -5.0);
    for (size_t k = 0; k < mesh.nodeX.size(); ++k) {
        double &x = mesh.nodeX[k], &y = mesh.nodeY[k];
        if (family == "graded") x = x*x;
        if (family == "skewed") x += 0.3 * std::sin(pi*x) * std::sin(pi*y);
        if (family == "rotated") { const auto a=x,b=y; x=1e6+.8*a-.6*b; y=6e6+.6*a+.8*b; }
    }
    const auto oldNodes = mesh.faceNodes;
    mesh.faceNodeOffsets = {0}; mesh.faceNodes.clear();
    auto face = [&](std::initializer_list<int64_t> nodes) {
        mesh.faceNodes.insert(mesh.faceNodes.end(),nodes);
        mesh.faceNodeOffsets.push_back(mesh.faceNodes.size());
    };
    for (int j=0;j<n;++j) for (int i=0;i<n;++i) {
        if (family=="island" && i>=n/3 && i<2*n/3 && j>=n/3 && j<2*n/3) continue;
        if (family=="corridor" && i>=n/2 && j>=n/2) continue;
        const size_t k=4*(j*n+i);
        const auto a=oldNodes[k], b=oldNodes[k+1], c=oldNodes[k+2], d=oldNodes[k+3];
        if (family=="triangles" || (family!="quads" && (i+j)%2==0)) {
            face({a,b,c}); face({a,c,d});
        } else face({a,b,c,d});
    }
    if (family=="reordered") {
        const auto offsets=mesh.faceNodeOffsets, nodes=mesh.faceNodes;
        mesh.faceNodeOffsets={0}; mesh.faceNodes.clear();
        for (int64_t f=offsets.size()-2;f>=0;--f) {
            for (auto k=offsets[f+1]; k>offsets[f]; --k) mesh.faceNodes.push_back(nodes[k-1]);
            mesh.faceNodeOffsets.push_back(mesh.faceNodes.size());
        }
    }
    return mesh;
}
const std::vector<std::string> families{
    "triangles","quads","mixed","corridor","island","graded","skewed","rotated","reordered"};

MeshDefinition channel(int n,const std::string &family,double length=1000,double width=50) {
    auto mesh=FVQual::Testing::channelStrip(n,length,width);
    const auto nodes=mesh.faceNodes;
    mesh.faceNodes.clear(); mesh.faceNodeOffsets={0}; mesh.nodeZ.assign(mesh.nodeCount(),-20);
    auto face=[&](std::initializer_list<int64_t> f) {
        mesh.faceNodes.insert(mesh.faceNodes.end(),f); mesh.faceNodeOffsets.push_back(mesh.faceNodes.size());
    };
    for (int i=0;i<n;++i) {
        auto a=nodes[4*i],b=nodes[4*i+1],c=nodes[4*i+2],d=nodes[4*i+3];
        if (family=="triangles" || (family=="mixed" && i%2==0)) { face({a,b,c}); face({a,c,d}); }
        else face({a,b,c,d});
    }
    return mesh;
}
struct Model {
    FVQual::Mesh::HorizontalMesh horizontal;
    FVQual::Mesh::LayeredMesh layered;
    FVQual::State::SpeciesRegistry registry;
    FVQual::State::StateStore state;
    FVQual::Hydro::InternalModeSolver solver;
    int tracer=-1;
    void build(const MeshDefinition &definition, double diffusivity=0, double maxStep=.002,int layers=2) {
        std::string message;
        ASSERT_TRUE(horizontal.build(definition,message)) << message;
        FVQual::Mesh::SigmaSpec sigma; sigma.layerCount=layers; sigma.minLayerThickness=0;
        ASSERT_TRUE(layered.build(horizontal,sigma,0,message)) << message;
        ASSERT_GE(registry.addTemperature(10,message),0);
        FVQual::State::SpeciesDefinition species; species.id="tracer"; species.units="kg/m3";
        tracer=registry.add(species,message); ASSERT_GE(tracer,0);
        state.allocate(layered,registry);
        FVQual::Hydro::InternalModeOptions options;
        options.manningN=0; options.horizontalDiffusivity=diffusivity; options.maxTimeStep=maxStep;
        options.neutralViscosity=0; options.backgroundViscosity=0; options.backgroundDiffusivity=0;
        ASSERT_TRUE(solver.initialize(layered,state,registry,options,message)) << message;
    }
};
double inventory(const std::vector<double> &values, const FVQual::Mesh::LayeredMesh &mesh) {
    double result=0; for (size_t i=0;i<values.size();++i) result+=values[i]*mesh.cellVolume()[i];
    return result;
}

void writeReviewPackage(const HydroCouple::SDK::IO::MeshPackage &package,
                        const std::filesystem::path &directory,const std::string &name) {
    std::string message;
    std::filesystem::create_directories(directory);
    const auto path=directory/(name+".mesh.json"); auto nc=path; nc.replace_extension(".nc");
    std::filesystem::remove(path); std::filesystem::remove(nc);
    ASSERT_TRUE(HydroCouple::SDK::IO::writeMeshPackage(path,package,message)) << message;
    nlohmann::json composition{
        {"schema_version","1.1"},
        {"metadata",{{"id",name},{"caption",name+" — synthetic wet reservoir mesh qualification"}}},
        {"components",nlohmann::json::array({{
            {"id","reservoir"},{"info",{{"component_info_id","fvqual.component"}}},
            {"arguments",{{"geometry_counts",{{"values",{package.mesh.faceCount(),2}}}},
                {"geometry_values",{{"values",{1,1,5,10,.02,.04}}}},
                {"mesh",{{"@uri",name+".mesh.json"}}}}}
        }})},
        {"workflow",{{"strategy","time_stepped"}}},
        {"writers",nlohmann::json::array({{{"type","csv"},{"path",name+"_output/temperature.csv"}}})},
        {"run",{{"manifest",name+"_output/run.json"}}}
    };
    std::ofstream(directory/(name+".composition.json")) << composition.dump(2);
}
}

TEST(MeshQualification, GlobalDefectsCannotBecomeDisconnectedFluxDomains) {
    FVQual::Mesh::HorizontalMesh mesh;
    std::string message;
    auto base=FVQual::Testing::quadStrip(1);
    auto append=[&](const std::vector<double> &x,const std::vector<double> &y) {
        auto d=base; const auto start=d.nodeCount();
        d.nodeX.insert(d.nodeX.end(),x.begin(),x.end()); d.nodeY.insert(d.nodeY.end(),y.begin(),y.end());
        for (size_t k=0;k<x.size();++k) d.faceNodes.push_back(start+k);
        d.faceNodeOffsets.push_back(d.faceNodes.size()); return d;
    };
    // Crossing faces, contained face, unmatched midpoint and coincident nodes.
    for (const auto &d : {append({.5,1.5,1.5,.5},{.2,.2,.8,.8}),
                          append({.2,.8,.8,.2},{.2,.2,.8,.8}),
                          append({1,2,2},{.5,0,1}),
                          append({1,2,2,1},{0,0,1,1})}) {
        EXPECT_FALSE(mesh.build(d,message)); EXPECT_FALSE(message.empty());
    }
    // Disconnected valid islands must be counted, not silently joined.
    ASSERT_TRUE(mesh.build(append({2,3,3,2},{0,0,1,1}),message)) << message;
    EXPECT_EQ(mesh.geometryReport().connectedDomains,2);
    EXPECT_EQ(mesh.geometryReport().boundaryLoops,2);
}

TEST(MeshQualification, ConvexPolygonsAndIrregularBedRemainValidControlVolumes) {
    FVQual::Mesh::HorizontalMesh mesh;
    MeshDefinition hexagon; hexagon.faceNodeOffsets={0,6}; hexagon.faceNodes={0,1,2,3,4,5};
    for (int i=0;i<6;++i) { hexagon.nodeX.push_back(std::cos(i*pi/3)); hexagon.nodeY.push_back(std::sin(i*pi/3)); }
    std::string message;
    ASSERT_TRUE(mesh.build(hexagon,message)) << message;
    EXPECT_EQ(mesh.geometryReport().polygons,1);
    for (const auto &family:families) {
        SCOPED_TRACE(family);
        ASSERT_TRUE(mesh.build(corpus(8,family),message)) << message;
        for (int64_t c=0;c<mesh.faceCount();++c)
            mesh.faceBed()[c]=-5+std::sin(3*mesh.faceCentroidX()[c])+std::cos(2*mesh.faceCentroidY()[c]);
        FVQual::Hydro::BarotropicSolver solver;
        FVQual::Hydro::BarotropicOptions options; options.frictionEnabled=false;
        ASSERT_TRUE(solver.initialize(mesh,options,message)); solver.setUniformSurface(0);
        const double before=solver.totalVolume();
        ASSERT_TRUE(solver.advance(.1,message)) << message;
        EXPECT_NEAR(solver.totalVolume(),before,1e-11);
        for (double u:solver.velocityX()) EXPECT_NEAR(u,0,1e-11);
        for (double v:solver.velocityY()) EXPECT_NEAR(v,0,1e-11);
    }
}

TEST(MeshQualification, SurfaceSeicheMatchesTheoryOnNativeTriQuadChannels) {
    nlohmann::json results=nlohmann::json::array();
    const double period=2000/std::sqrt(9.80665*20);
    for (const std::string family:{"triangles","quads","mixed"}) {
        SCOPED_TRACE(family);
        FVQual::Mesh::HorizontalMesh mesh; std::string message;
        ASSERT_TRUE(mesh.build(channel(80,family),message));
        FVQual::Hydro::BarotropicSolver solver;
        FVQual::Hydro::BarotropicOptions options; options.frictionEnabled=false; options.cflNumber=.45;
        ASSERT_TRUE(solver.initialize(mesh,options,message));
        std::vector<double> surface(mesh.faceCount());
        for (int64_t c=0;c<mesh.faceCount();++c) surface[c]=.02*std::cos(pi*mesh.faceCentroidX()[c]/1000);
        ASSERT_TRUE(solver.setSurface(surface,message));
        std::vector<double> crossings; double previous=surface[0],before=0;
        while (solver.time()<2.5*period) {
            ASSERT_TRUE(solver.advance(.5,message)) << message;
            const double value=solver.surfaceElevation()[0];
            if ((previous>0)!=(value>0)) crossings.push_back(before+.5*previous/(previous-value));
            previous=value; before=solver.time();
        }
        ASSERT_GE(crossings.size(),4u);
        const double measured=2*(crossings.back()-crossings.front())/(crossings.size()-1);
        EXPECT_NEAR(measured,period,.06*period);
        results.push_back({{"family",family},{"period_seconds",measured},{"theory_seconds",period}});
    }
    std::ofstream(std::string(FVQual_TEST_OUTPUT_DIR)+"/mesh_surface_seiche.json") << results.dump(2);
}

TEST(MeshQualification, OpenBoundaryClosesWaterAndSpeciesBudgetsAcrossCellTypes) {
    for (const std::string family:{"triangles","quads","mixed"}) {
        SCOPED_TRACE(family);
        Model model; ASSERT_NO_FATAL_FAILURE(model.build(channel(16,family),0,1));
        const auto &h=model.horizontal;
        std::vector<int64_t> inletEdges,outletEdges;
        for (int64_t e=0;e<h.edgeCount();++e) if (h.isBoundaryEdge(e)) {
            if (h.edgeMidX()[e]<1e-8) inletEdges.push_back(e);
            if (h.edgeMidX()[e]>1000-1e-8) outletEdges.push_back(e);
        }
        FVQual::Hydro::BarotropicSolver::BoundaryCondition inflow,level;
        inflow.type=FVQual::Hydro::BarotropicSolver::BoundaryType::Inflow; inflow.specificDischarge=.05;
        level.type=FVQual::Hydro::BarotropicSolver::BoundaryType::Stage; level.stage=0;
        std::string message;
        const int inlet=model.solver.external().addBoundary(inflow,inletEdges,message);
        const int outlet=model.solver.external().addBoundary(level,outletEdges,message);
        ASSERT_GE(inlet,0); ASSERT_GE(outlet,0);
        ASSERT_TRUE(model.solver.bindOpenBoundaries(FVQual::Forcing::Placement::Uniform,message)) << message;
        std::fill(model.state.speciesBlock(model.tracer).begin(),model.state.speciesBlock(model.tracer).end(),.25);
        const double water=model.layered.totalVolume(),mass=model.solver.totalSpeciesMass(model.tracer);
        double exchanged=0;
        for (int i=0;i<20;++i) {
            ASSERT_TRUE(model.solver.step(.5,message)) << message;
            exchanged+=model.solver.boundaryVolume(inlet)+model.solver.boundaryVolume(outlet);
            EXPECT_LT(std::abs(model.solver.diagnostics().surfaceFluxResidual),1e-10);
        }
        EXPECT_GT(model.solver.boundaryVolume(inlet),0);
        EXPECT_NEAR(model.layered.totalVolume()-water,exchanged,1e-7);
        EXPECT_NEAR(model.solver.totalSpeciesMass(model.tracer)-mass,.25*exchanged,1e-7);
        for (double value:model.state.speciesBlock(model.tracer)) EXPECT_NEAR(value,.25,1e-10);
    }
}

TEST(MeshQualification, SmoothAdvectionConvergesAndImprovesOnUpwindAcrossCellTypes) {
    nlohmann::json records=nlohmann::json::array();
    for (const std::string family:{"triangles","quads","mixed"}) {
        double previous=0;
        for (int n:{24,48,96}) {
            Model model; ASSERT_NO_FATAL_FAILURE(model.build(channel(n,family,1,.1)));
            const auto &h=model.horizontal; const auto &m=model.layered;
            std::vector<double> flux(h.edgeCount()*2,0),vertical(m.interfaceSlotCount(),0);
            for (int64_t e=0;e<h.edgeCount();++e) if (!h.isBoundaryEdge(e))
                for (int k=0;k<2;++k) flux[2*e+k]=.5*h.edgeNormalX()[e]*h.edgeLength()[e]*10;
            auto gaussian=[](double x,double center){ return std::exp(-.5*std::pow((x-center)/.04,2)); };
            double upwind=0,muscl=0;
            for (auto scheme:{FVQual::Transport::AdvectionScheme::Upwind,FVQual::Transport::AdvectionScheme::MusclFct}) {
                FVQual::Transport::EulerianTransport engine; engine.initialize(m);
                std::vector<double> values(m.cellCount());
                for (int64_t c=0;c<h.faceCount();++c) for (int k=0;k<2;++k)
                    values[m.cell(c,k)]=gaussian(h.faceCentroidX()[c],.3);
                const double mass=inventory(values,m),dt=.2/n;
                for (double t=0;t<.3-1e-14;t+=dt) engine.advect(m,m.cellVolume(),flux,vertical,std::min(dt,.3-t),scheme,values);
                EXPECT_NEAR(inventory(values,m),mass,1e-10);
                double error=0;
                for (int64_t c=0;c<h.faceCount();++c) {
                    const double delta=values[m.cell(c,0)]-gaussian(h.faceCentroidX()[c],.45);
                    error+=h.faceArea()[c]*delta*delta;
                }
                error=std::sqrt(error/h.totalArea());
                (scheme==FVQual::Transport::AdvectionScheme::Upwind?upwind:muscl)=error;
            }
            SCOPED_TRACE(family+"/"+std::to_string(n));
            EXPECT_LT(muscl,upwind);
            if (previous>0) EXPECT_LT(muscl,.85*previous);
            records.push_back({{"family",family},{"resolution",n},{"muscl_l2",muscl},{"upwind_l2",upwind}});
            previous=muscl;
        }
    }
    std::ofstream(std::string(FVQual_TEST_OUTPUT_DIR)+"/mesh_advection_convergence.json") << records.dump(2);
}

TEST(MeshQualification, InternalSeicheRemainsPhysicalAcrossCellTypes) {
    const double reduced=9.80665*(FVQual::Hydro::EquationOfState::density(5)-FVQual::Hydro::EquationOfState::density(25))/1000;
    const double period=2000/std::sqrt(reduced*5*15/20);
    nlohmann::json records=nlohmann::json::array();
    for (const std::string family:{"triangles","quads","mixed"}) {
        SCOPED_TRACE(family);
        Model model; ASSERT_NO_FATAL_FAILURE(model.build(channel(40,family),0,20,20));
        const int temperature=model.registry.temperatureIndex();
        for (int64_t c=0;c<model.layered.columnCount();++c) for (int k=0;k<20;++k) {
            const double z=.5*(model.layered.z(c,k)+model.layered.z(c,k+1));
            const double tilt=-5+1.5*(model.horizontal.faceCentroidX()[c]/1000-.5)*2;
            model.state.speciesAt(temperature,model.layered.cell(c,k))=z>tilt?25:5;
        }
        auto interface=[&]() {
            for (int k=0;k<20;++k) if (model.state.speciesAt(temperature,model.layered.cell(0,k))<15)
                return .5*(model.layered.z(0,k)+model.layered.z(0,k+1));
            return model.layered.z(0,20);
        };
        std::vector<double> crossings; double previous=interface()+5,before=0;
        std::string message;
        while (model.solver.time()<2.2*period) {
            ASSERT_TRUE(model.solver.advance(20,message)) << message;
            const double value=interface()+5;
            if ((previous>0 && value<=0)||(previous<0 && value>=0))
                crossings.push_back(before+20*previous/(previous-value));
            previous=value; before=model.solver.time();
        }
        ASSERT_GE(crossings.size(),2u);
        const double measured=2*(crossings.back()-crossings.front())/(crossings.size()-1);
        EXPECT_NEAR(measured,period,.15*period);
        EXPECT_LT(std::abs(model.solver.diagnostics().surfaceFluxResidual),1e-10);
        records.push_back({{"family",family},{"period_seconds",measured},{"theory_seconds",period}});
    }
    std::ofstream(std::string(FVQual_TEST_OUTPUT_DIR)+"/mesh_internal_seiche.json") << records.dump(2);
}

TEST(MeshQualification, CorpusClosesGeometryAndConservesStillWaterAndUniformTracer) {
    nlohmann::json records=nlohmann::json::array();
    for (const auto &family:families) {
        SCOPED_TRACE(family);
        Model model; ASSERT_NO_FATAL_FAILURE(model.build(corpus(8,family)));
        const auto report=model.horizontal.geometryReport();
        EXPECT_EQ(report.connectedDomains,1);
        EXPECT_EQ(report.boundaryLoops,family=="island" ? 2:1);
        EXPECT_LT(report.maximumRelativeClosureError,1e-12);
        EXPECT_GT(report.minimumNormalDistance,0);
        if (family=="triangles") EXPECT_EQ(report.quads,0);
        if (family=="quads") EXPECT_EQ(report.triangles,0);
        auto values=model.state.speciesBlock(model.tracer);
        std::fill(values.begin(),values.end(),.25);
        const double water=model.layered.totalVolume(), mass=model.solver.totalSpeciesMass(model.tracer);
        std::string message;
        ASSERT_TRUE(model.solver.advance(.03,message)) << message;
        EXPECT_NEAR(model.layered.totalVolume(),water,1e-11);
        EXPECT_NEAR(model.solver.totalSpeciesMass(model.tracer),mass,1e-11);
        for (double value:values) EXPECT_NEAR(value,.25,1e-12);
        for (double velocity:model.state.u) EXPECT_NEAR(velocity,0,1e-11);
        for (double velocity:model.state.v) EXPECT_NEAR(velocity,0,1e-11);
        records.push_back({{"family",family},{"triangles",report.triangles},{"quads",report.quads},
            {"boundary_loops",report.boundaryLoops},{"closure_error",report.maximumRelativeClosureError},
            {"nonorthogonality_degrees",report.maximumNonOrthogonalityDegrees},
            {"volume_residual",model.layered.totalVolume()-water},
            {"mass_residual",model.solver.totalSpeciesMass(model.tracer)-mass}});
#ifdef FVQual_HAVE_NETCDF
        HydroCouple::SDK::IO::MeshPackage package;
        package.mesh=corpus(8,family); package.horizontalCrs="LOCAL";
        package.verticalDatum="synthetic-flat-bed"; package.elevationKind="bed-elevation";
        for (int64_t n=0;n<package.mesh.nodeCount();++n) package.nodeIds.push_back(std::to_string(n));
        for (int64_t f=0;f<package.mesh.faceCount();++f) package.faceIds.push_back(std::to_string(f));
        package.regionTags.assign(package.mesh.faceCount(),family);
        package.provenance={{"generator","FVQual verification corpus"},{"family",family},{"quad_conversion","none"}};
        ASSERT_NO_FATAL_FAILURE(writeReviewPackage(package,std::filesystem::path(FVQual_TEST_OUTPUT_DIR)/"mesh_review",family));
#endif
    }
    std::ofstream(std::string(FVQual_TEST_OUTPUT_DIR)+"/mesh_corpus.json") << records.dump(2);
}

TEST(MeshQualification, DiffusionPreservesBoundsAndReplaysAfterLargeStep) {
    for (const std::string family:{"triangles","mixed","graded","skewed"}) {
        SCOPED_TRACE(family);
        Model model; ASSERT_NO_FATAL_FAILURE(model.build(corpus(8,family),.1,1));
        for (int64_t c=0;c<model.layered.columnCount();++c) for (int k=0;k<2;++k)
            model.state.speciesAt(model.tracer,model.layered.cell(c,k))=model.horizontal.faceCentroidX()[c]<.5?1:0;
        const double mass=model.solver.totalSpeciesMass(model.tracer);
        std::string message; auto checkpoint=model.solver.checkpoint(message); ASSERT_NE(checkpoint,nullptr);
        ASSERT_TRUE(model.solver.advance(.05,message)) << message;
        const auto result=model.state.species;
        for (double value:model.state.speciesBlock(model.tracer)) { EXPECT_GE(value,-1e-12); EXPECT_LE(value,1+1e-12); }
        EXPECT_NEAR(model.solver.totalSpeciesMass(model.tracer),mass,1e-10);
        ASSERT_TRUE(model.solver.restoreCheckpoint(*checkpoint,message));
        ASSERT_TRUE(model.solver.advance(.05,message));
        EXPECT_EQ(model.state.species,result);
    }
}

#ifdef FVQual_HAVE_NETCDF
TEST(MeshQualification, BellingeRetainsNativeTerrainAndHasAnExplicitSyntheticWetReview) {
    const auto source=std::getenv("FVQUAL_BELLINGE_MESH");
    if (!source) GTEST_SKIP() << "Set FVQUAL_BELLINGE_MESH for the local full-size qualification";
    HydroCouple::SDK::IO::MeshPackage package, restored;
    HydroCouple::SDK::IO::MeshImportOptions options;
    options.horizontalCrs="EPSG:25832"; options.verticalDatum="Bellinge-source-terrain-unverified";
    options.elevationKind="terrain-elevation";
    std::string message;
    ASSERT_TRUE(HydroCouple::SDK::IO::importSurfaceMesh(source,options,package,message)) << message;
    // The GUI's editable Bellinge mesh is regenerated. Verify its actual raw
    // records independently instead of pinning a superseded mesher's counts.
    std::ifstream raw(source);ASSERT_TRUE(raw);std::string line,section;
    int64_t triangles=0,quads=0,boundaries=0;std::vector<double> sourceX,sourceY,sourceZ;
    while(std::getline(raw,line)) {
        const auto start=line.find_first_not_of(" \t\r");if(start==std::string::npos || line[start]==';')continue;
        line=line.substr(start);
        if(line.front()=='['){section=line.substr(0,line.find(']')+1);continue;}
        if(section=="[2D_VERTICES]") {
            double x,y,z;std::istringstream record(line);ASSERT_TRUE(bool(record>>x>>y>>z));
            sourceX.push_back(x);sourceY.push_back(y);sourceZ.push_back(z);
        } else if(section=="[2D_TRIANGLES]")++triangles;
        else if(section=="[2D_QUADS]" || section=="[2D_QUADRILATERALS]")++quads;
        else if(section=="[2D_BOUNDARY_CONDITIONS]")++boundaries;
    }
    ASSERT_GT(sourceX.size(),1000u);ASSERT_EQ(package.mesh.nodeCount(),sourceX.size());
    ASSERT_EQ(package.mesh.faceCount(),triangles+quads);
    EXPECT_EQ(package.mesh.nodeX,sourceX);EXPECT_EQ(package.mesh.nodeY,sourceY);EXPECT_EQ(package.mesh.nodeZ,sourceZ);
    FVQual::Mesh::HorizontalMesh horizontal;
    ASSERT_TRUE(horizontal.build(package.mesh,message)) << message;
    auto report=horizontal.geometryReport(); EXPECT_EQ(report.triangles,triangles); EXPECT_EQ(report.quads,quads);
    EXPECT_EQ(report.boundaryLoops,1); EXPECT_EQ(report.connectedDomains,1);
    EXPECT_EQ(package.boundaries.size(),boundaries); EXPECT_LT(report.maximumRelativeClosureError,1e-12);
    const auto originalZ=package.mesh.nodeZ;
    package.provenance["qualification_original_elevation_kind"]="terrain-elevation";
    package.provenance["qualification_original_vertical_datum"]=package.verticalDatum;
    package.provenance["qualification_bed_override_m"]=-5.0;
    package.provenance["qualification_note"]="Synthetic bed for mesh/solver inspection only; not Bellinge flood predictions or reservoir bathymetry";
    package.sourceSections["FVQUAL_REVIEW_ORIGINAL_TERRAIN_Z"]=originalZ;
    package.mesh.nodeZ.assign(package.mesh.nodeCount(),-5);
    package.initialDepth.clear(); package.elevationKind="bed-elevation";
    package.verticalDatum="synthetic-flat-bed";
    // Native open boundaries and infiltration remain provenance. This test
    // advances a closed footprint without activating SWMM ownership twice.
    package.sourceSections["FVQUAL_REVIEW_NATIVE_BOUNDARIES"]=nlohmann::json::array();
    for (auto &boundary:package.boundaries) {
        package.sourceSections["FVQUAL_REVIEW_NATIVE_BOUNDARIES"].push_back({{"nodes",{boundary.node0,boundary.node1}},
            {"type",boundary.type},{"parameters",boundary.parameters},{"group",boundary.group}});
        boundary.type="WALL"; boundary.parameters.clear();
    }
    const auto directory=std::filesystem::path(FVQual_TEST_OUTPUT_DIR)/
        ("bellinge_mesh_current_"+std::to_string(sourceX.size())+"_"+std::to_string(triangles+quads));
    ASSERT_NO_FATAL_FAILURE(writeReviewPackage(package,directory,"bellinge.synthetic"));
    ASSERT_TRUE(HydroCouple::SDK::IO::readMeshPackage(directory/"bellinge.synthetic.mesh.json",restored,message)) << message;
    EXPECT_EQ(restored.mesh.faceNodes,package.mesh.faceNodes); EXPECT_EQ(restored.faceIds,package.faceIds);
    EXPECT_EQ(restored.nodeIds,package.nodeIds); EXPECT_EQ(restored.sourceSections,package.sourceSections);
    Model model; ASSERT_NO_FATAL_FAILURE(model.build(restored.mesh));
    std::fill(model.state.speciesBlock(model.tracer).begin(),model.state.speciesBlock(model.tracer).end(),.25);
    const double volume=model.layered.totalVolume(), mass=model.solver.totalSpeciesMass(model.tracer);
    ASSERT_TRUE(model.solver.advance(.004,message)) << message;
    EXPECT_NEAR(model.layered.totalVolume(),volume,volume*1e-12);
    EXPECT_NEAR(model.solver.totalSpeciesMass(model.tracer),mass,mass*1e-12);
    HydroCouple::SDK::IO::MeshPackageReference reference(restored);
    EXPECT_NE(reference.srText().find("Transverse_Mercator"),std::string::npos);
    EXPECT_NE(reference.verticalSrText().find("synthetic-flat-bed"),std::string::npos);
}
#endif

TEST(MeshQualification, MixedSkewedAndGradedCirculationIsBoundedAndConservative) {
    for (const std::string family:{"triangles","quads","mixed","graded","skewed","reordered"}) {
        SCOPED_TRACE(family);
        Model model; ASSERT_NO_FATAL_FAILURE(model.build(corpus(8,family)));
        const auto &h=model.horizontal; const auto &m=model.layered;
        std::vector<double> flux(h.edgeCount()*2,0), vertical(m.interfaceSlotCount(),0);
        auto psi=[&](int64_t n){ return .02*std::sin(pi*h.nodeX()[n])*std::sin(pi*h.nodeY()[n]); };
        for (int64_t e=0;e<h.edgeCount();++e) if (!h.isBoundaryEdge(e))
            for (int k=0;k<2;++k) flux[2*e+k]=2.5*(psi(h.edgeNodes()[2*e+1])-psi(h.edgeNodes()[2*e]));
        FVQual::Transport::EulerianTransport engine; engine.initialize(m);
        for (auto scheme:{FVQual::Transport::AdvectionScheme::Upwind,FVQual::Transport::AdvectionScheme::MusclFct}) {
            std::vector<double> values(m.cellCount(),.25);
            engine.advect(m,m.cellVolume(),flux,vertical,.001,scheme,values);
            for (double value:values) EXPECT_NEAR(value,.25,1e-12);
            for (int64_t c=0;c<h.faceCount();++c) for (int k=0;k<2;++k)
                values[m.cell(c,k)]=std::hypot(h.faceCentroidX()[c]-.5,h.faceCentroidY()[c]-.5)<.23?1:0;
            const auto mass=inventory(values,m);
            for (int step=0;step<100;++step) engine.advect(m,m.cellVolume(),flux,vertical,.001,scheme,values);
            EXPECT_NEAR(inventory(values,m),mass,1e-11);
            for (double value:values) { EXPECT_GE(value,-1e-12); EXPECT_LE(value,1+1e-12); }
        }
    }
}

TEST(MeshQualification, ProductionHorizontalDiffusionConvergesOnTriQuadFamilies) {
    nlohmann::json records=nlohmann::json::array();
    constexpr double D=.02, duration=.1;
    for (const std::string family:{"quads","triangles","mixed"}) {
        double previous=0;
        for (int n:{8,16,32}) {
            SCOPED_TRACE(family+"/"+std::to_string(n));
            Model model; ASSERT_NO_FATAL_FAILURE(model.build(corpus(n,family),D,.02/(n*n*D)));
            auto exact=[&](int64_t c,double time){ return 1+.25*std::cos(pi*model.horizontal.faceCentroidX()[c])*
                std::cos(pi*model.horizontal.faceCentroidY()[c])*std::exp(-2*pi*pi*D*time); };
            for (int64_t c=0;c<model.layered.columnCount();++c) for (int k=0;k<2;++k)
                model.state.speciesAt(model.tracer,model.layered.cell(c,k))=exact(c,0);
            const double mass=model.solver.totalSpeciesMass(model.tracer);
            std::string message; ASSERT_TRUE(model.solver.advance(duration,message)) << message;
            double error=0;
            for (int64_t c=0;c<model.layered.columnCount();++c) {
                const double delta=model.state.speciesAt(model.tracer,model.layered.cell(c,0))-exact(c,duration);
                error+=model.horizontal.faceArea()[c]*delta*delta;
            }
            error=std::sqrt(error/model.horizontal.totalArea());
            EXPECT_NEAR(model.solver.totalSpeciesMass(model.tracer),mass,1e-10);
            if (previous>0) EXPECT_LT(error,.8*previous);
            records.push_back({{"family",family},{"resolution",n},{"l2_error",error},
                {"refinement_ratio",previous>0?error/previous:0}});
            previous=error;
        }
        EXPECT_LT(previous,2e-4);
    }
    std::ofstream(std::string(FVQual_TEST_OUTPUT_DIR)+"/mesh_diffusion_convergence.json") << records.dump(2);
}
