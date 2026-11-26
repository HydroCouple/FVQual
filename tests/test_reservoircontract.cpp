// SPDX-License-Identifier: MIT
#include "fvqual/component/fvqualcomponent.h"
#include "hydrocouplesdk/io/csvwriter.h"
#include "hydrocouplesdk/io/runmanifest.h"
#include "hydrocouplesdk/component/resultsmodelcomponent.h"
#include "hydrocouplesdk/io/snapshot.h"
#include <gtest/gtest.h>
#include "hydrocouplespatiotemporal.h"
#include <filesystem>
#include <fstream>
#include <limits>
using namespace HydroCouple;
using namespace HydroCouple::SDK;
using namespace HydroCouple::SDK::IO;
using FVQual::Component::FVQualComponent;
namespace {
IArgument *argument(FVQualComponent &c,const std::string &id) {
    for (auto *a:c.arguments()) if(a->id()==id) return a;
    return nullptr;
}
void assign(FVQualComponent &c,const std::string &id,const nlohmann::json &j) {
    auto *a=argument(c,id); ASSERT_NE(a,nullptr); std::string message;
    ASSERT_TRUE(a->initialize(j.dump(),IArgument::ArgumentInputType::JSON,message))<<message;
}
IOutput *output(IModelComponent &c,const std::string &id) {
    for(auto *o:c.outputs()) if(o->id()==id) return o;
    return nullptr;
}
std::vector<double> values(const IComponentDataItem &item,std::vector<int64_t> start={},std::vector<int64_t> count={}) {
    if(count.empty())count=item.shape();if(start.empty())start.assign(count.size(),0);
    int64_t n=1;for(auto size:count)n*=size;std::vector<double> v(size_t(n),0.0);
    BufferDescriptor b;b.data=v.data();b.kind=DataKind::Float64;b.rank=int(count.size());b.shape=count.data();b.space=MemorySpace::Host;
    std::string m;EXPECT_TRUE(item.getValuesInto(b,start,count,&m))<<m;return v;
}
void configure(FVQualComponent &c) {
    assign(c,"geometry_counts",{2,2});assign(c,"geometry_values",{20.,10.,5.,10.,1.,2.});
    assign(c,"sigma_interfaces",{0.,.2,1.});
    assign(c,"initial_surface_elevation_m",{1.,1.});
    assign(c,"initial_temperature_c",{22.,12.,24.,14.});
    assign(c,"conservative_species",{{"tracer",.002}});
}
}
TEST(ReservoirContract, NamedOverridesProfilesAndInventoriesMatchSolver) {
    FVQualComponent c("reservoir");configure(c);
    assign(c,"duration_s",{4.});assign(c,"coupling_interval_s",{2.});
    c.initialize();ASSERT_EQ(c.status(),IModelComponent::ComponentStatus::Initialized);
    c.prepare();ASSERT_NE(c.status(),IModelComponent::ComponentStatus::Failed);
    EXPECT_EQ(values(*output(c,"temperature")),(std::vector<double>{22,12,24,14}));
    EXPECT_EQ(values(*output(c,"eta")),(std::vector<double>{1,1}));
    EXPECT_EQ(values(*output(c,"water_depth")),(std::vector<double>{6,6}));
    EXPECT_NEAR(c.mesh().z(0,1),-.2,1e-12);
    EXPECT_NEAR(values(*output(c,"water_volume"))[0],1200,1e-10);
    EXPECT_NEAR(values(*output(c,"inventory_tracer"))[0],2.4,1e-12);
    auto *q=dynamic_cast<IQuantity *>(output(c,"tracer")->valueDefinition());ASSERT_NE(q,nullptr);
    EXPECT_DOUBLE_EQ(q->unit()->dimensions()->power(IUnitDimensions::FundamentalUnitDimension::Mass),1);
    EXPECT_DOUBLE_EQ(q->unit()->dimensions()->power(IUnitDimensions::FundamentalUnitDimension::Length),-3);
    c.update();EXPECT_NEAR(c.currentTime(),2,1e-12);c.update();EXPECT_EQ(c.status(),IModelComponent::ComponentStatus::Done);
    EXPECT_NEAR(values(*output(c,"water_volume"))[0],1200,1e-8);
    EXPECT_NEAR(values(*output(c,"inventory_tracer"))[0],2.4,1e-10);
}
TEST(ReservoirContract, RejectsIncompatibleProfilesSigmaAndSpecies) {
    FVQualComponent badProfile("p");configure(badProfile);assign(badProfile,"initial_temperature_c",{1,2,3});
    badProfile.initialize();EXPECT_EQ(badProfile.status(),IModelComponent::ComponentStatus::Failed);
    FVQualComponent badSigma("s");configure(badSigma);assign(badSigma,"layer_count",{3});
    badSigma.initialize();EXPECT_EQ(badSigma.status(),IModelComponent::ComponentStatus::Failed);
    FVQualComponent badSpecies("b");configure(badSpecies);assign(badSpecies,"conservative_species",{{"temperature",1.}});
    badSpecies.initialize();EXPECT_EQ(badSpecies.status(),IModelComponent::ComponentStatus::Failed);
    FVQualComponent badSurface("d");configure(badSurface);assign(badSurface,"initial_surface_elevation_m",{-5.});
    badSurface.initialize();badSurface.prepare();EXPECT_EQ(badSurface.status(),IModelComponent::ComponentStatus::Failed);
    FVQualComponent nonfinite("n");auto *legacy=dynamic_cast<Argument1DDouble *>(argument(nonfinite,"geometry_values"));
    (*legacy)[0]=std::numeric_limits<double>::infinity();nonfinite.initialize();EXPECT_EQ(nonfinite.status(),IModelComponent::ComponentStatus::Failed);
}
TEST(ReservoirContract, RecordedCsvReopensWithGeometryUnitsLayersAndEveryCell) {
    const auto dir=std::filesystem::path(FVQual_TEST_OUTPUT_DIR)/"reservoir_contract";
    std::filesystem::create_directories(dir);const auto path=dir/"run.json";
    std::vector<double> last;
    {
        FVQualComponent c("reservoir");configure(c);c.initialize();c.prepare();
        ASSERT_NE(c.status(),IModelComponent::ComponentStatus::Failed);
        auto writer=std::make_shared<CSVWriter>(dir/"fields.csv");ASSERT_EQ(writer->initialize(),0);
        RunRecorder recorder("reservoir-contract");recorder.addComponent(&c);recorder.addWriter(writer);
        for (int step=0;step<2;++step) {
            c.update();ASSERT_NE(c.status(),IModelComponent::ComponentStatus::Failed);
            Snapshot snapshot(c.currentDateTime()->julianDay());std::string message;
            for (auto *o:c.outputs()) {if(o->dataKind()==DataKind::String)continue;ASSERT_TRUE(snapshot.capture(*o,message))<<message;ASSERT_TRUE(snapshot.captureGeometry(*o,message))<<message;}
            ASSERT_EQ(writer->write(snapshot),0);last=values(*output(c,"temperature"));
        }
        c.finish();ASSERT_EQ(writer->finalize(),0);std::string message;
        ASSERT_TRUE(recorder.finishToFile(path,message))<<message;
    }
    std::string message;auto c=ResultsModelComponent::open(path,"reservoir",message);ASSERT_NE(c,nullptr)<<message;
    c->initialize();ASSERT_EQ(c->status(),IModelComponent::ComponentStatus::Finished);
    auto *t=output(*c,"temperature");ASSERT_NE(t,nullptr);EXPECT_EQ(t->shape(),(std::vector<int64_t>{2,2,2}));
    auto *mesh=dynamic_cast<HydroCouple::Spatial::IPolyhedralSurfaceComponentDataItem *>(t);ASSERT_NE(mesh,nullptr);
    EXPECT_EQ(mesh->entityDimension(),t->dimensions()[1]);EXPECT_EQ(mesh->polyhedralSurface()->meshView()->faceCount(),2);
    ASSERT_NE(mesh->polyhedralSurface()->spatialReferenceSystem(),nullptr);
    EXPECT_NE(mesh->polyhedralSurface()->spatialReferenceSystem()->verticalSrText().find("synthetic"),std::string::npos);
    auto *layered=dynamic_cast<HydroCouple::SpatioTemporal::ITimeLayeredMeshComponentDataItem *>(t);ASSERT_NE(layered,nullptr);
    EXPECT_EQ(layered->layerDimension(),t->dimensions()[2]);EXPECT_EQ(layered->verticalCoordinate()->layerCount(),2);
    EXPECT_EQ(values(*t,{1,0,0},{1,2,2}),last);
    auto *q=dynamic_cast<IQuantity *>(t->valueDefinition());ASSERT_NE(q,nullptr);
    EXPECT_DOUBLE_EQ(q->unit()->offsetToSI(),273.15);
    EXPECT_DOUBLE_EQ(q->unit()->dimensions()->power(IUnitDimensions::FundamentalUnitDimension::Temperature),1);
    RunManifest manifest;ASSERT_TRUE(RunManifest::read(path,manifest,message));
    EXPECT_EQ(manifest.results.size(),59u); // spatial fields, inventories, receipts and independent process ledger (transfer IDs remain control outputs)
    EXPECT_TRUE(manifest.entry("reservoir","temperature")->geometry.contains("vertical"));
    auto g=manifest.entry("reservoir","temperature")->geometry.at("artifact").get<std::string>();
    const auto geometry=dir/g;const auto backup=dir/(g+".backup");std::filesystem::rename(geometry,backup);
    auto missing=ResultsModelComponent::open(path,"reservoir",message);ASSERT_NE(missing,nullptr);missing->initialize();
    EXPECT_EQ(missing->status(),IModelComponent::ComponentStatus::Failed);
    std::filesystem::rename(backup,geometry);
}

TEST(ReservoirContract, TransferBudgetsCheckpointAndRecordedElevationsFollowTheSelectedTime) {
    const auto dir=std::filesystem::path(FVQual_TEST_OUTPUT_DIR)/"reservoir_transfers";std::filesystem::create_directories(dir);
    const auto source=dir/"inflows.json";
    nlohmann::json data={{"schema","fvqual.inflow-series/1.0"},
        {"units",{{"time","elapsed_s"},{"water","m3"},{"heat","J_relative_to_0_degC"},{"species","kg"}}},
        {"species",{{"tracer",{{"units","kg/m3"},{"initial",.002}}}}}, {"frames",nlohmann::json::array()}};
    for(int i=0;i<2;++i) data["frames"].push_back({{"id",std::to_string(i)},{"start_s",i},{"end_s",i+1},
        {"inflows",nlohmann::json::array({{{"column",0},{"volume_m3",1.},{"heat_j",41860000.},
            {"species_mass_kg",{{"tracer",.1}}}}})}});
    std::ofstream(source)<<data.dump();
    FVQualComponent c("transfers");assign(c,"geometry_counts",{2,2});assign(c,"geometry_values",{20.,10.,5.,10.,1.,2.});
    assign(c,"recorded_inflows",{{"@uri",source.string()}});c.initialize();c.prepare();ASSERT_NE(c.status(),IModelComponent::ComponentStatus::Failed);
    const double initialWater=values(*output(c,"water_volume"))[0],initialHeat=values(*output(c,"heat_inventory"))[0];
    auto writer=std::make_shared<CSVWriter>(dir/"fields.csv");writer->setTimeUnits("elapsed_s");ASSERT_EQ(writer->initialize(),0);
    RunRecorder recorder("transfers");recorder.addComponent(&c);recorder.addWriter(writer);
    std::string message;std::vector<double> elevations;
    c.update();ASSERT_NE(c.status(),IModelComponent::ComponentStatus::Failed);
    EXPECT_NEAR(values(*output(c,"water_volume"))[0]-initialWater,1.,1e-9);
    EXPECT_DOUBLE_EQ(values(*output(c,"accepted_water"))[0],1.);
    EXPECT_NEAR(values(*output(c,"heat_inventory"))[0]-initialHeat,41860000.,1e-3);
    EXPECT_DOUBLE_EQ(values(*output(c,"accepted_tracer"))[0],.1);
    std::string token;ASSERT_TRUE(c.saveState(token,message));
    for(int step=0;step<2;++step) {
        Snapshot snapshot(c.currentTime());for(auto *o:c.outputs()){if(o->dataKind()==DataKind::String)continue;ASSERT_TRUE(snapshot.capture(*o,message));ASSERT_TRUE(snapshot.captureGeometry(*o,message));}
        ASSERT_EQ(writer->write(snapshot),0);elevations.push_back(c.mesh().z(0,0));if(step==0)c.update();
    }
    ASSERT_TRUE(c.restoreState(token,message));EXPECT_DOUBLE_EQ(values(*output(c,"accepted_water"))[0],1.);
    EXPECT_DOUBLE_EQ(values(*output(c,"accepted_tracer"))[0],.1);c.update();EXPECT_DOUBLE_EQ(values(*output(c,"accepted_water"))[0],2.);
    c.finish();ASSERT_EQ(writer->finalize(),0);ASSERT_TRUE(recorder.finishToFile(dir/"run.json",message))<<message;
    auto reopened=ResultsModelComponent::open(dir/"run.json","transfers",message);ASSERT_NE(reopened,nullptr);reopened->initialize();
    ASSERT_EQ(reopened->status(),IModelComponent::ComponentStatus::Finished);
    auto *t=output(*reopened,"temperature");auto *layer=dynamic_cast<HydroCouple::Spatial::ILayering *>(t);ASSERT_NE(layer,nullptr);
    ASSERT_GT(std::abs(elevations[0]-elevations[1]),1e-6);
    for(int step=0;step<2;++step) {
        values(*t,{step,0,0},{1,2,2});EXPECT_DOUBLE_EQ(layer->verticalCoordinate()->interfaceElevation(0,0),elevations[step]);
    }
    auto *times=dynamic_cast<HydroCouple::Temporal::ITimeSeriesComponentDataItem *>(t);ASSERT_NE(times,nullptr);
    EXPECT_NEAR(times->times()[0],2440587.5+1./86400,1e-9);
}

TEST(ReservoirContract, RelativeKineticsDeckPublishesDeclaredConcentrationAndKgInventory) {
    const auto dir=std::filesystem::path(FVQual_TEST_OUTPUT_DIR)/"reservoir_units";std::filesystem::create_directories(dir);
    std::ofstream(dir/"species.rxn")<<"[species]\nX mg/L 2.0\n";
    FVQualComponent c("units");c.setReferenceDirectory(dir.string());
    assign(c,"geometry_counts",{2,2});assign(c,"geometry_values",{20.,10.,5.,10.,1.,2.});
    assign(c,"kinetics_file",{"species.rxn"});c.initialize();c.prepare();ASSERT_NE(c.status(),IModelComponent::ComponentStatus::Failed);
    auto *q=dynamic_cast<IQuantity *>(output(c,"X")->valueDefinition());ASSERT_NE(q,nullptr);
    EXPECT_DOUBLE_EQ(q->unit()->conversionFactorToSI(),.001);
    EXPECT_DOUBLE_EQ(q->unit()->dimensions()->power(IUnitDimensions::FundamentalUnitDimension::Length),-3);
    EXPECT_NEAR(values(*output(c,"inventory_X"))[0],2.,1e-12);c.update();
    EXPECT_NEAR(values(*output(c,"inventory_X"))[0],2.,1e-10);
}
