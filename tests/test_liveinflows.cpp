// SPDX-License-Identifier: MIT
#include "fvqual/component/fvqualcomponent.h"
#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
using namespace HydroCouple;
using namespace HydroCouple::SDK;
using FVQual::Component::FVQualComponent;
namespace {
using Status=IModelComponent::ComponentStatus;
using Json=nlohmann::json;
class PacketOutput : public Output1DString {
public:
    using Output1DString::Output1DString;
    bool refuse=false;
    bool getValuesInto(const BufferDescriptor &b,std::span<const int64_t> s,std::span<const int64_t> n,std::string *m) const override {
        if(refuse) {if(m)*m="donor read failed";return false;}
        return Output1DString::getValuesInto(b,s,n,m);
    }
};
class PacketDonor : public AbstractModelComponent, public HydroCouple::Temporal::ITimeModelComponent {
public:
    PacketDonor() : AbstractModelComponent("donor"),definition("packet",typeid(std::string)),packet("outfall_transfers",&dimension,1,&definition,this) {}
    std::vector<std::string> validate() override {return {};}
    void prepare() override {setStatus(Status::Preparing);setStatus(Status::Updated);}
    void update(const std::vector<IOutput *> & ={}) override {}
    void finish() override {setStatus(Status::Finishing);setStatus(Status::Finished);}
    HydroCouple::Temporal::IDateTime *currentDateTime() const override {return &clock;}
    HydroCouple::Temporal::ITimeSpan *simulationPeriod() const override {return &period;}
    double nextDateTimeJulianDay() const override {return clock.julianDay()+.7/86400.;}
    void publish(const Json &j) {packet[0]=j.dump();clock.setJulianDay(period.julianDay()+j.at("end_s").get<double>()/86400.);}
    mutable SDK::Temporal::TimeData clock{"clock",2440587.5};
    mutable SDK::Temporal::TimeSpan period{"period",2440587.5,2./86400.};
    Dimension dimension{"packet"};ValueDefinition definition;PacketOutput packet;
protected:
    void createArguments() override {}
    bool initializeArguments(std::string &) override {return true;}
    void createInputs() override {}
    void createOutputs() override {addOutput(&packet);}
    void initializeFailureCleanUp() override {}
};
IArgument *arg(FVQualComponent &c,std::string_view id) {for(auto *a:c.arguments())if(a->id()==id)return a;return nullptr;}
IOutput *out(FVQualComponent &c,std::string_view id) {for(auto *o:c.outputs())if(o->id()==id)return o;return nullptr;}
void set(FVQualComponent &c,std::string_view id,const Json &j) {std::string m;ASSERT_TRUE(arg(c,id)->initialize(j.dump(),IArgument::ArgumentInputType::JSON,m))<<m;}
void configure(FVQualComponent &c) {set(c,"geometry_counts",{2,2});set(c,"geometry_values",{20.,10.,5.,10.,1.,2.});set(c,"conservative_species",{{"tracer",.002}});c.initialize();}
IInput *wire(FVQualComponent &c,PacketDonor &d) {
    IInput *i=nullptr;for(auto *p:c.inputs())if(p->id()=="integrated_inflows")i=p;
    EXPECT_NE(i,nullptr);std::string m;EXPECT_TRUE(i->canConsume(&d.packet,m))<<m;
    EXPECT_TRUE(i->setProvider(&d.packet));return i;
}
Json frame(int id,double start,double end) {return {{"id","native/"+std::to_string(id)},{"start_s",start},{"end_s",end},
    {"inflows",Json::array({{{"column",0},{"volume_m3",1.25},{"heat_j",1.25*1000*4186*18},
        {"placement","surface"},{"species_mass_kg",{{"tracer",.125}}}}})}};}
Json packet(Json j) {j["schema"]="hydrocouple.integrated-inflow/1.0";j["epoch_julian_day"]=2440587.5;j["ownership"]="already_debited";
    j["units"]={{"time","elapsed_s"},{"water","m3"},{"heat","J_relative_to_0_degC"},{"species","kg"}};return j;}
double scalar(FVQualComponent &c,std::string_view id) {double v=0;int64_t n=1,z=0;BufferDescriptor b;b.data=&v;b.kind=DataKind::Float64;b.rank=1;b.shape=&n;std::string m;
    EXPECT_TRUE(out(c,id)->getValuesInto(b,{&z,1},{&n,1},&m))<<m;return v;}
std::string receipt(FVQualComponent &c) {std::string v,m;int64_t n=1,z=0;BufferDescriptor b;b.data=&v;b.kind=DataKind::String;b.rank=1;b.shape=&n;b.itemSizeBytes=sizeof(std::string);
    EXPECT_TRUE(out(c,"accepted_transfer_id")->getValuesInto(b,{&z,1},{&n,1},&m))<<m;return v;}
}
TEST(LiveInflows, NonuniformFramesConserveAndExactlyMatchRecordedReplay) {
    const auto path=std::filesystem::path(FVQual_TEST_OUTPUT_DIR)/"live_inflows.json";
    auto first=frame(0,0,.7),second=frame(1,.7,2.);
    auto series=packet(first);series["schema"]="fvqual.inflow-series/1.0";series["frames"]={first,second};series["species"]={{"tracer",{{"units","kg/m3"},{"initial",.002}}}};
    std::ofstream(path)<<series.dump();
    FVQualComponent live("live"),recorded("recorded");configure(live);configure(recorded);
    set(recorded,"conservative_species",Json::object());
    set(recorded,"recorded_inflows",{{"@uri",path.string()}});recorded.initialize();
    PacketDonor donor;donor.initialize();donor.prepare();wire(live,donor);
    live.prepare();recorded.prepare();ASSERT_EQ(live.status(),Status::Updated);ASSERT_EQ(recorded.status(),Status::Updated);
    const auto water=scalar(live,"water_volume"),heat=scalar(live,"heat_inventory"),mass=scalar(live,"inventory_tracer");
    std::string token,m;
    for(const auto &f:{first,second}) {
        donor.publish(packet(f));live.update();recorded.update();ASSERT_NE(live.status(),Status::Failed);
        EXPECT_DOUBLE_EQ(live.currentTime(),f.at("end_s").get<double>());
        EXPECT_EQ(live.state().species,recorded.state().species);EXPECT_EQ(live.state().eta,recorded.state().eta);
        EXPECT_EQ(receipt(live),f.at("id").get<std::string>());
        if(f==first) {ASSERT_TRUE(live.saveState(token,m));live.update();EXPECT_DOUBLE_EQ(scalar(live,"accepted_water"),1.25);}
    }
    EXPECT_EQ(live.status(),Status::Done);
    EXPECT_NEAR(scalar(live,"water_volume")-water,2.5,1e-8);
    EXPECT_NEAR(scalar(live,"heat_inventory")-heat,2.5*1000*4186*18,1e-3);
    EXPECT_NEAR(scalar(live,"inventory_tracer")-mass,.25,1e-10);
    const auto final=live.state().species;
    ASSERT_TRUE(live.restoreState(token,m));EXPECT_EQ(receipt(live),first.at("id").get<std::string>());
    live.update();EXPECT_EQ(live.state().species,final);EXPECT_DOUBLE_EQ(scalar(live,"accepted_water"),2.5);
}
TEST(LiveInflows, InvalidPacketsAndReadFailuresLeaveTheReceiverUnchanged) {
    for(int fault=0;fault<10;++fault) {
        FVQualComponent c("invalid");configure(c);PacketDonor d;d.initialize();d.prepare();wire(c,d);c.prepare();
        const auto state=c.state().species,eta=c.state().eta;auto p=packet(frame(0,0,.7));
        if(fault==0)p["epoch_julian_day"]=2440588.5;
        if(fault==1)p["units"]["water"]="ft3";
        if(fault==2)p["ownership"]="withdraw_again";
        if(fault==3)p["inflows"][0]["column"]=2;
        if(fault==4)p["start_s"]=.1;
        if(fault==5)p["end_s"]=3.;
        if(fault==9)p["inflows"][0]["column"]=.5;
        if(fault==6)p["inflows"][0]["species_mass_kg"]=Json::object();
        d.publish(p);if(fault==7)d.clock.setJulianDay(2440587.5);if(fault==8)d.packet.refuse=true;
        c.update();EXPECT_EQ(c.status(),Status::Failed)<<fault;EXPECT_DOUBLE_EQ(c.currentTime(),0);
        EXPECT_EQ(c.state().species,state);EXPECT_EQ(c.state().eta,eta);EXPECT_TRUE(receipt(c).empty());
    }
}
TEST(LiveInflows, PreflightRejectsEpochAndHorizonMismatch) {
    for(bool epoch:{true,false}) {
        FVQualComponent c("preflight");configure(c);PacketDonor d;d.initialize();
        if(epoch)d.period.setJulianDay(2440588.5);else d.period.setDuration(3./86400.);
        wire(c,d);c.prepare();EXPECT_EQ(c.status(),Status::Failed);
    }
}
