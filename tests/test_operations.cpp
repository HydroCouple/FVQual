// SPDX-License-Identifier: MIT
#include "fvqual/component/fvqualcomponent.h"
#include "fvqual/capi/fvqual_configuration.h"
#include "fvqual/capi/fvqual_coupling.h"
#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
using J = nlohmann::json;
using FVQual::Component::FVQualComponent;
using AT = HydroCouple::IArgument::ArgumentInputType;
using CS = HydroCouple::IModelComponent::ComponentStatus;
namespace {
HydroCouple::IArgument *argument(FVQualComponent &c, const std::string &id) { for (auto *a : c.arguments()) if (a->id() == id) return a; return nullptr; }
void set(FVQualComponent &c, const std::string &id, const J &value) { std::string message; ASSERT_TRUE(argument(c,id)->initialize(value.dump(), AT::JSON, message)) << message; }
J wrap(const std::string &section, const J &value) { return {{"schema","fvqual.operations/1.0"},{"section",section},{"value",value}}; }
double ledger(FVQualComponent &c, const std::string &name) { return FVQual::Forcing::processBudgetValues(c.solver(), c.registry()).at(name); }
void setup(FVQualComponent &c) { set(c,"geometry_counts",{1,1}); set(c,"geometry_values",{10.,10.,5.,10.,1.,2.}); set(c,"conservative_species",{{"tracer",.002}}); }
J inAndOut() { return J::array({{{"id","in"},{"column",0},{"discharge_m3_s",1.},{"temperature_c",20.},{"concentrations",{{"tracer",.01}}},
    {"schedule",J::array({{{"time_s",0},{"discharge_m3_s",1.}},{{"time_s",.4},{"discharge_m3_s",2.}}})}},
    {{"id","out"},{"column",0},{"discharge_m3_s",-1.},{"placement","uniform"}}}); }
}
TEST(ReservoirOperations, CollocatedInletWithdrawalAndIrregularKnotsConserveIndependentLoads) {
    FVQualComponent c("operations"); setup(c); set(c,"sources",wrap("sources",inAndOut()));
    c.initialize(); ASSERT_EQ(c.status(),CS::Initialized); c.prepare(); ASSERT_NE(c.status(),CS::Failed);
    c.update(); ASSERT_NE(c.status(),CS::Failed);
    // Inlet integrates 1*.4 + 2*.6 = 1.6; outlet integrates exactly 1.
    EXPECT_NEAR(ledger(c,"budget_water_source_in"),1.6,1e-12);
    EXPECT_NEAR(ledger(c,"budget_water_source_out"),1,1e-12);
    EXPECT_NEAR(c.mesh().totalVolume(),500.6,1e-9);
    EXPECT_NEAR(ledger(c,"budget_tracer_source_in"),.016,1e-14);
    // .4 s first withdrawal at .002, .6 s at concentration after first load.
    const double expectedOut = .4*.002 + .6*(1 + .4*(.01-.002))/500;
    EXPECT_NEAR(ledger(c,"budget_tracer_source_out"),expectedOut,1e-13);
    EXPECT_NEAR(c.state().totalMass(1,c.mesh()),1+.016-expectedOut,1e-12);
    EXPECT_NEAR(ledger(c,"budget_water_residual"),0,1e-9);
    EXPECT_NEAR(ledger(c,"budget_tracer_residual"),0,1e-12);
    EXPECT_NEAR(ledger(c,"budget_heat_residual"),0,1e-3);
    std::string token,message; ASSERT_TRUE(c.saveState(token,message))<<message;
    c.update(); const auto final = FVQual::Forcing::processBudgetValues(c.solver(),c.registry());
    ASSERT_TRUE(c.restoreState(token,message))<<message; c.update();
    EXPECT_EQ(FVQual::Forcing::processBudgetValues(c.solver(),c.registry()),final);
}
TEST(ReservoirOperations, InvalidEditIsAtomicAndUnknownKeysAreRefused) {
    FVQualComponent c("operations"); setup(c); c.initialize(); ASSERT_EQ(c.status(),CS::Initialized); auto *a=argument(c,"sources"); std::string before,after,message;
    ASSERT_TRUE(a->serialize(AT::JSON,before,message));
    for (const auto &rows : {J::array({{{"id","bad"},{"column",99},{"discharge_m3_s",1}}}),
        J::array({{{"id","bad"},{"column",0},{"discharge_m3_s",1},{"typo",2}}}),
        J::array({{{"id","bad"},{"column",0},{"discharge_m3_s",1},{"schedule",J::array({{{"time_s",1},{"discharge_m3_s",1}}})}}})}) {
        EXPECT_FALSE(a->initialize(wrap("sources",rows).dump(),AT::JSON,message)); EXPECT_FALSE(message.empty());
        ASSERT_TRUE(a->serialize(AT::JSON,after,message)); EXPECT_EQ(after,before);
    }
}
TEST(ReservoirOperations, PrescribedHeatHasAnIndependentEnergyOracle) {
    FVQual::Mesh::HorizontalMesh horizontal; std::string message;
    HydroCouple::SDK::IO::MeshDefinition definition;
    definition.nodeX={0,10,10,0}; definition.nodeY={0,0,10,10}; definition.nodeZ={-5,-5,-5,-5};
    definition.faceNodeOffsets={0,4}; definition.faceNodes={0,1,2,3};
    ASSERT_TRUE(horizontal.build(definition,message)) << message;
    FVQual::Mesh::LayeredMesh mesh; FVQual::Mesh::SigmaSpec sigma; sigma.layerCount=1;
    ASSERT_TRUE(mesh.build(horizontal,sigma,std::vector<double>{0},message));
    FVQual::State::SpeciesRegistry registry; ASSERT_GE(registry.addTemperature(10,message),0);
    FVQual::State::StateStore state; state.allocate(mesh,registry);
    FVQual::Hydro::InternalModeSolver solver; FVQual::Hydro::InternalModeOptions options; options.surfaceHeatFlux=100;
    ASSERT_TRUE(solver.initialize(mesh,state,registry,options,message)); ASSERT_TRUE(solver.advance(10,message));
    EXPECT_NEAR(solver.budget().surfaceHeat[5],100000,1e-9);
    EXPECT_NEAR((solver.totalSpeciesMass(0)-5000)*1000*4186,100000,1e-3);
    EXPECT_NEAR(solver.scalarBudgetResidual(0)*1000*4186,0,1e-3);
}
TEST(ReservoirOperations, CApiAndComponentShareConfigurationAndBudgetAndRollback) {
    FVQualComponent component("shared"); setup(component); set(component,"sources",wrap("sources",inAndOut()));
    component.initialize(); ASSERT_EQ(component.status(),CS::Initialized); component.prepare(); ASSERT_NE(component.status(),CS::Failed); component.update();
    FVQ_Mesh *mesh=nullptr; FVQ_Simulation *sim=nullptr; FVQ_SigmaSpec sigma; fvq_sigma_default(&sigma); sigma.layerCount=1; FVQ_Options options; fvq_options_default(&options);
    ASSERT_EQ(fvq_mesh_create_channel(1,10,10,-5,&mesh),FVQ_OK);
    ASSERT_EQ(fvq_sim_create(mesh,&sigma,&options,&sim),FVQ_OK); fvq_mesh_destroy(mesh);
    ASSERT_EQ(fvq_sim_add_temperature(sim,10,nullptr),FVQ_OK); ASSERT_EQ(fvq_sim_add_species(sim,"tracer","kg/m3",.002,nullptr),FVQ_OK);
    const J config={{"sources",inAndOut()}};
    ASSERT_EQ(fvq_sim_configure_json(sim,config.dump().c_str()),FVQ_OK)<<fvq_last_error();
    const J invalid={{"sources",J::array({{{"id","bad"},{"column",2},{"discharge_m3_s",1}}})}};
    EXPECT_NE(fvq_sim_configure_json(sim,invalid.dump().c_str()),FVQ_OK);
    ASSERT_EQ(fvq_sim_initialize(sim),FVQ_OK)<<fvq_last_error();
    ASSERT_EQ(fvq_sim_advance(sim,1),FVQ_OK);
    // Named: a range-for over items() of a temporary json would iterate freed memory (C++20).
    const J budget=FVQual::Forcing::processBudgetValues(component.solver(),component.registry());
    for(const auto &[id,value]:budget.items()) {
        double amount=0;ASSERT_EQ(fvq_sim_budget_value(sim,id.c_str(),&amount),FVQ_OK)<<id; EXPECT_NEAR(amount,value.get<double>(),id.starts_with("budget_heat") ? 1e-3 : 1e-9)<<id;
    }
    uint64_t token;ASSERT_EQ(fvq_sim_save_state(sim,&token),FVQ_OK);ASSERT_EQ(fvq_sim_advance(sim,1),FVQ_OK);double before,after;
    ASSERT_EQ(fvq_sim_budget_value(sim,"budget_tracer_source_out",&before),FVQ_OK);ASSERT_EQ(fvq_sim_restore_state(sim,token),FVQ_OK);ASSERT_EQ(fvq_sim_advance(sim,1),FVQ_OK);ASSERT_EQ(fvq_sim_budget_value(sim,"budget_tracer_source_out",&after),FVQ_OK);EXPECT_DOUBLE_EQ(before,after);fvq_sim_destroy(sim);
}
TEST(ReservoirOperations, HeatTermsAndSettlingBedBasisCloseWithScheduledWeatherAndReactions) {
    FVQualComponent c("heat"); setup(c);
    set(c,"reaction_deck",wrap("reaction_deck",std::string("[species]\nreactive kg/m3 0.001\nparticles kg/m3 0.002\nsediment kg/m2 0 benthic\n[rates]\nreactive WATER 0 - 0.00001 * reactive\n[settling]\nparticles 0.001 sediment\n")));
    set(c,"meteorology",wrap("meteorology",{{"values",{{"air_temperature_c",19},{"wind_speed_m_s",0},{"shortwave_w_m2",100}}},{"records",J::array({{{"time_s",0},{"shortwave_w_m2",100}},{{"time_s",.3},{"shortwave_w_m2",200}}})}}));
    c.initialize();ASSERT_EQ(c.status(),CS::Initialized);c.prepare();ASSERT_NE(c.status(),CS::Failed);c.update();ASSERT_NE(c.status(),CS::Failed);
    EXPECT_DOUBLE_EQ(c.solver().meteorology().airTemperature.at(0),19);
    // Domain area = 100 m2. Default albedo=.06. Independently integrated solar input.
    EXPECT_NEAR(ledger(c,"budget_heat_shortwave"),100*.94*(100*.3+200*.7),1e-8);
    EXPECT_LT(ledger(c,"budget_heat_back_radiation"),0);EXPECT_LT(ledger(c,"budget_heat_latent"),0);
    EXPECT_NEAR(ledger(c,"budget_heat_residual"),0,1e-3);EXPECT_NEAR(ledger(c,"budget_water_residual"),0,1e-9);
    EXPECT_LT(ledger(c,"budget_reactive_reaction"),0);EXPECT_NEAR(ledger(c,"budget_reactive_residual"),0,1e-12);
    EXPECT_GT(ledger(c,"budget_sediment_settling"),0);
    EXPECT_NEAR(ledger(c,"budget_sediment_settling")+ledger(c,"budget_particles_settling"),0,1e-12);
    EXPECT_NEAR(ledger(c,"budget_sediment_residual"),0,1e-12);EXPECT_NEAR(ledger(c,"budget_particles_residual"),0,1e-12);
}
TEST(ReservoirOperations, StructuresGroupsBoundariesAndInternalTransfersHaveSeparateBudgets) {
    FVQualComponent c("all");set(c,"geometry_counts",{2,1});set(c,"geometry_values",{20.,10.,5.,10.,1.,2.});set(c,"conservative_species",{{"tracer",.002}});
    set(c,"sources",wrap("sources",J::array({{{"id","tributary"},{"column",0},{"discharge_m3_s",.01},{"temperature_c",20},{"concentrations",{{"tracer",.01}}}}})));
    set(c,"structures",wrap("structures",J::array({{{"id","port"},{"column",1},{"type","pump"},{"capacity_m3_s",-.001},{"placement","uniform"}}})));
    set(c,"outlet_groups",wrap("outlet_groups",J::array({{{"id","release"},{"ports",{"port"}},{"demand_m3_s",.0005},{"blend",false}}})));
    set(c,"internal_transfers",wrap("internal_transfers",J::array({{{"id","circulation"},{"from_column",0},{"to_column",1},{"discharge_m3_s",.003},{"from_placement","uniform"},{"to_placement","uniform"}}})));
    set(c,"open_boundaries",wrap("open_boundaries",J::array({{{"id","edge"},{"type","inflow"},{"node_pairs",J::array({{0,1}})},{"specific_discharge_m2_s",.001},{"temperature_c",20},{"concentrations",{{"tracer",.01}}}}})));
    c.initialize();ASSERT_EQ(c.status(),CS::Initialized);c.prepare();ASSERT_NE(c.status(),CS::Failed);c.update();ASSERT_NE(c.status(),CS::Failed);
    EXPECT_NEAR(ledger(c,"budget_water_source_in"),.01,1e-12);EXPECT_NEAR(ledger(c,"budget_water_structure_out"),.0005,1e-12);
    EXPECT_NEAR(ledger(c,"budget_water_transfer_in"),.003,1e-12);EXPECT_NEAR(ledger(c,"budget_water_transfer_out"),.003,1e-12);
    EXPECT_NEAR(ledger(c,"budget_tracer_transfer_in"),ledger(c,"budget_tracer_transfer_out"),1e-14);
    EXPECT_NEAR(ledger(c,"budget_water_boundary_in"),.01,1e-9);
    EXPECT_NEAR(ledger(c,"budget_water_residual"),0,1e-8);EXPECT_NEAR(ledger(c,"budget_heat_residual"),0,1e-2);EXPECT_NEAR(ledger(c,"budget_tracer_residual"),0,1e-10);
}
