// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
//
// Gates for the forcing arguments: a whole-series met argument drives air
// temperature when no exchange input is wired (and enables meteorology by
// itself), and a kinetics deck registers its species and actually reacts.

#include "fvqual/component/fvqualcomponent.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <span>
#include <string>
#include <vector>

using FVQual::Component::FVQualComponent;
using ComponentStatus = HydroCouple::IModelComponent::ComponentStatus;

namespace {

//! Runs the default basin to Done and returns the final temperature block.
std::vector<double> runToDone(FVQualComponent &component)
{
    component.initialize();
    EXPECT_EQ(component.status(), ComponentStatus::Initialized);
    EXPECT_TRUE(component.validate().empty());
    component.prepare();
    EXPECT_NE(component.status(), ComponentStatus::Failed);
    int steps = 0;
    while (component.status() != ComponentStatus::Done &&
           component.status() != ComponentStatus::Failed && steps < 100)
    {
        component.update();
        ++steps;
    }
    EXPECT_EQ(component.status(), ComponentStatus::Done);
    const auto block =
        component.state().speciesBlock(component.temperatureIndex());
    return {block.begin(), block.end()};
}

HydroCouple::SDK::Temporal::TimeSeriesArgumentDouble *metArgumentOf(
    FVQualComponent &component)
{
    for (HydroCouple::IArgument *argument : component.arguments())
    {
        if (argument->id() == "met_air_temperature")
            return dynamic_cast<
                HydroCouple::SDK::Temporal::TimeSeriesArgumentDouble *>(
                argument);
    }
    return nullptr;
}

void setKineticsFile(FVQualComponent &component, const std::string &path)
{
    for (HydroCouple::IArgument *argument : component.arguments())
    {
        if (argument->id() != "kinetics_file")
            continue;
        std::string message;
        ASSERT_TRUE(argument->initialize(
            nlohmann::json{{"values", {path}}}.dump(),
            HydroCouple::IArgument::ArgumentInputType::JSON, message))
            << message;
    }
}

} // namespace

TEST(FVQualForcingArguments, AMetSeriesCoolsTheWaterWhenNoInputIsWired)
{
    FVQualComponent still("fvqual.nomet");
    const std::vector<double> unforced = runToDone(still);
    ASSERT_FALSE(unforced.empty());
    // No forcing: the uniform initial temperature just sits there.
    EXPECT_NEAR(unforced[0], 10.0, 1.0e-9);

    // COLD air, deliberately: MeteorologyField's default air temperature
    // is 20 degC, warmer than the 10-degree water — so "it warmed" cannot
    // tell an applied series from an enabled-but-unfed default. Cooling
    // can only come from the series actually being applied.
    FVQualComponent chilled("fvqual.met");
    auto *met = metArgumentOf(chilled);
    ASSERT_NE(met, nullptr);
    std::string message;
    ASSERT_TRUE(met->initialize(
        R"({"times": [0.0, 10.0], "values": [-5.0, -5.0]})",
        HydroCouple::IArgument::ArgumentInputType::JSON, message))
        << message;
    const std::vector<double> forced = runToDone(chilled);
    ASSERT_EQ(forced.size(), unforced.size());

    const double coldest = *std::min_element(forced.begin(), forced.end());
    EXPECT_LT(coldest, unforced[0] - 1.0e-6);
}

TEST(FVQualForcingArguments, AKineticsDeckRegistersItsSpeciesAndReacts)
{
    const std::string deck =
        (std::filesystem::path(FVQual_TEST_OUTPUT_DIR) / "decay.rxn")
            .string();
    ASSERT_TRUE(std::filesystem::exists(deck)) << deck;

    FVQualComponent component("fvqual.kinetics");
    setKineticsFile(component, deck);
    component.initialize();
    ASSERT_EQ(component.status(), ComponentStatus::Initialized);
    component.prepare();
    ASSERT_NE(component.status(), ComponentStatus::Failed);

    // Temperature plus the deck's X.
    ASSERT_EQ(component.registry().count(), 2);

    int steps = 0;
    while (component.status() != ComponentStatus::Done &&
           component.status() != ComponentStatus::Failed && steps < 100)
    {
        component.update();
        ++steps;
    }
    ASSERT_EQ(component.status(), ComponentStatus::Done);

    // decay.rxn: X starts at 7 g/m3 with dX/dt = -k X, k = 3/day, run for
    // one day: 7 e^-3 ~ 0.349. A configure() that installed nothing would
    // leave 7 everywhere.
    const auto block = component.state().speciesBlock(1);
    ASSERT_FALSE(block.empty());
    EXPECT_NEAR(block[0], 7.0 * std::exp(-3.0), 0.05);
}

TEST(FVQualForcingArguments, AMalformedKineticsDeckFailsByName)
{
    const std::string bad =
        (std::filesystem::path(FVQual_TEST_OUTPUT_DIR) / "malformed.rxn")
            .string();
    {
        std::ofstream file(bad);
        file << "[SPECIES]\n   X\n"; // a species missing units and value
    }

    FVQualComponent component("fvqual.badkinetics");
    setKineticsFile(component, bad);
    component.initialize();
    EXPECT_NE(component.status(), ComponentStatus::Initialized);
    std::string joined;
    for (const auto &entry : component.errors())
        joined += entry.message + "\n";
    EXPECT_NE(joined.find("malformed.rxn"), std::string::npos) << joined;
}
