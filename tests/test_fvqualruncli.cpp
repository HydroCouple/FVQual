// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
//
// The standalone pathway's gate: `fvqual run` executes a composition
// document as a subprocess (the property "the shipped binary runs a
// composition" is not one a library test can show), and the run it leaves
// behind is consumable — the manifest reopens as a ResultsModelComponent
// whose recorded temperature is readable, which is exactly what a later
// composition's execution mode "open" does.

#include <gtest/gtest.h>

#include "hydrocouplesdk/component/resultsmodelcomponent.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>

namespace {

std::filesystem::path outputDirectory()
{
    return std::filesystem::path(FVQual_TEST_OUTPUT_DIR) / "standalone";
}

} // namespace

TEST(FVQualRunCli, RunsACompositionAndItsManifestReopensAsAComponent)
{
    std::filesystem::remove_all(outputDirectory());

    const std::string command =
        std::string("\"") + FVQUAL_CLI_PATH + "\" run \"" +
        FVQUAL_COMPOSITIONS_DIR "/standalone_basin.json\"";
    ASSERT_EQ(std::system(command.c_str()), 0) << command;

    // The artifacts are real: rows in the CSV, a manifest naming the run.
    std::ifstream csv(outputDirectory() / "results.csv");
    ASSERT_TRUE(csv.good());
    int lines = 0;
    for (std::string line; std::getline(csv, line);)
        ++lines;
    EXPECT_GT(lines, 6) << "a header plus six coupling steps of rows";

    // Pathway (b)'s point: the run reopens with no model binary in sight.
    std::string message;
    auto reopened = HydroCouple::SDK::ResultsModelComponent::open(
        (outputDirectory() / "run.json").string(), "basin", message);
    ASSERT_NE(reopened, nullptr) << message;
    reopened->initialize();
    std::string errors;
    for (const auto &entry : reopened->errors())
        errors += entry.message + "\n";
    EXPECT_EQ(reopened->status(),
              HydroCouple::IModelComponent::ComponentStatus::Finished)
        << errors;
    ASSERT_FALSE(reopened->outputs().empty());
    EXPECT_EQ(reopened->outputs()[0]->id(), "temperature");
}
