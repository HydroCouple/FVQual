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
#include "hydrocouplesdk/io/runmanifest.h"
#include "hydrocouplesdk/io/meshpackage.h"
#include "hydrocouplesdk/differentiation/differentials.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>

namespace {

std::filesystem::path outputDirectory()
{
    return std::filesystem::path(FVQual_TEST_OUTPUT_DIR) / "standalone";
}

// Runs a command line through the shell. cmd.exe drops the first and last
// quote of a command that holds more than two, so on Windows the whole line
// goes in one more pair.
int run(std::string command)
{
#ifdef _WIN32
    command = "\"" + command + "\"";
#endif
    return std::system(command.c_str());
}

// The composition reads a mesh package, which needs NetCDF.
#ifdef FVQual_HAVE_NETCDF
TEST(FVQualRunCli, RecordedMixedMeshCompositionResolvesInputsAndRecordsElapsedTime)
{
    namespace IO = HydroCouple::SDK::IO;
    const auto directory = std::filesystem::path(FVQual_TEST_OUTPUT_DIR) / "recorded_composition";
    std::filesystem::remove_all(directory);
    std::filesystem::create_directories(directory);
    const auto example = std::filesystem::path(FVQUAL_FIXTURES_DIR).parent_path().parent_path() /
        "examples/swmm_reservoir";
    IO::MeshPackage package;
    IO::MeshImportOptions options;
    options.horizontalCrs = "LOCAL";
    options.verticalDatum = "synthetic-bed";
    options.elevationKind = "bed-elevation";
    std::string message;
    ASSERT_TRUE(IO::importSurfaceMesh(example / "reservoir.2dm", options, package, message)) << message;
    ASSERT_TRUE(IO::writeMeshPackage(directory / "reservoir.mesh.json", package, message)) << message;
    std::filesystem::copy_file(example / "inflow_series.json", directory / "inflow_series.json");
    std::filesystem::copy_file(example / "composition.json", directory / "composition.json");
    const std::string command = std::string("\"") + FVQUAL_CLI_PATH + "\" run \"" +
        (directory / "composition.json").string() + "\"";
    ASSERT_EQ(run(command), 0) << command;
    IO::RunManifest manifest;
    ASSERT_TRUE(IO::RunManifest::read(directory / "output/run.json", manifest, message)) << message;
    EXPECT_EQ(manifest.status, IO::RunStatus::Completed);
    ASSERT_EQ(manifest.components.size(), 1u);
    EXPECT_TRUE(manifest.components[0].errors.empty());
    ASSERT_GE(manifest.results.size(), 13u);
    ASSERT_TRUE(manifest.results[0].time.has_value());
    EXPECT_EQ(manifest.results[0].time->count, 4);
    EXPECT_DOUBLE_EQ(manifest.results[0].time->start, 1);
    EXPECT_DOUBLE_EQ(manifest.results[0].time->end, 4);
    EXPECT_EQ(manifest.results[0].time->units, "seconds_since_start");
}
#endif

} // namespace

TEST(FVQualRunCli, RelativeReactionDeckRunsAndReopensWithConservedConversion) {
    namespace IO=HydroCouple::SDK::IO;
    const auto dir=std::filesystem::path(FVQual_TEST_OUTPUT_DIR)/"reactive_composition";
    std::filesystem::remove_all(dir);std::filesystem::create_directories(dir);
    const auto example=std::filesystem::path(FVQUAL_FIXTURES_DIR).parent_path().parent_path()/"examples/reactive_transport";
    for(const char *file:{"composition.json","reaction.rxn"})std::filesystem::copy_file(example/file,dir/file);
    const std::string command=std::string("\"")+FVQUAL_CLI_PATH+"\" run \""+(dir/"composition.json").string()+"\"";
    ASSERT_EQ(run(command),0)<<command;std::string message;IO::RunManifest manifest;
    ASSERT_TRUE(IO::RunManifest::read(dir/"output/run.json",manifest,message))<<message;
    ASSERT_NE(manifest.entry("reservoir","A"),nullptr);EXPECT_EQ(manifest.entry("reservoir","A")->time->count,20);
    auto reopened=HydroCouple::SDK::ResultsModelComponent::open(dir/"output/run.json","reservoir",message);ASSERT_NE(reopened,nullptr)<<message;reopened->initialize();
    ASSERT_EQ(reopened->status(),HydroCouple::IModelComponent::ComponentStatus::Finished);
    double mass=0;int inventories=0;
    for(auto *o:reopened->outputs())if(o->id()=="inventory_A" || o->id()=="inventory_B") {
        const auto shape=o->shape();std::vector<double> data(HydroCouple::SDK::Differentiation::elementCount(shape));
        ASSERT_TRUE(o->getValuesInto(HydroCouple::SDK::Differentiation::hostBuffer(data,shape),std::vector<int64_t>(shape.size(),0),shape,&message))<<message;
        mass+=data.back();++inventories;
    }
    EXPECT_EQ(inventories,2);EXPECT_NEAR(mass,.25+.0216,1e-10);
}

TEST(FVQualRunCli, RunsACompositionAndItsManifestReopensAsAComponent)
{
    std::filesystem::remove_all(outputDirectory());

    const std::string command =
        std::string("\"") + FVQUAL_CLI_PATH + "\" run \"" +
        FVQUAL_COMPOSITIONS_DIR "/standalone_basin.json\"";
    ASSERT_EQ(run(command), 0) << command;

    // The artifacts are real: rows in the CSV, a manifest naming the run.
    std::ifstream csv(outputDirectory() / "results.csv");
    ASSERT_TRUE(csv.good());
    int lines = 0;
    double firstTime = -1.0, lastTime = -1.0;
    for (std::string line; std::getline(csv, line);)
    {
        if (lines > 0)
        {
            const double time = std::stod(line.substr(0, line.find(',')));
            if (firstTime < 0.0) firstTime = time;
            lastTime = time;
        }
        ++lines;
    }
    EXPECT_GT(lines, 6) << "a header plus six coupling steps of rows";
    EXPECT_DOUBLE_EQ(firstTime, 3600.0);
    EXPECT_DOUBLE_EQ(lastTime, 21600.0);

    // Pathway (b)'s point: the run reopens with no model binary in sight.
    std::string message;
    HydroCouple::SDK::IO::RunManifest manifest;
    ASSERT_TRUE(HydroCouple::SDK::IO::RunManifest::read(
        outputDirectory() / "run.json", manifest, message)) << message;
    ASSERT_FALSE(manifest.results.empty());
    ASSERT_TRUE(manifest.results[0].time.has_value());
    EXPECT_EQ(manifest.results[0].time->units, "seconds_since_start");
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
