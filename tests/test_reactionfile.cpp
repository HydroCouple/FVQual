// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_reactionfile.cpp
 * \brief Verification of text-configured reaction networks (F6.3).
 *
 * The claim a configuration file makes is equivalence: a model written in text
 * must behave exactly like the same model written in code. So the central gate
 * here is not that the file parses — it is that a file expressing first-order
 * decay produces bit-for-bit the same concentrations as the template that
 * expresses it in C++.
 *
 * The other half is refusal. A configuration file is the one part of a model a
 * user edits by hand, so every malformed case has to be rejected with a line
 * number rather than half-applied or silently ignored. Files here are written
 * under the test output directory so a failed run leaves something to look at.
 */

#include "fvqual/hydro/internalmode.h"
#include "fvqual/kinetics/reactionfile.h"
#include "fvqual/kinetics/templates.h"
#include "fvqual/mesh/horizontalmesh.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/mesh/sigmaspec.h"
#include "fvqual/state/speciesregistry.h"
#include "fvqual/state/statestore.h"
#include "testmeshes.h"

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <sstream>
#include <string>

using FVQual::Hydro::InternalModeOptions;
using FVQual::Hydro::InternalModeSolver;
using FVQual::Kinetics::kSecondsPerDay;
using FVQual::Kinetics::ReactionFile;
using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;
using FVQual::State::SpeciesKind;
using FVQual::State::SpeciesRegistry;
using FVQual::State::StateStore;
namespace Testing = FVQual::Testing;

namespace {

/*!
 * \brief Parses \p text, expecting success.
 */
bool parseText(ReactionFile &file, const std::string &text, std::string &message)
{
    std::istringstream input(text);

    return file.parse(input, "<test>", message);
}

/*!
 * \brief Writes \p text to the reviewable test output directory and returns it.
 */
std::string writeFile(const std::string &name, const std::string &text)
{
    const std::string path = std::string(FVQual_TEST_OUTPUT_DIR) + "/" + name;

    std::ofstream out(path);
    out << text;
    out.close();

    return path;
}

/*!
 * \brief A still basin, as elsewhere: rate laws are the only thing moving.
 */
struct Basin
{
    HorizontalMesh     horizontal;
    LayeredMesh        layered;
    SpeciesRegistry    registry;
    StateStore         state;
    InternalModeSolver solver;

    void build(int layers, double totalDepth)
    {
        std::string message;

        ASSERT_TRUE(horizontal.build(Testing::channelStrip(3, 300.0, 40.0), message))
            << message;

        for (int64_t f = 0; f < horizontal.faceCount(); ++f)
            horizontal.faceBed()[static_cast<size_t>(f)] = -totalDepth;

        SigmaSpec spec;
        spec.layerCount        = layers;
        spec.minLayerThickness = 0.0;

        ASSERT_TRUE(layered.build(horizontal, spec, 0.0, message)) << message;

        state.allocate(layered, registry);

        InternalModeOptions options;
        options.gravity               = 9.80665;
        options.manningN              = 0.0;
        options.cflNumber             = 0.5;
        options.surfaceHeatFlux       = 0.0;
        options.neutralViscosity      = 0.0;
        options.backgroundViscosity   = 0.0;
        options.backgroundDiffusivity = 0.0;
        options.convectiveMixing      = 0.0;

        ASSERT_TRUE(solver.initialize(layered, state, registry, options, message))
            << message;
    }

    void setUniform(int species, double value)
    {
        for (int64_t cell = 0; cell < layered.cellCount(); ++cell)
            state.speciesAt(species, cell) = value;
    }

    [[nodiscard]] double at(int species, int k) const
    {
        return state.speciesAt(species, layered.cell(0, k));
    }

    void run(double dt, int steps)
    {
        std::string message;

        for (int i = 0; i < steps; ++i)
            ASSERT_TRUE(solver.step(dt, message)) << message;
    }
};

} // namespace

// ── Parsing ───────────────────────────────────────────────────────────────

TEST(ReactionFile, ACompleteFileParses)
{
    const std::string text = R"(
; A comment, and a blank line above.

[species]
   DO       g/m3   9.0
   CBOD     g/m3   5.0
   TSS      g/m3   1.0
   TSS_BED  g/m2   0.0    benthic

[Constants]
   kd   0.23/86400        ; per day, converted in place

[RATES]
   DO     WATER       0 - kd * CBOD
   DO     SURFACE     2.3e-5 * (9.0 - DO)
   DO     bed         0 - 1.16e-5
   CBOD   WATER       0 - kd * CBOD

[SETTLING]
   TSS    1.0/86400       TSS_BED

[LIGHT]
   TSS    0.05
)";

    ReactionFile file;
    std::string  message;

    ASSERT_TRUE(parseText(file, text, message)) << message;

    EXPECT_EQ(file.speciesCount(), 4);
    EXPECT_EQ(file.constantCount(), 1);
    EXPECT_EQ(file.rateCount(), 4);

    // Section names, scope keywords and the BENTHIC qualifier are all
    // case-insensitive; species ids are not, because they must match a registry
    // that treats them as distinct.
    SpeciesRegistry registry;
    ASSERT_TRUE(file.registerSpecies(registry, message)) << message;

    ASSERT_EQ(registry.count(), 4);
    EXPECT_EQ(registry.definition(registry.indexOf("DO")).units, "g/m3");
    EXPECT_NEAR(registry.definition(registry.indexOf("DO")).initialValue, 9.0, 1.0e-12);

    const int bed = registry.indexOf("TSS_BED");
    ASSERT_GE(bed, 0);
    EXPECT_EQ(registry.definition(bed).kind, SpeciesKind::MsxBenthic);
    EXPECT_FALSE(registry.definition(bed).isTransported());
}

TEST(ReactionFile, NumericFieldsMayBeArithmetic)
{
    // The one place the file undoes the per-day trap without guessing: a
    // coefficient stays readable as the division a modeller actually means.
    const std::string text = R"(
[CONSTANTS]
   a    0.23/86400
   b    2.662037037037037e-6
   c    1/24/3600
)";

    ReactionFile file;
    std::string  message;

    ASSERT_TRUE(parseText(file, text, message)) << message;
    EXPECT_EQ(file.constantCount(), 3);

    // Verified through a model, since the values are private: a rate built from
    // `a` must decay at 0.23/day.
    const std::string model = R"(
[SPECIES]
   X    g/m3   4.0
[CONSTANTS]
   k    0.23/86400
[RATES]
   X    WATER   0 - k * X
)";

    ReactionFile configured;
    ASSERT_TRUE(parseText(configured, model, message)) << message;

    Basin basin;
    ASSERT_TRUE(configured.registerSpecies(basin.registry, message)) << message;
    basin.build(3, 9.0);

    const int x = basin.solver.speciesIndex("X");
    ASSERT_GE(x, 0);
    basin.setUniform(x, 4.0);

    ASSERT_TRUE(configured.configure(basin.solver, message)) << message;
    ASSERT_TRUE(basin.solver.compileReactions(message)) << message;

    basin.run(864.0, 10); // a tenth of a day

    EXPECT_NEAR(basin.at(x, 1), 4.0 * std::exp(-0.23 * 0.1), 1.0e-8);
}

// ── Equivalence with the code path ────────────────────────────────────────

TEST(ReactionFile, AFileAndATemplateProduceTheSameModel)
{
    // The claim that matters. Text configuration is only worth having if it is
    // the same engine underneath, so the two paths are run side by side and
    // required to agree to round-off — not approximately, and not to a
    // tolerance chosen to make them agree.
    constexpr double perDay = 3.0;
    constexpr double c0     = 7.0;
    constexpr double dt     = 864.0;
    constexpr int    steps  = 10;

    // The code path.
    double fromTemplate = 0.0;
    {
        FVQual::Kinetics::FirstOrderDecayTemplate decay;
        decay.id               = "X";
        decay.decayRatePerDay  = perDay;
        decay.temperatureTheta = 1.0;

        Basin basin;
        std::string message;

        ASSERT_TRUE(decay.registerSpecies(basin.registry, message)) << message;
        basin.build(4, 12.0);

        const int x = basin.solver.speciesIndex("X");
        basin.setUniform(x, c0);

        ASSERT_TRUE(decay.configure(basin.solver, message)) << message;
        ASSERT_TRUE(basin.solver.compileReactions(message)) << message;

        basin.run(dt, steps);

        fromTemplate = basin.at(x, 2);
    }

    // The same model, written out.
    const std::string text = R"(
[SPECIES]
   X    g/m3   7.0
[CONSTANTS]
   k    3.0/86400
[RATES]
   X    WATER   0 - k * X
)";

    const std::string path = writeFile("decay.rxn", text);

    double fromFile = 0.0;
    {
        ReactionFile file;
        std::string  message;

        ASSERT_TRUE(file.read(path, message)) << message;

        Basin basin;
        ASSERT_TRUE(file.registerSpecies(basin.registry, message)) << message;
        basin.build(4, 12.0);

        const int x = basin.solver.speciesIndex("X");
        basin.setUniform(x, c0);

        ASSERT_TRUE(file.configure(basin.solver, message)) << message;
        ASSERT_TRUE(basin.solver.compileReactions(message)) << message;

        basin.run(dt, steps);

        fromFile = basin.at(x, 2);
    }

    EXPECT_EQ(fromFile, fromTemplate)
        << "file " << fromFile << " against template " << fromTemplate;

    // And both are the right answer, so agreeing on a wrong one would not pass.
    EXPECT_NEAR(fromFile, c0 * std::exp(-perDay * 0.1), 1.0e-8);
}

TEST(ReactionFile, ScopesAndSettlingSurviveTheRoundTrip)
{
    // Everything the templates need the format to express, in one model:
    // scoped laws at both ends of the column, settling with a benthic receiver,
    // and self-shading.
    const std::string text = R"(
[SPECIES]
   DO       g/m3   8.0
   ALGAE    g/m3   2.0
   ALGAE_BED g/m2  0.0   BENTHIC

[RATES]
   DO   SURFACE   2.0/86400 * (9.0 - DO)
   DO   BED       0 - 0.5/86400

[SETTLING]
   ALGAE   0.5/86400     ALGAE_BED

[LIGHT]
   ALGAE   0.02
)";

    const std::string path = writeFile("scoped.rxn", text);

    ReactionFile file;
    std::string  message;

    ASSERT_TRUE(file.read(path, message)) << message;

    Basin basin;
    ASSERT_TRUE(file.registerSpecies(basin.registry, message)) << message;
    basin.build(5, 15.0);

    const int oxygen = basin.solver.speciesIndex("DO");
    const int algae  = basin.solver.speciesIndex("ALGAE");
    const int bed    = basin.solver.speciesIndex("ALGAE_BED");

    basin.setUniform(oxygen, 8.0);
    basin.setUniform(algae, 2.0);

    ASSERT_TRUE(file.configure(basin.solver, message)) << message;
    ASSERT_TRUE(basin.solver.compileReactions(message)) << message;

    EXPECT_NEAR(basin.solver.settling().velocity(algae), 0.5 / kSecondsPerDay, 1.0e-18);
    EXPECT_EQ(basin.solver.settling().benthicReceiver(algae), bed);
    EXPECT_NEAR(basin.solver.light().specificAttenuation(algae), 0.02, 1.0e-15);

    basin.run(864.0, 10);

    // Reaeration raised the surface, sediment demand lowered the bed, and the
    // interior — which has neither law — did not move.
    EXPECT_GT(basin.at(oxygen, 0), 8.0);
    EXPECT_LT(basin.at(oxygen, 4), 8.0);
    EXPECT_NEAR(basin.at(oxygen, 2), 8.0, 1.0e-9);

    EXPECT_GT(basin.solver.settling().deliveredToBed(algae), 0.0);
}

// ── Refusal ───────────────────────────────────────────────────────────────

TEST(ReactionFile, MalformedInputIsRefusedWithALineNumber)
{
    struct Case
    {
        const char *text;
        const char *line;   //!< Line the failure should be reported on.
        const char *needle; //!< Something the message must mention.
    };

    const Case cases[] = {
        {"[SPECIES]\n   X g/m3 1.0\n[NOPE]\n", "3", "NOPE"},
        {"[SPECIES\n", "1", "closed"},
        {"   X g/m3 1.0\n", "1", "before any section"},
        {"[SPECIES]\n   X g/m3\n", "2", "initial value"},
        {"[SPECIES]\n   X g/m3 1.0\n   X g/m3 2.0\n", "3", "twice"},
        {"[SPECIES]\n   X g/m3 nonsense\n", "2", "nonsense"},
        {"[CONSTANTS]\n   k\n", "2", "name and a value"},
        {"[CONSTANTS]\n   k 1/\n", "2", "expression"},
        {"[CONSTANTS]\n   k other*2\n", "2", "variables"},
        {"[RATES]\n   X SIDEWAYS 0 - X\n", "2", "scope"},
        {"[RATES]\n   X\n", "2", "scope"},
        {"[RATES]\n   X 0 - X\n", "2", "scope"},
        {"[SETTLING]\n   X\n", "2", "velocity"},
        {"[LIGHT]\n   X\n", "2", "attenuation"},
    };

    for (const Case &c : cases)
    {
        ReactionFile file;
        std::string  message;

        EXPECT_FALSE(parseText(file, c.text, message))
            << "accepted malformed input: " << c.text;

        EXPECT_NE(message.find(std::string("<test>:") + c.line), std::string::npos)
            << "wrong or missing line number in: " << message;

        EXPECT_NE(message.find(c.needle), std::string::npos)
            << "message did not explain the problem: " << message;
    }
}

TEST(ReactionFile, AMalformedFileConfiguresNothing)
{
    // Parsing and applying are separate so that a bad file cannot leave a model
    // half-built. Verified by the state of the network afterwards, not by the
    // return value alone.
    const std::string text = R"(
[SPECIES]
   X    g/m3   1.0
[RATES]
   X    0 - 0.1 * X
   X    SIDEWAYS   0 - 0.2 * X
)";

    ReactionFile file;
    std::string  message;

    EXPECT_FALSE(parseText(file, text, message));

    EXPECT_EQ(file.speciesCount(), 0) << "a rejected file left species behind";
    EXPECT_EQ(file.rateCount(), 0) << "a rejected file left rate laws behind";
}

TEST(ReactionFile, UnknownReferencesAreRefusedAtConfiguration)
{
    // The parser cannot know whether a species exists — the registry does. So
    // these are caught when the file is applied, and they must be caught rather
    // than quietly skipped.
    const std::string text = R"(
[SPECIES]
   X    g/m3   1.0
[SETTLING]
   GHOST   1.0e-5
)";

    ReactionFile file;
    std::string  message;

    ASSERT_TRUE(parseText(file, text, message)) << message;

    Basin basin;
    ASSERT_TRUE(file.registerSpecies(basin.registry, message)) << message;
    basin.build(3, 9.0);

    EXPECT_FALSE(file.configure(basin.solver, message));
    EXPECT_NE(message.find("GHOST"), std::string::npos) << message;
}

TEST(ReactionFile, AMissingFileIsReported)
{
    ReactionFile file;
    std::string  message;

    EXPECT_FALSE(file.read(std::string(FVQual_TEST_OUTPUT_DIR) + "/absent.rxn", message));
    EXPECT_NE(message.find("absent.rxn"), std::string::npos) << message;
}
