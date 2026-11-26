// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  runcommand.cpp
 * \brief `fvqual run <composition>` — the standalone pathway.
 *
 * Runs a Composition Specification document whose components are FVQual
 * instances, recording outputs through the SDK writers and finishing with
 * a run manifest — so a mesh generator's or FVQual's own results are
 * consumable downstream via execution mode "open", and an @from-bound
 * composition runs headless exactly as it runs in a GUI host.
 *
 * The recipe is the SDK's own (examples/serial_coupling): ModelInitializer
 * applies the document (bindings and all), TimeSteppedWorkflow derives the
 * schedule, IOThread carries per-step snapshots to the writers, and
 * RunRecorder writes the manifest once the writers have finalized.
 */

#include "fvqual/component/fvqualcomponent.h"
#include "fvqual/component/fvqualcomponentinfo.h"
#include "fvqual/version.h"

#include "hydrocouplesdk/component/timesteppedworkflow.h"
#include "hydrocouplesdk/io/csvwriter.h"
#include "hydrocouplesdk/io/iothread.h"
#include "hydrocouplesdk/io/modelinitializer.h"
#include "hydrocouplesdk/io/runmanifest.h"
#include "hydrocouplesdk/io/snapshot.h"

#include <filesystem>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace {

using namespace HydroCouple;
using namespace HydroCouple::SDK;

} // namespace

int runCompositionCommand(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cerr << "usage: fvqual run <composition.json|.yaml>\n";
        return EXIT_FAILURE;
    }
    const std::filesystem::path documentPath = argv[2];
    const std::filesystem::path outputBase =
        documentPath.has_parent_path() ? documentPath.parent_path()
                                       : std::filesystem::path(".");

    // Every component id resolves to a fresh FVQual instance CARRYING THE
    // DOCUMENT'S id — the run manifest catalogs components by their own
    // id() and a later composition's "open" block looks them up by the
    // id it declared, so the two must be the same string. The info object
    // outlives them all, as the loading contract requires.
    FVQual::Component::FVQualComponentInfo info;
    std::map<std::string, std::unique_ptr<IModelComponent>> instances;
    IO::ModelInitializer initializer(
        [&](const std::string &id) -> IModelComponent * {
            auto found = instances.find(id);
            if (found == instances.end())
                found = instances
                            .emplace(id, std::make_unique<
                                             FVQual::Component::FVQualComponent>(
                                             id, &info))
                            .first;
            return found->second.get();
        });

    std::string message;
    if (!initializer.initializeFromFile(documentPath, message))
    {
        std::cerr << "composition failed: " << message << '\n';
        return EXIT_FAILURE;
    }
    const IO::CompositionSpec &spec = initializer.spec();
    if (spec.workflow.strategy == IO::WorkflowStrategy::PullDriven)
    {
        std::cerr << "fvqual run drives time-stepped compositions; "
                     "pull-driven ones need a pulling host\n";
        return EXIT_FAILURE;
    }

    // ── Writers behind the IO thread ─────────────────────────────────────
    IO::IOThread ioThread(/*capacity=*/32);
    IO::RunRecorder recorder;
    recorder.setCaption("fvqual run " +
                        std::string(FVQual_FULL_VERSION));
    recorder.setComposition(documentPath, spec.toJson(), true);

    for (const IO::WriterSpec &writer : spec.writers)
    {
        if (writer.type != "csv")
        {
            std::cerr << "writer type \"" << writer.type
                      << "\" is not available in fvqual run; csv is\n";
            return EXIT_FAILURE;
        }
        std::filesystem::path path = writer.path;
        if (path.is_relative())
            path = outputBase / path;
        std::filesystem::create_directories(path.parent_path());
        auto csv = std::make_shared<IO::CSVWriter>(path);
        csv->setTimeUnits("seconds_since_start");
        ioThread.addWriter(csv);
        recorder.addWriter(csv);
    }
    if (!ioThread.start())
    {
        std::cerr << "IO thread failed to start (code "
                  << ioThread.errorCode() << ")\n";
        return EXIT_FAILURE;
    }

    // ── Workflow ─────────────────────────────────────────────────────────
    TimeSteppedWorkflow workflow("fvqual-run");
    std::vector<IModelComponent *> components;
    for (const IO::ComponentSpec &block : spec.components)
    {
        IModelComponent *component = initializer.component(block.id);
        if (component == nullptr)
            continue;
        components.push_back(component);
        recorder.addComponent(component);
        if (!workflow.addModelComponent(component, nullptr, &message))
        {
            std::cerr << "workflow wiring failed: " << message << '\n';
            return EXIT_FAILURE;
        }
    }

    workflow.initialize();
    const std::vector<std::string> problems = workflow.validate();
    if (!problems.empty())
    {
        for (const std::string &problem : problems)
            std::cerr << "invalid composition: " << problem << '\n';
        return EXIT_FAILURE;
    }
    workflow.prepare();

    // ── Run, snapshotting every capturable output every step ─────────────
    using WorkflowStatus = IWorkflowComponent::WorkflowStatus;
    int step = 0;
    const int maximumSteps =
        spec.workflow.maxSteps > 0 ? spec.workflow.maxSteps : 1000000;
    while (workflow.status() != WorkflowStatus::Done &&
           workflow.status() != WorkflowStatus::Failed &&
           step < maximumSteps)
    {
        workflow.update();
        ++step;
        if (workflow.status() == WorkflowStatus::Failed)
        {
            for (const ErrorEntry &entry : workflow.errors())
                std::cerr << "workflow error: " << entry.message << '\n';
            return EXIT_FAILURE;
        }
        double elapsed = -1.0;
        for (IModelComponent *component : components)
        {
            auto *fvqual = dynamic_cast<FVQual::Component::FVQualComponent *>(component);
            if (!fvqual)
            {
                std::cerr << "fvqual run requires live FVQual components for time-stamped output\n";
                return EXIT_FAILURE;
            }
            const double current = fvqual->currentTime();
            if (elapsed >= 0.0 && std::abs(elapsed - current) > 1e-9)
            {
                std::cerr << "components reached different times; output requires a shared time coordinate\n";
                return EXIT_FAILURE;
            }
            elapsed = current;
        }
        IO::Snapshot snapshot(elapsed);
        bool captured = false;
        for (IModelComponent *component : components)
            for (IOutput *output : component->outputs())
            {
                if (output->dataKind()==DataKind::String || output->dataKind()==DataKind::Opaque) continue;
                std::string ignored;
                if (!snapshot.capture(*output, ignored) || !snapshot.captureGeometry(*output, ignored)) {
                    std::cerr << "recording output: " << ignored << '\n'; return EXIT_FAILURE;
                }
                captured = true;
            }
        if (captured && !ioThread.enqueue(std::move(snapshot)))
        {
            std::cerr << "writer error (code " << ioThread.errorCode()
                      << ")\n";
            return EXIT_FAILURE;
        }
    }

    if (workflow.status() != WorkflowStatus::Done)
    {
        std::cerr << "composition did not finish within its maximum steps\n";
        return EXIT_FAILURE;
    }
    workflow.finish();
    ioThread.finish();
    if (ioThread.errorCode() != 0)
    {
        std::cerr << "IO finished with error code " << ioThread.errorCode()
                  << '\n';
        return EXIT_FAILURE;
    }

    // ── The manifest, once the writer catalogs are final ─────────────────
    std::filesystem::path manifestPath =
        spec.runManifest.empty() ? std::filesystem::path("run.json")
                                 : std::filesystem::path(spec.runManifest);
    if (manifestPath.is_relative())
        manifestPath = outputBase / manifestPath;
    std::filesystem::create_directories(manifestPath.parent_path());
    if (!recorder.finishToFile(manifestPath, message,
                               IO::RunStatus::Completed))
    {
        std::cerr << "run manifest failed: " << message << '\n';
        return EXIT_FAILURE;
    }

    std::cout << "completed " << step << " steps; manifest at "
              << manifestPath.string() << '\n';
    return EXIT_SUCCESS;
}
