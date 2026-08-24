// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  fvqual_bench.cpp
 * \brief Performance baseline for the internal mode (F7.2).
 *
 * Phase 7 exists to make the model fast. It cannot start from intuition about
 * which loops look expensive: a kernel that dominates the page can take two
 * percent of the run, and the one that actually costs is usually the one nobody
 * reads. So the first deliverable of the phase is a measurement, and the port
 * order follows from it.
 *
 * The scenario is a stratified basin stepped with everything switched on —
 * baroclinic pressure, mode splitting, transport, vertical mixing — because a
 * profile of a stripped-down configuration would answer a question nobody
 * asked. Sizes sweep so the *scaling* is visible too: a phase whose share grows
 * with the mesh is a different problem from one with a fixed overhead.
 */

#include "fvqual/hydro/internalmode.h"
#include "fvqual/mesh/horizontalmesh.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/mesh/sigmaspec.h"
#include "fvqual/state/speciesregistry.h"
#include "fvqual/state/statestore.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

using namespace FVQual;

/*!
 * \brief A rectangular basin of \p columnsX by \p columnsY cells.
 */
HydroCouple::SDK::IO::MeshDefinition basin(int columnsX, int columnsY, double length,
                                           double width)
{
    HydroCouple::SDK::IO::MeshDefinition mesh;
    mesh.meshName = "mesh2d";

    const double dx = length / columnsX;
    const double dy = width / columnsY;

    for (int j = 0; j <= columnsY; ++j)
    {
        for (int i = 0; i <= columnsX; ++i)
        {
            mesh.nodeX.push_back(i * dx);
            mesh.nodeY.push_back(j * dy);
        }
    }

    const auto node = [&](int i, int j) { return j * (columnsX + 1) + i; };

    mesh.faceNodeOffsets.push_back(0);

    for (int j = 0; j < columnsY; ++j)
    {
        for (int i = 0; i < columnsX; ++i)
        {
            mesh.faceNodes.push_back(node(i, j));
            mesh.faceNodes.push_back(node(i + 1, j));
            mesh.faceNodes.push_back(node(i + 1, j + 1));
            mesh.faceNodes.push_back(node(i, j + 1));
            mesh.faceNodeOffsets.push_back(
                static_cast<int64_t>(mesh.faceNodes.size()));
        }
    }

    return mesh;
}

struct Result
{
    int64_t columns    = 0;
    int64_t substeps   = 0;
    double  dt         = 0.0;
    int64_t cells      = 0;
    int64_t steps      = 0;
    double  seconds    = 0.0;
    Hydro::StepProfile profile;
};

Result run(int columnsX, int columnsY, int layers, int steps)
{
    constexpr double kLength = 20000.0;
    constexpr double kWidth  = 4000.0;
    constexpr double kDepth  = 40.0;

    Mesh::HorizontalMesh horizontal;
    std::string          message;

    if (!horizontal.build(basin(columnsX, columnsY, kLength, kWidth), message))
    {
        std::fprintf(stderr, "mesh: %s\n", message.c_str());
        std::exit(1);
    }

    for (int64_t f = 0; f < horizontal.faceCount(); ++f)
        horizontal.faceBed()[static_cast<size_t>(f)] = -kDepth;

    Mesh::SigmaSpec spec;
    spec.layerCount        = layers;
    spec.minLayerThickness = 0.0;

    Mesh::LayeredMesh layered;

    if (!layered.build(horizontal, spec, 0.0, message))
    {
        std::fprintf(stderr, "layers: %s\n", message.c_str());
        std::exit(1);
    }

    State::SpeciesRegistry registry;
    const int temperature = registry.addTemperature(20.0, message);

    State::SpeciesDefinition tracer;
    tracer.id           = "tracer";
    tracer.units        = "g/m3";
    tracer.initialValue = 1.0;
    registry.add(tracer, message);

    State::StateStore state;
    state.allocate(layered, registry);

    // Summer stratification, so the baroclinic terms are doing real work.
    for (int64_t c = 0; c < layered.columnCount(); ++c)
    {
        for (int k = 0; k < layers; ++k)
        {
            const int64_t cell = layered.cell(c, k);
            const double  z    = 0.5 * (layered.z(c, k) + layered.z(c, k + 1));

            state.speciesAt(temperature, cell) = 8.0 + 14.0 * std::exp(z / 8.0);
        }
    }

    Hydro::InternalModeOptions options;
    options.windStressX = 0.05;
    options.cflNumber   = 0.5;

    Hydro::InternalModeSolver solver;

    if (!solver.initialize(layered, state, registry, options, message))
    {
        std::fprintf(stderr, "solver: %s\n", message.c_str());
        std::exit(1);
    }

    // One untimed step so first-touch allocation and cold caches do not land
    // in the measurement.
    if (!solver.step(solver.stableTimeStep(), message))
        std::fprintf(stderr, "warmup: %s\n", message.c_str());

    solver.resetProfile();

    const auto start = std::chrono::steady_clock::now();

    for (int i = 0; i < steps; ++i)
    {
        if (!solver.step(solver.stableTimeStep(), message))
        {
            std::fprintf(stderr, "step: %s\n", message.c_str());
            std::exit(1);
        }
    }

    const auto finish = std::chrono::steady_clock::now();

    Result result;
    result.columns  = layered.columnCount();
    result.substeps = solver.diagnostics().externalSteps;
    result.dt       = solver.diagnostics().timeStep;
    result.cells   = layered.cellCount();
    result.steps   = steps;
    result.seconds = std::chrono::duration<double>(finish - start).count();
    result.profile = solver.profile();

    return result;
}

void report(const Result &r)
{
    const auto &p = r.profile;

    const double total = static_cast<double>(p.total);

    const auto share = [&](int64_t part) {
        return total > 0.0 ? 100.0 * static_cast<double>(part) / total : 0.0;
    };

    const double perStep = r.seconds / static_cast<double>(r.steps);

    std::printf("\n%lld cells, %lld steps, %.3f s  (%.3f ms/step, %.2f Mcell-steps/s)\n",
                static_cast<long long>(r.cells), static_cast<long long>(r.steps),
                r.seconds, 1000.0 * perStep,
                static_cast<double>(r.cells) * r.steps / r.seconds / 1.0e6);

    std::printf("  %lld columns, internal dt %.3f s, %lld external substeps per step\n",
                static_cast<long long>(r.columns), r.dt,
                static_cast<long long>(r.substeps));
    std::printf("  external cost per column-substep: %.1f ns\n",
                r.substeps > 0 ? static_cast<double>(p.external) / r.steps /
                                     (r.columns * r.substeps)
                               : 0.0);
    std::printf("  %-12s %8s %7s\n", "phase", "ms/step", "share");

    const struct
    {
        const char *name;
        int64_t     value;
    } phases[] = {
        {"density", p.density},   {"pressure", p.pressure}, {"forcing", p.forcing},
        {"external", p.external}, {"momentum", p.momentum}, {"geometry", p.geometry},
        {"fluxes", p.fluxes},     {"transport", p.transport},
    };

    int64_t accounted = 0;

    for (const auto &phase : phases)
    {
        accounted += phase.value;
        std::printf("  %-12s %8.4f %6.1f%%\n", phase.name,
                    static_cast<double>(phase.value) / r.steps / 1.0e6,
                    share(phase.value));
    }

    std::printf("  %-12s %8.4f %6.1f%%\n", "unaccounted",
                static_cast<double>(p.total - accounted) / r.steps / 1.0e6,
                share(p.total - accounted));
}

} // namespace

int main(int argc, char **argv)
{
    const int steps = argc > 1 ? std::atoi(argv[1]) : 20;

    std::printf("FVQual internal-mode baseline (%d steps per size)\n", steps);

    // A sweep, so the scaling of each phase is visible and not just its share
    // at one arbitrary size.
    report(run(40, 8, 10, steps));
    report(run(80, 16, 15, steps));
    report(run(160, 32, 20, steps));

    return 0;
}
