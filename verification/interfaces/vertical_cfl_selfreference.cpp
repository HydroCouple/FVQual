// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  vertical_cfl_selfreference.cpp
 * \brief Does the diagnosed vertical velocity depend on the step that produced
 *        it? (task 81 — measured, and the answer is no)
 *
 * **The suspicion, as recorded.** `diagnoseVerticalFlux(dt)` closes layer
 * continuity with `volumeRate = (V_now - V_old) / dt`, and `stableTimeStep()`
 * then reads the resulting `m_verticalFlux` to choose the *next* step. If
 * `V_now - V_old` were dominated by cancellation noise, dividing it by a
 * vanishing `dt` would inflate the diagnosed velocity, which would demand a
 * smaller step, which would inflate it further — a collapse that feeds itself.
 * The task list carried this as an open defect.
 *
 * **Measured, it does not happen.**
 *
 * | | at rest | tilted surface, 600 s of seiche |
 * |---|---|---|
 * | dt range probed | 100 s down to 1e-4 s | cap 10 s down to 0.25 s |
 * | max abs w | 1.2e-16 falling to 7.7e-39 m/s | 1.89e-6 to 1.78e-6 m/s |
 * | behaviour | falls with dt | flat to within **1.06x** |
 *
 * At rest the volume change is *exactly* zero rather than noisy — the stage does
 * not move, so the two volumes are bitwise equal — and there is nothing to
 * amplify. With real motion the ratio converges to a genuine rate: a fortyfold
 * change in step size moves it by six percent. The lag is real (the limit uses
 * the previous step's flux) but benign, which is the ordinary situation in an
 * explicit code.
 *
 * **Two things this file got wrong before it got them right**, recorded because
 * both are easy to repeat:
 *
 * 1. The first tilted experiment took a *single* step from rest. Nothing has
 *    accelerated after one step, so every value was round-off, and comparing
 *    them produced a 16000x "spread" that was noise over noise — and a verdict
 *    of "the self-reference is real" that was pure artefact. It now runs the
 *    seiche for 600 s, and refuses to interpret any value below 1e-6 m/s.
 * 2. The cap has to actually bind, or four runs take identical steps and the
 *    table compares one calculation with itself. The achieved internal step is
 *    printed next to each cap so that can be checked rather than assumed.
 *
 * **Not covered:** a flat bed throughout. The hypsographic case, where interface
 * area falls to a fraction of plan area, is the configuration that produced the
 * F82 blow-up; the evidence there is DeGray itself, which runs without the
 * vertical limit binding spuriously.
 *
 * Build (from the repository root, against a configured build tree):
 *
 *     c++ -std=c++20 -I include -I <build>/generated -I <sdk>/include \
 *         -I <hydrocouple>/include \
 *         verification/interfaces/vertical_cfl_selfreference.cpp \
 *         <build>/src/libfvqual-core.a -o vertical_cfl_selfreference
 */

#include "fvqual/hydro/internalmode.h"
#include "fvqual/mesh/horizontalmesh.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/mesh/sigmaspec.h"
#include "fvqual/state/speciesregistry.h"
#include "fvqual/state/statestore.h"

#include "hydrocouplesdk/io/meshdefinition.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

constexpr int    kColumns = 8;
constexpr int    kLayers  = 12;
constexpr double kDepth   = 30.0;
constexpr double kProbeDuration = 600.0; //!< Long enough for a seiche to move.

HydroCouple::SDK::IO::MeshDefinition strip(int columns, double length, double width)
{
    HydroCouple::SDK::IO::MeshDefinition mesh;
    mesh.meshName = "mesh2d";

    const double dx = length / static_cast<double>(columns);

    for (int j = 0; j <= 1; ++j)
        for (int i = 0; i <= columns; ++i)
        {
            mesh.nodeX.push_back(static_cast<double>(i) * dx);
            mesh.nodeY.push_back(static_cast<double>(j) * width);
        }

    mesh.faceNodeOffsets.push_back(0);

    const auto node = [&](int i, int j) { return j * (columns + 1) + i; };

    for (int i = 0; i < columns; ++i)
    {
        mesh.faceNodes.push_back(node(i, 0));
        mesh.faceNodes.push_back(node(i + 1, 0));
        mesh.faceNodes.push_back(node(i + 1, 1));
        mesh.faceNodes.push_back(node(i, 1));
        mesh.faceNodeOffsets.push_back(static_cast<int64_t>(mesh.faceNodes.size()));
    }

    return mesh;
}

//! One step at the given dt; returns the diagnosed maximum vertical velocity.
//! With \p tilt non-zero the surface starts sloped, so there is genuine motion
//! and a genuine volume change — the case where ΔV/dt has to converge to a real
//! rate rather than merely being zero.
double maxVerticalFluxAfterOneStep(double dt, double tilt = 0.0)
{
    std::string message;

    FVQual::Mesh::HorizontalMesh horizontal;

    if (!horizontal.build(strip(kColumns, 2000.0, 400.0), message))
    {
        std::printf("mesh build failed: %s\n", message.c_str());
        return -1.0;
    }

    for (int64_t f = 0; f < horizontal.faceCount(); ++f)
        horizontal.faceBed()[static_cast<size_t>(f)] = -kDepth;

    FVQual::Mesh::SigmaSpec spec;
    spec.layerCount = kLayers;

    FVQual::Mesh::LayeredMesh layered;

    std::vector<double> surface(static_cast<size_t>(horizontal.faceCount()), 0.0);

    for (size_t c = 0; c < surface.size(); ++c)
        surface[c] = tilt * (static_cast<double>(c) /
                                 static_cast<double>(surface.size() - 1) - 0.5);

    if (!layered.build(horizontal, spec, surface, message))
    {
        std::printf("layered build failed: %s\n", message.c_str());
        return -1.0;
    }

    FVQual::State::SpeciesRegistry registry;

    if (registry.addTemperature(10.0, message) < 0)
    {
        std::printf("registry failed: %s\n", message.c_str());
        return -1.0;
    }

    FVQual::State::StateStore state;
    state.allocate(layered, registry);

    FVQual::Hydro::InternalModeSolver  solver;
    FVQual::Hydro::InternalModeOptions options;

    // The floor must not clamp the very steps being probed, and with a tilt the
    // cap is the experiment's independent variable.
    options.minTimeStep = 1.0e-9;

    if (tilt != 0.0)
        options.maxTimeStep = dt;

    if (!solver.initialize(layered, state, registry, options, message))
    {
        std::printf("solver init failed: %s\n", message.c_str());
        return -1.0;
    }

    // At rest when tilt is zero: no forcing, no slope, uniform density, so any
    // diagnosed vertical velocity is numerical.
    if (tilt == 0.0)
    {
        if (!solver.step(dt, message))
        {
            std::printf("step(%g) failed: %s\n", dt, message.c_str());
            return -1.0;
        }

        return solver.diagnostics().maxVerticalFlux;
    }

    // With a tilt, one step is not an experiment: the flow starts from rest, so
    // after a single step nothing has accelerated and the diagnosed velocity is
    // still round-off. An earlier version of this file stopped there, compared
    // one round-off value against another, and reported a 16000x "spread" that
    // was noise over noise. Let the seiche actually develop, capping the step at
    // dt so the solver takes many of them, and report the largest genuine
    // vertical velocity over a fixed physical duration.
    double largest = 0.0;

    for (int i = 0; i < 60; ++i)
    {
        if (!solver.advance(kProbeDuration / 60.0, message))
        {
            std::printf("advance failed at dt cap %g: %s\n", dt, message.c_str());
            return -1.0;
        }

        largest = std::max(largest, solver.diagnostics().maxVerticalFlux);
    }

    // Report the step actually taken, so the experiment can be checked for
    // being an experiment: if the cap never binds, every run took the same
    // steps and the table below compares four identical calculations.
    std::printf("       (cap %8g s -> last internal step %10.4g s)\n", dt,
                solver.diagnostics().timeStep);

    return largest;
}

} // namespace

int main()
{
    std::printf("Vertical velocity diagnosed in a lake at rest, one step (task 81)\n\n");
    std::printf("  %12s  %18s  %14s\n", "dt (s)", "max |w| (m/s)", "w * dt (m)");
    std::printf("  %12s  %18s  %14s\n", "------", "-------------", "----------");

    const double steps[] = {100.0, 10.0, 1.0, 0.1, 0.01, 1.0e-3, 1.0e-4};

    double previousProduct = 0.0;
    bool   scalesLikeInverseDt = true;
    bool   any                 = false;

    for (double dt : steps)
    {
        const double w = maxVerticalFluxAfterOneStep(dt);

        if (w < 0.0)
            return 1;

        const double product = w * dt;
        std::printf("  %12g  %18.6e  %14.6e\n", dt, w, product);

        // If w ~ C/dt then w*dt is constant. If w is a genuine velocity it
        // stays put while w*dt falls with dt.
        if (any && previousProduct > 0.0 && product > 0.0 &&
            std::abs(product / previousProduct - 1.0) > 0.5)
            scalesLikeInverseDt = false;

        if (product > 0.0)
        {
            previousProduct = product;
            any             = true;
        }
    }

    std::printf("\n");

    if (!any)
    {
        std::printf("VERDICT: the diagnosed velocity is identically zero at every step.\n"
                    "         There is no 1/dt amplification to remove here.\n");
        return 0;
    }

    std::printf("At rest: %s\n\n",
                scalesLikeInverseDt
                    ? "w * dt is constant — the velocity scales like 1/dt, so the"
                      " cancellation noise IS being amplified as suspected."
                    : "w * dt falls with dt — no 1/dt amplification.");

    // The case that actually matters. A lake at rest has an exactly zero volume
    // change, so there is nothing for a 1/dt to amplify and the first table
    // cannot settle the question on its own. Tilt the surface and there is a
    // real volume change: ΔV/dt must now converge to a genuine rate as dt
    // shrinks. If instead it diverges, the diagnosed velocity is an artefact of
    // the step that produced it, and feeding it to the next step's Courant
    // limit is the self-reference the task suspected.
    std::printf("Vertical velocity with a tilted surface (real motion)\n\n");
    std::printf("  %12s  %18s\n", "dt (s)", "max |w| (m/s)");
    std::printf("  %12s  %18s\n", "------", "-------------");

    double smallest = 0.0;
    double largest  = 0.0;
    bool   first    = true;

    // Coarser caps than the at-rest table: each of these advances a full 600 s
    // of physical time, and a 1e-4 cap would be six million steps for no extra
    // information.
    const double caps[] = {10.0, 2.0, 1.0, 0.25};

    for (double dt : caps)
    {
        const double w = maxVerticalFluxAfterOneStep(dt, /*tilt=*/0.5);

        if (w < 0.0)
            return 1;

        std::printf("  %12g  %18.6e\n", dt, w);

        if (first)
        {
            smallest = largest = w;
            first    = false;
        }
        else
        {
            smallest = std::min(smallest, w);
            largest  = std::max(largest, w);
        }
    }

    const double spread = smallest > 0.0 ? largest / smallest : 0.0;

    // Refuse to read ratios of round-off. Anything below a micrometre per second
    // is not motion in a reservoir, and comparing two such numbers says nothing
    // about the code that produced them.
    constexpr double kMeaningful = 1.0e-6;

    std::printf("\nTilted: largest %.3e m/s, smallest %.3e m/s, spread %.2fx\n",
                largest, smallest, spread);

    if (largest < kMeaningful)
    {
        std::printf("VERDICT: INCONCLUSIVE — every value is at round-off (< %g m/s),\n"
                    "         so the spread is noise over noise and means nothing.\n",
                    kMeaningful);
        return 0;
    }

    std::printf("VERDICT: %s\n",
                spread > 10.0
                    ? "the diagnosed velocity depends strongly on dt — it is an\n"
                      "         artefact of the step, and the self-reference is real."
                    : "the diagnosed velocity is essentially independent of dt — it\n"
                      "         converges to a genuine rate, so there is no 1/dt\n"
                      "         amplification to remove.");
    return 0;
}
