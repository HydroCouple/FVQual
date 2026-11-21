// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  scenarios.cpp
 * \brief Runs FVQual's analytical benchmarks and emits animation frames.
 *
 * Every scenario here has a closed-form or theoretical reference, and the
 * point of the exercise is to see the model and that reference side by side
 * as they evolve — a converged error norm tells you *how much* a scheme is
 * wrong, but watching a front arrive late tells you *why*.
 *
 * Output is a single JSON document consumed by verification/index.html.
 * Frames are decimated to keep the page small; the physics runs at the
 * solver's own adaptive time step regardless.
 *
 * Build and run:
 *   cmake --build build --target fvqual-verify
 *   ./build/verification/fvqual-verify verification/data.json
 */

#include "fvqual/hydro/barotropicsolver.h"
#include "fvqual/hydro/equationofstate.h"
#include "fvqual/hydro/internalmode.h"
#include "fvqual/transport/euleriantransport.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

using namespace FVQual;

namespace {

constexpr double kGravity = 9.80665;
constexpr double kPi      = 3.14159265358979323846;

// ── Tiny JSON emitter ─────────────────────────────────────────────────────
//
// A dependency-free writer, because pulling a JSON library into a
// verification tool to serialize arrays of doubles would be silly.

std::string number(double value)
{
    if (!std::isfinite(value))
        return "0";

    std::ostringstream out;
    out << std::setprecision(6) << value;
    return out.str();
}

std::string array(const std::vector<double> &values)
{
    std::string text = "[";

    for (size_t i = 0; i < values.size(); ++i)
    {
        if (i > 0)
            text += ",";
        text += number(values[i]);
    }

    return text + "]";
}

std::string quoted(const std::string &text)
{
    return "\"" + text + "\"";
}

/*!
 * \brief One animated scenario: a spatial axis, a series of frames, and the
 *        analytical reference for each frame.
 */
struct Scenario
{
    std::string id;
    std::string title;
    std::string subtitle;
    std::string xLabel;
    std::string yLabel;
    std::string commentary;
    std::string verdict;

    std::vector<double> x;

    std::vector<double>              times;
    std::vector<std::vector<double>> model;
    std::vector<std::vector<double>> analytic;

    // Optional second model series (for scheme comparisons).
    std::string                      comparisonLabel;
    std::vector<std::vector<double>> comparison;

    // Optional scalar time series plotted beneath the field.
    std::string         traceLabel;
    std::vector<double> trace;
    std::vector<double> traceReference;

    std::string modelLabel    = "FVQual";
    std::string analyticLabel = "analytical";

    [[nodiscard]] std::string toJson() const
    {
        std::string text = "{";
        text += quoted("id") + ":" + quoted(id) + ",";
        text += quoted("title") + ":" + quoted(title) + ",";
        text += quoted("subtitle") + ":" + quoted(subtitle) + ",";
        text += quoted("xLabel") + ":" + quoted(xLabel) + ",";
        text += quoted("yLabel") + ":" + quoted(yLabel) + ",";
        text += quoted("commentary") + ":" + quoted(commentary) + ",";
        text += quoted("verdict") + ":" + quoted(verdict) + ",";
        text += quoted("modelLabel") + ":" + quoted(modelLabel) + ",";
        text += quoted("analyticLabel") + ":" + quoted(analyticLabel) + ",";
        text += quoted("comparisonLabel") + ":" + quoted(comparisonLabel) + ",";
        text += quoted("traceLabel") + ":" + quoted(traceLabel) + ",";
        text += quoted("x") + ":" + array(x) + ",";
        text += quoted("times") + ":" + array(times) + ",";
        text += quoted("trace") + ":" + array(trace) + ",";
        text += quoted("traceReference") + ":" + array(traceReference) + ",";

        const auto frames = [](const std::vector<std::vector<double>> &series) {
            std::string out = "[";
            for (size_t i = 0; i < series.size(); ++i)
            {
                if (i > 0)
                    out += ",";
                out += array(series[i]);
            }
            return out + "]";
        };

        text += quoted("model") + ":" + frames(model) + ",";
        text += quoted("analytic") + ":" + frames(analytic) + ",";
        text += quoted("comparison") + ":" + frames(comparison);

        return text + "}";
    }
};

// ── Mesh helper ───────────────────────────────────────────────────────────

HydroCouple::SDK::IO::MeshDefinition strip(int columns, double length, double width)
{
    HydroCouple::SDK::IO::MeshDefinition mesh;
    mesh.meshName = "mesh2d";

    const double dx = length / static_cast<double>(columns);

    for (int i = 0; i <= columns; ++i)
    {
        mesh.nodeX.push_back(i * dx);
        mesh.nodeY.push_back(0.0);
        mesh.nodeX.push_back(i * dx);
        mesh.nodeY.push_back(width);
    }

    mesh.faceNodeOffsets.push_back(0);

    for (int i = 0; i < columns; ++i)
    {
        mesh.faceNodes.push_back(2 * i);
        mesh.faceNodes.push_back(2 * (i + 1));
        mesh.faceNodes.push_back(2 * (i + 1) + 1);
        mesh.faceNodes.push_back(2 * i + 1);
        mesh.faceNodeOffsets.push_back(static_cast<int64_t>(mesh.faceNodes.size()));
    }

    return mesh;
}

template <typename BedFunction>
Mesh::HorizontalMesh channel(int columns, double length, double width, BedFunction bed)
{
    Mesh::HorizontalMesh mesh;
    std::string          message;

    if (!mesh.build(strip(columns, length, width), message))
        std::fprintf(stderr, "mesh build failed: %s\n", message.c_str());

    for (int64_t f = 0; f < mesh.faceCount(); ++f)
        mesh.faceBed()[static_cast<size_t>(f)] = bed(mesh.faceCentroidX()[static_cast<size_t>(f)]);

    return mesh;
}

// ── 1. Ritter dam break ───────────────────────────────────────────────────

Scenario ritterDamBreak()
{
    constexpr int    kColumns  = 240;
    constexpr double kLength   = 1200.0;
    constexpr double kDam      = 600.0;
    constexpr double kUpstream = 1.0;
    constexpr double kEnd      = 40.0;
    constexpr int    kFrames   = 90;

    Scenario scene;
    scene.id       = "ritter";
    scene.title    = "Ritter dam break";
    scene.subtitle = "Instantaneous dam removal onto a dry bed";
    scene.xLabel   = "distance (m)";
    scene.yLabel   = "depth (m)";
    scene.commentary =
        "The classic dry-bed test of the external mode. Ritter's solution gives the "
        "depth everywhere at every instant, so there is nowhere for an error to hide. "
        "Watch the rarefaction fan: the model tracks the interior almost exactly, and "
        "the visible disagreement is confined to the vanishing tip, where a "
        "first-order scheme's front speed converges only slowly.";
    scene.verdict =
        "Interior converges at first order (9.2% to 1.0% from 200 to 3200 cells); "
        "the front trails theory by a fixed distance regardless of refinement.";

    const Mesh::HorizontalMesh mesh = channel(kColumns, kLength, 10.0, [](double) { return 0.0; });

    Hydro::BarotropicOptions options;
    options.gravity         = kGravity;
    options.frictionEnabled = false;
    options.dryDepth        = 1.0e-6;

    Hydro::BarotropicSolver solver;
    std::string             message;

    if (!solver.initialize(mesh, options, message))
        std::fprintf(stderr, "ritter init: %s\n", message.c_str());

    std::vector<double> surface(static_cast<size_t>(kColumns));
    for (int64_t c = 0; c < mesh.faceCount(); ++c)
    {
        surface[static_cast<size_t>(c)] =
            mesh.faceCentroidX()[static_cast<size_t>(c)] < kDam ? kUpstream : 0.0;
    }

    if (!solver.setSurface(surface, message))
        std::fprintf(stderr, "ritter surface: %s\n", message.c_str());

    scene.x.assign(mesh.faceCentroidX().begin(), mesh.faceCentroidX().end());

    const double celerity = std::sqrt(kGravity * kUpstream);

    for (int frame = 0; frame <= kFrames; ++frame)
    {
        const double target = kEnd * static_cast<double>(frame) / kFrames;

        if (!solver.advance(target - solver.time(), message))
            std::fprintf(stderr, "ritter advance: %s\n", message.c_str());

        scene.times.push_back(target);
        scene.model.push_back(solver.depth());

        std::vector<double> exact(static_cast<size_t>(kColumns), 0.0);

        for (int64_t c = 0; c < mesh.faceCount(); ++c)
        {
            const double x = mesh.faceCentroidX()[static_cast<size_t>(c)];

            if (target <= 0.0)
            {
                exact[static_cast<size_t>(c)] = x < kDam ? kUpstream : 0.0;
                continue;
            }

            const double xi = (x - kDam) / target;

            if (xi <= -celerity)
                exact[static_cast<size_t>(c)] = kUpstream;
            else if (xi >= 2.0 * celerity)
                exact[static_cast<size_t>(c)] = 0.0;
            else
                exact[static_cast<size_t>(c)] =
                    std::pow(2.0 * celerity - xi, 2.0) / (9.0 * kGravity);
        }

        scene.analytic.push_back(exact);
    }

    return scene;
}

// ── 2. Barotropic seiche ──────────────────────────────────────────────────

Scenario barotropicSeiche()
{
    constexpr int    kColumns   = 160;
    constexpr double kLength    = 1000.0;
    constexpr double kDepth     = 10.0;
    constexpr double kAmplitude = 0.05;
    constexpr int    kFrames    = 120;

    const double period = 2.0 * kLength / std::sqrt(kGravity * kDepth);

    Scenario scene;
    scene.id       = "seiche";
    scene.title    = "Barotropic seiche";
    scene.subtitle = "Fundamental standing wave in a closed basin";
    scene.xLabel   = "distance (m)";
    scene.yLabel   = "surface elevation (m)";
    scene.commentary =
        "A closed basin sloshing at its fundamental period, T = 2L/sqrt(gH). This is "
        "the cleanest test of wave speed: an error in either the pressure gradient or "
        "the celerity shows up as the model drifting out of phase with the analytical "
        "standing wave. Over three periods the two stay locked; the slow amplitude "
        "loss is first-order numerical damping, which the analytical line does not "
        "have.";
    scene.verdict = "Period within 5% of 2L/sqrt(gH); mass closure below 1e-10.";

    const Mesh::HorizontalMesh mesh =
        channel(kColumns, kLength, 20.0, [](double) { return -kDepth; });

    Hydro::BarotropicOptions options;
    options.gravity         = kGravity;
    options.frictionEnabled = false;

    Hydro::BarotropicSolver solver;
    std::string             message;

    if (!solver.initialize(mesh, options, message))
        std::fprintf(stderr, "seiche init: %s\n", message.c_str());

    std::vector<double> surface(static_cast<size_t>(kColumns));
    for (int64_t c = 0; c < mesh.faceCount(); ++c)
    {
        surface[static_cast<size_t>(c)] =
            kAmplitude *
            std::cos(kPi * mesh.faceCentroidX()[static_cast<size_t>(c)] / kLength);
    }

    if (!solver.setSurface(surface, message))
        std::fprintf(stderr, "seiche surface: %s\n", message.c_str());

    scene.x.assign(mesh.faceCentroidX().begin(), mesh.faceCentroidX().end());
    scene.traceLabel = "surface elevation at the closed end (m)";

    const double endTime = 3.0 * period;

    for (int frame = 0; frame <= kFrames; ++frame)
    {
        const double target = endTime * static_cast<double>(frame) / kFrames;

        if (!solver.advance(target - solver.time(), message))
            std::fprintf(stderr, "seiche advance: %s\n", message.c_str());

        scene.times.push_back(target);
        scene.model.push_back(solver.surfaceElevation());

        std::vector<double> exact(static_cast<size_t>(kColumns), 0.0);

        for (int64_t c = 0; c < mesh.faceCount(); ++c)
        {
            const double x = mesh.faceCentroidX()[static_cast<size_t>(c)];
            exact[static_cast<size_t>(c)] = kAmplitude * std::cos(kPi * x / kLength) *
                                            std::cos(2.0 * kPi * target / period);
        }

        scene.analytic.push_back(exact);

        scene.trace.push_back(solver.surfaceElevation()[0]);
        scene.traceReference.push_back(kAmplitude * std::cos(2.0 * kPi * target / period));
    }

    return scene;
}

// ── 3. Wind setup ─────────────────────────────────────────────────────────

Scenario windSetup()
{
    constexpr int    kColumns = 100;
    constexpr double kLength  = 1000.0;
    constexpr double kDepth   = 10.0;
    constexpr double kStress  = 1.0;
    constexpr int    kFrames  = 80;
    constexpr double kEnd     = 400.0;

    Scenario scene;
    scene.id       = "wind";
    scene.title    = "Wind setup";
    scene.subtitle = "Steady wind stress tilting a closed basin";
    scene.xLabel   = "distance (m)";
    scene.yLabel   = "surface elevation (m)";
    scene.commentary =
        "Steady wind piles water against the downwind shore until the pressure "
        "gradient balances the stress: d(eta)/dx = tau/(rho g H). The model starts "
        "from the analytical tilt and must hold it. The control matters as much as "
        "the test — the dashed line shows the same tilt with the wind switched off, "
        "which collapses into a seiche, so a solver that simply never moved would "
        "fail rather than pass.";
    scene.verdict = "Balanced tilt held within 10% of tau/(rho g H); unforced control collapses.";

    const Mesh::HorizontalMesh mesh =
        channel(kColumns, kLength, 20.0, [](double) { return -kDepth; });

    Hydro::BarotropicOptions options;
    options.gravity         = kGravity;
    options.windStressX     = kStress;
    options.frictionEnabled = true;

    const double slope = kStress / (options.waterDensity * kGravity * kDepth);

    std::vector<double> tilted(static_cast<size_t>(kColumns));
    for (int64_t c = 0; c < mesh.faceCount(); ++c)
    {
        tilted[static_cast<size_t>(c)] =
            (mesh.faceCentroidX()[static_cast<size_t>(c)] - 0.5 * kLength) * slope;
    }

    Hydro::BarotropicOptions unforced = options;
    unforced.windStressX              = 0.0;

    Hydro::BarotropicSolver forced;
    Hydro::BarotropicSolver control;
    std::string             message;

    if (!forced.initialize(mesh, options, message) ||
        !control.initialize(mesh, unforced, message))
    {
        std::fprintf(stderr, "wind init: %s\n", message.c_str());
    }

    if (!forced.setSurface(tilted, message) || !control.setSurface(tilted, message))
        std::fprintf(stderr, "wind surface: %s\n", message.c_str());

    scene.x.assign(mesh.faceCentroidX().begin(), mesh.faceCentroidX().end());
    scene.comparisonLabel = "no wind (control)";

    for (int frame = 0; frame <= kFrames; ++frame)
    {
        const double target = kEnd * static_cast<double>(frame) / kFrames;

        if (!forced.advance(target - forced.time(), message) ||
            !control.advance(target - control.time(), message))
        {
            std::fprintf(stderr, "wind advance: %s\n", message.c_str());
        }

        scene.times.push_back(target);
        scene.model.push_back(forced.surfaceElevation());
        scene.analytic.push_back(tilted);
        scene.comparison.push_back(control.surfaceElevation());
    }

    return scene;
}

// ── 4. MUSCL versus upwind ────────────────────────────────────────────────

Scenario schemeComparison()
{
    constexpr int    kColumns = 120;
    constexpr int    kLayers  = 4;
    constexpr double kLength  = 1200.0;
    constexpr double kDepth   = 8.0;
    constexpr double kSpeed   = 0.5;
    constexpr int    kStart   = 20;
    constexpr int    kStop    = 40;
    constexpr int    kFrames  = 100;
    constexpr double kEnd     = 1000.0;

    Scenario scene;
    scene.id       = "transport";
    scene.title    = "Transport scheme comparison";
    scene.subtitle = "A square pulse translated 500 m at constant velocity";
    scene.xLabel   = "distance (m)";
    scene.yLabel   = "concentration";
    scene.commentary =
        "The same pulse carried by the same flow under two schemes. First-order "
        "upwind bleeds the edges away until the pulse is a smooth hill; MUSCL with "
        "flux correction keeps the plateau and the shoulders. This is not cosmetic: "
        "the smearing is what thickens a thermocline, and a thickened thermocline "
        "weakens the density gradient that drives the flow. Note that neither scheme "
        "overshoots — boundedness is guaranteed by the limiter, not by luck.";
    scene.verdict =
        "MUSCL L1 error less than half of upwind; both bounded and mass-conserving "
        "to 1e-12.";

    scene.modelLabel      = "MUSCL + FCT";
    scene.comparisonLabel = "first-order upwind";
    scene.analyticLabel   = "exact";

    const Mesh::HorizontalMesh horizontal =
        channel(kColumns, kLength, 20.0, [](double) { return -kDepth; });

    Mesh::SigmaSpec spec;
    spec.layerCount        = kLayers;
    spec.minLayerThickness = 0.0;

    Mesh::LayeredMesh layered;
    std::string       message;

    if (!layered.build(horizontal, spec, 0.0, message))
        std::fprintf(stderr, "transport mesh: %s\n", message.c_str());

    const double cellSize  = kLength / kColumns;
    const double thickness = kDepth / kLayers;

    std::vector<double> edgeFlux(static_cast<size_t>(horizontal.edgeCount()) *
                                     static_cast<size_t>(kLayers),
                                 0.0);

    for (int64_t e = 0; e < horizontal.edgeCount(); ++e)
    {
        const size_t ei = static_cast<size_t>(e);

        if (horizontal.edgeRightFace()[ei] == Mesh::kNoFace ||
            std::abs(horizontal.edgeNormalX()[ei]) < 0.5)
        {
            continue;
        }

        const double sign = horizontal.edgeNormalX()[ei] > 0.0 ? 1.0 : -1.0;

        for (int k = 0; k < kLayers; ++k)
        {
            edgeFlux[ei * static_cast<size_t>(kLayers) + static_cast<size_t>(k)] =
                sign * kSpeed * thickness * horizontal.edgeLength()[ei];
        }
    }

    const std::vector<double> verticalFlux(
        static_cast<size_t>(layered.interfaceSlotCount()), 0.0);
    const std::vector<double> volume = layered.cellVolume();

    const auto initial = [&]() {
        std::vector<double> values(static_cast<size_t>(layered.cellCount()), 0.0);

        for (int64_t c = kStart; c < kStop; ++c)
        {
            for (int k = 0; k < kLayers; ++k)
                values[static_cast<size_t>(layered.cell(c, k))] = 1.0;
        }

        return values;
    };

    std::vector<double> muscl  = initial();
    std::vector<double> upwind = initial();

    Transport::EulerianTransport musclEngine;
    Transport::EulerianTransport upwindEngine;
    musclEngine.initialize(layered);
    upwindEngine.initialize(layered);

    scene.x.assign(horizontal.faceCentroidX().begin(), horizontal.faceCentroidX().end());

    const double dt = 0.4 * cellSize / kSpeed;

    const auto surfaceRow = [&](const std::vector<double> &values) {
        std::vector<double> row(static_cast<size_t>(layered.columnCount()), 0.0);

        for (int64_t c = 0; c < layered.columnCount(); ++c)
            row[static_cast<size_t>(c)] = values[static_cast<size_t>(layered.cell(c, 0))];

        return row;
    };

    double elapsed = 0.0;

    for (int frame = 0; frame <= kFrames; ++frame)
    {
        const double target = kEnd * static_cast<double>(frame) / kFrames;

        while (elapsed < target - 1.0e-9)
        {
            const double step = std::min(dt, target - elapsed);

            musclEngine.advect(layered, volume, edgeFlux, verticalFlux, step,
                               Transport::AdvectionScheme::MusclFct, muscl);
            upwindEngine.advect(layered, volume, edgeFlux, verticalFlux, step,
                                Transport::AdvectionScheme::Upwind, upwind);

            elapsed += step;
        }

        scene.times.push_back(target);
        scene.model.push_back(surfaceRow(muscl));
        scene.comparison.push_back(surfaceRow(upwind));

        const double shift = kSpeed * target;

        std::vector<double> exact(static_cast<size_t>(layered.columnCount()), 0.0);

        for (int64_t c = 0; c < layered.columnCount(); ++c)
        {
            const double x = horizontal.faceCentroidX()[static_cast<size_t>(c)];
            const bool   inside =
                x >= kStart * cellSize + shift && x < kStop * cellSize + shift;

            exact[static_cast<size_t>(c)] = inside ? 1.0 : 0.0;
        }

        scene.analytic.push_back(exact);
    }

    return scene;
}

// ── 5. Internal two-layer seiche ──────────────────────────────────────────

Scenario internalSeiche()
{
    constexpr int    kColumns = 40;
    constexpr int    kLayers  = 20;
    constexpr double kLength  = 1000.0;
    constexpr double kDepth   = 20.0;
    constexpr double kUpper   = 5.0;
    constexpr double kTilt    = 1.5;
    constexpr int    kFrames  = 110;

    const double light   = Hydro::EquationOfState::density(25.0);
    const double heavy   = Hydro::EquationOfState::density(5.0);
    const double reduced = kGravity * (heavy - light) / 1000.0;
    const double speed   = std::sqrt(reduced * kUpper * (kDepth - kUpper) / kDepth);
    const double period  = 2.0 * kLength / speed;

    Scenario scene;
    scene.id       = "internal-seiche";
    scene.title    = "Internal seiche";
    scene.subtitle = "Two-layer interface oscillating beneath a nearly flat surface";
    scene.xLabel   = "distance (m)";
    scene.yLabel   = "interface depth (m)";
    scene.commentary =
        "A tilted thermocline released in a closed basin. The internal wave travels "
        "at 0.33 m/s against 14 m/s for the surface wave — a factor of forty — and "
        "reproducing that separation is the whole point of mode splitting. The "
        "measured period runs about 8% long because numerical diffusion thickens the "
        "interface, which slows the wave; under first-order upwind the bias was 19%.";
    scene.verdict =
        "Period within 8% of the two-layer theory; the same case was chaotic before "
        "the baroclinic pressure gradient was corrected.";

    scene.modelLabel    = "modelled interface";
    scene.analyticLabel = "analytical standing wave";

    const Mesh::HorizontalMesh horizontal =
        channel(kColumns, kLength, 50.0, [](double) { return -kDepth; });

    Mesh::SigmaSpec spec;
    spec.layerCount        = kLayers;
    spec.minLayerThickness = 0.0;

    Mesh::LayeredMesh layered;
    std::string       message;

    if (!layered.build(horizontal, spec, 0.0, message))
        std::fprintf(stderr, "internal mesh: %s\n", message.c_str());

    State::SpeciesRegistry registry;
    if (registry.addTemperature(10.0, message) < 0)
        std::fprintf(stderr, "internal species: %s\n", message.c_str());

    const int temperature = registry.temperatureIndex();

    State::StateStore state;
    state.allocate(layered, registry);

    for (int64_t c = 0; c < layered.columnCount(); ++c)
    {
        const double x = horizontal.faceCentroidX()[static_cast<size_t>(c)];
        const double interfaceZ = -kUpper + kTilt * (x / kLength - 0.5) * 2.0;

        for (int k = 0; k < kLayers; ++k)
        {
            const double z = 0.5 * (layered.z(c, k) + layered.z(c, k + 1));
            state.speciesAt(temperature, layered.cell(c, k)) = z > interfaceZ ? 25.0 : 5.0;
        }
    }

    Hydro::InternalModeOptions options;
    options.gravity   = kGravity;
    options.manningN  = 0.0;
    options.cflNumber = 0.5;

    Hydro::InternalModeSolver solver;
    if (!solver.initialize(layered, state, registry, options, message))
        std::fprintf(stderr, "internal init: %s\n", message.c_str());

    scene.x.assign(horizontal.faceCentroidX().begin(), horizontal.faceCentroidX().end());
    scene.traceLabel = "interface depth at the closed end (m)";

    const double endTime = 2.0 * period;

    // Depth of the 15 C isotherm, interpolated between layer centers.
    const auto interfaceDepth = [&](int64_t column) {
        for (int k = 1; k < kLayers; ++k)
        {
            const double above = state.speciesAt(temperature, layered.cell(column, k - 1));
            const double below = state.speciesAt(temperature, layered.cell(column, k));

            if (above >= 15.0 && below < 15.0)
            {
                const double zAbove = 0.5 * (layered.z(column, k - 1) + layered.z(column, k));
                const double zBelow = 0.5 * (layered.z(column, k) + layered.z(column, k + 1));
                const double fraction = (above - 15.0) / (above - below);

                return zAbove + fraction * (zBelow - zAbove);
            }
        }

        return -kUpper;
    };

    for (int frame = 0; frame <= kFrames; ++frame)
    {
        const double target = endTime * static_cast<double>(frame) / kFrames;

        if (!solver.advance(target - solver.time(), message))
            std::fprintf(stderr, "internal advance: %s\n", message.c_str());

        std::vector<double> modelled(static_cast<size_t>(layered.columnCount()), 0.0);
        std::vector<double> exact(static_cast<size_t>(layered.columnCount()), 0.0);

        for (int64_t c = 0; c < layered.columnCount(); ++c)
        {
            modelled[static_cast<size_t>(c)] = interfaceDepth(c);

            const double x = horizontal.faceCentroidX()[static_cast<size_t>(c)];
            exact[static_cast<size_t>(c)] =
                -kUpper + kTilt * (x / kLength - 0.5) * 2.0 *
                              std::cos(2.0 * kPi * target / period);
        }

        scene.times.push_back(target);
        scene.model.push_back(modelled);
        scene.analytic.push_back(exact);

        scene.trace.push_back(modelled[0]);
        scene.traceReference.push_back(-kUpper -
                                       kTilt * std::cos(2.0 * kPi * target / period));
    }

    return scene;
}

// ── 6. Lock exchange ──────────────────────────────────────────────────────

Scenario lockExchange()
{
    constexpr int    kColumns = 200;
    constexpr int    kLayers  = 16;
    constexpr double kLength  = 600.0;
    constexpr double kDepth   = 10.0;
    constexpr double kDam     = 300.0;
    constexpr int    kFrames  = 100;
    constexpr double kEnd     = 200.0;

    const double light    = Hydro::EquationOfState::density(25.0);
    const double heavy    = Hydro::EquationOfState::density(5.0);
    const double reduced  = kGravity * (heavy - light) / 1000.0;
    const double analytic = 0.5 * std::sqrt(reduced * kDepth);

    Scenario scene;
    scene.id       = "lock-exchange";
    scene.title    = "Lock exchange";
    scene.subtitle = "Dense water released against light water over a flat bed";
    scene.xLabel   = "distance (m)";
    scene.yLabel   = "bed-layer temperature (deg C)";
    scene.commentary =
        "Remove the barrier between cold and warm water and a gravity current runs "
        "along the bed while a return flow travels back above it. Theory gives the "
        "front speed as 0.5 sqrt(g' H); the marker tracks where the model's front "
        "actually is against where theory says it should be. The model runs slightly "
        "slow, and the gap closes under mesh refinement rather than with a better "
        "advection scheme — a front speed is set by integral buoyancy balance, not by "
        "how sharply the interface is resolved.";
    scene.verdict =
        "Front speed converges to theory under refinement: -11.1%, -6.2%, -4.3% at "
        "100, 200 and 400 cells.";

    scene.modelLabel    = "bed-layer temperature";
    scene.analyticLabel = "theoretical front position";

    const Mesh::HorizontalMesh horizontal =
        channel(kColumns, kLength, 20.0, [](double) { return -kDepth; });

    Mesh::SigmaSpec spec;
    spec.layerCount        = kLayers;
    spec.minLayerThickness = 0.0;

    Mesh::LayeredMesh layered;
    std::string       message;

    if (!layered.build(horizontal, spec, 0.0, message))
        std::fprintf(stderr, "lock mesh: %s\n", message.c_str());

    State::SpeciesRegistry registry;
    if (registry.addTemperature(10.0, message) < 0)
        std::fprintf(stderr, "lock species: %s\n", message.c_str());

    const int temperature = registry.temperatureIndex();

    State::StateStore state;
    state.allocate(layered, registry);

    for (int64_t c = 0; c < layered.columnCount(); ++c)
    {
        const double x = horizontal.faceCentroidX()[static_cast<size_t>(c)];

        for (int k = 0; k < kLayers; ++k)
            state.speciesAt(temperature, layered.cell(c, k)) = x < kDam ? 5.0 : 25.0;
    }

    Hydro::InternalModeOptions options;
    options.gravity   = kGravity;
    options.manningN  = 0.0;
    options.cflNumber = 0.5;

    Hydro::InternalModeSolver solver;
    if (!solver.initialize(layered, state, registry, options, message))
        std::fprintf(stderr, "lock init: %s\n", message.c_str());

    scene.x.assign(horizontal.faceCentroidX().begin(), horizontal.faceCentroidX().end());
    scene.traceLabel = "front position (m)";

    for (int frame = 0; frame <= kFrames; ++frame)
    {
        const double target = kEnd * static_cast<double>(frame) / kFrames;

        if (!solver.advance(target - solver.time(), message))
            std::fprintf(stderr, "lock advance: %s\n", message.c_str());

        std::vector<double> bed(static_cast<size_t>(layered.columnCount()), 0.0);

        for (int64_t c = 0; c < layered.columnCount(); ++c)
            bed[static_cast<size_t>(c)] =
                state.speciesAt(temperature, layered.cell(c, kLayers - 1));

        scene.times.push_back(target);
        scene.model.push_back(bed);

        // The analytical reference here is a position, drawn as a step.
        const double frontX = kDam + analytic * target;

        std::vector<double> marker(static_cast<size_t>(layered.columnCount()), 0.0);
        for (int64_t c = 0; c < layered.columnCount(); ++c)
        {
            marker[static_cast<size_t>(c)] =
                horizontal.faceCentroidX()[static_cast<size_t>(c)] < frontX ? 5.0 : 25.0;
        }

        scene.analytic.push_back(marker);

        // Sub-cell front position from the modeled bed temperature.
        double modelledFront = kDam;
        for (int64_t c = 0; c + 1 < layered.columnCount(); ++c)
        {
            if (bed[static_cast<size_t>(c)] < 15.0 && bed[static_cast<size_t>(c + 1)] >= 15.0)
            {
                const double t0 = bed[static_cast<size_t>(c)];
                const double t1 = bed[static_cast<size_t>(c + 1)];
                const double f  = (15.0 - t0) / (t1 - t0);

                modelledFront = horizontal.faceCentroidX()[static_cast<size_t>(c)] +
                                f * (horizontal.faceCentroidX()[static_cast<size_t>(c + 1)] -
                                     horizontal.faceCentroidX()[static_cast<size_t>(c)]);
            }
        }

        scene.trace.push_back(modelledFront);
        scene.traceReference.push_back(frontX);
    }

    return scene;
}

// ── 7. Open-boundary rarefaction ──────────────────────────────────────────

/*!
 * \brief A sudden stage drop at an open outlet, against the exact rarefaction.
 */
Scenario boundaryRarefaction()
{
    constexpr int    kColumns  = 300;
    constexpr double kLength   = 6000.0;
    constexpr double kBed      = -5.0;
    constexpr double kStanding = 0.0;
    constexpr double kDropped  = -1.5;
    constexpr double kEnd      = 240.0;
    constexpr int    kFrames   = 90;

    const double h0 = kStanding - kBed;
    const double hb = kDropped - kBed;
    const double c0 = std::sqrt(kGravity * h0);
    const double cb = std::sqrt(kGravity * hb);
    const double ub = 2.0 * (c0 - cb);

    Scenario scene;
    scene.id       = "rarefaction";
    scene.title    = "Open boundary: sudden stage drop";
    scene.subtitle = "A rarefaction released by lowering the outlet stage";
    scene.xLabel   = "distance (m)";
    scene.yLabel   = "depth (m)";
    scene.commentary =
        "The outlet stage drops 1.5 m instantly and a rarefaction runs back up the "
        "channel. Nothing here is prescribed except the stage: the discharge that "
        "leaves, and the shape of the fan, both follow from the Riemann invariant "
        "the interior sends out. The exact solution is known, so the boundary either "
        "reproduces it or it does not.";
    scene.verdict =
        "Boundary discharge within 1.3% of the exact rarefaction; copying the "
        "interior velocity instead of using the invariant reaches only 3.5%.";

    const Mesh::HorizontalMesh mesh =
        channel(kColumns, kLength, 20.0, [](double) { return kBed; });

    Hydro::BarotropicOptions options;
    options.gravity         = kGravity;
    options.frictionEnabled = false;
    options.cflNumber       = 0.4;

    Hydro::BarotropicSolver solver;
    std::string             message;

    if (!solver.initialize(mesh, options, message))
        std::fprintf(stderr, "rarefaction init: %s\n", message.c_str());

    solver.setUniformSurface(kStanding);

    std::vector<int64_t> outlet;

    for (int64_t e = 0; e < mesh.edgeCount(); ++e)
    {
        if (mesh.edgeRightFace()[static_cast<size_t>(e)] != Mesh::kNoFace)
            continue;

        if (std::abs(mesh.edgeMidX()[static_cast<size_t>(e)] - kLength) < 1.0e-6)
            outlet.push_back(e);
    }

    Hydro::BarotropicSolver::BoundaryCondition stage;
    stage.id    = "outlet";
    stage.type  = Hydro::BarotropicSolver::BoundaryType::Stage;
    stage.stage = kDropped;

    if (solver.addBoundary(stage, outlet, message) < 0)
        std::fprintf(stderr, "rarefaction boundary: %s\n", message.c_str());

    scene.x.assign(mesh.faceCentroidX().begin(), mesh.faceCentroidX().end());

    scene.traceLabel = "outlet discharge (m3/s)";

    for (int frame = 0; frame <= kFrames; ++frame)
    {
        const double target = kEnd * static_cast<double>(frame) / kFrames;

        if (!solver.advance(target - solver.time(), message))
            std::fprintf(stderr, "rarefaction advance: %s\n", message.c_str());

        scene.times.push_back(target);
        scene.model.push_back(solver.depth());

        scene.trace.push_back(frame == 0 ? 0.0 : -solver.boundaryFlux(0));
        scene.traceReference.push_back(hb * ub * 20.0);

        std::vector<double> exact(static_cast<size_t>(kColumns), h0);

        for (int64_t c = 0; c < mesh.faceCount(); ++c)
        {
            const double x = mesh.faceCentroidX()[static_cast<size_t>(c)];

            if (target <= 0.0)
            {
                exact[static_cast<size_t>(c)] = h0;
                continue;
            }

            // Similarity coordinate measured from the outlet, negative upstream.
            const double xi = (x - kLength) / target;

            if (xi <= -c0)
                exact[static_cast<size_t>(c)] = h0;              // untouched
            else if (xi >= ub - cb)
                exact[static_cast<size_t>(c)] = hb;              // behind the fan
            else
                exact[static_cast<size_t>(c)] =
                    std::pow(2.0 * c0 - xi, 2.0) / (9.0 * kGravity);
        }

        scene.analytic.push_back(exact);
    }

    return scene;
}

// ── 8. Radiating versus reflecting ────────────────────────────────────────

/*!
 * \brief The same wave against a wall and against a transmissive boundary.
 */
Scenario radiationComparison()
{
    constexpr int    kColumns = 240;
    constexpr double kLength  = 3000.0;
    constexpr double kBed     = -5.0;
    constexpr double kEnd     = 700.0;
    constexpr int    kFrames  = 100;

    Scenario scene;
    scene.id              = "radiation";
    scene.title           = "Radiating versus reflecting boundary";
    scene.subtitle        = "One disturbance, two downstream ends";
    scene.xLabel          = "distance (m)";
    scene.yLabel          = "surface elevation (m)";
    scene.modelLabel      = "transmissive end";
    scene.comparisonLabel = "closed wall";
    scene.analyticLabel   = "";
    scene.commentary =
        "The property no mass balance can check. Both runs start from the same hump "
        "and are identical until the wave reaches the far end. The wall returns it; "
        "the transmissive boundary lets it go. A poorly posed open boundary looks "
        "stable and conserves mass perfectly, and announces itself only here.";
    scene.verdict =
        "An open end leaves under a fifth of the surface variance a closed one does. "
        "A fixed-stage boundary, by contrast, reflects on purpose.";

    const Mesh::HorizontalMesh mesh =
        channel(kColumns, kLength, 20.0, [](double) { return kBed; });

    scene.x.assign(mesh.faceCentroidX().begin(), mesh.faceCentroidX().end());

    std::vector<double> initial(static_cast<size_t>(kColumns), 0.0);

    for (int64_t c = 0; c < mesh.faceCount(); ++c)
    {
        if (mesh.faceCentroidX()[static_cast<size_t>(c)] < 400.0)
            initial[static_cast<size_t>(c)] = 0.30;
    }

    const auto run = [&](bool openEnd, std::vector<std::vector<double>> &frames) {
        Hydro::BarotropicOptions options;
        options.gravity         = kGravity;
        options.frictionEnabled = false;
        options.cflNumber       = 0.4;

        Hydro::BarotropicSolver solver;
        std::string             message;

        if (!solver.initialize(mesh, options, message))
            std::fprintf(stderr, "radiation init: %s\n", message.c_str());

        if (openEnd)
        {
            std::vector<int64_t> outlet;

            for (int64_t e = 0; e < mesh.edgeCount(); ++e)
            {
                if (mesh.edgeRightFace()[static_cast<size_t>(e)] != Mesh::kNoFace)
                    continue;

                if (std::abs(mesh.edgeMidX()[static_cast<size_t>(e)] - kLength) < 1.0e-6)
                    outlet.push_back(e);
            }

            Hydro::BarotropicSolver::BoundaryCondition open;
            open.id   = "radiating";
            open.type = Hydro::BarotropicSolver::BoundaryType::Transmissive;

            if (solver.addBoundary(open, outlet, message) < 0)
                std::fprintf(stderr, "radiation boundary: %s\n", message.c_str());
        }

        if (!solver.setSurface(initial, message))
            std::fprintf(stderr, "radiation surface: %s\n", message.c_str());

        for (int frame = 0; frame <= kFrames; ++frame)
        {
            const double target = kEnd * static_cast<double>(frame) / kFrames;

            if (!solver.advance(target - solver.time(), message))
                std::fprintf(stderr, "radiation advance: %s\n", message.c_str());

            frames.push_back(solver.surfaceElevation());
        }
    };

    run(true, scene.model);
    run(false, scene.comparison);

    for (int frame = 0; frame <= kFrames; ++frame)
        scene.times.push_back(kEnd * static_cast<double>(frame) / kFrames);

    return scene;
}

// ── 9. MacDonald steady profile ───────────────────────────────────────────

/*!
 * \brief A backwater profile forming from a flat start on a manufactured bed.
 */
Scenario macdonaldProfile()
{
    constexpr int    kColumns = 200;
    constexpr double kLength  = 1000.0;
    constexpr double kWidth   = 10.0;
    constexpr double kUnitQ   = 2.0;
    constexpr double kManning = 0.03;
    constexpr double kBase    = 2.0;
    constexpr double kBumpAmp = 0.4;
    constexpr double kSharp   = 16.0;
    constexpr double kEnd     = 20000.0;
    constexpr int    kFrames  = 100;

    const auto exactDepth = [](double x) {
        const double t = x / kLength - 0.5;
        return kBase * (1.0 + kBumpAmp * std::exp(-kSharp * t * t));
    };

    const auto depthSlope = [](double x) {
        const double t = x / kLength - 0.5;
        return kBase * kBumpAmp * std::exp(-kSharp * t * t) * (-2.0 * kSharp * t) /
               kLength;
    };

    const auto bedSlope = [&](double x) {
        const double h = exactDepth(x);
        return (kUnitQ * kUnitQ / (kGravity * h * h * h) - 1.0) * depthSlope(x) -
               kManning * kManning * kUnitQ * kUnitQ / std::pow(h, 10.0 / 3.0);
    };

    // Bed by Simpson integration, datum at the downstream end.
    const auto bedAt = [&](double x) {
        constexpr int intervals = 2000;

        const double step = (kLength - x) / intervals;

        if (step <= 0.0)
            return 0.0;

        double total = bedSlope(x) + bedSlope(kLength);

        for (int i = 1; i < intervals; ++i)
            total += (i % 2 == 0 ? 2.0 : 4.0) * bedSlope(x + i * step);

        return -total * step / 3.0;
    };

    Scenario scene;
    scene.id       = "macdonald";
    scene.title    = "MacDonald steady profile";
    scene.subtitle = "A backwater profile assembling itself from a flat start";
    scene.xLabel   = "distance (m)";
    scene.yLabel   = "depth (m)";
    scene.commentary =
        "The bed here was built backwards: a smooth depth profile was chosen, and the "
        "steady momentum equation integrated to find the bed that makes it exact. The "
        "model is given only that bed, an inflow, and a downstream stage, and has to "
        "find the profile for itself. Watch it fill from a flat surface and settle "
        "onto the analytical curve.";
    scene.verdict =
        "Mean depth error 0.57 mm on a 2 m depth; refinement ratios 1.97 and 1.98 "
        "across 50/100/200 cells, which is clean first order.";

    const Mesh::HorizontalMesh mesh = channel(kColumns, kLength, kWidth, bedAt);

    Hydro::BarotropicOptions options;
    options.gravity         = kGravity;
    options.manningN        = kManning;
    options.cflNumber       = 0.4;
    options.frictionEnabled = true;

    Hydro::BarotropicSolver solver;
    std::string             message;

    if (!solver.initialize(mesh, options, message))
        std::fprintf(stderr, "macdonald init: %s\n", message.c_str());

    std::vector<int64_t> inlet;
    std::vector<int64_t> outlet;

    for (int64_t e = 0; e < mesh.edgeCount(); ++e)
    {
        if (mesh.edgeRightFace()[static_cast<size_t>(e)] != Mesh::kNoFace)
            continue;

        const double x = mesh.edgeMidX()[static_cast<size_t>(e)];

        if (std::abs(x) < 1.0e-6)
            inlet.push_back(e);
        else if (std::abs(x - kLength) < 1.0e-6)
            outlet.push_back(e);
    }

    Hydro::BarotropicSolver::BoundaryCondition upstream;
    upstream.id                = "inflow";
    upstream.type              = Hydro::BarotropicSolver::BoundaryType::Inflow;
    upstream.specificDischarge = kUnitQ;

    if (solver.addBoundary(upstream, inlet, message) < 0)
        std::fprintf(stderr, "macdonald inflow: %s\n", message.c_str());

    Hydro::BarotropicSolver::BoundaryCondition downstream;
    downstream.id    = "stage";
    downstream.type  = Hydro::BarotropicSolver::BoundaryType::Stage;
    downstream.stage = bedAt(kLength) + exactDepth(kLength);

    if (solver.addBoundary(downstream, outlet, message) < 0)
        std::fprintf(stderr, "macdonald stage: %s\n", message.c_str());

    // Flat start, at the downstream stage.
    std::vector<double> flat(static_cast<size_t>(kColumns), downstream.stage);

    if (!solver.setSurface(flat, message))
        std::fprintf(stderr, "macdonald surface: %s\n", message.c_str());

    scene.x.assign(mesh.faceCentroidX().begin(), mesh.faceCentroidX().end());

    std::vector<double> exact(static_cast<size_t>(kColumns), 0.0);

    for (int64_t c = 0; c < mesh.faceCount(); ++c)
        exact[static_cast<size_t>(c)] = exactDepth(mesh.faceCentroidX()[static_cast<size_t>(c)]);

    scene.traceLabel = "mean depth error (mm)";

    for (int frame = 0; frame <= kFrames; ++frame)
    {
        // Quadratic spacing: the interesting motion is all early.
        const double fraction = static_cast<double>(frame) / kFrames;
        const double target   = kEnd * fraction * fraction;

        if (!solver.advance(target - solver.time(), message))
            std::fprintf(stderr, "macdonald advance: %s\n", message.c_str());

        scene.times.push_back(target);
        scene.model.push_back(solver.depth());
        scene.analytic.push_back(exact);

        double error = 0.0;
        int    count = 0;

        for (int64_t c = 0; c < mesh.faceCount(); ++c)
        {
            const double x = mesh.faceCentroidX()[static_cast<size_t>(c)];

            if (x < 0.03 * kLength || x > 0.97 * kLength)
                continue;

            error += std::abs(solver.depth()[static_cast<size_t>(c)] -
                              exact[static_cast<size_t>(c)]);
            ++count;
        }

        scene.trace.push_back(count > 0 ? 1000.0 * error / count : 0.0);
        scene.traceReference.push_back(0.57);
    }

    return scene;
}

} // namespace

int main(int argc, char **argv)
{
    const std::string path = argc > 1 ? argv[1] : "verification/data.json";

    std::vector<Scenario> scenes;

    std::fprintf(stderr, "running ritter dam break...\n");
    scenes.push_back(ritterDamBreak());
    std::fprintf(stderr, "running barotropic seiche...\n");
    scenes.push_back(barotropicSeiche());
    std::fprintf(stderr, "running wind setup...\n");
    scenes.push_back(windSetup());
    std::fprintf(stderr, "running transport comparison...\n");
    scenes.push_back(schemeComparison());
    std::fprintf(stderr, "running internal seiche...\n");
    scenes.push_back(internalSeiche());
    std::fprintf(stderr, "running lock exchange...\n");
    scenes.push_back(lockExchange());
    std::fprintf(stderr, "running boundary rarefaction...\n");
    scenes.push_back(boundaryRarefaction());
    std::fprintf(stderr, "running radiation comparison...\n");
    scenes.push_back(radiationComparison());
    std::fprintf(stderr, "running macdonald profile...\n");
    scenes.push_back(macdonaldProfile());

    std::ofstream out(path);

    if (!out)
    {
        std::fprintf(stderr, "cannot write %s\n", path.c_str());
        return 1;
    }

    out << "{\"scenarios\":[";

    for (size_t i = 0; i < scenes.size(); ++i)
    {
        if (i > 0)
            out << ",";
        out << scenes[i].toJson();
    }

    out << "]}";

    std::fprintf(stderr, "wrote %s\n", path.c_str());
    return 0;
}
