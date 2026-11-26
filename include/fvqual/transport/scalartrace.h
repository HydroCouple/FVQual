// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
#ifndef FVQUAL_TRANSPORT_SCALARTRACE_H
#define FVQUAL_TRANSPORT_SCALARTRACE_H
#include "fvqual/fvqual_export.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/forcing/sources.h"
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace FVQual::Transport {
//! Discrete passive-scalar derivatives on the executed hydraulic trajectory.
//! Geometry, flows, placement, mixing coefficients and substep choices are
//! prescribed. This is NOT a temperature/hydrodynamic or SWMM adjoint.
struct ScalarEntry { size_t row, column; double value; };
struct FVQual_EXPORT ScalarStep {
    double time = 0, dt = 0;
    int layers = 0;
    std::vector<ScalarEntry> advection, inflow;
    std::vector<std::string> sourceIds;
    std::vector<double> sourceValues, lower, diagonal, upper, before, after;
    [[nodiscard]] std::vector<double> apply(std::span<const double> values,std::span<const double> sources) const;
    [[nodiscard]] std::pair<std::vector<double>,std::vector<double>> transpose(std::span<const double> cotangents) const;
    static ScalarStep capture(const Mesh::LayeredMesh &, std::span<const double> oldVolume,
        std::span<const double> edgeFlux, std::span<const double> verticalFlux,
        const Forcing::SourceSet &, int species, std::span<const double> values,
        double time, double dt);
};
//! Immutable persistent history: saving a solver checkpoint shares completed
//! records instead of copying all prior steps. Each record remains O(cells+faces).
struct ScalarTraceNode {
    ScalarStep step;
    std::shared_ptr<const ScalarTraceNode> previous;
};
struct ScalarPullback {
    std::vector<double> initial, sources;
    double diffusionMultiplier = 0;
};
class FVQual_EXPORT ScalarTrajectory {
public:
    std::vector<double> initial;
    double startTime = 0;
    size_t count = 0;
    std::shared_ptr<const ScalarTraceNode> head;
    [[nodiscard]] std::vector<const ScalarStep *> steps() const;
    [[nodiscard]] std::vector<double> sourceValues() const;
    //! All states including initial, step-major. Positive multiplier scales
    //! only the recorded vertical diffusion operator; hydraulics stay frozen.
    [[nodiscard]] std::vector<double> evaluate(std::span<const double> initialValues,
        std::span<const double> sources, double diffusionMultiplier = 1) const;
    [[nodiscard]] std::vector<double> jvp(std::span<const double> initialDirection,
        std::span<const double> sourceDirection, double diffusionDirection = 0) const;
    //! Cotangents for ALL states, including initial, accumulated in reverse order.
    [[nodiscard]] ScalarPullback vjp(std::span<const double> stateCotangents) const;
    [[nodiscard]] size_t bytes() const;
};
}
#endif
