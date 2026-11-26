// SPDX-License-Identifier: MIT
#ifndef FVQUAL_TRANSPORT_REACTIVETRACE_H
#define FVQUAL_TRANSPORT_REACTIVETRACE_H
#include "fvqual/transport/scalartrace.h"
#include "fvqual/kinetics/reactiontrace.h"
namespace FVQual::Transport {
struct ReactiveStep {
    std::vector<ScalarStep> transport;
    Kinetics::ReactionStep reaction;
    std::vector<double> after;
};
struct ReactiveTraceNode {ReactiveStep step;std::shared_ptr<const ReactiveTraceNode> previous;};
struct ReactivePullback {std::vector<double> initial,constants;size_t peakTapeNodes=0;};
//! Multispecies native upwind/diffusion/reaction trajectory. Hydraulic state,
//! mixing, source flows/composition/placement and geometry are prescribed.
//! Initial concentrations and kinetic constants are differentiated. State
//! blocks are selected-species-major within each recorded time.
struct FVQual_EXPORT ReactiveTrajectory {
    std::vector<int> species;
    std::vector<std::string> constantNames;
    std::vector<double> initial,constants;
    double startTime=0;
    size_t cells=0,count=0;
    std::shared_ptr<const ReactiveTraceNode> head;
    [[nodiscard]] std::vector<const ReactiveStep *> steps() const;
    [[nodiscard]] std::vector<double> evaluate(std::span<const double> initialValues,std::span<const double> constantValues) const;
    [[nodiscard]] std::vector<double> jvp(std::span<const double> initialDirection,std::span<const double> constantDirection) const;
    [[nodiscard]] ReactivePullback vjp(std::span<const double> stateCotangents) const;
};
}
#endif
