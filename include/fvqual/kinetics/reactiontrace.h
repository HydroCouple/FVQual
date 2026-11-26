// SPDX-License-Identifier: MIT
#ifndef FVQUAL_KINETICS_REACTIONTRACE_H
#define FVQUAL_KINETICS_REACTIONTRACE_H
#include "fvqual/kinetics/reactionnetwork.h"
#include <memory>
namespace FVQual::Kinetics {
struct ReactionPullback {
    std::vector<double> species, constants;
    // light, depth, thickness, surface-area/volume and bed-area/volume, field-major.
    std::vector<double> geometry;
    size_t peakTapeNodes=0;
};
//! Immutable native reaction operator. Derivatives include the smooth executed
//! adaptive Cash-Karp controller, not a frozen internal time schedule. Branch
//! boundaries, invalid expression domains and unconverged solves are refused.
class FVQual_EXPORT ReactionStep {
public:
    ReactionStep()=default;
    static ReactionStep capture(const ReactionNetwork &,const State::SpeciesRegistry &,
        const State::StateStore &,double dt,const DerivedFields &,
        std::span<const unsigned char> position);
    [[nodiscard]] std::vector<double> evaluate(std::span<const double> species,
        std::span<const double> constants) const;
    [[nodiscard]] std::vector<double> jvp(std::span<const double> speciesDirection,
        std::span<const double> constantDirection,std::span<const double> geometryDirection={}) const;
    [[nodiscard]] ReactionPullback vjp(std::span<const double> cotangents) const;
    [[nodiscard]] const std::vector<double> &before() const;
    [[nodiscard]] const std::vector<double> &constants() const;
private:
    struct Data;
    std::shared_ptr<const Data> m_data;
};
}
#endif
