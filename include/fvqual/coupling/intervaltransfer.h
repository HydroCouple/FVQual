// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
#ifndef FVQUAL_COUPLING_INTERVALTRANSFER_H
#define FVQUAL_COUPLING_INTERVALTRANSFER_H
#include "fvqual/hydro/internalmode.h"
#include <map>

namespace FVQual::Coupling {

//! Accepted donor exports, already debited by SWMM. Times are elapsed seconds.
//! Heat is rho*cp*volume*temperature relative to 0 Celsius; species mass is kg.
//! Only positive one-way transfers are supported here. Reversals need a donor
//! extraction/capacity protocol and are explicitly refused.
struct InflowTransfer {
    int64_t column = 0;
    double volumeM3 = 0.0, heatJ = 0.0;
    Forcing::Placement placement = Forcing::Placement::Uniform;
    double elevationM = 0.0;
    std::map<std::string, double> speciesMassKg;
};
struct TransferFrame {
    std::string id;
    double startSeconds = 0.0, endSeconds = 0.0;
    std::vector<InflowTransfer> inflows;
};
struct TransferReceipt {
    std::string id;
    double startSeconds = 0.0, endSeconds = 0.0;
    double requestedVolumeM3 = 0.0, acceptedVolumeM3 = 0.0;
    double requestedHeatJ = 0.0, acceptedHeatJ = 0.0;
    std::map<std::string, double> acceptedSpeciesKg;
    bool replayed = false;
};

//! Executes one synchronized, conservative interval with rollback on failure.
//! Transfer IDs are idempotent only for identical frames. Receipts and original
//! frames are retained for the session. Not safe for concurrent calls.
//! Native species concentration units must be kg/m3 for this initial contract.
class FVQual_EXPORT IntervalTransferReceiver {
public:
    IntervalTransferReceiver(Hydro::InternalModeSolver &solver,
                             const State::SpeciesRegistry &registry,
                             const Mesh::LayeredMesh &mesh);
    [[nodiscard]] bool apply(const TransferFrame &frame, TransferReceipt &receipt,
                             std::string &message);
    [[nodiscard]] const std::map<std::string, TransferReceipt> &receipts() const noexcept {
        return m_receipts;
    }
private:
    Hydro::InternalModeSolver *m_solver;
    const State::SpeciesRegistry *m_registry;
    const Mesh::LayeredMesh *m_mesh;
    std::map<std::string, TransferFrame> m_frames;
    std::map<std::string, TransferReceipt> m_receipts;
};

} // namespace FVQual::Coupling
#endif
