// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
#include "capi_internal.hpp"
#include "fvqual/capi/fvqual_coupling.h"
#include <atomic>
using FVQual::CApi::fail;
using FVQual::CApi::guard;
using FVQual::CApi::succeed;
namespace {
FVQ_ErrorCode ready(const FVQ_Simulation *sim) {
    if (!sim) return fail(FVQ_ERR_NULL_ARGUMENT, "simulation is null");
    if (!sim->initialized) return fail(FVQ_ERR_LIFECYCLE, "coupling/replay requires initialized simulation");
    return FVQ_OK;
}
void receiver(FVQ_Simulation *sim) {
    if (!sim->receiver)
        sim->receiver = std::make_unique<FVQual::Coupling::IntervalTransferReceiver>(
            sim->solver, sim->registry, sim->layered);
}
}
extern "C" {
FVQ_ErrorCode fvq_sim_apply_inflow_frame(FVQ_Simulation *sim, const char *id,
    double start, double end, const FVQ_InflowTransfer *transfers,
    int64_t count, FVQ_TransferReceipt *receipt)
{
    return guard([&]() -> FVQ_ErrorCode {
        if (auto error = ready(sim); error != FVQ_OK) return error;
        if (!id || !receipt || (count > 0 && !transfers))
            return fail(FVQ_ERR_NULL_ARGUMENT, "transfer ID, receipt and transfer buffer are required");
        if (count < 0) return fail(FVQ_ERR_INVALID_ARGUMENT, "negative transfer count");
        FVQual::Coupling::TransferFrame frame;
        frame.id = id; frame.startSeconds = start; frame.endSeconds = end;
        for (int64_t i = 0; i < count; ++i) {
            const auto &t = transfers[i];
            if (t.speciesCount < 0 || (t.speciesCount > 0 && !t.species))
                return fail(FVQ_ERR_INVALID_ARGUMENT, "invalid species count/buffer");
            if (t.placement < FVQ_INFLOW_UNIFORM || t.placement > FVQ_INFLOW_DENSITY_MATCHED)
                return fail(FVQ_ERR_INVALID_ARGUMENT, "invalid inflow placement");
            FVQual::Coupling::InflowTransfer inflow;
            inflow.column = t.column; inflow.volumeM3 = t.volumeM3; inflow.heatJ = t.heatJ;
            inflow.placement = static_cast<FVQual::Forcing::Placement>(t.placement);
            inflow.elevationM = t.elevationM;
            for (int64_t s = 0; s < t.speciesCount; ++s) {
                if (!t.species[s].id) return fail(FVQ_ERR_NULL_ARGUMENT, "species ID is null");
                if (!inflow.speciesMassKg.emplace(t.species[s].id, t.species[s].massKg).second)
                    return fail(FVQ_ERR_SPECIES, "duplicate transfer species ID");
            }
            frame.inflows.push_back(std::move(inflow));
        }
        receiver(sim);
        FVQual::Coupling::TransferReceipt accepted;
        std::string message;
        if (!sim->receiver->apply(frame, accepted, message))
            return fail(FVQ_ERR_INVALID_ARGUMENT, message);
        receipt->requestedVolumeM3 = accepted.requestedVolumeM3;
        receipt->acceptedVolumeM3 = accepted.acceptedVolumeM3;
        receipt->requestedHeatJ = accepted.requestedHeatJ;
        receipt->acceptedHeatJ = accepted.acceptedHeatJ;
        receipt->replayed = accepted.replayed;
        return succeed();
    });
}
FVQ_ErrorCode fvq_sim_transfer_species_receipt(const FVQ_Simulation *sim,
    const char *transferId, const char *speciesId, double *acceptedKg)
{
    return guard([&]() -> FVQ_ErrorCode {
        if (auto error = ready(sim); error != FVQ_OK) return error;
        if (!transferId || !speciesId || !acceptedKg)
            return fail(FVQ_ERR_NULL_ARGUMENT, "receipt identity and output are required");
        if (!sim->receiver) return fail(FVQ_ERR_INVALID_ARGUMENT, "no transfer receipts");
        const auto &receipts = sim->receiver->receipts();
        auto found = receipts.find(transferId);
        if (found == receipts.end() || !found->second.acceptedSpeciesKg.contains(speciesId))
            return fail(FVQ_ERR_SPECIES, "transfer/species receipt not found");
        *acceptedKg = found->second.acceptedSpeciesKg.at(speciesId);
        return succeed();
    });
}
FVQ_ErrorCode fvq_sim_save_state(FVQ_Simulation *sim, uint64_t *token)
{
    return guard([&]() -> FVQ_ErrorCode {
        if (auto error = ready(sim); error != FVQ_OK) return error;
        if (!token) return fail(FVQ_ERR_NULL_ARGUMENT, "checkpoint output is null");
        FVQ_Simulation_::Checkpoint saved;
        std::string message;
        saved.solver = sim->solver.checkpoint(message);
        if (!saved.solver) return fail(FVQ_ERR_UNSUPPORTED, message);
        if (sim->receiver)
            saved.receiver = std::make_unique<FVQual::Coupling::IntervalTransferReceiver>(*sim->receiver);
        static std::atomic<uint64_t> next{1};
        const auto id = next.fetch_add(1);
        sim->checkpoints.emplace(id, std::move(saved));
        *token = id;
        return succeed();
    });
}
FVQ_ErrorCode fvq_sim_restore_state(FVQ_Simulation *sim, uint64_t token)
{
    return guard([&]() -> FVQ_ErrorCode {
        if (auto error = ready(sim); error != FVQ_OK) return error;
        const auto found = sim->checkpoints.find(token);
        if (found == sim->checkpoints.end()) return fail(FVQ_ERR_INVALID_ARGUMENT, "invalid checkpoint token");
        std::unique_ptr<FVQual::Coupling::IntervalTransferReceiver> restored;
        if (found->second.receiver)
            restored = std::make_unique<FVQual::Coupling::IntervalTransferReceiver>(*found->second.receiver);
        std::string message;
        if (!sim->solver.restoreCheckpoint(*found->second.solver, message))
            return fail(FVQ_ERR_INVALID_ARGUMENT, message);
        sim->receiver = std::move(restored);
        return succeed();
    });
}
FVQ_ErrorCode fvq_sim_release_state(FVQ_Simulation *sim, uint64_t token)
{
    return guard([&]() -> FVQ_ErrorCode {
        if (auto error = ready(sim); error != FVQ_OK) return error;
        if (sim->checkpoints.erase(token) != 1)
            return fail(FVQ_ERR_INVALID_ARGUMENT, "invalid or released checkpoint token");
        return succeed();
    });
}
}
