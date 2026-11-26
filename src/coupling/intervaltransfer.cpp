// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
#include "fvqual/coupling/intervaltransfer.h"
#include <cmath>
#include <stdexcept>

namespace FVQual::Coupling {
namespace {
bool identical(const TransferFrame &a, const TransferFrame &b) {
    if (a.startSeconds != b.startSeconds || a.endSeconds != b.endSeconds ||
        a.inflows.size() != b.inflows.size()) return false;
    for (size_t i = 0; i < a.inflows.size(); ++i) {
        const auto &x = a.inflows[i], &y = b.inflows[i];
        if (x.column != y.column || x.volumeM3 != y.volumeM3 || x.heatJ != y.heatJ ||
            x.placement != y.placement || x.elevationM != y.elevationM ||
            x.speciesMassKg != y.speciesMassKg) return false;
    }
    return true;
}
}

IntervalTransferReceiver::IntervalTransferReceiver(Hydro::InternalModeSolver &solver,
    const State::SpeciesRegistry &registry, const Mesh::LayeredMesh &mesh)
    : m_solver(&solver), m_registry(&registry), m_mesh(&mesh) {}

bool IntervalTransferReceiver::apply(const TransferFrame &frame,
                                    TransferReceipt &receipt, std::string &message)
{
    try {
        auto found = m_frames.find(frame.id);
        if (found != m_frames.end()) {
            if (!identical(frame, found->second))
                throw std::invalid_argument("transfer ID reused with different contents");
            receipt = m_receipts.at(frame.id);
            receipt.replayed = true;
            message.clear();
            return true;
        }
        const double dt = frame.endSeconds - frame.startSeconds;
        if (frame.id.empty() || !std::isfinite(dt) || !std::isfinite(frame.startSeconds) ||
            !std::isfinite(frame.endSeconds) || !(dt > 0.0) ||
            std::abs(frame.startSeconds - m_solver->time()) > 1.0e-8)
            throw std::invalid_argument("transfer needs unique ID and interval starting at receiver time");
        auto originalSources = m_solver->sources();
        auto sources = originalSources;
        TransferReceipt candidate;
        candidate.id = frame.id;
        candidate.startSeconds = frame.startSeconds;
        candidate.endSeconds = frame.endSeconds;
        const double heatPerVolumeDegree = m_solver->options().referenceDensity *
                                           m_solver->options().heatCapacity;
        for (const auto &inflow : frame.inflows) {
            if (inflow.column < 0 || inflow.column >= m_mesh->columnCount() ||
                !std::isfinite(inflow.volumeM3) || inflow.volumeM3 < 0.0 ||
                !std::isfinite(inflow.heatJ) || !std::isfinite(inflow.elevationM) ||
                !std::isfinite(heatPerVolumeDegree) || !(heatPerVolumeDegree > 0.0))
                throw std::invalid_argument("invalid inflow column, amount or heat convention");
            if (inflow.placement != Forcing::Placement::Uniform &&
                inflow.placement != Forcing::Placement::Surface &&
                inflow.placement != Forcing::Placement::Bed &&
                inflow.placement != Forcing::Placement::Elevation &&
                inflow.placement != Forcing::Placement::DensityMatched)
                throw std::invalid_argument("unsupported inflow placement");
            Forcing::Source source;
            source.id = frame.id + ":" + std::to_string(sources.count()-originalSources.count());
            source.origin = Forcing::SourceOrigin::Coupled;
            source.column = inflow.column;
            source.discharge = inflow.volumeM3 / dt;
            source.placement = inflow.placement;
            source.elevation = inflow.elevationM;
            source.concentration.assign(static_cast<size_t>(m_registry->count()), 0.0);
            if (inflow.volumeM3 == 0.0 && inflow.heatJ != 0.0)
                throw std::invalid_argument("water-associated heat requires positive water volume");
            if (m_registry->temperatureIndex() < 0 && inflow.heatJ != 0.0)
                throw std::invalid_argument("receiver has no temperature species");
            source.temperature = inflow.volumeM3 > 0.0
                ? (inflow.heatJ / inflow.volumeM3) / heatPerVolumeDegree : 0.0;
            if (!std::isfinite(source.discharge) || !std::isfinite(source.temperature))
                throw std::invalid_argument("inflow discharge or temperature overflows");
            for (const auto &[id, amount] : inflow.speciesMassKg) {
                const int index = m_registry->indexOf(id);
                if (index < 0 || !std::isfinite(amount) || amount < 0.0 ||
                    (inflow.volumeM3 == 0.0 && amount != 0.0))
                    throw std::invalid_argument("invalid species identity or mass");
                const auto &definition = m_registry->definition(index);
                if (!definition.isTransported() ||
                    definition.kind == State::SpeciesKind::Temperature ||
                    definition.kind == State::SpeciesKind::Age ||
                    (definition.units != "kg/m3" && definition.units != "kg m-3"))
                    throw std::invalid_argument("species transfer requires transported kg/m3 species");
                source.concentration[static_cast<size_t>(index)] =
                    inflow.volumeM3 > 0.0 ? amount / inflow.volumeM3 : 0.0;
                if (!std::isfinite(source.concentration[static_cast<size_t>(index)]))
                    throw std::invalid_argument("inflow concentration overflows");
                candidate.acceptedSpeciesKg[id] += amount;
            }
            for (const auto &definition : m_registry->definitions())
                if (definition.isTransported() &&
                    definition.kind != State::SpeciesKind::Temperature &&
                    definition.kind != State::SpeciesKind::Age &&
                    !inflow.speciesMassKg.contains(definition.id))
                    throw std::invalid_argument("declare mass, including zero, for every incoming species");
            if (sources.add(source, message) < 0)
                throw std::invalid_argument(message);
            candidate.requestedVolumeM3 += inflow.volumeM3;
            candidate.requestedHeatJ += inflow.heatJ;
        }
        candidate.acceptedVolumeM3 = candidate.requestedVolumeM3;
        candidate.acceptedHeatJ = candidate.requestedHeatJ;
        if (!std::isfinite(candidate.acceptedVolumeM3) ||
            !std::isfinite(candidate.acceptedHeatJ))
            throw std::invalid_argument("transfer totals overflow");
        for (const auto &[id, mass] : candidate.acceptedSpeciesKg)
            if (!std::isfinite(mass))
                throw std::invalid_argument("species transfer total overflows: " + id);
        TransferReceipt output = candidate;
        auto checkpoint = m_solver->checkpoint(message);
        if (!checkpoint) return false;
        // Allocate ledger entries before touching numerical state. On failure
        // remove them and restore the complete solver, not just concentrations.
        m_frames.emplace(frame.id, frame);
        try { m_receipts.emplace(frame.id, candidate); }
        catch (...) { m_frames.erase(frame.id); throw; }
        m_solver->sources() = std::move(sources);
        bool advanced = false;
        try { advanced = m_solver->advance(dt, message); }
        catch (const std::exception &error) { message = error.what(); }
        if (!advanced) {
            const std::string failure = message;
            std::string restoreMessage;
            const bool restored = m_solver->restoreCheckpoint(*checkpoint, restoreMessage);
            m_frames.erase(frame.id);
            m_receipts.erase(frame.id);
            message = failure + (restored ? "; interval rolled back" : "; rollback failed: " + restoreMessage);
            return false;
        }
        // The checkpoint's source configuration is needed again. The temporary
        // interval sources expire exactly here, independent of reporting time.
        m_solver->sources() = std::move(originalSources);
        receipt = std::move(output);
        message.clear();
        return true;
    } catch (const std::exception &error) {
        message = error.what();
        return false;
    }
}
} // namespace FVQual::Coupling
