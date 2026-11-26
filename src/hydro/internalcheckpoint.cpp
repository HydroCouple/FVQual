// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
#include "fvqual/hydro/internalmode.h"
#include <algorithm>
#include <stdexcept>

namespace FVQual::Hydro {
void InternalModeSolver::reset() noexcept
{
    const uint64_t next = m_checkpointGeneration + 1;
    *this = InternalModeSolver{};
    m_checkpointGeneration = next;
}

struct InternalModeCheckpoint
{
    const InternalModeSolver *owner = nullptr;
    uint64_t generation = 0;
    State::StateStore state;
    Mesh::LayeredMesh mesh;
    std::unique_ptr<InternalModeSolver> solver;
};

// Keep this list complete when adding private solver state. The pressure
// closure is copied, never shared; barotropic callbacks remain bound to owner.
InternalModeSolver::InternalModeSolver(const InternalModeSolver &other, CheckpointCopy)
    :
    m_checkpointGeneration(other.m_checkpointGeneration),
    m_checkpointReady(other.m_checkpointReady),
    m_traceSpecies(other.m_traceSpecies),
    m_traceLimit(other.m_traceLimit),
    m_scalarTrace(other.m_scalarTrace),
    m_recordReactive(other.m_recordReactive),
    m_reactiveLimit(other.m_reactiveLimit),
    m_reactiveTrace(other.m_reactiveTrace),
    m_boundaryLinks(other.m_boundaryLinks),
    m_boundaryValue(other.m_boundaryValue),
    m_mesh(other.m_mesh),
    m_state(other.m_state),
    m_species(other.m_species),
    m_options(other.m_options),
    m_momentumAudit(other.m_momentumAudit),
    m_budget(other.m_budget),
    m_heatAudit(other.m_heatAudit),
    m_advectionAudit(other.m_advectionAudit),
    m_heatAuditWork(other.m_heatAuditWork),
    m_budgetStarted(other.m_budgetStarted),
    m_operationControls(other.m_operationControls),
    m_manning(other.m_manning),
    m_diagnostics(other.m_diagnostics),
    m_external(other.m_external),
    m_transport(other.m_transport),
    m_sources(other.m_sources),
    m_structures(other.m_structures),
    m_groups(other.m_groups),
    m_transfers(other.m_transfers),
    m_settling(other.m_settling),
    m_reactions(other.m_reactions),
    m_light(other.m_light),
    m_solar(other.m_solar),
    m_turbulence(other.m_turbulence),
    m_meteorology(other.m_meteorology),
    m_useMeteorology(other.m_useMeteorology),
    m_solarAttenuation(other.m_solarAttenuation),
    m_incidentShortwave(other.m_incidentShortwave),
    m_windStressX(other.m_windStressX),
    m_windStressY(other.m_windStressY),
    m_shearSquared(other.m_shearSquared),
    m_buoyancyFrequency(other.m_buoyancyFrequency),
    m_lightOptions(other.m_lightOptions),
    m_profile(other.m_profile),
    m_pressure(other.m_pressure ? other.m_pressure->clone() : nullptr),
    m_time(other.m_time),
    m_temperatureIndex(other.m_temperatureIndex),
    m_salinityIndex(other.m_salinityIndex),
    m_accelerationX(other.m_accelerationX),
    m_accelerationY(other.m_accelerationY),
    m_columnAccelerationX(other.m_columnAccelerationX),
    m_columnAccelerationY(other.m_columnAccelerationY),
    m_km(other.m_km),
    m_kz(other.m_kz),
    m_cellDepth(other.m_cellDepth),
    m_cellPosition(other.m_cellPosition),
    m_edgeLayerFlux(other.m_edgeLayerFlux),
    m_verticalFlux(other.m_verticalFlux),
    m_oldVolume(other.m_oldVolume),
    m_workU(other.m_workU),
    m_workV(other.m_workV),
    m_scalarWork(other.m_scalarWork),
    m_scalarLoad(other.m_scalarLoad),
    m_diffusionGradX(other.m_scalarWork.size(), 0.0),
    m_diffusionGradY(other.m_scalarWork.size(), 0.0),
    m_diffusionLow(other.m_scalarWork.size(), 0.0),
    m_diffusionPlus(other.m_scalarWork.size(), 0.0),
    m_diffusionMinus(other.m_scalarWork.size(), 0.0),
    m_cellSource(other.m_cellSource),
    m_columnSource(other.m_columnSource),
    m_structureVolume(other.m_structureVolume),
    m_sourceValue(other.m_sourceValue),
    m_lower(other.m_lower),
    m_diagonal(other.m_diagonal),
    m_upper(other.m_upper),
    m_rhs(other.m_rhs),
    m_solution(other.m_solution),
    m_sweep(other.m_sweep)
{
    if (other.m_pressure && !m_pressure)
        throw std::runtime_error("pressure closure does not support checkpoint replay");
}

std::shared_ptr<const InternalModeCheckpoint> InternalModeSolver::checkpoint(std::string &message) const
{
    try {
        if (!m_checkpointReady || !m_mesh || !m_state)
            throw std::invalid_argument("checkpoint requires an initialized solver");
        auto saved = std::make_shared<InternalModeCheckpoint>();
        saved->owner = this;
        saved->generation = m_checkpointGeneration;
        saved->mesh = *m_mesh;
        saved->state = *m_state;
        saved->solver.reset(new InternalModeSolver(*this, CheckpointCopy{}));
        message.clear();
        return saved;
    } catch (const std::exception &error) {
        message = error.what();
        return {};
    }
}

bool InternalModeSolver::restoreCheckpoint(const InternalModeCheckpoint &saved, std::string &message)
{
    try {
        if (!m_checkpointReady || saved.owner != this ||
            saved.generation != m_checkpointGeneration || !saved.solver)
            throw std::invalid_argument("checkpoint belongs to another solver configuration");
        State::StateStore state = saved.state;
        Mesh::LayeredMesh mesh = saved.mesh;
        InternalModeSolver solver(*saved.solver, CheckpointCopy{});
        // Allocate everything before committing, and keep C API state views
        // allocation-stable. Configuration dimensions were fixed at initialize.
        for (auto pair : {std::pair{&m_state->eta, &state.eta},
                          std::pair{&m_state->ubar, &state.ubar},
                          std::pair{&m_state->vbar, &state.vbar},
                          std::pair{&m_state->u, &state.u},
                          std::pair{&m_state->v, &state.v},
                          std::pair{&m_state->density, &state.density},
                          std::pair{&m_state->w, &state.w},
                          std::pair{&m_state->kz, &state.kz},
                          std::pair{&m_state->km, &state.km},
                          std::pair{&m_state->species, &state.species}})
            if (pair.first->size() != pair.second->size())
                throw std::invalid_argument("checkpoint state dimensions changed");
        mesh.setGeometryEpochForRestore(m_mesh->geometryEpoch() + 1);
        *m_mesh = std::move(mesh);
        for (auto pair : {std::pair{&m_state->eta, &state.eta},
                          std::pair{&m_state->ubar, &state.ubar},
                          std::pair{&m_state->vbar, &state.vbar},
                          std::pair{&m_state->u, &state.u},
                          std::pair{&m_state->v, &state.v},
                          std::pair{&m_state->density, &state.density},
                          std::pair{&m_state->w, &state.w},
                          std::pair{&m_state->kz, &state.kz},
                          std::pair{&m_state->km, &state.km},
                          std::pair{&m_state->species, &state.species}})
            std::copy(pair.second->begin(), pair.second->end(), pair.first->begin());
        *this = std::move(solver);
        message.clear();
        return true;
    } catch (const std::exception &error) {
        message = error.what();
        return false;
    }
}
} // namespace FVQual::Hydro
