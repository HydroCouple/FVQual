// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "capi_internal.hpp"

#include "fvqual/capi/fvqual_sim.h"

#include <memory>
#include <string>

using FVQual::CApi::fail;
using FVQual::CApi::guard;
using FVQual::CApi::succeed;

namespace {

/*!
 * \brief Rejects calls that arrive before or after initialize(), as required.
 */
FVQ_ErrorCode requireState(const FVQ_Simulation *sim, bool wantInitialized,
                           const char *what)
{
    if (sim == nullptr)
        return fail(FVQ_ERR_NULL_ARGUMENT, "simulation handle is null");

    if (sim->initialized != wantInitialized)
    {
        return fail(FVQ_ERR_LIFECYCLE,
                    std::string(what) + (wantInitialized
                                             ? " requires an initialized simulation"
                                             : " must be called before initialize"));
    }

    return FVQ_OK;
}

} // namespace

extern "C" {

FVQ_ErrorCode fvq_sim_create(const FVQ_Mesh *mesh, const FVQ_SigmaSpec *sigma,
                             const FVQ_Options *options, FVQ_Simulation **out)
{
    return guard([&]() -> FVQ_ErrorCode {
        if (mesh == nullptr || sigma == nullptr || out == nullptr)
            return fail(FVQ_ERR_NULL_ARGUMENT, "mesh, sigma, and output are required");

        std::string message;

        auto handle = std::make_unique<FVQ_Simulation_>();

        // A private copy: the caller may destroy their mesh immediately, and
        // the layered mesh holds a pointer into this one for the whole run.
        handle->horizontal = mesh->mesh;

        if (!FVQual::CApi::toNative(*sigma, handle->sigma, message))
            return fail(FVQ_ERR_SIGMA, message);

        FVQ_Options resolved;
        if (options != nullptr)
            resolved = *options;
        else
            fvq_options_default(&resolved);

        FVQual::CApi::toNative(resolved, handle->options);

        handle->initialSurface.assign(
            static_cast<size_t>(handle->horizontal.faceCount()), 0.0);

        *out = handle.release();
        return succeed();
    });
}

void fvq_sim_destroy(FVQ_Simulation *sim)
{
    delete sim;
}

FVQ_ErrorCode fvq_sim_add_species(FVQ_Simulation *sim, const char *id, const char *units,
                                  double initialValue, int *index)
{
    return guard([&]() -> FVQ_ErrorCode {
        const FVQ_ErrorCode state = requireState(sim, false, "adding a species");
        if (state != FVQ_OK)
            return state;

        if (id == nullptr)
            return fail(FVQ_ERR_NULL_ARGUMENT, "species id is null");

        FVQual::State::SpeciesDefinition definition;
        definition.id           = id;
        definition.units        = units == nullptr ? "" : units;
        definition.initialValue = initialValue;

        std::string message;
        const int   assigned = sim->registry.add(definition, message);

        if (assigned < 0)
            return fail(FVQ_ERR_SPECIES, message);

        if (index != nullptr)
            *index = assigned;

        return succeed();
    });
}

FVQ_ErrorCode fvq_sim_add_temperature(FVQ_Simulation *sim, double initialCelsius,
                                      int *index)
{
    return guard([&]() -> FVQ_ErrorCode {
        const FVQ_ErrorCode state = requireState(sim, false, "adding temperature");
        if (state != FVQ_OK)
            return state;

        std::string message;
        const int   assigned = sim->registry.addTemperature(initialCelsius, message);

        if (assigned < 0)
            return fail(FVQ_ERR_SPECIES, message);

        if (index != nullptr)
            *index = assigned;

        return succeed();
    });
}

FVQ_ErrorCode fvq_sim_add_age(FVQ_Simulation *sim, int *index)
{
    return guard([&]() -> FVQ_ErrorCode {
        const FVQ_ErrorCode state = requireState(sim, false, "adding water age");
        if (state != FVQ_OK)
            return state;

        std::string message;
        const int   assigned = sim->registry.addAge(message);

        if (assigned < 0)
            return fail(FVQ_ERR_SPECIES, message);

        if (index != nullptr)
            *index = assigned;

        return succeed();
    });
}

FVQ_ErrorCode fvq_sim_set_initial_surface(FVQ_Simulation *sim, const double *elevation,
                                          int64_t count)
{
    return guard([&]() -> FVQ_ErrorCode {
        const FVQ_ErrorCode state =
            requireState(sim, false, "setting the initial surface");
        if (state != FVQ_OK)
            return state;

        if (elevation == nullptr)
            return fail(FVQ_ERR_NULL_ARGUMENT, "elevation array is null");

        if (count != sim->horizontal.faceCount())
            return fail(FVQ_ERR_INVALID_ARGUMENT,
                        "surface length must equal the column count");

        sim->initialSurface.assign(elevation, elevation + count);
        return succeed();
    });
}

FVQ_ErrorCode fvq_sim_initialize(FVQ_Simulation *sim)
{
    return guard([&]() -> FVQ_ErrorCode {
        const FVQ_ErrorCode state = requireState(sim, false, "initialize");
        if (state != FVQ_OK)
            return state;

        std::string message;

        if (!sim->layered.build(sim->horizontal, sim->sigma, sim->initialSurface, message))
        {
            // A surface at or below the bed is the common cause, and it is
            // worth its own code so a binding can say "dry column" rather than
            // "mesh error".
            const bool dry = message.find("bed") != std::string::npos;
            return fail(dry ? FVQ_ERR_DRY_COLUMN : FVQ_ERR_MESH, message);
        }

        sim->state.allocate(sim->layered, sim->registry);

        if (!sim->solver.initialize(sim->layered, sim->state, sim->registry, sim->options,
                                    message))
        {
            const bool dry = message.find("dry") != std::string::npos;
            return fail(dry ? FVQ_ERR_DRY_COLUMN : FVQ_ERR_INVALID_ARGUMENT, message);
        }

        sim->initialized = true;
        return succeed();
    });
}

double fvq_sim_stable_time_step(const FVQ_Simulation *sim)
{
    if (sim == nullptr || !sim->initialized)
        return 0.0;

    return sim->solver.stableTimeStep();
}

FVQ_ErrorCode fvq_sim_step(FVQ_Simulation *sim, double dt)
{
    return guard([&]() -> FVQ_ErrorCode {
        const FVQ_ErrorCode state = requireState(sim, true, "step");
        if (state != FVQ_OK)
            return state;

        std::string message;

        if (!sim->solver.step(dt, message))
            return fail(FVQ_ERR_INVALID_ARGUMENT, message);

        return succeed();
    });
}

FVQ_ErrorCode fvq_sim_advance(FVQ_Simulation *sim, double duration)
{
    return guard([&]() -> FVQ_ErrorCode {
        const FVQ_ErrorCode state = requireState(sim, true, "advance");
        if (state != FVQ_OK)
            return state;

        std::string message;

        if (!sim->solver.advance(duration, message))
            return fail(FVQ_ERR_INVALID_ARGUMENT, message);

        return succeed();
    });
}

double fvq_sim_time(const FVQ_Simulation *sim)
{
    return sim == nullptr ? 0.0 : sim->solver.time();
}

int64_t fvq_sim_column_count(const FVQ_Simulation *sim)
{
    return sim == nullptr ? 0 : sim->horizontal.faceCount();
}

int fvq_sim_layer_count(const FVQ_Simulation *sim)
{
    return sim == nullptr ? 0 : sim->sigma.resolvedLayerCount();
}

int64_t fvq_sim_cell_count(const FVQ_Simulation *sim)
{
    if (sim == nullptr || !sim->initialized)
        return 0;

    return sim->layered.cellCount();
}

int fvq_sim_species_count(const FVQ_Simulation *sim)
{
    return sim == nullptr ? 0 : sim->registry.count();
}

FVQ_ErrorCode fvq_sim_diagnostics(const FVQ_Simulation *sim, FVQ_Diagnostics *out)
{
    if (sim == nullptr || out == nullptr)
        return fail(FVQ_ERR_NULL_ARGUMENT, "simulation and output must not be null");

    if (!sim->initialized)
        return fail(FVQ_ERR_LIFECYCLE, "diagnostics require an initialized simulation");

    const FVQual::Hydro::InternalModeDiagnostics &source = sim->solver.diagnostics();

    out->timeStep              = source.timeStep;
    out->maxHorizontalVelocity = source.maxHorizontalVelocity;
    out->maxVerticalFlux       = source.maxVerticalFlux;
    out->surfaceFluxResidual   = source.surfaceFluxResidual;
    out->totalVolume           = source.totalVolume;
    out->maxRichardson         = source.maxRichardson;
    out->externalSteps         = source.externalSteps;
    out->steps                 = source.steps;

    return succeed();
}

double fvq_sim_total_volume(const FVQ_Simulation *sim)
{
    if (sim == nullptr || !sim->initialized)
        return 0.0;

    return sim->layered.totalVolume();
}

FVQ_ErrorCode fvq_sim_species_mass(const FVQ_Simulation *sim, int species, double *out)
{
    if (sim == nullptr || out == nullptr)
        return fail(FVQ_ERR_NULL_ARGUMENT, "simulation and output must not be null");

    if (!sim->initialized)
        return fail(FVQ_ERR_LIFECYCLE, "species mass requires an initialized simulation");

    if (species < 0 || species >= sim->registry.count())
        return fail(FVQ_ERR_INVALID_ARGUMENT, "species index is out of range");

    *out = sim->solver.totalSpeciesMass(species);
    return succeed();
}

} // extern "C"
