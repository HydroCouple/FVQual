// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "capi_internal.hpp"

#include "fvqual/capi/fvqual_options.h"
#include "fvqual/capi/fvqual_sigma.h"
#include "fvqual/version.h"

#include <cstring>
#include <string>

namespace FVQual::CApi {

namespace {

/*!
 * \brief The calling thread's most recent failure detail.
 */
std::string &detailBuffer()
{
    static thread_local std::string buffer;
    return buffer;
}

} // namespace

FVQ_ErrorCode fail(FVQ_ErrorCode code, const std::string &detail)
{
    detailBuffer() = detail;
    return code;
}

FVQ_ErrorCode succeed()
{
    detailBuffer().clear();
    return FVQ_OK;
}

const char *lastDetail()
{
    return detailBuffer().c_str();
}

bool toNative(const FVQ_SigmaSpec &spec, FVQual::Mesh::SigmaSpec &out,
              std::string &message)
{
    switch (spec.distribution)
    {
    case FVQ_SIGMA_UNIFORM:
        out.distribution = FVQual::Mesh::SigmaSpec::Distribution::Uniform;
        break;
    case FVQ_SIGMA_STRETCHED:
        out.distribution = FVQual::Mesh::SigmaSpec::Distribution::Stretched;
        break;
    case FVQ_SIGMA_EXPLICIT:
        out.distribution = FVQual::Mesh::SigmaSpec::Distribution::Explicit;
        break;
    default:
        message = "unrecognized sigma distribution";
        return false;
    }

    out.layerCount        = spec.layerCount;
    out.theta             = spec.theta;
    out.b                 = spec.blend;
    out.minLayerThickness = spec.minLayerThickness;
    out.explicitInterfaces.clear();

    if (spec.distribution == FVQ_SIGMA_EXPLICIT)
    {
        if (spec.interfaces == nullptr || spec.interfaceCount < 2)
        {
            message = "explicit layering needs at least two interfaces";
            return false;
        }

        out.explicitInterfaces.assign(spec.interfaces,
                                      spec.interfaces +
                                          static_cast<size_t>(spec.interfaceCount));
    }

    return true;
}

void toNative(const FVQ_Options &options, FVQual::Hydro::InternalModeOptions &out)
{
    out.gravity               = options.gravity;
    out.referenceDensity      = options.referenceDensity;
    out.coriolis              = options.coriolis;
    out.neutralViscosity      = options.neutralViscosity;
    out.backgroundViscosity   = options.backgroundViscosity;
    out.backgroundDiffusivity = options.backgroundDiffusivity;
    out.convectiveMixing      = options.convectiveMixing;
    out.horizontalViscosity   = options.horizontalViscosity;
    out.horizontalDiffusivity = options.horizontalDiffusivity;
    out.windStressX           = options.windStressX;
    out.windStressY           = options.windStressY;
    out.surfaceHeatFlux       = options.surfaceHeatFlux;
    out.heatCapacity          = options.heatCapacity;
    out.manningN              = options.manningN;
    out.cflNumber             = options.cflNumber;
    out.minTimeStep           = options.minTimeStep;
    out.maxTimeStep           = options.maxTimeStep;

    out.advectionScheme = options.advectionScheme == FVQ_ADVECTION_UPWIND
                              ? FVQual::Transport::AdvectionScheme::Upwind
                              : FVQual::Transport::AdvectionScheme::MusclFct;
}

} // namespace FVQual::CApi

using FVQual::CApi::fail;

extern "C" {

const char *fvq_error_message(FVQ_ErrorCode code)
{
    switch (code)
    {
    case FVQ_OK:                   return "success";
    case FVQ_ERR_NULL_ARGUMENT:    return "a required argument was null";
    case FVQ_ERR_INVALID_ARGUMENT: return "an argument was out of range or inconsistent";
    case FVQ_ERR_MESH:             return "mesh construction or validation failed";
    case FVQ_ERR_SIGMA:            return "the vertical layering specification is invalid";
    case FVQ_ERR_LIFECYCLE:        return "the call was made out of order";
    case FVQ_ERR_DRY_COLUMN:       return "a column is dry; the 3-D mode requires every column wet";
    case FVQ_ERR_SPECIES:          return "species registration failed";
    case FVQ_ERR_FILE:             return "a file could not be read or written";
    case FVQ_ERR_UNSUPPORTED:      return "this build lacks the required feature";
    case FVQ_ERR_ALLOCATION:       return "out of memory";
    case FVQ_ERR_INTERNAL:         return "an internal error occurred";
    }

    return "unrecognized error code";
}

const char *fvq_last_error(void)
{
    return FVQual::CApi::lastDetail();
}

void fvq_clear_error(void)
{
    (void)FVQual::CApi::succeed();
}

const char *fvq_version(void)
{
    return FVQual_FULL_VERSION;
}

void fvq_version_numbers(int *major, int *minor, int *patch)
{
    if (major != nullptr)
        *major = FVQual_VERSION_MAJOR;
    if (minor != nullptr)
        *minor = FVQual_VERSION_MINOR;
    if (patch != nullptr)
        *patch = FVQual_VERSION_PATCH;
}

/* ── Sigma ─────────────────────────────────────────────────────────────── */

void fvq_sigma_default(FVQ_SigmaSpec *spec)
{
    if (spec == nullptr)
        return;

    spec->distribution      = FVQ_SIGMA_UNIFORM;
    spec->layerCount        = 10;
    spec->theta             = 0.0;
    spec->blend             = 0.0;
    spec->minLayerThickness = 0.05;
    spec->interfaces        = nullptr;
    spec->interfaceCount    = 0;
}

FVQ_ErrorCode fvq_sigma_validate(const FVQ_SigmaSpec *spec)
{
    if (spec == nullptr)
        return fail(FVQ_ERR_NULL_ARGUMENT, "sigma specification is null");

    FVQual::Mesh::SigmaSpec native;
    std::string             message;

    if (!FVQual::CApi::toNative(*spec, native, message))
        return fail(FVQ_ERR_SIGMA, message);

    std::vector<double> interfaces;

    if (!native.resolve(interfaces, message))
        return fail(FVQ_ERR_SIGMA, message);

    return FVQual::CApi::succeed();
}

/* ── Options ───────────────────────────────────────────────────────────── */

void fvq_options_default(FVQ_Options *options)
{
    if (options == nullptr)
        return;

    const FVQual::Hydro::InternalModeOptions defaults;

    options->gravity               = defaults.gravity;
    options->referenceDensity      = defaults.referenceDensity;
    options->coriolis              = defaults.coriolis;
    options->neutralViscosity      = defaults.neutralViscosity;
    options->backgroundViscosity   = defaults.backgroundViscosity;
    options->backgroundDiffusivity = defaults.backgroundDiffusivity;
    options->convectiveMixing      = defaults.convectiveMixing;
    options->horizontalViscosity   = defaults.horizontalViscosity;
    options->horizontalDiffusivity = defaults.horizontalDiffusivity;
    options->advectionScheme       = FVQ_ADVECTION_MUSCL_FCT;
    options->windStressX           = defaults.windStressX;
    options->windStressY           = defaults.windStressY;
    options->surfaceHeatFlux       = defaults.surfaceHeatFlux;
    options->heatCapacity          = defaults.heatCapacity;
    options->manningN              = defaults.manningN;
    options->cflNumber             = defaults.cflNumber;
    options->minTimeStep           = defaults.minTimeStep;
    options->maxTimeStep           = defaults.maxTimeStep;
}

} // extern "C"
