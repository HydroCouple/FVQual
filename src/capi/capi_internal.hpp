// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  capi_internal.hpp
 * \brief Handle definitions and error plumbing shared by the C API sources.
 *
 * Not installed: this is the private side of the ABI firewall. Callers see
 * incomplete types; only these translation units know what a handle contains.
 */

#ifndef FVQual_CAPI_INTERNAL_HPP
#define FVQual_CAPI_INTERNAL_HPP

#include "fvqual/capi/fvqual_engine.h"
#include "fvqual/capi/fvqual_options.h"
#include "fvqual/capi/fvqual_sigma.h"
#include "fvqual/hydro/internalmode.h"
#include "fvqual/mesh/horizontalmesh.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/state/speciesregistry.h"
#include "fvqual/state/statestore.h"

#include <string>
#include <map>
#include "fvqual/coupling/intervaltransfer.h"
#include "hydrocouplesdk/io/meshpackage.h"
#include <optional>
#include "fvqual/kinetics/reactionfile.h"

namespace FVQual::CApi {

/*!
 * \brief Records the detail message for the calling thread and returns \p code.
 *
 * Thread-local rather than per-handle because failures happen before a handle
 * exists — a rejected mesh has nowhere else to put its explanation.
 */
FVQ_ErrorCode fail(FVQ_ErrorCode code, const std::string &detail);

/*!
 * \brief Clears the thread-local detail and returns FVQ_OK.
 */
FVQ_ErrorCode succeed();

/*!
 * \brief The calling thread's detail message, valid until its next API call.
 */
const char *lastDetail();

/*!
 * \brief Translates a C sigma specification into the engine's own.
 * \return false with \p message set when the specification is malformed.
 */
bool toNative(const FVQ_SigmaSpec &spec, FVQual::Mesh::SigmaSpec &out,
              std::string &message);

/*!
 * \brief Translates C options into the engine's own.
 */
void toNative(const FVQ_Options &options, FVQual::Hydro::InternalModeOptions &out);

/*!
 * \brief Guards a C API body against exceptions escaping into C.
 *
 * An exception crossing the C boundary is undefined behaviour, so every entry
 * point that can allocate runs its body through this.
 */
template <typename Body>
FVQ_ErrorCode guard(Body &&body)
{
    try
    {
        return body();
    }
    catch (const std::bad_alloc &)
    {
        return fail(FVQ_ERR_ALLOCATION, "out of memory");
    }
    catch (const std::exception &error)
    {
        return fail(FVQ_ERR_INTERNAL, error.what());
    }
    catch (...)
    {
        return fail(FVQ_ERR_INTERNAL, "unknown error");
    }
}

} // namespace FVQual::CApi

/*!
 * \brief A horizontal mesh, owned independently of any simulation.
 */
struct FVQ_Mesh_
{
    FVQual::Mesh::HorizontalMesh mesh;
    std::optional<HydroCouple::SDK::IO::MeshPackage> package;
};

/*!
 * \brief A simulation and everything it owns.
 *
 * The horizontal mesh is a *copy* of the one passed at creation: LayeredMesh
 * holds a pointer into it, and coupling that lifetime to a caller-owned handle
 * would make destroying a mesh mid-run a use-after-free. Copying costs a few
 * arrays once and removes the failure mode entirely.
 *
 * Never move or reallocate this object after initialize(): the solver holds
 * pointers to the mesh and state members.
 */
struct FVQ_Simulation_
{
    FVQual::Mesh::HorizontalMesh    horizontal;
    FVQual::Mesh::SigmaSpec         sigma;
    FVQual::Mesh::LayeredMesh       layered;
    FVQual::State::SpeciesRegistry  registry;
    FVQual::State::StateStore       state;
    FVQual::Hydro::InternalModeSolver solver;
    FVQual::Hydro::InternalModeOptions options;

    std::vector<double> initialSurface;
    std::vector<double> manning;
    bool                initialized = false;
    nlohmann::json operations = nlohmann::json::object();
    std::optional<FVQual::Kinetics::ReactionFile> reactionDeck;
    std::unique_ptr<FVQual::Coupling::IntervalTransferReceiver> receiver;
    struct Checkpoint {
        std::shared_ptr<const FVQual::Hydro::InternalModeCheckpoint> solver;
        std::unique_ptr<FVQual::Coupling::IntervalTransferReceiver> receiver;
    };
    std::map<uint64_t, Checkpoint> checkpoints;
};

#endif // FVQual_CAPI_INTERNAL_HPP
