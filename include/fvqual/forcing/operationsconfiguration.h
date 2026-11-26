// SPDX-License-Identifier: MIT
#ifndef FVQUAL_OPERATIONS_CONFIGURATION_H
#define FVQUAL_OPERATIONS_CONFIGURATION_H
#include "fvqual/fvqual_export.h"
#include "fvqual/forcing/forcingfield.h"
#include <nlohmann/json.hpp>
#include <string>
#include <vector>
namespace FVQual::Hydro { class InternalModeSolver; }
namespace FVQual::Mesh { class HorizontalMesh; }
namespace FVQual::State { class SpeciesRegistry; }
namespace FVQual::Forcing {
//! Piecewise constant controls on elapsed seconds. The solver splits at knots.
struct OperationControl {
    enum class Kind { Meteorology, SourceDischarge, SourceValue, StructureOpening,
        StructureTailwater, StructureDemand, GroupDemand, GroupTemperature,
        TransferDischarge, BoundaryDischarge, BoundaryStage, BoundaryValue };
    Kind kind{};
    int index = 0, species = -1;
    std::string speciesId;
    std::vector<double> times, values;
    std::vector<MeteorologyField> weather;
};
//! Validates the entire deck before installation; shared by component and C API.
FVQual_EXPORT bool configureOperations(Hydro::InternalModeSolver &solver,
    const Mesh::HorizontalMesh &mesh, const State::SpeciesRegistry &species,
    const nlohmann::json &configuration, std::string &message);
//! Published budget names/values. Known concentration units are converted to kg;
//! other species explicitly retain their declared concentration * m3 basis.
FVQual_EXPORT nlohmann::json processBudgetValues(const Hydro::InternalModeSolver &solver,
    const State::SpeciesRegistry &species);
}
#endif
