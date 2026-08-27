// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   fvqualcomponent.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  FVQualComponent — FVQual as a HydroCouple model component (F8.1b).
 *
 * Wraps the mesh, state and internal-mode solver behind
 * `HydroCouple::SDK::AbstractModelComponent`, so a stratified reservoir can be
 * coupled to a watershed, a network model, or a water-quality engine through
 * the framework rather than through a file.
 *
 * **What it publishes.** Temperature on `{column, layer}`, as a
 * `LayeredMeshOutput`, carrying the vertical coordinate so a partner can find
 * out where those layers actually are. That is the whole point of the F8.2
 * interfaces: without the coordinate a consumer has an array it cannot place in
 * space, and under a sigma grid layer 7 moves every step.
 *
 * **What it consumes.** Air temperature per column, as a 1-D input driving the
 * surface heat budget. One input is enough to exercise the consuming side
 * honestly; the remaining meteorological variables follow the same pattern and
 * are not written until something asks for them.
 *
 * **Scope, stated plainly.** The mesh is built from geometry arguments, not read
 * from a bathymetry file. Reading W2 or UGRID bathymetry through an argument is
 * a natural extension and deliberately not done here: it would make every gate
 * depend on an external dataset, and six FVQual gates have already been caught
 * silently skipping for exactly that reason. `W2Bathymetry` and `UgridReader`
 * remain the supported paths into a real mesh, through the applications.
 *
 * \threadsafety Not thread-safe; the framework drives one component serially.
 */

#ifndef FVQual_COMPONENT_FVQUALCOMPONENT_H
#define FVQual_COMPONENT_FVQUALCOMPONENT_H

#include "fvqual/component/sigmaverticalcoordinate.h"
#include "fvqual/fvqual_export.h"
#include "fvqual/hydro/internalmode.h"
#include "fvqual/mesh/horizontalmesh.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/mesh/sigmaspec.h"
#include "fvqual/state/speciesregistry.h"
#include "fvqual/state/statestore.h"

#include "hydrocouplesdk/component/abstractmodelcomponent.h"
#include "hydrocouplesdk/core/valuedefinition.h"
#include "hydrocouplesdk/data/argument1d.h"
#include "hydrocouplesdk/data/exchangeitems1d.h"
#include "hydrocouplesdk/spatial/layeredmeshexchangeitems.h"

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace FVQual::Component {

/*!
 * \brief FVQual behind the HydroCouple component contract.
 */
class FVQual_EXPORT FVQualComponent : public HydroCouple::SDK::AbstractModelComponent
{
    using ComponentStatus = HydroCouple::IModelComponent::ComponentStatus;

public:
    /*!
     * \brief Construct.
     * \param[in] id   Component identifier.
     * \param[in] info Owning component info (borrowed, may be null).
     */
    explicit FVQualComponent(
        std::string_view id,
        HydroCouple::SDK::AbstractModelComponentInfo *info = nullptr);

    ~FVQualComponent() override = default;

    // ── IModelComponent ───────────────────────────────────────────────────

    /*!
     * \brief Reports every configuration problem, not just the first.
     *
     * A caller fixing a configuration one message per run is a caller who runs
     * five times to learn five things it could have been told once.
     */
    std::vector<std::string> validate() override;

    //! Builds the mesh, allocates state, and initialises the solver.
    void prepare() override;

    //! Advances one coupling interval and republishes the outputs.
    void update(const std::vector<HydroCouple::IOutput *> &requiredOutputs = {}) override;

    void finish() override;

    // ── Read-only accessors, for gates and for embedding ──────────────────

    [[nodiscard]] double currentTime() const noexcept { return m_solver.time(); }
    [[nodiscard]] const FVQual::Mesh::LayeredMesh &mesh() const noexcept
    {
        return m_layered;
    }

    /*!
     * \brief The model's own state, which is what the outputs are published
     *        from.
     *
     * Mutable, so an embedder can set an initial condition richer than the
     * single uniform temperature the arguments express, and so a gate can
     * compare what the component publishes against the source it publishes
     * from. That comparison is only worth anything on a field that varies:
     * against a uniform one, a publish loop that scrambled its indices reads
     * the same number everywhere and passes.
     */
    [[nodiscard]] FVQual::State::StateStore &state() noexcept { return m_state; }
    [[nodiscard]] const FVQual::State::StateStore &state() const noexcept
    {
        return m_state;
    }

    //! Index of temperature in the species registry, or -1 before prepare().
    [[nodiscard]] int temperatureIndex() const noexcept { return m_temperatureIndex; }

protected:
    void createArguments() override;
    bool initializeArguments(std::string &message) override;
    void createInputs() override;
    void createOutputs() override;

    /*!
     * \brief Returns the component to a state where `initialize()` can be run
     *        again rather than leaving a half-built model behind.
     */
    void initializeFailureCleanUp() override;

private:
    //! Reads the geometry arguments into the plain fields the rest uses.
    [[nodiscard]] bool readGeometryArguments(std::string &message);

    //! Copies the solver's temperature block into the layered output.
    void publishTemperature();

    // ── Configuration, from arguments ─────────────────────────────────────

    int    m_columnCount   = 10;
    int    m_layerCount    = 10;
    double m_length        = 1000.0;
    double m_width         = 200.0;
    double m_depth         = 20.0;
    double m_initialTemperature = 10.0;
    double m_couplingInterval   = 3600.0;
    double m_endTime            = 86400.0;

    // ── The model ─────────────────────────────────────────────────────────

    FVQual::Mesh::HorizontalMesh      m_horizontal;
    FVQual::Mesh::SigmaSpec           m_sigma;
    FVQual::Mesh::LayeredMesh         m_layered;
    FVQual::State::SpeciesRegistry    m_registry;
    FVQual::State::StateStore         m_state;
    FVQual::Hydro::InternalModeSolver m_solver;
    FVQual::Hydro::InternalModeOptions m_options;

    std::optional<SigmaVerticalCoordinate> m_vertical;

    int m_temperatureIndex = -1;

    // ── Framework objects. Raw, because the base stores raw pointers and
    //    re-runs createInputs()/createOutputs() on re-initialize; the unique_ptrs
    //    below own them across that cycle.
    HydroCouple::SDK::Argument1DInt    *m_geometryInts    = nullptr;
    HydroCouple::SDK::Argument1DDouble *m_geometryDoubles = nullptr;
    HydroCouple::SDK::Input1DDouble    *m_airTemperature  = nullptr;
    HydroCouple::SDK::Spatial::LayeredMeshOutputDouble *m_temperature = nullptr;

    std::unique_ptr<HydroCouple::SDK::Quantity>  m_countQuantity;
    std::unique_ptr<HydroCouple::SDK::Quantity>  m_lengthQuantity;
    std::unique_ptr<HydroCouple::SDK::Quantity>  m_temperatureQuantity;
    std::unique_ptr<HydroCouple::SDK::Argument1DInt>    m_ownedGeometryInts;
    std::unique_ptr<HydroCouple::SDK::Argument1DDouble> m_ownedGeometryDoubles;
    std::unique_ptr<HydroCouple::SDK::Input1DDouble>    m_ownedAirTemperature;
    std::unique_ptr<HydroCouple::SDK::Spatial::LayeredMeshOutputDouble> m_ownedTemperature;

    HydroCouple::SDK::Dimension m_columnDimension{"columns", "Column index"};
    HydroCouple::SDK::Dimension m_settingDimension{"settings", "Setting index"};

    bool m_prepared = false;
};

} // namespace FVQual::Component

#endif // FVQual_COMPONENT_FVQUALCOMPONENT_H
