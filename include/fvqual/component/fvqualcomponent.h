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
 * **What it publishes.** Layered temperature, velocity, cell volume and
 * configured species, face water level/depth, and domain water/heat/species
 * inventories with cumulative accepted inflow totals. Layered fields share
 * the live sigma coordinate; heat is relative to 0 Celsius.
 *
 * **What it consumes.** Air temperature per column and integrated_inflows,
 * an acknowledged interval packet of water, heat and species already exported
 * by its donor. Live and recorded inflows are mutually exclusive; the donor
 * must be scheduled before this receiver.
 *
 * Mesh packages retain mixed topology, initial depths and roughness. A
 * recorded_inflows argument consumes integrated donor exports through a
 * conservative receiver. Temperature output carries its horizontal surface
 * and live vertical coordinates. Checkpoints include all solver history and
 * the transfer ledger; they are same-instance memory tokens.
 *
 * \threadsafety Not thread-safe; the framework drives one component serially.
 */

#ifndef FVQual_COMPONENT_FVQUALCOMPONENT_H
#define FVQual_COMPONENT_FVQUALCOMPONENT_H

#include "fvqual/component/sigmaverticalcoordinate.h"
#include "fvqual/component/recordedinflowargument.h"
#include "fvqual/component/operationsargument.h"
#include "hydrocouplesdk/spatial/polyhedralsurfaceargument.h"
#include "hydrocouplesdk/temporal/timeseriesargument.h"
#include "hydrocouplesdk/temporal/timedata.h"
#include "fvqual/kinetics/reactionfile.h"
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
#include "fvqual/component/reservoirargument.h"
#include "hydrocouplesdk/data/idbasedargument.h"
#include "hydrocouplesdk/data/exchangeitems1d.h"
#include "hydrocouplesdk/spatial/layeredmeshexchangeitems.h"
#include "hydrocouplesdk/spatial/meshadapters.h"

#include <memory>
#include <map>
#include "hydrocouplesdk/core/unit.h"
#include <optional>
#include <string>
#include <vector>

namespace FVQual::Component {

/*!
 * \brief FVQual behind the HydroCouple component contract.
 */
class FVQual_EXPORT FVQualComponent : public HydroCouple::SDK::AbstractModelComponent,
                                     public virtual HydroCouple::ICheckpointableModelComponent,
                                     public virtual HydroCouple::IDifferentiableModelComponent,
                                     public virtual HydroCouple::Temporal::ITimeModelComponent
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

    [[nodiscard]] std::set<HydroCouple::Capability> capabilities() const override;
    [[nodiscard]] bool saveState(std::string &token, std::string &message) override;
    [[nodiscard]] bool restoreState(const std::string &token, std::string &message) override;
    [[nodiscard]] bool releaseState(const std::string &token, std::string &message) override;
    [[nodiscard]] HydroCouple::Temporal::IDateTime *currentDateTime() const override;
    [[nodiscard]] HydroCouple::Temporal::ITimeSpan *simulationPeriod() const override;
    [[nodiscard]] double nextDateTimeJulianDay() const override;
    //! Opt-in conditional species derivatives; hydraulic/thermal inputs stay
    //! prescribed. Only explicitly listed species and kinetic constants qualify.
    [[nodiscard]] std::vector<HydroCouple::IInput *> differentiableInputs() const override { return {}; }
    [[nodiscard]] std::vector<HydroCouple::IArgument *> differentiableArguments() const override;
    [[nodiscard]] std::vector<HydroCouple::IOutput *> differentiableOutputs() const override;
    [[nodiscard]] std::vector<HydroCouple::IComponentDataItem *> differentiableStates() const override;
    [[nodiscard]] std::vector<HydroCouple::IComponentDataItem *> states() const override { return differentiableStates(); }
    [[nodiscard]] bool jvp(HydroCouple::DifferentialSet seeds,HydroCouple::DifferentialSet results,std::string *message=nullptr) override;
    [[nodiscard]] bool vjp(HydroCouple::DifferentialSet seeds,HydroCouple::DifferentialSet results,std::string *message=nullptr) override;

    // ── Read-only accessors, for gates and for embedding ──────────────────

    [[nodiscard]] double currentTime() const noexcept { return m_solver.time(); }
    [[nodiscard]] const Hydro::InternalModeSolver &solver() const noexcept { return m_solver; }
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

    //! The species registry, for gates and for embedding (as state() is).
    [[nodiscard]] const FVQual::State::SpeciesRegistry &registry() const noexcept
    {
        return m_registry;
    }

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
    [[nodiscard]] std::vector<std::string> configurationProblems() const;
    //! Reads the geometry arguments into the plain fields the rest uses.
    [[nodiscard]] bool readGeometryArguments(std::string &message);

    //! Copies the solver's temperature block into the layered output.
    void publishTemperature();
    bool applyLiveInflow(std::string &message);
    bool readReservoirArguments(std::string &message);
    bool registerConfiguredSpecies(std::string &message);
    bool configureDifferentiation(std::string &message);
    bool applyReactionParameters(std::string &message);
    bool differential(bool reverse,HydroCouple::DifferentialSet,HydroCouple::DifferentialSet,std::string *);

    // ── Configuration, from arguments ─────────────────────────────────────

    int    m_columnCount   = 10;
    int    m_layerCount    = 10;
    double m_length        = 1000.0;
    double m_width         = 200.0;
    double m_depth         = 20.0;
    double m_initialTemperature = 10.0;
    double m_couplingInterval   = 3600.0;
    double m_endTime            = 86400.0;
    double m_epochJulianDay     = 2440587.5; // Synthetic default: 1970-01-01 UTC.
    mutable HydroCouple::SDK::Temporal::TimeData m_clock{"CurrentTime", 2440587.5};
    mutable HydroCouple::SDK::Temporal::TimeSpan m_period{"SimulationPeriod", 2440587.5, 1};

    std::map<std::string, std::unique_ptr<ReservoirArgument<double>>> m_settings;
    std::unique_ptr<ReservoirArgument<int>> m_layerArgument;
    std::unique_ptr<ReservoirArgument<int>> m_differentiationArgument;
    std::unique_ptr<HydroCouple::SDK::IdBasedArgumentDouble> m_reactionParameters;
    std::unique_ptr<HydroCouple::SDK::ValueDefinition> m_reactionParameterDefinition;
    bool m_differentiationReady=false,m_linearized=false;
    std::vector<HydroCouple::IOutput *> m_differentiableFields;
    std::unique_ptr<HydroCouple::SDK::IdBasedArgumentDouble> m_speciesArgument;
    std::vector<std::unique_ptr<HydroCouple::SDK::Unit>> m_units;
    std::vector<std::unique_ptr<HydroCouple::SDK::Quantity>> m_quantities;
    std::vector<std::unique_ptr<HydroCouple::SDK::ValueDefinition>> m_speciesValueDefinitions;
    std::vector<std::unique_ptr<HydroCouple::SDK::Spatial::LayeredMeshOutputDouble>> m_fields;
    std::vector<std::unique_ptr<ColumnOutput>> m_columnFields;
    std::vector<std::unique_ptr<HydroCouple::SDK::Output1DDouble>> m_budgetFields;
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
    HydroCouple::SDK::Spatial::PolyhedralSurfaceArgument *m_meshArgument =
        nullptr;
    HydroCouple::SDK::Temporal::TimeSeriesArgumentDouble *m_metArgument =
        nullptr;
    HydroCouple::SDK::Argument1DString *m_kineticsFile = nullptr;

    FVQual::Kinetics::ReactionFile m_reactionFile;
    bool m_hasKinetics = false;
    HydroCouple::SDK::Input1DDouble    *m_airTemperature  = nullptr;
    HydroCouple::SDK::Spatial::LayeredMeshOutputDouble *m_temperature = nullptr;

    std::unique_ptr<HydroCouple::SDK::Quantity>  m_countQuantity;
    std::unique_ptr<HydroCouple::SDK::Quantity>  m_lengthQuantity;
    std::unique_ptr<HydroCouple::SDK::ValueDefinition> m_fileDefinition;
    std::unique_ptr<HydroCouple::SDK::Unit> m_temperatureUnit;
    std::unique_ptr<HydroCouple::SDK::Quantity>  m_temperatureQuantity;
    std::unique_ptr<HydroCouple::SDK::Argument1DInt>    m_ownedGeometryInts;
    std::unique_ptr<HydroCouple::SDK::Argument1DDouble> m_ownedGeometryDoubles;
    std::unique_ptr<HydroCouple::SDK::Spatial::PolyhedralSurfaceArgument>
        m_ownedMesh;
    std::unique_ptr<HydroCouple::SDK::Temporal::TimeSeriesArgumentDouble>
        m_ownedMet;
    std::unique_ptr<HydroCouple::SDK::Argument1DString> m_ownedKineticsFile;
    std::unique_ptr<RecordedInflowArgument> m_inflowArgument;
    std::vector<std::unique_ptr<OperationsArgument>> m_operations;
    bool m_validationProbe = false;
    nlohmann::json operationsConfiguration() const;
    bool validateOperationsEdit(const nlohmann::json &payload, std::string &message) const;
    std::unique_ptr<HydroCouple::SDK::Unit> m_epochUnit;
    std::unique_ptr<HydroCouple::SDK::Quantity> m_epochQuantity;
    std::unique_ptr<HydroCouple::SDK::Argument1DDouble> m_epochArgument;
    std::unique_ptr<Coupling::IntervalTransferReceiver> m_receiver;
    size_t m_nextFrame = 0;
    std::unique_ptr<HydroCouple::SDK::Input1DString> m_liveInflow;
    std::unique_ptr<HydroCouple::SDK::Output1DString> m_acceptedTransfer;
    std::string m_lastAcceptedTransfer;
    std::unique_ptr<HydroCouple::SDK::Input1DDouble>    m_ownedAirTemperature;
    std::unique_ptr<HydroCouple::SDK::Spatial::LayeredMeshOutputDouble> m_ownedTemperature;
    std::unique_ptr<HydroCouple::SDK::IO::MeshPackageReference> m_outputReference;
    std::unique_ptr<HydroCouple::SDK::Spatial::PolyhedralSurfaceAdapter> m_outputSurface;

    HydroCouple::SDK::Dimension m_columnDimension{"columns", "Column index"};
    HydroCouple::SDK::Dimension m_settingDimension{"settings", "Setting index"};
    HydroCouple::SDK::Dimension m_timeDimension{"time", "Time"};

    struct SavedState {
        std::shared_ptr<const FVQual::Hydro::InternalModeCheckpoint> solver;
        std::vector<double> airTemperature;
        std::unique_ptr<Coupling::IntervalTransferReceiver> receiver;
        size_t nextFrame = 0;
        std::string acceptedTransfer;
    };
    std::map<std::string, SavedState> m_checkpoints;
    bool m_prepared = false;
};

} // namespace FVQual::Component

#endif // FVQual_COMPONENT_FVQUALCOMPONENT_H
