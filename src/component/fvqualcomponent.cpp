// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/component/fvqualcomponent.h"

#include "hydrocouplesdk/io/meshdefinition.h"

#include <algorithm>
#include <cmath>

namespace FVQual::Component {

using HydroCouple::SDK::Argument1DDouble;
using HydroCouple::SDK::Argument1DInt;
using HydroCouple::SDK::Input1DDouble;
using HydroCouple::SDK::Quantity;
using HydroCouple::SDK::Spatial::LayeredMeshOutputDouble;

namespace {

// Argument slots. Named constants rather than bare subscripts, because an
// off-by-one into an argument array reads as a plausible configuration.
enum GeometryInt
{
    kColumns = 0,
    kLayers,
    kGeometryIntCount
};

enum GeometryDouble
{
    kLength = 0,
    kWidth,
    kDepth,
    kInitialTemperature,
    kCouplingInterval,
    kEndTime,
    kGeometryDoubleCount
};

//! A rectangular basin, one row of columns along its length.
HydroCouple::SDK::IO::MeshDefinition rectangularBasin(int columns, double length,
                                                      double width)
{
    HydroCouple::SDK::IO::MeshDefinition mesh;
    mesh.meshName = "mesh2d";

    const double dx = length / static_cast<double>(columns);

    for (int j = 0; j <= 1; ++j)
        for (int i = 0; i <= columns; ++i)
        {
            mesh.nodeX.push_back(static_cast<double>(i) * dx);
            mesh.nodeY.push_back(static_cast<double>(j) * width);
        }

    mesh.faceNodeOffsets.push_back(0);

    const auto node = [&](int i, int j) { return j * (columns + 1) + i; };

    for (int i = 0; i < columns; ++i)
    {
        mesh.faceNodes.push_back(node(i, 0));
        mesh.faceNodes.push_back(node(i + 1, 0));
        mesh.faceNodes.push_back(node(i + 1, 1));
        mesh.faceNodes.push_back(node(i, 1));
        mesh.faceNodeOffsets.push_back(static_cast<int64_t>(mesh.faceNodes.size()));
    }

    return mesh;
}

} // namespace

FVQualComponent::FVQualComponent(std::string_view id,
                                 HydroCouple::SDK::AbstractModelComponentInfo *info)
    : AbstractModelComponent(id, "FVQual", info)
{
}

// ── Arguments ─────────────────────────────────────────────────────────────

void FVQualComponent::createArguments()
{
    m_countQuantity.reset(Quantity::unitLess("Count"));
    m_lengthQuantity.reset(Quantity::unitLess("Geometry"));

    m_ownedGeometryInts = std::make_unique<Argument1DInt>(
        "geometry_counts", &m_settingDimension, kGeometryIntCount,
        m_countQuantity.get(), this);
    m_geometryInts = m_ownedGeometryInts.get();

    (*m_geometryInts)[kColumns] = m_columnCount;
    (*m_geometryInts)[kLayers]  = m_layerCount;
    addArgument(m_geometryInts);

    m_ownedGeometryDoubles = std::make_unique<Argument1DDouble>(
        "geometry_values", &m_settingDimension, kGeometryDoubleCount,
        m_lengthQuantity.get(), this);
    m_geometryDoubles = m_ownedGeometryDoubles.get();

    (*m_geometryDoubles)[kLength]             = m_length;
    (*m_geometryDoubles)[kWidth]              = m_width;
    (*m_geometryDoubles)[kDepth]              = m_depth;
    (*m_geometryDoubles)[kInitialTemperature] = m_initialTemperature;
    (*m_geometryDoubles)[kCouplingInterval]   = m_couplingInterval;
    (*m_geometryDoubles)[kEndTime]            = m_endTime;
    addArgument(m_geometryDoubles);

    // The mesh, when someone has one: a UGRID file, or a mesh generator's
    // live output through an @from binding. Optional on purpose — the
    // rectangular basin the geometry arguments describe remains the
    // default, so an unconfigured component is still a valid one.
    m_ownedMesh =
        std::make_unique<HydroCouple::SDK::Spatial::PolyhedralSurfaceArgument>(
            "mesh", this);
    m_meshArgument = m_ownedMesh.get();
    addArgument(m_meshArgument);

    // Whole-series meteorology, for a standalone run or an @from binding
    // from a timeseries provider. Times are DAYS on the model's own clock
    // (the W2 JDAY convention this lineage already speaks); a wired
    // air_temperature exchange input takes precedence, being the livelier
    // data.
    m_ownedMet = std::make_unique<HydroCouple::SDK::Temporal::TimeSeriesArgumentDouble>(
        "met_air_temperature", "Air temperature series", this,
        &m_timeDimension);
    m_metArgument = m_ownedMet.get();
    m_metArgument->setIsOptional(true);
    addArgument(m_metArgument);

    // The kinetics deck: species, constants, rate laws, settling, light.
    m_ownedKineticsFile = std::make_unique<HydroCouple::SDK::Argument1DString>(
        "kinetics_file", &m_settingDimension, 1, m_lengthQuantity.get(),
        this);
    m_kineticsFile = m_ownedKineticsFile.get();
    (*m_kineticsFile)[0] = std::string();
    m_kineticsFile->addFileFilter("Reaction Decks (*.rxn)");
    m_kineticsFile->setIsOptional(true);
    addArgument(m_kineticsFile);
}

bool FVQualComponent::readGeometryArguments(std::string &message)
{
    if (m_geometryInts == nullptr || m_geometryDoubles == nullptr)
    {
        message = "geometry arguments were not created";
        return false;
    }

    m_columnCount = (*m_geometryInts)[kColumns];
    m_layerCount  = (*m_geometryInts)[kLayers];

    m_length             = (*m_geometryDoubles)[kLength];
    m_width              = (*m_geometryDoubles)[kWidth];
    m_depth              = (*m_geometryDoubles)[kDepth];
    m_initialTemperature = (*m_geometryDoubles)[kInitialTemperature];
    m_couplingInterval   = (*m_geometryDoubles)[kCouplingInterval];
    m_endTime            = (*m_geometryDoubles)[kEndTime];

    return true;
}

bool FVQualComponent::initializeArguments(std::string &message)
{
    message.clear();

    if (!readGeometryArguments(message))
        return false;

    // A supplied mesh is the geometry: its faces are the columns, and the
    // rectangle the geometry arguments describe stands down.
    if (m_meshArgument != nullptr && m_meshArgument->hasMesh())
        m_columnCount =
            static_cast<int>(m_meshArgument->mesh().faceCount());

    // The kinetics deck is parsed here so a malformed file fails the
    // initialize with its line number, before anything is half-built.
    m_hasKinetics = false;
    const std::string kineticsPath =
        m_kineticsFile != nullptr ? (*m_kineticsFile)[0] : std::string();
    if (!kineticsPath.empty())
    {
        if (!m_reactionFile.read(kineticsPath, message))
            return false;
        m_hasKinetics = true;
    }

    // Every problem, not the first one. `validate()` reports the same set; this
    // is the gate that stops `initialize()` building a mesh it cannot use.
    const std::vector<std::string> problems = validate();

    if (!problems.empty())
    {
        message = problems.front();
        for (size_t i = 1; i < problems.size(); ++i)
            message += "; " + problems[i];
        return false;
    }

    return true;
}

// ── Validation ────────────────────────────────────────────────────────────

std::vector<std::string> FVQualComponent::validate()
{
    setStatus(ComponentStatus::Validating);

    std::vector<std::string> problems;

    if (m_columnCount < 1)
        problems.emplace_back("column count must be at least 1");
    if (m_layerCount < 1)
        problems.emplace_back("layer count must be at least 1");
    if (!(m_length > 0.0))
        problems.emplace_back("basin length must be positive");
    if (!(m_width > 0.0))
        problems.emplace_back("basin width must be positive");
    if (!(m_depth > 0.0))
        problems.emplace_back("basin depth must be positive");
    if (!(m_couplingInterval > 0.0))
        problems.emplace_back("coupling interval must be positive");
    if (!(m_endTime > 0.0))
        problems.emplace_back("end time must be positive");
    if (m_endTime > 0.0 && m_couplingInterval > m_endTime)
        problems.emplace_back("coupling interval exceeds the end time, so the run "
                              "would finish before its first exchange");

    setStatus(problems.empty() ? ComponentStatus::Valid : ComponentStatus::Invalid);
    return problems;
}

// ── Exchange items ────────────────────────────────────────────────────────

void FVQualComponent::createInputs()
{
    m_temperatureQuantity.reset(Quantity::unitLess("Temperature"));

    m_ownedAirTemperature = std::make_unique<Input1DDouble>(
        "air_temperature", &m_columnDimension, std::max(m_columnCount, 1),
        m_temperatureQuantity.get(), this);
    m_airTemperature = m_ownedAirTemperature.get();
    addInput(m_airTemperature);
}

void FVQualComponent::createOutputs()
{
    if (!m_temperatureQuantity)
        m_temperatureQuantity.reset(Quantity::unitLess("Temperature"));

    // The vertical coordinate has to exist before the output that publishes it,
    // and the mesh it reads is not built until prepare(). Binding it here to a
    // mesh that is still empty is safe because SigmaVerticalCoordinate holds a
    // reference and reads through it on every call — the same liveness the
    // layered item is designed around.
    if (!m_vertical)
        m_vertical.emplace(m_layered);

    m_ownedTemperature = std::make_unique<LayeredMeshOutputDouble>(
        "temperature", /*surface=*/nullptr, &*m_vertical,
        HydroCouple::Spatial::MeshDataObjectType::Cell,
        HydroCouple::Spatial::SpatialDataType::Scalar, std::max(m_columnCount, 1),
        std::max(m_layerCount, 1), m_temperatureQuantity.get(), this);
    m_temperature = m_ownedTemperature.get();
    addOutput(m_temperature);
}

// ── Lifecycle ─────────────────────────────────────────────────────────────

void FVQualComponent::prepare()
{
    setStatus(ComponentStatus::Preparing);

    std::string message;

    const bool suppliedMesh =
        m_meshArgument != nullptr && m_meshArgument->hasMesh();
    const HydroCouple::SDK::IO::MeshDefinition &definition =
        suppliedMesh ? m_meshArgument->mesh()
                     : rectangularBasin(m_columnCount, m_length, m_width);
    if (!m_horizontal.build(definition, message))
    {
        setStatus(ComponentStatus::Failed, "building the horizontal mesh: " + message);
        return;
    }

    // Bed elevation: a supplied mesh that carries node z speaks for
    // itself (per-face mean); everything else sits at the uniform -depth
    // the arguments describe.
    if (suppliedMesh && !definition.nodeZ.empty())
    {
        for (int64_t f = 0; f < m_horizontal.faceCount(); ++f)
        {
            const int64_t begin = definition.faceNodeOffsets[static_cast<size_t>(f)];
            const int64_t end   = definition.faceNodeOffsets[static_cast<size_t>(f) + 1];
            double sum = 0.0;
            for (int64_t k = begin; k < end; ++k)
                sum += definition.nodeZ[static_cast<size_t>(
                    definition.faceNodes[static_cast<size_t>(k)])];
            m_horizontal.faceBed()[static_cast<size_t>(f)] =
                sum / static_cast<double>(end - begin);
        }
    }
    else
    {
        for (int64_t f = 0; f < m_horizontal.faceCount(); ++f)
            m_horizontal.faceBed()[static_cast<size_t>(f)] = -m_depth;
    }

    m_sigma.layerCount = m_layerCount;

    if (!m_layered.build(m_horizontal, m_sigma, 0.0, message))
    {
        setStatus(ComponentStatus::Failed, "building the layered mesh: " + message);
        return;
    }

    m_temperatureIndex = m_registry.addTemperature(m_initialTemperature, message);

    if (m_temperatureIndex < 0)
    {
        setStatus(ComponentStatus::Failed, "registering temperature: " + message);
        return;
    }

    if (m_hasKinetics && !m_reactionFile.registerSpecies(m_registry, message))
    {
        setStatus(ComponentStatus::Failed,
                  "registering kinetics species: " + message);
        return;
    }

    m_state.allocate(m_layered, m_registry);

    if (!m_solver.initialize(m_layered, m_state, m_registry, m_options, message))
    {
        setStatus(ComponentStatus::Failed, "initialising the solver: " + message);
        return;
    }

    if (m_hasKinetics)
    {
        if (!m_reactionFile.configure(m_solver, message))
        {
            setStatus(ComponentStatus::Failed,
                      "configuring kinetics: " + message);
            return;
        }
        if (!m_solver.compileReactions(message))
        {
            setStatus(ComponentStatus::Failed,
                      "compiling kinetics: " + message);
            return;
        }
    }

    // A met series makes meteorology mean something. Only the argument
    // enables it: a wired input keeps the embedder-controlled behaviour it
    // has always had, so no existing composition changes underneath.
    if (m_metArgument != nullptr && m_metArgument->timeSeries() &&
        m_metArgument->timeSeries()->count() > 0 &&
        !m_solver.meteorologyEnabled())
        m_solver.enableMeteorology();

    // No rebinding of m_vertical here, deliberately. It was bound in
    // createOutputs() to `m_layered`, and `build()` mutates that object in
    // place rather than replacing it, so the coordinate is already describing
    // the mesh just built. An earlier version re-emplaced it defensively; an
    // ablation that deleted that line changed nothing, which is what showed the
    // line was dead. It is the same liveness the layered item is designed
    // around, arriving one level down.
    m_prepared = true;
    publishTemperature();

    setStatus(ComponentStatus::Updated, "prepared");
}

void FVQualComponent::update(const std::vector<HydroCouple::IOutput *> &requiredOutputs)
{
    if (!m_prepared)
    {
        setStatus(ComponentStatus::Failed, "update() called before prepare()");
        return;
    }

    setStatus(ComponentStatus::Updating);

    applyInputValues();

    // Air temperature. A wired exchange input carries spatial variation and
    // wins; the whole-series argument stands behind it for standalone runs,
    // applied uniformly at the model's own day count.
    const bool airInputWired =
        m_airTemperature != nullptr && m_airTemperature->provider() != nullptr;
    if (!airInputWired && m_solver.meteorologyEnabled() &&
        m_metArgument != nullptr && m_metArgument->timeSeries() &&
        m_metArgument->timeSeries()->count() > 0)
    {
        const double air = m_metArgument->getInterpolatedValue(
            m_solver.time() / 86400.0);
        m_solver.meteorology().airTemperature.setUniform(air);
    }
    else if (m_airTemperature != nullptr && m_solver.meteorologyEnabled())
    {
        std::vector<double> air(static_cast<size_t>(m_columnCount));
        for (int c = 0; c < m_columnCount; ++c)
            air[static_cast<size_t>(c)] = (*m_airTemperature)[c];

        std::string message;
        if (!m_solver.meteorology().airTemperature.setPerColumn(
                air, m_layered.columnCount(), message))
            reportError(HydroCouple::ErrorEntry::Severity::Warning, 0,
                        "air temperature input ignored: " + message);
    }

    std::string message;
    const double target = std::min(m_couplingInterval, m_endTime - m_solver.time());

    if (target > 0.0 && !m_solver.advance(target, message))
    {
        setStatus(ComponentStatus::Failed, "advancing the solver: " + message);
        return;
    }

    publishTemperature();
    updateOutputValues(requiredOutputs);

    const bool done = m_solver.time() >= m_endTime - 1.0e-9;
    setStatus(done ? ComponentStatus::Done : ComponentStatus::Updated, "stepped");
}

void FVQualComponent::finish()
{
    setStatus(ComponentStatus::Finishing);
    m_prepared = false;
    setStatus(ComponentStatus::Finished);
}

void FVQualComponent::initializeFailureCleanUp()
{
    // Leave nothing half-built. The framework may call initialize() again, and
    // a solver left initialised against a mesh that failed to build is the kind
    // of state that produces a plausible second run.
    m_prepared         = false;
    m_temperatureIndex = -1;
    m_vertical.reset();
}

// ── Publishing ────────────────────────────────────────────────────────────

void FVQualComponent::publishTemperature()
{
    if (m_temperature == nullptr || m_temperatureIndex < 0)
        return;

    const std::span<const double> block = m_state.speciesBlock(m_temperatureIndex);

    for (int c = 0; c < m_columnCount; ++c)
        for (int k = 0; k < m_layerCount; ++k)
        {
            const int64_t cell = m_layered.cell(c, k);

            if (cell >= 0 && cell < static_cast<int64_t>(block.size()))
                (*m_temperature)(c, k) = block[static_cast<size_t>(cell)];
        }
}

} // namespace FVQual::Component
