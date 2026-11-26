// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/component/fvqualcomponent.h"

#include "hydrocouplesdk/io/meshdefinition.h"
#include "hydrocouplesdk/io/uriresolver.h"

#include <algorithm>
#include <cmath>
#include <atomic>
#include <limits>
#include <filesystem>
#include <set>
#include <sstream>

namespace FVQual::Component {

using HydroCouple::SDK::Argument1DDouble;
using HydroCouple::SDK::Argument1DInt;
using HydroCouple::SDK::Input1DDouble;
using HydroCouple::SDK::Quantity;
using HydroCouple::SDK::Spatial::LayeredMeshOutputDouble;

namespace {

// A deck reference names a reaction file, not a JSON array stored in that file.
// Keep it as @uri so composition loading resolves paths before initialize().
class ReactionDeckArgument final : public HydroCouple::SDK::Argument1DString {
public:
    using Base=HydroCouple::SDK::Argument1DString;
    using Base::Base;using Base::initialize;
    bool initialize(const std::string &value,HydroCouple::IArgument::ArgumentInputType type,std::string &message) override {
        using AT=HydroCouple::IArgument::ArgumentInputType;
        if(type==AT::File || type==AT::URL){(*this)[0]=value;setCurrentArgumentInputType(type);message.clear();return true;}
        if(type==AT::JSON) {
            try {
                const auto j=nlohmann::json::parse(value);std::string uri;
                if(HydroCouple::SDK::IO::jsonUriReference(value,uri)){(*this)[0]=uri;setCurrentArgumentInputType(type);message.clear();return true;}
                const auto &v=j.is_object()?j.at("values"):j;
                if(!v.is_array() || v.size()!=1 || !v[0].is_string()){message="kinetics_file requires one reaction-deck reference";return false;}
            }catch(const std::exception &e){message=e.what();return false;}
        }
        return Base::initialize(value,type,message);
    }
    bool serialize(HydroCouple::IArgument::ArgumentInputType type,std::string &text,std::string &message) const override {
        if(type==HydroCouple::IArgument::ArgumentInputType::JSON && !(*this)[0].empty()) {
            text=nlohmann::json{{"@uri",(*this)[0]}}.dump();message.clear();return true;
        }
        return Base::serialize(type,text,message);
    }
};

// Packets are read explicitly inside the receiver transaction. A provider read
// is observational: physical donor advancement belongs to the host schedule.
class IntegratedInflowInput final : public HydroCouple::SDK::Input1DString {
public:
    using HydroCouple::SDK::Input1DString::Input1DString;
    void retrieveValuesFromProvider() override {}
    bool canConsume(HydroCouple::IOutput *p, std::string &message) const override {
        if (!p || p->dataKind()!=HydroCouple::DataKind::String || p->shape()!=std::vector<int64_t>{1} ||
            (!dynamic_cast<HydroCouple::Temporal::ITimeModelComponent *>(p->modelComponent()) || dynamic_cast<HydroCouple::IAdaptedOutput *>(p))) {
            message="Integrated inflows need one String packet from a component declaring physical time";
            return false;
        }
        message.clear(); return true;
    }
};

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
    m_fileDefinition = std::make_unique<HydroCouple::SDK::ValueDefinition>("File", typeid(std::string));
    m_temperatureUnit = std::make_unique<HydroCouple::SDK::Unit>("degC", 1.0, 273.15);
    m_temperatureUnit->dimensionsInternal()->setPower(
        HydroCouple::IUnitDimensions::FundamentalUnitDimension::Temperature, 1.0);
    m_temperatureQuantity = std::make_unique<Quantity>(
        "Temperature", typeid(double), m_temperatureUnit.get());
    m_temperatureQuantity->setValueKind(HydroCouple::ValueKind::Intensive);

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
        &m_timeDimension, m_temperatureQuantity.get());
    m_metArgument = m_ownedMet.get();
    m_metArgument->setIsOptional(true);
    addArgument(m_metArgument);

    // The kinetics deck: species, constants, rate laws, settling, light.
    m_ownedKineticsFile = std::make_unique<ReactionDeckArgument>(
        "kinetics_file", &m_settingDimension, 1, m_fileDefinition.get(),
        this);
    m_kineticsFile = m_ownedKineticsFile.get();
    (*m_kineticsFile)[0] = std::string();
    m_kineticsFile->addFileFilter("Reaction Decks (*.rxn)");
    m_kineticsFile->setIsOptional(true);
    addArgument(m_kineticsFile);
    m_epochUnit = std::make_unique<HydroCouple::SDK::Unit>("day", 86400.0, 0.0);
    m_epochUnit->dimensionsInternal()->setPower(
        HydroCouple::IUnitDimensions::FundamentalUnitDimension::Time, 1.0);
    m_epochQuantity = std::make_unique<Quantity>("UTC epoch as Julian day", typeid(double), m_epochUnit.get());
    m_epochArgument = std::make_unique<Argument1DDouble>(
        "epoch_julian_day", &m_settingDimension, 1, m_epochQuantity.get(), this);
    (*m_epochArgument)[0] = m_epochJulianDay;
    m_epochArgument->setIsOptional(true);
    addArgument(m_epochArgument.get());
    using D = HydroCouple::IUnitDimensions::FundamentalUnitDimension;
    auto setting = [&](const char *id, const char *caption, double initial,
                       const char *units, double length, double time, double temperature = 0) {
        auto unit = std::make_unique<HydroCouple::SDK::Unit>(units,1.0,temperature ? 273.15 : 0.0);
        unit->dimensionsInternal()->setPower(D::Length,length);
        unit->dimensionsInternal()->setPower(D::Time,time);
        unit->dimensionsInternal()->setPower(D::Temperature,temperature);
        auto q = std::make_unique<Quantity>(caption,typeid(double),unit.get());
        auto a = std::make_unique<ReservoirArgument<double>>(id,&m_settingDimension,1,q.get(),this);
        a->setCaption(caption); a->setIsOptional(true); (*a)[0] = initial;
        addArgument(a.get()); m_settings.emplace(id,std::move(a));
        m_units.push_back(std::move(unit)); m_quantities.push_back(std::move(q));
    };
    m_layerArgument = std::make_unique<ReservoirArgument<int>>(
        "layer_count",&m_settingDimension,1,m_countQuantity.get(),this);
    m_layerArgument->setCaption("Number of vertical layers");
    m_layerArgument->setIsOptional(true); (*m_layerArgument)[0] = m_layerCount;
    addArgument(m_layerArgument.get());
    m_differentiationArgument=std::make_unique<ReservoirArgument<int>>(
        "conditional_species_derivatives",&m_settingDimension,1,m_countQuantity.get(),this);
    m_differentiationArgument->setCaption("Conditional species derivatives (0 off, 1 upwind)");
    m_differentiationArgument->setDescription("Advanced: 1 selects upwind transport and enables species/kinetic derivatives conditional on prescribed hydraulics, mixing, geometry and inflows. Thermal/momentum and SWMM gradients are excluded.");
    m_differentiationArgument->setIsOptional(true);(*m_differentiationArgument)[0]=0;addArgument(m_differentiationArgument.get());
    m_reactionParameterDefinition=std::make_unique<HydroCouple::SDK::ValueDefinition>("Kinetic constants in reaction-deck units",typeid(double));
    m_reactionParameters=std::make_unique<HydroCouple::SDK::IdBasedArgumentDouble>(
        "reaction_constants",std::vector<std::string>{},&m_settingDimension,m_reactionParameterDefinition.get(),this);
    m_reactionParameters->setCaption("Kinetic constant overrides (reaction-deck units)");
    m_reactionParameters->setIsOptional(true);addArgument(m_reactionParameters.get());
    setting("duration_s","Run duration",m_endTime,"s",0,1);
    setting("coupling_interval_s","Exchange interval",m_couplingInterval,"s",0,1);
    setting("initial_surface_elevation_m","Initial water surface: one value or one per face",0,"m",1,0);
    setting("initial_temperature_c","Initial temperature: uniform, surface-to-bed profile or face-major cells",10,"degC",0,0,1);
    setting("initial_velocity_x_m_s","Initial x velocity: uniform or face-major cells",0,"m/s",1,-1);
    setting("initial_velocity_y_m_s","Initial y velocity: uniform or face-major cells",0,"m/s",1,-1);
    setting("horizontal_diffusivity_m2_s","Horizontal scalar diffusivity",0,"m2/s",2,-1);
    setting("maximum_time_step_s","Maximum internal time step",300,"s",0,1);
    setting("sigma_interfaces","Optional normalized interfaces: surface 0 to bed 1",0,"1",0,0);
    auto massUnit = std::make_unique<HydroCouple::SDK::Unit>("kg/m3",1,0);
    massUnit->dimensionsInternal()->setPower(D::Mass,1);
    massUnit->dimensionsInternal()->setPower(D::Length,-3);
    auto massQuantity = std::make_unique<Quantity>("Conservative species initial concentrations",typeid(double),massUnit.get());
    massQuantity->setValueKind(HydroCouple::ValueKind::Intensive);
    m_speciesArgument = std::make_unique<HydroCouple::SDK::IdBasedArgumentDouble>(
        "conservative_species",std::vector<std::string>{},&m_settingDimension,massQuantity.get(),this);
    m_speciesArgument->setCaption("Conservative species: identifier and initial concentration");
    m_speciesArgument->setIsOptional(true); addArgument(m_speciesArgument.get());
    m_units.push_back(std::move(massUnit)); m_quantities.push_back(std::move(massQuantity));
    m_inflowArgument = std::make_unique<RecordedInflowArgument>(
        &m_settingDimension, m_fileDefinition.get(), this);
    addArgument(m_inflowArgument.get());
    for (const std::string section : {"meteorology", "sources", "structures", "outlet_groups", "open_boundaries", "internal_transfers", "reaction_deck"}) {
        auto arg = std::make_unique<OperationsArgument>(section, &m_settingDimension, m_fileDefinition.get(), this);
        arg->validator = [this](const nlohmann::json &payload, std::string &message) {
            return m_validationProbe || status() == ComponentStatus::Created || validateOperationsEdit(payload, message);
        };
        addArgument(arg.get()); m_operations.push_back(std::move(arg));
    }

}

bool FVQualComponent::readGeometryArguments(std::string &message)
{
    if (m_geometryInts == nullptr || m_geometryDoubles == nullptr)
    {
        message = "geometry arguments were not created";
        return false;
    }

    if (m_geometryInts->shape()!=std::vector<int64_t>{kGeometryIntCount} ||
        m_geometryDoubles->shape()!=std::vector<int64_t>{kGeometryDoubleCount}) {
        message="geometry_counts needs two values and geometry_values needs six";return false;
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
    if(m_differentiationArgument->shape()!=std::vector<int64_t>{1} ||
        ((*m_differentiationArgument)[0]!=0 && (*m_differentiationArgument)[0]!=1)) {
        message="conditional_species_derivatives needs one value: 0 or 1";return false;
    }
    if (m_epochArgument->shape()!=std::vector<int64_t>{1}) {message="epoch_julian_day requires one value";return false;}
    m_epochJulianDay = (*m_epochArgument)[0];
    if (m_inflowArgument && !m_inflowArgument->load(message)) return false;

    // A supplied mesh is the geometry: its faces are the columns, and the
    // rectangle the geometry arguments describe stands down.
    if (m_meshArgument != nullptr && m_meshArgument->hasMesh())
    {
        if (m_meshArgument->mesh().faceCount() > std::numeric_limits<int>::max()) {
            message = "mesh exceeds FVQual component column index range";
            return false;
        }
        m_columnCount = static_cast<int>(m_meshArgument->mesh().faceCount());
    }

    if (m_meshArgument && m_meshArgument->meshPackage()) {
        const auto &package = *m_meshArgument->meshPackage();
        if (package.elevationKind != "bed-elevation") {
            message = "FVQual reservoir mesh requires declared bed elevations, not terrain";
            return false;
        }
        const auto operations = operationsConfiguration();
        for (const auto &boundary : package.boundaries) {
            if (boundary.type == "WALL") continue;
            bool mapped = false;
            for (const auto &row : operations.at("open_boundaries"))
                for (const auto &pair : row.at("node_pairs"))
                    if (pair == nlohmann::json::array({boundary.node0, boundary.node1}) ||
                        pair == nlohmann::json::array({boundary.node1, boundary.node0})) mapped = true;
            if (!mapped) { message = "imported boundary '" + boundary.type + "' needs an explicit FVQual boundary mapping"; return false; }
        }
        // Native infiltration/coupling rows describe SWMM ownership; preserving
        // them is provenance, not permission to activate the same physics twice.
    }

    // The kinetics deck is parsed here so a malformed file fails the
    // initialize with its line number, before anything is half-built.
    m_hasKinetics = false;
    std::string kineticsPath =
        m_kineticsFile != nullptr ? (*m_kineticsFile)[0] : std::string();
    if (!kineticsPath.empty())
    {
        if (HydroCouple::SDK::IO::uriScheme(kineticsPath).empty() && std::filesystem::path(kineticsPath).is_relative())
            kineticsPath = (std::filesystem::path(referenceDirectory()) / kineticsPath).string();
        std::string text;
        if(!HydroCouple::SDK::IO::readReference(kineticsPath,text,message))return false;
        std::istringstream stream(text);
        if (!m_reactionFile.parse(stream,kineticsPath, message))
            return false;
        m_hasKinetics = true;
    }

    for (const auto &arg : m_operations) if (arg->id() == "reaction_deck" && !arg->value().get<std::string>().empty()) {
        if (m_hasKinetics) { message = "choose an inline reaction deck or kinetics_file, not both"; return false; }
        std::istringstream input(arg->value().get<std::string>());
        if (!m_reactionFile.parse(input, "reaction_deck", message)) return false;
        m_hasKinetics = true;
    }

    if (!readReservoirArguments(message) || !registerConfiguredSpecies(message)) return false;

    // Every problem, not the first one. `validate()` reports the same set; this
    // is the gate that stops `initialize()` building a mesh it cannot use.
    const std::vector<std::string> problems = configurationProblems();

    if (!problems.empty())
    {
        message = problems.front();
        for (size_t i = 1; i < problems.size(); ++i)
            message += "; " + problems[i];
        return false;
    }

    if (!m_validationProbe) for (const auto &arg : m_operations) if (!arg->value().empty()) {
        std::string payload;
        if (!arg->serialize(HydroCouple::IArgument::ArgumentInputType::JSON, payload, message) ||
            !validateOperationsEdit(nlohmann::json::parse(payload), message)) return false;
        break;
    }
    return true;
}

nlohmann::json FVQualComponent::operationsConfiguration() const
{
    nlohmann::json result = nlohmann::json::object();
    for (const auto &arg : m_operations) if (arg->id() != "reaction_deck") result[arg->id()] = arg->value();
    return result;
}

bool FVQualComponent::validateOperationsEdit(const nlohmann::json &payload, std::string &message) const
{
    FVQualComponent probe("operation-validation"); probe.m_validationProbe = true;
    using AT = HydroCouple::IArgument::ArgumentInputType;
    probe.setReferenceDirectory(referenceDirectory());
    const auto probeArguments = probe.arguments();
    for (auto *arg : arguments()) {
        if (auto *override = dynamic_cast<ReservoirArgument<double> *>(arg); override && !override->supplied) continue;
        if (auto *override = dynamic_cast<ReservoirArgument<int> *>(arg); override && !override->supplied) continue;
        auto found = std::find_if(probeArguments.begin(), probeArguments.end(), [&](auto *candidate) { return candidate->id() == arg->id(); });
        if (found == probeArguments.end()) { message = "missing probe argument"; return false; }
        std::string data;
        if (!arg->serialize(AT::JSON, data, message)) {
            if (arg == m_meshArgument) {
                if (!m_meshArgument->hasMesh()) continue;
                if (!(*found)->initialize(*arg, message)) return false;
                continue;
            }
            return false;
        }
        if (arg->id() == payload.at("section").get<std::string>()) data = payload.dump();
        if (!(*found)->initialize(data, AT::JSON, message)) return false;
    }
    probe.setReferenceDirectory(referenceDirectory());
    probe.initialize();
    if (probe.status() == ComponentStatus::Failed) { message = probe.errors().empty() ? "operation configuration failed" : probe.errors().back().message; return false; }
    probe.prepare();
    if (probe.status() == ComponentStatus::Failed) { message = probe.errors().empty() ? "operation configuration failed" : probe.errors().back().message; return false; }
    message.clear(); return true;
}

bool FVQualComponent::readReservoirArguments(std::string &message)
{
    if (m_layerArgument->supplied) {
        if (m_layerArgument->shape() != std::vector<int64_t>{1}) { message="layer_count requires one value"; return false; }
        m_layerCount=(*m_layerArgument)[0];
    } else (*m_layerArgument)[0]=m_layerCount;
    auto scalar = [&](const char *id,double &value) {
        auto &a=*m_settings.at(id);
        if (!a.supplied) { a[0]=value; return true; }
        if (a.shape()!=std::vector<int64_t>{1} || !std::isfinite(a[0])) { message=std::string(id)+" requires one finite value"; return false; }
        value=a[0]; return true;
    };
    if (!scalar("duration_s",m_endTime) || !scalar("coupling_interval_s",m_couplingInterval) ||
        !scalar("horizontal_diffusivity_m2_s",m_options.horizontalDiffusivity) ||
        !scalar("maximum_time_step_s",m_options.maxTimeStep)) return false;
    auto &sigma=*m_settings.at("sigma_interfaces");
    m_sigma = Mesh::SigmaSpec{}; m_sigma.layerCount=m_layerCount;
    if (sigma.supplied) {
        m_sigma.distribution=Mesh::SigmaSpec::Distribution::Explicit;
        for (int64_t i=0;i<sigma.shape()[0];++i) m_sigma.explicitInterfaces.push_back(sigma[i]);
        std::vector<double> resolved;
        if (!m_sigma.resolve(resolved,message)) return false;
        if (m_layerArgument->supplied && m_layerCount!=m_sigma.resolvedLayerCount()) { message="layer_count disagrees with sigma_interfaces"; return false; }
        m_layerCount=m_sigma.resolvedLayerCount();
    }
    for (const char *id : {"initial_surface_elevation_m","initial_temperature_c","initial_velocity_x_m_s","initial_velocity_y_m_s"}) {
        const auto &a=*m_settings.at(id);
        if (!a.supplied) continue;
        const auto n=a.shape()[0];
        const bool surface=std::string(id)=="initial_surface_elevation_m";
        if (n!=1 && n!=(surface ? m_columnCount : int64_t(m_columnCount)*m_layerCount) &&
            !(std::string(id)=="initial_temperature_c" && n==m_layerCount)) {
            message=std::string(id)+" has an incompatible shape"; return false;
        }
        for (int64_t i=0;i<n;++i) if (!std::isfinite(a[i])) { message=std::string(id)+" must be finite"; return false; }
    }
    if (!m_settings.at("initial_temperature_c")->supplied) (*m_settings.at("initial_temperature_c"))[0]=m_initialTemperature;
    if (m_options.horizontalDiffusivity<0 || m_options.maxTimeStep<m_options.minTimeStep) {
        message="diffusivity must be nonnegative and maximum time step must exceed the minimum"; return false;
    }
    return true;
}

bool FVQualComponent::registerConfiguredSpecies(std::string &message)
{
    m_registry=State::SpeciesRegistry{};
    m_temperatureIndex=m_registry.addTemperature(m_initialTemperature,message);
    if (m_temperatureIndex<0) return false;
    if (m_hasKinetics && !m_reactionFile.registerSpecies(m_registry,message)) return false;
    if (m_inflowArgument) for (const auto &d:m_inflowArgument->species())
        if (m_registry.add(d,message)<0) return false;
    const auto ids=m_speciesArgument->identifiers();
    const std::set<std::string> outputIds={"temperature","eta","water_depth","velocity_x","velocity_y","cell_volume","water_volume","heat_inventory","accepted_water","accepted_heat"};
    for (size_t i=0;i<ids.size();++i) {
        if (outputIds.count(ids[i]) || ids[i].starts_with("__") || (ids[i].starts_with("inventory_") || ids[i].starts_with("accepted_") || ids[i].starts_with("budget_")) ||
            ids[i].find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")!=std::string::npos ||
            !std::isfinite((*m_speciesArgument)[i]) || (*m_speciesArgument)[i]<0) {
            message="invalid conservative species identifier or initial concentration: "+ids[i]; return false;
        }
        State::SpeciesDefinition d; d.id=ids[i]; d.units="kg/m3"; d.initialValue=(*m_speciesArgument)[i];
        if (m_registry.add(d,message)<0) return false;
    }
    for (int i=1;i<m_registry.count();++i) {
        const auto &id=m_registry.definition(i).id;
        if (outputIds.count(id) || id.starts_with("inventory_") || id.starts_with("accepted_") || id.starts_with("budget_") || id.starts_with("__hc_") ||
            id.find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")!=std::string::npos) {
            message="species identifier is unsafe or collides with a reservoir output: "+id;return false;
        }
    }
    return true;
}

// ── Validation ────────────────────────────────────────────────────────────

std::vector<std::string> FVQualComponent::validate()
{
    setStatus(ComponentStatus::Validating);
    auto problems = configurationProblems();
    setStatus(problems.empty() ? ComponentStatus::Valid : ComponentStatus::Invalid);
    return problems;
}

std::vector<std::string> FVQualComponent::configurationProblems() const
{
    std::vector<std::string> problems;
    if (!std::isfinite(m_epochJulianDay) ||
        !std::isfinite(m_epochJulianDay + m_endTime / 86400.0))
        problems.emplace_back("model epoch and horizon must be finite Julian days");

    if (!std::isfinite(m_initialTemperature)) problems.emplace_back("initial temperature must be finite");
    if (m_columnCount < 1)
        problems.emplace_back("column count must be at least 1");
    if (m_layerCount < 1)
        problems.emplace_back("layer count must be at least 1");
    if (!std::isfinite(m_length) || !(m_length > 0.0))
        problems.emplace_back("basin length must be positive");
    if (!std::isfinite(m_width) || !(m_width > 0.0))
        problems.emplace_back("basin width must be positive");
    if (!std::isfinite(m_depth) || !(m_depth > 0.0))
        problems.emplace_back("basin depth must be positive");
    if (!std::isfinite(m_couplingInterval) || !(m_couplingInterval > 0.0))
        problems.emplace_back("coupling interval must be positive");
    if (!std::isfinite(m_endTime) || !(m_endTime > 0.0))
        problems.emplace_back("end time must be positive");
    if (m_endTime > 0.0 && m_couplingInterval > m_endTime)
        problems.emplace_back("coupling interval exceeds the end time, so the run "
                              "would finish before its first exchange");
    if (m_inflowArgument && !m_inflowArgument->frames().empty()) {
        if (std::abs(m_inflowArgument->frames().back().endSeconds - m_endTime) > 1.0e-8)
            problems.emplace_back("inflow series horizon must equal the configured end time");
        for (const auto &frame : m_inflowArgument->frames())
            for (const auto &inflow : frame.inflows)
                if (inflow.column >= m_columnCount)
                    problems.emplace_back("inflow column is outside the configured mesh");
    }

    if (m_liveInflow && m_liveInflow->provider()) {
        if (m_inflowArgument && !m_inflowArgument->frames().empty())
            problems.emplace_back("Choose live integrated inflows or recorded inflows, not both");
        auto *donor=dynamic_cast<HydroCouple::Temporal::ITimeModelComponent *>(m_liveInflow->provider()->modelComponent());
        if (!donor || std::abs(donor->simulationPeriod()->julianDay()-m_epochJulianDay)>1.e-10 ||
            std::abs(donor->simulationPeriod()->duration()*86400.0-m_endTime)>1.e-6)
            problems.emplace_back("Donor and receiver must have the same physical epoch and exchange horizon");
    }
    return problems;
}

// ── Exchange items ────────────────────────────────────────────────────────

void FVQualComponent::createInputs()
{
    m_ownedAirTemperature = std::make_unique<Input1DDouble>(
        "air_temperature", &m_columnDimension, std::max(m_columnCount, 1),
        m_temperatureQuantity.get(), this);
    m_airTemperature = m_ownedAirTemperature.get();
    addInput(m_airTemperature);
    m_liveInflow=std::make_unique<IntegratedInflowInput>("integrated_inflows", &m_settingDimension, 1, m_fileDefinition.get(), this);
    m_liveInflow->setCaption("Live water, heat and species transfers");
    m_liveInflow->setDescription("hydrocouple.integrated-inflow/1.0; already-debited amounts, elapsed seconds, m3, J relative to 0 Celsius and kg. Schedule donor first.");
    addInput(m_liveInflow.get());
}

void FVQualComponent::createOutputs()
{
    // createInputs() has established the physical Celsius quantity.

    // The vertical coordinate has to exist before the output that publishes it,
    // and the mesh it reads is not built until prepare(). Binding it here to a
    // mesh that is still empty is safe because SigmaVerticalCoordinate holds a
    // reference and reads through it on every call — the same liveness the
    // layered item is designed around.
    if (!m_vertical)
        m_vertical.emplace(m_layered);

    auto definition = m_meshArgument && m_meshArgument->hasMesh()
        ? m_meshArgument->mesh() : rectangularBasin(m_columnCount, m_length, m_width);
    if (definition.nodeZ.empty()) definition.nodeZ.assign(definition.nodeX.size(), -m_depth);
    m_outputReference.reset();
    if (m_meshArgument && m_meshArgument->meshPackage())
        m_outputReference = std::make_unique<HydroCouple::SDK::IO::MeshPackageReference>(
            *m_meshArgument->meshPackage());
    else if (!m_meshArgument || !m_meshArgument->hasMesh()) {
        HydroCouple::SDK::IO::MeshPackage reference;
        reference.horizontalCrs="LOCAL"; reference.verticalDatum="FVQual synthetic rectangular bed";
        reference.elevationKind="bed-elevation";
        m_outputReference = std::make_unique<HydroCouple::SDK::IO::MeshPackageReference>(reference);
    }
    m_outputSurface = std::make_unique<HydroCouple::SDK::Spatial::PolyhedralSurfaceAdapter>(
        HydroCouple::SDK::Spatial::polyhedralSurfaceFromMesh(definition), definition, m_outputReference.get());
    m_ownedTemperature = std::make_unique<LayeredMeshOutputDouble>(
        "temperature", m_outputSurface.get(), &*m_vertical,
        HydroCouple::Spatial::MeshLocation::Face,
        HydroCouple::Spatial::SpatialDataType::Scalar, std::max(m_columnCount, 1),
        std::max(m_layerCount, 1), m_temperatureQuantity.get(), this);
    m_temperature = m_ownedTemperature.get();
    addOutput(m_temperature);
    m_fields.clear(); m_columnFields.clear(); m_budgetFields.clear();
    using D=HydroCouple::IUnitDimensions::FundamentalUnitDimension;
    auto quantity=[&](const std::string &caption,const std::string &units,double length,double mass,double time,HydroCouple::ValueKind kind,double scale=1.0) {
        auto u=std::make_unique<HydroCouple::SDK::Unit>(units,scale,0);
        u->dimensionsInternal()->setPower(D::Length,length); u->dimensionsInternal()->setPower(D::Mass,mass); u->dimensionsInternal()->setPower(D::Time,time);
        auto q=std::make_unique<Quantity>(caption,typeid(double),u.get()); q->setValueKind(kind);
        auto *result=q.get(); m_units.push_back(std::move(u)); m_quantities.push_back(std::move(q)); return result;
    };
    auto *length=quantity("Elevation / water depth","m",1,0,0,HydroCouple::ValueKind::Intensive);
    auto *volume=quantity("Water volume","m3",3,0,0,HydroCouple::ValueKind::Extensive);
    auto *speed=quantity("Velocity along mesh reference axes","m/s",1,0,-1,HydroCouple::ValueKind::Intensive);
    auto *heat=quantity("Heat inventory relative to 0 Celsius","J",2,1,-2,HydroCouple::ValueKind::Extensive);
    auto *mass=quantity("Species inventory","kg",0,1,0,HydroCouple::ValueKind::Extensive);
    auto field=[&](const std::string &id,HydroCouple::SDK::ValueDefinition *q) {
        auto o=std::make_unique<LayeredMeshOutputDouble>(id,m_outputSurface.get(),&*m_vertical,
            HydroCouple::Spatial::MeshLocation::Face,HydroCouple::Spatial::SpatialDataType::Scalar,m_columnCount,m_layerCount,q,this);
        addOutput(o.get()); m_fields.push_back(std::move(o));
    };
    auto budget=[&](const std::string &id,HydroCouple::SDK::ValueDefinition *q) {
        auto o=std::make_unique<HydroCouple::SDK::Output1DDouble>(id,&m_settingDimension,1,q,this);
        addOutput(o.get()); m_budgetFields.push_back(std::move(o));
    };
    for (const char *id:{"eta","water_depth"}) {
        auto o=std::make_unique<ColumnOutput>(id,&m_columnDimension,m_columnCount,length,this,m_outputSurface.get());
        addOutput(o.get()); m_columnFields.push_back(std::move(o));
    }
    field("velocity_x",speed); field("velocity_y",speed); field("cell_volume",volume);
    for (int i=1;i<m_registry.count();++i) {
        const auto &d=m_registry.definition(i);
        const double concentrationScale=d.units=="kg/m3" ? 1.0 :
            (d.units=="mg/L" || d.units=="g/m3") ? .001 :
            (d.units=="ug/L" || d.units=="mg/m3") ? .000001 : 0.0;
        HydroCouple::SDK::ValueDefinition *q=nullptr;
        if (concentrationScale>0) {
            q=quantity(d.id,d.units,-3,1,0,HydroCouple::ValueKind::Intensive,concentrationScale);
            budget("inventory_"+d.id,mass);budget("accepted_"+d.id,mass);
        } else if (d.kind==State::SpeciesKind::Age && d.units=="days")
            q=quantity(d.id,"days",0,0,1,HydroCouple::ValueKind::Intensive,86400);
        else {
            // Preserve a deck's label without inventing a dimensional conversion.
            auto definition=std::make_unique<HydroCouple::SDK::ValueDefinition>(d.id+" ("+d.units+"; dimensions undeclared)",typeid(double));
            q=definition.get();m_speciesValueDefinitions.push_back(std::move(definition));
        }
        field(d.id,q);
        if (!concentrationScale && d.kind!=State::SpeciesKind::Age)
            m_fields.back()->setCaption(d.id+" ("+d.units+"; dimensions undeclared)");
    }
    budget("water_volume",volume); budget("heat_inventory",heat);
    budget("accepted_water",volume); budget("accepted_heat",heat);
    for (const std::string family : {"water", "heat"}) {
        auto *q = family == "water" ? volume : heat;
        for (const std::string term : {"initial", "residual"}) budget("budget_" + family + "_" + term, q);
        for (const std::string origin : {"source", "structure", "boundary", "coupled", "transfer"})
            for (const std::string direction : {"in", "out"}) budget("budget_" + family + "_" + origin + "_" + direction, q);
    }
    for (const std::string term : {"settling", "reaction", "shortwave", "longwave_in", "back_radiation", "latent", "sensible", "prescribed"})
        budget("budget_heat_" + term, heat);
    for (int i = 1; i < m_registry.count(); ++i) {
        const auto &d = m_registry.definition(i);
        const bool known = d.units == "kg/m3" || d.units == "kg/m2" || d.units == "mg/L" || d.units == "g/m3" || d.units == "g/m2" || d.units == "ug/L" || d.units == "mg/m3";
        HydroCouple::SDK::ValueDefinition *q = mass;
        if (!known) {
            auto definition = std::make_unique<HydroCouple::SDK::ValueDefinition>(d.id + " inventory (" + d.units + (d.kind == State::SpeciesKind::MsxBenthic ? " * m2)" : " * m3)"), typeid(double));
            q = definition.get(); m_speciesValueDefinitions.push_back(std::move(definition));
        }
        for (const std::string term : {"initial", "residual", "settling", "reaction"}) budget("budget_" + d.id + "_" + term, q);
        for (const std::string origin : {"source", "structure", "boundary", "coupled", "transfer"})
            for (const std::string direction : {"in", "out"}) budget("budget_" + d.id + "_" + origin + "_" + direction, q);
    }

    m_acceptedTransfer=std::make_unique<HydroCouple::SDK::Output1DString>("accepted_transfer_id", &m_settingDimension, 1, m_fileDefinition.get(), this);
    m_acceptedTransfer->setCaption("Last accepted donor transfer receipt");
    addOutput(m_acceptedTransfer.get());
}

// ── Lifecycle ─────────────────────────────────────────────────────────────

void FVQualComponent::prepare()
{
    setStatus(ComponentStatus::Preparing);
    m_prepared = false;
    m_differentiationReady=false;m_linearized=false;m_differentiableFields.clear();
    m_checkpoints.clear();
    m_solver.reset();
    m_receiver.reset();
    m_nextFrame = 0;
    m_lastAcceptedTransfer.clear();

    std::string message;
    if (m_liveInflow && m_liveInflow->provider()) {
        const auto problems=configurationProblems();
        if (!problems.empty()) { setStatus(ComponentStatus::Failed, problems.front()); return; }
    }

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

    std::vector<double> initialSurface(static_cast<size_t>(m_columnCount), 0.0);
    if (m_meshArgument && m_meshArgument->meshPackage() &&
        !m_meshArgument->meshPackage()->initialDepth.empty()) {
        const auto &depth = m_meshArgument->meshPackage()->initialDepth;
        for (int c = 0; c < m_columnCount; ++c)
            initialSurface[static_cast<size_t>(c)] =
                m_horizontal.faceBed()[static_cast<size_t>(c)] + depth[static_cast<size_t>(c)];
    }
    const auto &surface=*m_settings.at("initial_surface_elevation_m");
    if (surface.supplied) for (int c=0;c<m_columnCount;++c) initialSurface[c]=surface[surface.shape()[0]==1 ? 0 : c];
    for (int c=0;c<m_columnCount;++c) if (!(initialSurface[c]>m_horizontal.faceBed()[c])) {
        setStatus(ComponentStatus::Failed,"initial water surface must lie above bed in every face"); return;
    }
    if (!m_layered.build(m_horizontal, m_sigma, initialSurface, message))
    {
        setStatus(ComponentStatus::Failed, "building the layered mesh: " + message);
        return;
    }

    if (!registerConfiguredSpecies(message)) { setStatus(ComponentStatus::Failed,message); return; }
    m_state.allocate(m_layered, m_registry);
    for (int c=0;c<m_columnCount;++c) for (int k=0;k<m_layerCount;++k) {
        const auto cell=m_layered.cell(c,k);
        auto value=[&](const char *id,double fallback,bool profile=false) {
            const auto &a=*m_settings.at(id); if (!a.supplied) return fallback;
            const auto n=a.shape()[0]; return a[n==1 ? 0 : profile && n==m_layerCount ? k : cell];
        };
        m_state.speciesAt(m_temperatureIndex,cell)=value("initial_temperature_c",m_initialTemperature,true);
        m_state.u[cell]=value("initial_velocity_x_m_s",0); m_state.v[cell]=value("initial_velocity_y_m_s",0);
    }

    auto solverOptions=m_options;
    if((*m_differentiationArgument)[0])solverOptions.advectionScheme=Transport::AdvectionScheme::Upwind;
    if (!m_solver.initialize(m_layered, m_state, m_registry, solverOptions, message))
    {
        setStatus(ComponentStatus::Failed, "initialising the solver: " + message);
        return;
    }

    const auto operations = operationsConfiguration();
    if (!operations.at("meteorology").empty() &&
        ((m_metArgument && m_metArgument->timeSeries() && m_metArgument->timeSeries()->count() > 0) ||
         (m_airTemperature && m_airTemperature->provider()))) {
        setStatus(ComponentStatus::Failed, "Use either meteorology operations or the air-temperature series/input, so the weather source is unambiguous."); return;
    }
    if (!Forcing::configureOperations(m_solver, m_horizontal, m_registry, operations, message)) {
        setStatus(ComponentStatus::Failed, "configuring reservoir operations: " + message); return;
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

    if (m_meshArgument && m_meshArgument->meshPackage() &&
        !m_meshArgument->meshPackage()->roughness.empty() &&
        !m_solver.setManning(m_meshArgument->meshPackage()->roughness, message)) {
        setStatus(ComponentStatus::Failed, "applying mesh roughness: " + message);
        return;
    }
    if(!configureDifferentiation(message)){setStatus(ComponentStatus::Failed,"configuring derivatives: "+message);return;}

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
    m_receiver = std::make_unique<Coupling::IntervalTransferReceiver>(m_solver, m_registry, m_layered);
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
    m_linearized=false;
    std::string derivativeMessage;
    if(!applyReactionParameters(derivativeMessage) ||
       (m_differentiationReady && !m_solver.beginReactiveTrace(10000,derivativeMessage))) {
        setStatus(ComponentStatus::Failed,derivativeMessage);return;
    }

    applyInputValues();

    // Air temperature. A wired exchange input carries spatial variation and
    // wins; the whole-series argument stands behind it for standalone runs,
    // applied uniformly at the model's own day count.
    const bool airInputWired =
        m_airTemperature != nullptr && m_airTemperature->provider() != nullptr;
    const bool weatherOperations = !operationsConfiguration().at("meteorology").empty();
    if (!weatherOperations && !airInputWired && m_solver.meteorologyEnabled() &&
        m_metArgument != nullptr && m_metArgument->timeSeries() &&
        m_metArgument->timeSeries()->count() > 0)
    {
        const double air = m_metArgument->getInterpolatedValue(
            m_solver.time() / 86400.0);
        m_solver.meteorology().airTemperature.setUniform(air);
    }
    else if (!weatherOperations && m_airTemperature != nullptr && m_solver.meteorologyEnabled())
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

    bool advanced = true;
    if (m_liveInflow && m_liveInflow->provider()) {
        advanced=applyLiveInflow(message);
    } else if (m_inflowArgument && !m_inflowArgument->frames().empty()) {
        if (m_nextFrame < m_inflowArgument->frames().size()) {
            Coupling::TransferReceipt receipt;
            advanced = m_receiver->apply(m_inflowArgument->frames()[m_nextFrame], receipt, message);
            if (advanced) { m_lastAcceptedTransfer=m_inflowArgument->frames()[m_nextFrame].id; ++m_nextFrame; }
        }
    } else if (target > 0.0) advanced = m_solver.advance(target, message);
    if (!advanced)
    {
        m_solver.endReactiveTrace();
        setStatus(ComponentStatus::Failed, "advancing the solver: " + message);
        return;
    }
    if(m_differentiationReady){m_solver.endReactiveTrace();m_linearized=m_solver.reactiveTrace().count>0;}

    publishTemperature();
    updateOutputValues(requiredOutputs);

    const bool done = m_solver.time() >= m_endTime - 1.0e-9;
    setStatus(done ? ComponentStatus::Done : ComponentStatus::Updated, "stepped");
}

void FVQualComponent::finish()
{
    setStatus(ComponentStatus::Finishing);
    m_prepared = false;
    m_differentiationReady=false;m_linearized=false;
    setStatus(ComponentStatus::Finished);
}

void FVQualComponent::initializeFailureCleanUp()
{
    m_differentiationReady=false;m_linearized=false;m_differentiableFields.clear();
    m_checkpoints.clear();
    // Leave nothing half-built. The framework may call initialize() again, and
    // a solver left initialised against a mesh that failed to build is the kind
    // of state that produces a plausible second run.
    m_prepared         = false;
    m_temperatureIndex = -1;
    m_vertical.reset();
}

std::set<HydroCouple::Capability> FVQualComponent::capabilities() const
{
    std::set<HydroCouple::Capability> result{HydroCouple::Capability::Checkpointing};
    if(m_differentiationReady)result.insert(HydroCouple::Capability::Differentiable);return result;
}

HydroCouple::Temporal::IDateTime *FVQualComponent::currentDateTime() const
{
    const double day = m_epochJulianDay + currentTime() / 86400.0;
    if (m_clock.julianDay() != day) m_clock.setJulianDay(day);
    return &m_clock;
}

HydroCouple::Temporal::ITimeSpan *FVQualComponent::simulationPeriod() const
{
    if (m_period.julianDay() != m_epochJulianDay) m_period.setJulianDay(m_epochJulianDay);
    const double duration = m_endTime / 86400.0;
    if (m_period.duration() != duration) m_period.setDuration(duration);
    return &m_period;
}

double FVQualComponent::nextDateTimeJulianDay() const
{
    const double current = currentDateTime()->julianDay();
    if (!m_prepared || status() == ComponentStatus::Done || status() == ComponentStatus::Finished)
        return current;
    if (m_liveInflow && m_liveInflow->provider()) {
        const auto *donor=dynamic_cast<HydroCouple::Temporal::ITimeModelComponent *>(m_liveInflow->provider()->modelComponent());
        if (!donor) return current;
        const double published=donor->currentDateTime()->julianDay();
        return published>current ? published : donor->nextDateTimeJulianDay();
    }
    const double next = m_inflowArgument && !m_inflowArgument->frames().empty()
        ? (m_nextFrame < m_inflowArgument->frames().size()
            ? m_inflowArgument->frames()[m_nextFrame].endSeconds : currentTime())
        : std::min(m_endTime, currentTime() + m_couplingInterval);
    return m_epochJulianDay + next / 86400.0;
}

bool FVQualComponent::saveState(std::string &token, std::string &message)
{
    const auto prior = status();
    if (!m_prepared || (prior != ComponentStatus::Updated && prior != ComponentStatus::Done)) {
        message = "checkpoint requires a prepared Updated or Done component";
        setStatus(ComponentStatus::Failed, message);
        return false;
    }
    setStatus(ComponentStatus::Checkpointing);
    auto saved = m_solver.checkpoint(message);
    if (!saved) {
        setStatus(ComponentStatus::Failed, message);
        return false;
    }
    SavedState state;
    state.solver = std::move(saved);
    if (m_receiver) state.receiver = std::make_unique<Coupling::IntervalTransferReceiver>(*m_receiver);
    state.nextFrame = m_nextFrame;
    state.acceptedTransfer = m_lastAcceptedTransfer;
    for (int c = 0; c < m_columnCount; ++c) state.airTemperature.push_back((*m_airTemperature)[c]);
    static std::atomic<uint64_t> next{1};
    token = "fvqual-memory-v1:" + std::to_string(next.fetch_add(1));
    m_checkpoints.emplace(token, std::move(state));
    message.clear();
    setStatus(prior);
    return true;
}

bool FVQualComponent::restoreState(const std::string &token, std::string &message)
{
    const auto prior = status();
    auto found = m_checkpoints.find(token);
    if (!m_prepared || (prior != ComponentStatus::Updated && prior != ComponentStatus::Done) ||
        found == m_checkpoints.end()) {
        message = "checkpoint token is invalid, released or belongs to another component";
        setStatus(ComponentStatus::Failed, message);
        return false;
    }
    setStatus(ComponentStatus::Checkpointing);
    std::unique_ptr<Coupling::IntervalTransferReceiver> receiver;
    if (found->second.receiver) receiver = std::make_unique<Coupling::IntervalTransferReceiver>(*found->second.receiver);
    if (!m_solver.restoreCheckpoint(*found->second.solver, message)) {
        setStatus(ComponentStatus::Failed, message);
        return false;
    }
    for (int c = 0; c < m_columnCount; ++c) (*m_airTemperature)[c] = found->second.airTemperature[c];
    m_receiver = std::move(receiver);
    m_nextFrame = found->second.nextFrame;
    m_lastAcceptedTransfer = found->second.acceptedTransfer;
    m_linearized=false;
    publishTemperature();
    setStatus(ComponentStatus::Updated);
    return true;
}

bool FVQualComponent::releaseState(const std::string &token, std::string &message)
{
    if (m_checkpoints.erase(token) != 1) {
        message = "checkpoint token is invalid or already released";
        return false;
    }
    message.clear();
    return true;
}

bool FVQualComponent::applyLiveInflow(std::string &message)
{
    try {
        auto *provider=m_liveInflow->provider();
        const auto *donor=dynamic_cast<HydroCouple::Temporal::ITimeModelComponent *>(provider->modelComponent());
        if (!donor || donor->status()==ComponentStatus::Failed)
            throw std::invalid_argument("Live inflow donor is unavailable or failed");
        std::string text;
        const int64_t shape[1]={1},start[1]={0},count[1]={1};
        HydroCouple::BufferDescriptor buffer; buffer.data=&text; buffer.kind=HydroCouple::DataKind::String;
        buffer.rank=1; buffer.shape=shape; buffer.itemSizeBytes=sizeof(std::string);
        if (!provider->getValuesInto(buffer,start,count,&message)) return false;
        const auto packet=nlohmann::json::parse(text);
        if (packet.at("schema")!="hydrocouple.integrated-inflow/1.0" || packet.at("ownership")!="already_debited" ||
            packet.at("units")!=nlohmann::json{{"time","elapsed_s"},{"water","m3"},{"heat","J_relative_to_0_degC"},{"species","kg"}})
            throw std::invalid_argument("Unsupported live transfer schema, ownership or units");
        const double epoch=packet.at("epoch_julian_day").get<double>();
        if (!std::isfinite(epoch) || std::abs(epoch-m_epochJulianDay)>1.e-10)
            throw std::invalid_argument("Live transfer epoch differs from the receiver");
        Coupling::TransferFrame frame;
        frame.id=packet.at("id").get<std::string>();
        frame.startSeconds=packet.at("start_s"); frame.endSeconds=packet.at("end_s");
        if (!std::isfinite(frame.endSeconds) || frame.endSeconds>m_endTime+1.e-8 ||
            std::abs(donor->currentDateTime()->julianDay()-(epoch+frame.endSeconds/86400.0))>1.e-10)
            throw std::invalid_argument("Live transfer is outside the horizon or differs from the donor clock");
        const std::map<std::string,Forcing::Placement> placements{
            {"uniform",Forcing::Placement::Uniform},{"surface",Forcing::Placement::Surface},
            {"bed",Forcing::Placement::Bed},{"elevation",Forcing::Placement::Elevation},
            {"density-matched",Forcing::Placement::DensityMatched}};
        if (!packet.at("inflows").is_array()) throw std::invalid_argument("Live inflows must be an array");
        for (const auto &row:packet.at("inflows")) {
            Coupling::InflowTransfer flow;
            if (!row.at("column").is_number_integer()) throw std::invalid_argument("Live receiver column must be an integer index");
            flow.column=row.at("column"); flow.volumeM3=row.at("volume_m3"); flow.heatJ=row.at("heat_j");
            flow.placement=placements.at(row.value("placement",std::string("uniform")));
            flow.elevationM=row.value("elevation_m",0.0);
            flow.speciesMassKg=row.value("species_mass_kg",std::map<std::string,double>{});
            frame.inflows.push_back(std::move(flow));
        }
        Coupling::TransferReceipt receipt;
        if (!m_receiver->apply(frame,receipt,message)) return false;
        m_lastAcceptedTransfer=frame.id;
        message.clear(); return true;
    } catch (const std::exception &e) { message=e.what(); return false; }
}

// ── Publishing ────────────────────────────────────────────────────────────

void FVQualComponent::publishTemperature()
{
    if (m_acceptedTransfer) (*m_acceptedTransfer)[0]=m_lastAcceptedTransfer;
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
    for (int c=0;c<m_columnCount;++c) {
        (*m_columnFields[0])[c]=m_layered.z(c,0);
        (*m_columnFields[1])[c]=m_layered.z(c,0)-m_layered.z(c,m_layerCount);
        for (int k=0;k<m_layerCount;++k) {
            const auto cell=m_layered.cell(c,k);
            for (auto &o:m_fields) {
                double v=0;
                if (o->id()=="velocity_x") v=m_state.u[cell];
                else if (o->id()=="velocity_y") v=m_state.v[cell];
                else if (o->id()=="cell_volume") v=m_layered.cellVolume()[cell];
                else v=m_state.speciesAt(m_registry.indexOf(o->id()),cell);
                (*o)(c,k)=v;
            }
        }
    }
    double water=0,heat=0;
    if (m_receiver) for (const auto &[id,r]:m_receiver->receipts()) {water+=r.acceptedVolumeM3; heat+=r.acceptedHeatJ;}
    const auto ledger = Forcing::processBudgetValues(m_solver, m_registry);
    for (auto &o:m_budgetFields) {
        double v=0;
        if (o->id().starts_with("budget_")) v = ledger.at(o->id());
        else if (o->id()=="water_volume") for (double volume:m_layered.cellVolume()) v+=volume;
        else if (o->id()=="heat_inventory") v=m_state.totalMass(m_temperatureIndex,m_layered)*m_options.referenceDensity*m_options.heatCapacity;
        else if (o->id()=="accepted_water") v=water;
        else if (o->id()=="accepted_heat") v=heat;
        else if (o->id().starts_with("accepted_")) {
            if (m_receiver) for (const auto &[id,r]:m_receiver->receipts()) {
                const auto it=r.acceptedSpeciesKg.find(o->id().substr(9));
                if (it!=r.acceptedSpeciesKg.end()) v+=it->second;
            }
        }
        else {
            const int species=m_registry.indexOf(o->id().substr(10));
            const auto &units=m_registry.definition(species).units;
            const double scale=units=="kg/m3" ? 1.0 : (units=="mg/L" || units=="g/m3") ? .001 : .000001;
            v=m_state.totalMass(species,m_layered)*scale;
        }
        (*o)[0]=v;
    }
}

} // namespace FVQual::Component
