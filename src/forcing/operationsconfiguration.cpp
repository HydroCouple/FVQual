// SPDX-License-Identifier: MIT
#include "fvqual/forcing/operationsconfiguration.h"
#include "fvqual/hydro/internalmode.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>

namespace FVQual::Forcing {
namespace {
using J = nlohmann::json;
using K = OperationControl::Kind;
void require(bool condition, const std::string &text) { if (!condition) throw std::invalid_argument(text); }
void keys(const J &row, std::initializer_list<const char *> names) {
    require(row.is_object(), "operation must be an object");
    for (const auto &[name, value] : row.items()) {
        (void)value;
        require(std::any_of(names.begin(), names.end(), [&](const char *n) { return name == n; }), "unknown operation field: " + name);
    }
}
double number(const J &row, const char *key, double fallback) {
    if (!row.contains(key)) return fallback;
    require(row[key].is_number(), std::string(key) + " requires a number");
    double value = row[key]; require(std::isfinite(value), std::string(key) + " must be finite"); return value;
}
int64_t column(const J &row, const char *key, const Mesh::HorizontalMesh &mesh) {
    require(row.contains(key) && row[key].is_number_integer(), std::string(key) + " requires a zero-based column index");
    int64_t value = row[key]; require(value >= 0 && value < mesh.faceCount(), std::string(key) + " is outside the reservoir mesh"); return value;
}
Placement placement(const J &row, const char *key, Placement fallback) {
    if (!row.contains(key)) return fallback;
    static const std::map<std::string, Placement> choices{{"uniform", Placement::Uniform}, {"surface", Placement::Surface},
        {"bed", Placement::Bed}, {"elevation", Placement::Elevation}, {"density-matched", Placement::DensityMatched},
        {"selective-withdrawal", Placement::SelectiveWithdrawal}};
    const std::string name = row[key]; require(choices.contains(name), "unknown placement: " + name); return choices.at(name);
}
std::string identifier(const J &row, std::set<std::string> &used) {
    const std::string id = row.at("id"); require(!id.empty() && used.insert(id).second, "operation IDs must be nonempty and unique: " + id); return id;
}
std::vector<double> concentrations(const J &row, const State::SpeciesRegistry &registry) {
    std::vector<double> values(static_cast<size_t>(registry.count()), 0);
    auto data = row.value("concentrations", J::object()); require(data.is_object(), "concentrations require a species/value object in registered units");
    for (const auto &[id, value] : data.items()) {
        int s = registry.indexOf(id); require(s >= 0 && s != registry.temperatureIndex(), "unregistered concentration: " + id);
        require(registry.definition(s).isTransported(), "inflow concentrations must be transported species");
        require(value.is_number() && std::isfinite(value.get<double>()) && value.get<double>() >= 0, "concentrations must be finite and nonnegative"); values[s] = value;
    }
    return values;
}
void schedule(const J &row, std::initializer_list<std::pair<const char *, K>> fields,
              int index, std::vector<OperationControl> &controls,
              const State::SpeciesRegistry &registry, bool sourceValues = false) {
    if (!row.contains("schedule")) return;
    const auto &records = row.at("schedule"); require(records.is_array() && !records.empty(), "schedule must have records starting at time_s = 0");
    double prior = -1;
    std::set<std::string> fieldsUsed;
    for (const auto &r : records) {
        require(r.is_object(), "schedule record must be an object");
        for (const auto &[name, value] : r.items()) {
            (void)value;
            require(name == "time_s" || (sourceValues && (name == "temperature_c" || name == "concentrations")) ||
                std::any_of(fields.begin(), fields.end(), [&](auto f) { return name == f.first; }), "unsupported scheduled control: " + name);
            if (name != "time_s") fieldsUsed.insert(name);
        }
        double time = number(r, "time_s", -1);
        require(time >= 0 && time > prior && (prior >= 0 || time == 0), "schedule times must increase from zero in elapsed seconds"); prior = time;
    }
    auto add = [&](const std::string &name, K kind, int species = -1, const std::string &speciesId = "") {
        OperationControl c; c.kind = kind; c.index = index; c.species = species; c.speciesId = speciesId;
        for (const auto &r : records) {
            const auto &v = name == "concentrations" ? r.at(name).at(speciesId) : r.at(name);
            require(v.is_number() && std::isfinite(v.get<double>()), "scheduled values must be finite numbers");
            double value = v;
            if (kind == K::StructureOpening) require(value >= 0 && value <= 1, "opening must be in [0,1]");
            if (kind == K::StructureDemand || kind == K::GroupDemand || kind == K::TransferDischarge || (kind == K::SourceValue && species != registry.temperatureIndex()))
                require(value >= 0, "scheduled magnitude/concentration must be nonnegative");
            c.times.push_back(r.at("time_s")); c.values.push_back(value);
        }
        controls.push_back(std::move(c));
    };
    for (auto [name, kind] : fields) if (fieldsUsed.contains(name)) add(name, kind);
    if (sourceValues && fieldsUsed.contains("temperature_c")) add("temperature_c", K::SourceValue, registry.temperatureIndex(), "__TEMPERATURE__");
    if (sourceValues && fieldsUsed.contains("concentrations")) {
        const auto &first = records[0].at("concentrations"); require(first.is_object(), "schedule concentrations must be an object");
        for (const auto &[id, value] : first.items()) {
            (void)value; int species = registry.indexOf(id);
            require(species >= 0 && species != registry.temperatureIndex() && registry.definition(species).isTransported(), "unregistered scheduled species: " + id);
            add("concentrations", K::SourceValue, species, id);
        }
        for (const auto &r : records) require(r.at("concentrations").size() == first.size(), "each schedule row must contain the same species");
    }
}
void weatherFields(const J &row, MeteorologyField &weather, int64_t columns) {
    keys(row, {"time_s", "air_temperature_c", "dew_point_c", "wind_speed_m_s", "wind_towards_rad", "cloud_fraction", "shortwave_w_m2", "pressure_pa", "albedo"});
    for (auto [name, field] : std::initializer_list<std::pair<const char *, ForcingField *>>{
        {"air_temperature_c", &weather.airTemperature}, {"dew_point_c", &weather.dewPoint}, {"wind_speed_m_s", &weather.windSpeed},
        {"wind_towards_rad", &weather.windDirection}, {"cloud_fraction", &weather.cloudCover}, {"shortwave_w_m2", &weather.shortwave},
        {"pressure_pa", &weather.pressure}, {"albedo", &weather.albedo}}) {
        if (!row.contains(name)) continue;
        std::vector<double> values;
        if (row[name].is_array()) { require(row[name].size() == static_cast<size_t>(columns), "meteorology array must match reservoir column count"); values = row[name].get<std::vector<double>>(); }
        else values = {number(row, name, 0)};
        for (double value : values) {
            require(std::isfinite(value), "weather values must be finite");
            if (field == &weather.windSpeed || field == &weather.shortwave) require(value >= 0, "wind and shortwave must be nonnegative");
            if (field == &weather.cloudCover || field == &weather.albedo) require(value >= 0 && value <= 1, "cloud and albedo must be fractions");
            if (field == &weather.pressure) require(value > 0, "pressure must be positive");
            if (field == &weather.airTemperature || field == &weather.dewPoint) require(value > -273.15, "temperature must exceed absolute zero");
        }
        if (row[name].is_array()) { std::string message; require(field->setPerColumn(values, columns, message), message); }
        else field->setUniform(values[0]);
    }
}
}

bool configureOperations(Hydro::InternalModeSolver &solver, const Mesh::HorizontalMesh &mesh,
                         const State::SpeciesRegistry &registry, const J &data, std::string &message) {
    try {
        keys(data, {"meteorology", "sources", "structures", "outlet_groups", "open_boundaries", "internal_transfers"});
        SourceSet sources; StructureSet structures; OutletGroupSet groups; TransferSet transfers;
        std::vector<OperationControl> controls;
        std::set<std::string> ids;
        std::map<std::string, int> structureIds;
        double attenuation = .45;
        bool hasWeather = false;
        if (data.contains("meteorology") && !data["meteorology"].empty()) {
            const auto &met = data["meteorology"];
            keys(met, {"values", "records", "coefficients", "solar_attenuation_m_inv"});
            require(registry.temperatureIndex() >= 0, "meteorology requires temperature");
            attenuation = number(met, "solar_attenuation_m_inv", .45); require(attenuation >= 0, "solar attenuation must be nonnegative");
            MeteorologyField weather;
            const auto coefficients = met.value("coefficients", J::object());
            keys(coefficients, {"surface_absorption", "water_emissivity", "longwave_reflection", "evaporation_a", "evaporation_b", "wind_exponent", "wind_height_m", "evaporation_height_m", "wind_roughness_m", "psychrometric_pa_k", "air_density_kg_m3", "drag_coefficient"});
            for (auto [key, dest] : std::initializer_list<std::pair<const char *, double *>>{
                {"surface_absorption", &weather.base.surfaceAbsorption}, {"water_emissivity", &weather.base.waterEmissivity},
                {"longwave_reflection", &weather.base.longwaveReflection}, {"evaporation_a", &weather.base.evaporationA}, {"evaporation_b", &weather.base.evaporationB},
                {"wind_exponent", &weather.base.windExponent}, {"wind_height_m", &weather.base.windHeight}, {"evaporation_height_m", &weather.base.evaporationHeight},
                {"wind_roughness_m", &weather.base.windRoughness}, {"psychrometric_pa_k", &weather.base.psychrometric}, {"air_density_kg_m3", &weather.base.airDensity}, {"drag_coefficient", &weather.base.dragCoefficient}}) {
                *dest = number(coefficients, key, *dest); require(*dest >= 0, std::string(key) + " must be nonnegative");
            }
            require(weather.base.surfaceAbsorption <= 1 && weather.base.waterEmissivity <= 1 && weather.base.longwaveReflection <= 1, "radiation fractions must be in [0,1]");
            require(weather.base.windRoughness > 0 && weather.base.windHeight > weather.base.windRoughness && weather.base.evaporationHeight > weather.base.windRoughness && weather.base.psychrometric > 0 && weather.base.airDensity > 0, "invalid meteorology heights, density or psychrometric coefficient");
            weatherFields(met.value("values", J::object()), weather, mesh.faceCount());
            OperationControl c; c.kind = K::Meteorology;
            if (met.contains("records")) {
                const auto &records = met["records"]; require(records.is_array() && !records.empty(), "meteorology records must start at zero");
                double prior = -1;
                for (const auto &r : records) {
                    double t = number(r, "time_s", -1); require(t >= 0 && t > prior && (prior >= 0 || t == 0), "weather records must increase from time_s = 0"); prior = t;
                    weatherFields(r, weather, mesh.faceCount()); c.times.push_back(t); c.weather.push_back(weather);
                }
            } else { c.times = {0}; c.weather = {weather}; }
            controls.push_back(std::move(c)); hasWeather = true;
        }
        auto rows = [&](const char *key) { auto value = data.value(key, J::array()); require(value.is_array(), std::string(key) + " must be an array"); return value; };
        for (const auto &row : rows("sources")) {
            keys(row, {"id", "column", "discharge_m3_s", "placement", "elevation_m", "width_m", "withdrawal_angle_rad", "temperature_c", "concentrations", "schedule"});
            Source source; source.id = identifier(row, ids); source.column = column(row, "column", mesh);
            source.discharge = number(row, "discharge_m3_s", 0); source.placement = placement(row, "placement", Placement::Uniform);
            source.elevation = number(row, "elevation_m", 0); source.outletWidth = number(row, "width_m", 0);
            source.withdrawalAngle = number(row, "withdrawal_angle_rad", 3.141592653589793); require(source.outletWidth >= 0 && source.withdrawalAngle > 0, "invalid withdrawal width/angle");
            source.temperature = number(row, "temperature_c", 10); source.concentration = concentrations(row, registry);
            int i = sources.add(source, message); require(i >= 0, message);
            schedule(row, {{"discharge_m3_s", K::SourceDischarge}}, i, controls, registry, true);
            for (const auto &c : controls) if (c.kind == K::SourceDischarge && c.index == i)
                for (double q : c.values) require((source.placement != Placement::SelectiveWithdrawal || q <= 0) && (source.placement != Placement::DensityMatched || q >= 0), "scheduled discharge contradicts placement");
        }
        for (const auto &row : rows("structures")) {
            keys(row, {"id", "type", "column", "crest_elevation_m", "coefficient", "width_m", "area_m2", "opening", "tailwater_elevation_m", "capacity_m3_s", "scheduled_discharge_m3_s", "start_elevation_m", "stop_elevation_m", "placement", "outlet_width_m", "temperature_c", "concentrations", "schedule"});
            Structure s; s.id = identifier(row, ids); s.column = column(row, "column", mesh);
            const std::map<std::string, StructureType> types{{"weir", StructureType::Weir}, {"orifice", StructureType::Orifice}, {"pipe", StructureType::Pipe}, {"pump", StructureType::Pump}};
            const std::string type = row.at("type"); require(types.contains(type), "unknown structure type"); s.type = types.at(type);
            s.crestElevation = number(row, "crest_elevation_m", 0); s.coefficient = number(row, "coefficient", 0);
            s.width = number(row, "width_m", 0); s.area = number(row, "area_m2", 0); s.opening = number(row, "opening", 1);
            s.tailwaterElevation = number(row, "tailwater_elevation_m", -1e30); s.capacity = number(row, "capacity_m3_s", 0);
            s.scheduled = row.contains("scheduled_discharge_m3_s"); s.scheduledDischarge = number(row, "scheduled_discharge_m3_s", 0);
            s.controlled = row.contains("start_elevation_m") || row.contains("stop_elevation_m");
            require(!s.controlled || (row.contains("start_elevation_m") && row.contains("stop_elevation_m")), "both pump control elevations are required");
            s.startElevation = number(row, "start_elevation_m", 0); s.stopElevation = number(row, "stop_elevation_m", 0);
            s.placement = placement(row, "placement", Placement::SelectiveWithdrawal); s.outletWidth = number(row, "outlet_width_m", 0);
            require(s.outletWidth >= 0, "outlet width must be nonnegative");
            s.temperature = number(row, "temperature_c", 10); s.concentration = concentrations(row, registry);
            int i = structures.add(s, sources, message); require(i >= 0, message); structureIds[s.id] = i;
            schedule(row, {{"opening", K::StructureOpening}, {"tailwater_elevation_m", K::StructureTailwater}, {"scheduled_discharge_m3_s", K::StructureDemand}}, i, controls, registry);
        }
        std::set<int> assignedPorts;
        for (const auto &row : rows("outlet_groups")) {
            keys(row, {"id", "ports", "demand_m3_s", "target_temperature_c", "blend", "schedule"});
            OutletGroup group; group.id = identifier(row, ids); group.demand = number(row, "demand_m3_s", 0);
            group.targetTemperature = number(row, "target_temperature_c", 0); group.blend = row.value("blend", true);
            require(!group.blend || registry.temperatureIndex() >= 0, "temperature blending requires temperature");
            require(row.at("ports").is_array() && !row.at("ports").empty(), "outlet group needs named structure ports");
            for (const auto &port : row.at("ports")) {
                const std::string id = port; require(structureIds.contains(id), "unknown outlet port: " + id);
                int i = structureIds.at(id); require(assignedPorts.insert(i).second, "a port may belong to only one group"); group.ports.push_back(i);
            }
            int i = groups.add(group, structures, message); require(i >= 0, message);
            schedule(row, {{"demand_m3_s", K::GroupDemand}, {"target_temperature_c", K::GroupTemperature}}, i, controls, registry);
        }
        for (const auto &row : rows("internal_transfers")) {
            keys(row, {"id", "from_column", "to_column", "discharge_m3_s", "from_placement", "to_placement", "from_elevation_m", "to_elevation_m", "width_m", "schedule"});
            Transfer transfer; transfer.id = identifier(row, ids); transfer.fromColumn = column(row, "from_column", mesh); transfer.toColumn = column(row, "to_column", mesh);
            require(transfer.fromColumn != transfer.toColumn, "internal transfer must connect different columns");
            transfer.discharge = number(row, "discharge_m3_s", 0); transfer.fromPlacement = placement(row, "from_placement", Placement::SelectiveWithdrawal);
            transfer.toPlacement = placement(row, "to_placement", Placement::DensityMatched);
            transfer.fromElevation = number(row, "from_elevation_m", 0); transfer.toElevation = number(row, "to_elevation_m", 0); transfer.outletWidth = number(row, "width_m", 0);
            require(transfer.outletWidth >= 0 && transfer.fromPlacement != Placement::DensityMatched && transfer.toPlacement != Placement::SelectiveWithdrawal, "invalid transfer placement or width");
            int i = transfers.add(transfer, sources, message); require(i >= 0, message);
            schedule(row, {{"discharge_m3_s", K::TransferDischarge}}, i, controls, registry);
        }
        struct Boundary { Hydro::BarotropicSolver::BoundaryCondition condition; std::vector<int64_t> edges; J concentrations; };
        std::vector<Boundary> boundaries;
        std::set<int64_t> usedEdges;
        Hydro::BarotropicSolver checker; require(checker.initialize(mesh, {}, message), message);
        for (const auto &row : rows("open_boundaries")) {
            keys(row, {"id", "type", "node_pairs", "specific_discharge_m2_s", "stage_m", "supercritical_depth_m", "temperature_c", "concentrations", "schedule"});
            Boundary b; b.condition.id = identifier(row, ids);
            using BT = Hydro::BarotropicSolver::BoundaryType;
            const std::map<std::string, BT> types{{"inflow", BT::Inflow}, {"stage", BT::Stage}, {"transmissive", BT::Transmissive}, {"wall", BT::Wall}};
            std::string type = row.at("type"); require(types.contains(type), "unknown boundary type"); b.condition.type = types.at(type);
            b.condition.specificDischarge = number(row, "specific_discharge_m2_s", 0); b.condition.stage = number(row, "stage_m", 0); b.condition.depth = number(row, "supercritical_depth_m", 0);
            require(b.condition.depth >= 0, "boundary depth must be nonnegative");
            require(row.at("node_pairs").is_array() && !row.at("node_pairs").empty(), "boundary needs zero-based endpoint node pairs");
            for (const auto &pair : row.at("node_pairs")) {
                require(pair.is_array() && pair.size() == 2 && pair[0].is_number_integer() && pair[1].is_number_integer(), "boundary node pair needs two integer indices");
                int64_t a = pair[0], d = pair[1], edge = -1;
                for (int64_t e = 0; e < mesh.edgeCount(); ++e) {
                    const auto &n = mesh.edgeNodes();
                    if ((n[2*e] == a && n[2*e+1] == d) || (n[2*e] == d && n[2*e+1] == a)) { edge = e; break; }
                }
                require(edge >= 0 && mesh.edgeRightFace()[edge] < 0, "boundary pair is not an exterior reservoir edge");
                require(usedEdges.insert(edge).second, "boundary edge assigned more than once"); b.edges.push_back(edge);
            }
            int i = checker.addBoundary(b.condition, b.edges, message); require(i >= 0, message);
            b.concentrations = row.value("concentrations", J::object()); (void)concentrations(row, registry);
            if (row.contains("temperature_c")) { require(registry.temperatureIndex() >= 0, "boundary temperature requires temperature"); b.concentrations[registry.definition(registry.temperatureIndex()).id] = number(row, "temperature_c", 10); }
            const size_t first = controls.size();
            schedule(row, {{"specific_discharge_m2_s", K::BoundaryDischarge}, {"stage_m", K::BoundaryStage}}, i, controls, registry, true);
            for (size_t c = first; c < controls.size(); ++c) if (controls[c].kind == K::SourceValue) {
                controls[c].kind = K::BoundaryValue;
                controls[c].speciesId = registry.definition(controls[c].species).id;
            }
            boundaries.push_back(std::move(b));
        }
        // Configuration only: replace operation sets before the first step.
        require(solver.time() == 0 && solver.sources().count() == 0 && solver.external().boundaryCount() == 0, "operations must be installed once before advancing");
        solver.sources() = std::move(sources); solver.structures() = std::move(structures);
        solver.outletGroups() = std::move(groups); solver.transfers() = std::move(transfers);
        for (const auto &b : boundaries) require(solver.external().addBoundary(b.condition, b.edges, message) >= 0, message);
        if (!boundaries.empty()) {
            require(solver.bindOpenBoundaries(message), message);
            for (size_t i = 0; i < boundaries.size(); ++i)
                for (const auto &[id, value] : boundaries[i].concentrations.items()) require(solver.setBoundaryValue(static_cast<int>(i), id, value, message), message);
        }
        if (hasWeather) solver.enableMeteorology(attenuation);
        solver.setOperationControls(std::move(controls));
        message.clear(); return true;
    } catch (const std::exception &e) { message = e.what(); return false; }
}

nlohmann::json processBudgetValues(const Hydro::InternalModeSolver &solver, const State::SpeciesRegistry &registry) {
    J result = J::object(); const auto &b = solver.budget();
    result["budget_water_initial"] = b.initialWater;
    result["budget_water_residual"] = solver.waterBudgetResidual();
    static const std::array<const char *, 5> origins{"source", "structure", "boundary", "coupled", "transfer"};
    for (size_t i = 0; i < origins.size(); ++i) {
        result["budget_water_" + std::string(origins[i]) + "_in"] = b.flows[i].waterIn;
        result["budget_water_" + std::string(origins[i]) + "_out"] = b.flows[i].waterOut;
    }
    for (int s = 0; s < registry.count(); ++s) {
        const auto &d = registry.definition(s);
        double scale = (d.units == "kg/m3" || d.units == "kg/m2") ? 1 : (d.units == "mg/L" || d.units == "g/m3" || d.units == "g/m2") ? .001 : (d.units == "ug/L" || d.units == "mg/m3") ? .000001 : 1;
        const bool heat = s == registry.temperatureIndex();
        if (heat) scale = solver.options().referenceDensity * solver.options().heatCapacity;
        const std::string prefix = heat ? "budget_heat_" : "budget_" + d.id + "_";
        result[prefix + "initial"] = b.initialScalar.at(s) * scale;
        result[prefix + "residual"] = solver.scalarBudgetResidual(s) * scale;
        result[prefix + "settling"] = b.settling.at(s) * scale;
        result[prefix + "reaction"] = b.reactions.at(s) * scale;
        for (size_t i = 0; i < origins.size(); ++i) {
            result[prefix + origins[i] + "_in"] = b.flows[i].scalarIn.at(s) * scale;
            result[prefix + origins[i] + "_out"] = b.flows[i].scalarOut.at(s) * scale;
        }
    }
    if (registry.temperatureIndex() >= 0) {
        const std::array<const char *, 6> names{"shortwave", "longwave_in", "back_radiation", "latent", "sensible", "prescribed"};
        for (size_t i = 0; i < names.size(); ++i) result["budget_heat_" + std::string(names[i])] = b.surfaceHeat[i];
    }
    return result;
}
}
