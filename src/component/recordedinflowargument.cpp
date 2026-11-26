// SPDX-License-Identifier: MIT
#include "fvqual/component/recordedinflowargument.h"
#include "hydrocouplesdk/io/uriresolver.h"
#include <fstream>
#include <set>
#include <cmath>
#include <stdexcept>
namespace FVQual::Component {
using AT=HydroCouple::IArgument::ArgumentInputType;
using HydroCouple::SDK::Argument1DString;
RecordedInflowArgument::RecordedInflowArgument(HydroCouple::SDK::Dimension *dimension,
    HydroCouple::SDK::ValueDefinition *quantity, HydroCouple::IModelComponent *component)
    : Argument1DString("recorded_inflows",dimension,1,quantity,component) {
    setCaption("SWMM inflow series");
    setDescription("Accepted donor exports: elapsed seconds, m3 water, J relative to 0 Celsius and kg species. One-way only.");
    setIsOptional(true);
    setRole(HydroCouple::IArgument::ArgumentRole::Forcing);
    addFileFilter("Integrated inflow series (*.json)");
    (*this)[0]="";
}
bool RecordedInflowArgument::initialize(const std::string &value,AT type,std::string &message) {
    const std::string prior=(*this)[0];
    if (type==AT::JSON) {
        try {
            auto data=nlohmann::json::parse(value);
            if (!data.contains("@uri")) {
                const auto values=data.is_object() ? data.at("values") : data;
                if (!values.is_array() || values.size()!=1 || !values[0].is_string()) {
                    message="inflow argument needs exactly one file reference"; return false;
                }
            }
        } catch (const std::exception &e) { message=e.what(); return false; }
    }
    if (type==AT::File || type==AT::URL) {
        (*this)[0]=value;
    } else {
        std::string uri;
        if (type==AT::JSON && HydroCouple::SDK::IO::jsonUriReference(value,uri)) (*this)[0]=uri;
        else if (!Argument1DString::initialize(value,type,message)) return false;
    }
    if (!load(message)) { (*this)[0]=prior; return false; }
    setCurrentArgumentInputType(type);
    return true;
}
bool RecordedInflowArgument::initialize(const HydroCouple::IComponentDataItem &item,std::string &message) {
    const std::string prior=(*this)[0];
    if (!Argument1DString::initialize(item,message)) return false;
    if (!load(message)) { (*this)[0]=prior; return false; }
    return true;
}
bool RecordedInflowArgument::serialize(AT type,std::string &value,std::string &message) const {
    if (type==AT::JSON && !(*this)[0].empty()) {
        value=nlohmann::json{{"@uri",(*this)[0]}}.dump();
        message.clear(); return true;
    }
    return Argument1DString::serialize(type,value,message);
}
bool RecordedInflowArgument::load(std::string &message) {
    try {
        std::vector<Coupling::TransferFrame> frames;
        std::vector<State::SpeciesDefinition> species;
        if ((*this)[0].empty()) {
            m_frames.clear(); m_species.clear(); message.clear(); return true;
        }
        const auto scheme=HydroCouple::SDK::IO::uriScheme((*this)[0]);
        if (!scheme.empty() && scheme!="file")
            throw std::invalid_argument("recorded inflows require a local file reference");
        std::string text;
        if (!HydroCouple::SDK::IO::readReference((*this)[0],text,message)) return false;
        const auto data=nlohmann::json::parse(text);
        if (data.at("schema")!="fvqual.inflow-series/1.0" ||
            data.at("units")!=nlohmann::json{{"time","elapsed_s"},{"water","m3"},
                    {"heat","J_relative_to_0_degC"},{"species","kg"}})
            throw std::invalid_argument("unsupported inflow series schema or units");
        const auto speciesData = data.value("species",nlohmann::json::object());
        for (const auto &[id,properties] : speciesData.items()) {
            State::SpeciesDefinition definition;
            definition.id=id; definition.units=properties.at("units");
            definition.initialValue=properties.value("initial",0.0);
            if (definition.units!="kg/m3" || !std::isfinite(definition.initialValue) ||
                definition.initialValue<0 || id.empty())
                throw std::invalid_argument("inflow species require kg/m3 and finite nonnegative initial concentration");
            species.push_back(std::move(definition));
        }
        const std::map<std::string,Forcing::Placement> placements{
            {"uniform",Forcing::Placement::Uniform},{"surface",Forcing::Placement::Surface},
            {"bed",Forcing::Placement::Bed},{"elevation",Forcing::Placement::Elevation},
            {"density-matched",Forcing::Placement::DensityMatched}};
        double previous=0;
        std::set<std::string> ids;
        for (const auto &row : data.at("frames")) {
            Coupling::TransferFrame frame;
            frame.id=row.at("id");
            frame.startSeconds=row.at("start_s"); frame.endSeconds=row.at("end_s");
            if (frame.id.empty() || !ids.insert(frame.id).second ||
                frame.startSeconds!=previous || !std::isfinite(frame.endSeconds) || frame.endSeconds<=previous)
                throw std::invalid_argument("inflow intervals must be contiguous from zero with unique IDs");
            previous=frame.endSeconds;
            for (const auto &r : row.at("inflows")) {
                Coupling::InflowTransfer transfer;
                if (!r.at("column").is_number_integer()) throw std::invalid_argument("Recorded receiver column must be an integer index");
                transfer.column=r.at("column"); transfer.volumeM3=r.at("volume_m3"); transfer.heatJ=r.at("heat_j");
                transfer.elevationM=r.value("elevation_m",0.0);
                transfer.placement=placements.at(r.value("placement",std::string("uniform")));
                transfer.speciesMassKg=r.value("species_mass_kg",std::map<std::string,double>{});
                if (transfer.column<0 || !std::isfinite(transfer.volumeM3) || transfer.volumeM3<0 ||
                    !std::isfinite(transfer.heatJ) || !std::isfinite(transfer.elevationM))
                    throw std::invalid_argument("invalid inflow column/water/heat/elevation");
                for (const auto &[id,mass] : transfer.speciesMassKg)
                    if (!std::isfinite(mass) || mass<0 || (transfer.volumeM3==0 && mass!=0))
                        throw std::invalid_argument("invalid species mass");
                if (transfer.volumeM3==0 && transfer.heatJ!=0)
                    throw std::invalid_argument("water-associated heat requires water");
                frame.inflows.push_back(std::move(transfer));
            }
            frames.push_back(std::move(frame));
        }
        if (frames.empty()) throw std::invalid_argument("inflow series contains no intervals");
        m_frames=std::move(frames); m_species=std::move(species);
        message.clear(); return true;
    } catch (const std::exception &e) { message=e.what(); return false; }
}
}
