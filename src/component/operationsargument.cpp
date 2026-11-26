// SPDX-License-Identifier: MIT
#include "fvqual/component/operationsargument.h"
#include "hydrocouplesdk/io/uriresolver.h"
#include <stdexcept>
namespace FVQual::Component {
using AT = HydroCouple::IArgument::ArgumentInputType;
OperationsArgument::OperationsArgument(const std::string &section, HydroCouple::SDK::Dimension *d,
    HydroCouple::SDK::ValueDefinition *q, HydroCouple::IModelComponent *c)
    : HydroCouple::SDK::Argument1DString(section, d, 1, q, c), m_section(section) {
    setCaption(section); setDescription("Reservoir operations. Schedules hold values between knots on elapsed seconds. Columns and boundary nodes use zero-based mesh indices.");
    setIsOptional(true); setRole(HydroCouple::IArgument::ArgumentRole::Forcing);
    addFileFilter("Reservoir operations (*.json)");
    (*this)[0] = (section == "meteorology" ? nlohmann::json::object() : section == "reaction_deck" ? nlohmann::json("") : nlohmann::json::array()).dump();
}
nlohmann::json OperationsArgument::value() const { return nlohmann::json::parse((*this)[0]); }
bool OperationsArgument::initialize(const std::string &text, AT type, std::string &message) {
    try {
        std::string content = text;
        if (type == AT::File || type == AT::URL) {
            if (!HydroCouple::SDK::IO::readReference(text, content, message)) return false;
        } else if (type != AT::JSON) throw std::invalid_argument("operations require JSON or a local JSON file");
        auto data = nlohmann::json::parse(content);
        if (data.value("schema", std::string()) != "fvqual.operations/1.0" || data.at("section") != m_section || data.size() != 3)
            throw std::invalid_argument("expected fvqual.operations/1.0 with matching section and value");
        const auto &v = data.at("value");
        if (!(m_section == "meteorology" ? v.is_object() : m_section == "reaction_deck" ? v.is_string() : v.is_array()))
            throw std::invalid_argument("incorrect operation value shape");
        if (validator && !validator(data, message)) return false;
        (*this)[0] = v.dump(); setCurrentArgumentInputType(AT::JSON); message.clear(); return true;
    } catch (const std::exception &e) { message = e.what(); return false; }
}
bool OperationsArgument::serialize(AT type, std::string &text, std::string &message) const {
    if (type != AT::JSON) { message = "operations serialize as JSON"; return false; }
    text = nlohmann::json{{"schema", "fvqual.operations/1.0"}, {"section", m_section}, {"value", value()}}.dump(); message.clear(); return true;
}
}
