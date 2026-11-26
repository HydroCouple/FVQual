// SPDX-License-Identifier: MIT
#ifndef FVQUAL_OPERATIONSARGUMENT_H
#define FVQUAL_OPERATIONSARGUMENT_H
#include "hydrocouplesdk/data/argument1d.h"
#include "fvqual/fvqual_export.h"
#include <functional>
#include <nlohmann/json.hpp>
namespace FVQual::Component {
class FVQual_EXPORT OperationsArgument : public HydroCouple::SDK::Argument1DString {
public:
    OperationsArgument(const std::string &section, HydroCouple::SDK::Dimension *dimension,
        HydroCouple::SDK::ValueDefinition *quantity, HydroCouple::IModelComponent *component);
    bool initialize(const std::string &value, HydroCouple::IArgument::ArgumentInputType type, std::string &message) override;
    bool serialize(HydroCouple::IArgument::ArgumentInputType type, std::string &value, std::string &message) const override;
    nlohmann::json value() const;
    std::function<bool(const nlohmann::json &, std::string &)> validator;
private:
    std::string m_section;
};
}
#endif
