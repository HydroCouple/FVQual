// SPDX-License-Identifier: MIT
#ifndef FVQUAL_RECORDEDINFLOWARGUMENT_H
#define FVQUAL_RECORDEDINFLOWARGUMENT_H
#include "fvqual/coupling/intervaltransfer.h"
#include "hydrocouplesdk/data/argument1d.h"
namespace FVQual::Component {
class FVQual_EXPORT RecordedInflowArgument : public HydroCouple::SDK::Argument1DString {
public:
    RecordedInflowArgument(HydroCouple::SDK::Dimension *dimension,
                           HydroCouple::SDK::ValueDefinition *quantity,
                           HydroCouple::IModelComponent *component);
    bool initialize(const std::string &value,
                    HydroCouple::IArgument::ArgumentInputType type,
                    std::string &message) override;
    bool initialize(const HydroCouple::IComponentDataItem &item, std::string &message) override;
    bool serialize(HydroCouple::IArgument::ArgumentInputType type, std::string &value, std::string &message) const override;
    bool load(std::string &message);
    const std::vector<Coupling::TransferFrame> &frames() const { return m_frames; }
    const std::vector<State::SpeciesDefinition> &species() const { return m_species; }
private:
    std::vector<Coupling::TransferFrame> m_frames;
    std::vector<State::SpeciesDefinition> m_species;
};
}
#endif
