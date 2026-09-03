// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/component/fvqualcomponentinfo.h"

#include "fvqual/component/fvqualcomponent.h"
#include "fvqual/version.h"

namespace FVQual::Component {

FVQualComponentInfo::FVQualComponentInfo()
    : AbstractModelComponentInfo("fvqual.component")
{
    setCaption("FVQual");
    setDescription(
        "Finite-volume water-quality and hydrodynamics model on layered "
        "unstructured meshes (CE-QUAL-W2 lineage).");
    setDeveloper("Caleb Buahin");
    setVersion(FVQual_VERSION);
    setTags({"Model", "Water Quality", "Hydrodynamics"});
}

std::unique_ptr<HydroCouple::IModelComponent>
FVQualComponentInfo::createComponentInstance()
{
    return std::make_unique<FVQualComponent>("fvqual.instance", this);
}

} // namespace FVQual::Component
