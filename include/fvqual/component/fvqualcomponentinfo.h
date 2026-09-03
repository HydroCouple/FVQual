// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   fvqualcomponentinfo.h
 * \brief  FVQualComponentInfo — the component's metadata and factory.
 *
 * The half that makes FVQual loadable: a host's registry asks a library
 * for its IComponentInfo and creates instances through it. Kept in the
 * core library (not only the plugin) so an embedder can create components
 * the same way a host does.
 */

#ifndef FVQUAL_COMPONENT_FVQUALCOMPONENTINFO_H
#define FVQUAL_COMPONENT_FVQUALCOMPONENTINFO_H

#include "fvqual/fvqual_export.h"

#include "hydrocouplesdk/component/abstractmodelcomponentinfo.h"

#include <memory>

namespace FVQual::Component {

/*!
 * \brief Metadata and instance factory for FVQualComponent.
 */
class FVQual_EXPORT FVQualComponentInfo
    : public HydroCouple::SDK::AbstractModelComponentInfo
{
public:
    FVQualComponentInfo();

    ~FVQualComponentInfo() override = default;

    [[nodiscard]] std::unique_ptr<HydroCouple::IModelComponent>
    createComponentInstance() override;
};

} // namespace FVQual::Component

#endif // FVQUAL_COMPONENT_FVQUALCOMPONENTINFO_H
