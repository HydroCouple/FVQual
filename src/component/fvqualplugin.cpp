// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   fvqualplugin.cpp
 * \brief  The loadable half: FVQual published through the stamped ABI
 *         entry points, so HydroCoupleComposer's registry, the headless
 *         runner, and the Python loader all pick it up identically.
 *
 * This translation unit exists only in the plugin library — the core
 * stays a plain library an embedder links, and the entry points live
 * where a loader looks for them.
 */

#include "hydrocouplesdk/componentabi.h"

#include "fvqual/component/fvqualcomponentinfo.h"

HYDROCOUPLE_DECLARE_COMPONENT(FVQual::Component::FVQualComponentInfo)
