// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
//
// The plugin gate: FVQual loaded exactly as a host loads it — dlopen, the
// ABI stamp compared before any C++ symbol, the info entry point cast
// across the image boundary, and a created instance initialized far
// enough to show its ports and run its default basin to Finished.

#include <gtest/gtest.h>

#include "hydrocouplesdk/componentabi.h"

#include <dlfcn.h>

#include <filesystem>
#include <memory>
#include <string>

namespace {

std::filesystem::path pluginPath()
{
    std::filesystem::path path =
        std::filesystem::path(FVQUAL_COMPONENTS_DIR) /
        "libfvqual-component.so";
    if (!std::filesystem::exists(path))
        path.replace_extension(".dylib");
    return path;
}

} // namespace

TEST(FVQualPluginTest, LoadsThroughTheStampedConventionAndRunsItsDefaults)
{
    void *handle = dlopen(pluginPath().c_str(), RTLD_NOW | RTLD_LOCAL);
    ASSERT_NE(handle, nullptr) << dlerror();

    auto abi = reinterpret_cast<HydroCoupleComponentAbiFn>(
        dlsym(handle, HYDROCOUPLE_COMPONENT_ABI_SYMBOL));
    ASSERT_NE(abi, nullptr) << dlerror();
    EXPECT_STREQ(abi(), HYDROCOUPLE_COMPONENT_ABI_STAMP);

    auto info = reinterpret_cast<HydroCoupleComponentInfoFn>(
        dlsym(handle, HYDROCOUPLE_COMPONENT_INFO_SYMBOL));
    ASSERT_NE(info, nullptr) << dlerror();
    auto *modelInfo =
        dynamic_cast<HydroCouple::IModelComponentInfo *>(info());
    ASSERT_NE(modelInfo, nullptr)
        << "cross-image dynamic_cast failed: interface RTTI is not shared";
    EXPECT_EQ(modelInfo->id(), "fvqual.component");

    using Status = HydroCouple::IModelComponent::ComponentStatus;
    {
        std::unique_ptr<HydroCouple::IModelComponent> component =
            modelInfo->createComponentInstance();
        ASSERT_NE(component, nullptr);

        component->initialize();
        ASSERT_EQ(component->status(), Status::Initialized);

        // The ports a composition wires: the met input, the field output.
        ASSERT_EQ(component->inputs().size(), 1u);
        EXPECT_EQ(component->inputs()[0]->id(), "air_temperature");
        ASSERT_EQ(component->outputs().size(), 1u);
        EXPECT_EQ(component->outputs()[0]->id(), "temperature");

        // The default rectangular basin must run — the contract that lets
        // a host instantiate the component just to see it.
        EXPECT_TRUE(component->validate().empty());
        component->prepare();
        int steps = 0;
        while (component->status() != Status::Done &&
               component->status() != Status::Failed && steps < 100)
        {
            component->update();
            ++steps;
        }
        component->finish();
        ASSERT_EQ(component->status(), Status::Finished);
    } // the instance dies before its library does
    dlclose(handle);
}
