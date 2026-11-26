// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
//
// The plugin gate: FVQual loaded exactly as a host loads it — dlopen, the
// ABI stamp compared before any C++ symbol, the info entry point cast
// across the image boundary, and a created instance initialized far
// enough to show its ports and run its default basin to Finished.

#include <gtest/gtest.h>

#include "hydrocouplesdk/componentabi.h"
#include "hydrocouplesdk/differentiation/differentials.h"

#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#include <filesystem>
#include <memory>
#include <string>

namespace {

std::filesystem::path pluginPath()
{
    return std::filesystem::path(FVQUAL_PLUGIN_PATH);
}

// What a host does with a component library: open it, look up a symbol,
// close it. dlopen/dlsym elsewhere, LoadLibrary/GetProcAddress on Windows.
#ifdef _WIN32
void *openLibrary(const std::filesystem::path &path)
{
    return LoadLibraryW(path.c_str());
}
void *librarySymbol(void *library, const char *name)
{
    return reinterpret_cast<void *>(
        GetProcAddress(static_cast<HMODULE>(library), name));
}
void closeLibrary(void *library)
{
    FreeLibrary(static_cast<HMODULE>(library));
}
std::string libraryError()
{
    return "Windows error " + std::to_string(GetLastError());
}
#else
void *openLibrary(const std::filesystem::path &path)
{
    return dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
}
void *librarySymbol(void *library, const char *name)
{
    return dlsym(library, name);
}
void closeLibrary(void *library)
{
    dlclose(library);
}
std::string libraryError()
{
    const char *error = dlerror();
    return error ? error : "";
}
#endif

} // namespace

TEST(FVQualPluginTest, ConditionalDifferentiationCrossesThePluginBoundary) {
    void *handle=openLibrary(pluginPath());ASSERT_NE(handle,nullptr)<<libraryError();
    const auto info=reinterpret_cast<HydroCoupleComponentInfoFn>(librarySymbol(handle, HYDROCOUPLE_COMPONENT_INFO_SYMBOL));ASSERT_NE(info,nullptr);
    auto *modelInfo=dynamic_cast<HydroCouple::IModelComponentInfo *>(info());ASSERT_NE(modelInfo,nullptr);
    {
        auto component=modelInfo->createComponentInstance();std::string message;
        const auto assign=[&](const char *id,const std::string &text){
            for(auto *a:component->arguments())if(a->id()==id){EXPECT_TRUE(a->initialize(text,HydroCouple::IArgument::ArgumentInputType::JSON,message))<<message;return;}
            FAIL()<<"missing argument "<<id;
        };
        EXPECT_FALSE(component->capabilities().count(HydroCouple::Capability::Differentiable));
        assign("geometry_counts","[2,2]");assign("geometry_values","[20,10,5,12,1,2]");
        assign("conservative_species","{\"tracer\":0.002}");assign("conditional_species_derivatives","[1]");
        component->initialize();component->prepare();ASSERT_NE(component->status(),HydroCouple::IModelComponent::ComponentStatus::Failed);
        auto *d=dynamic_cast<HydroCouple::IDifferentiableModelComponent *>(component.get());ASSERT_NE(d,nullptr);
        EXPECT_TRUE(component->capabilities().count(HydroCouple::Capability::Differentiable));ASSERT_EQ(d->differentiableOutputs().size(),1u);
        EXPECT_TRUE(d->differentiableArguments().empty());component->update();auto *item=d->differentiableOutputs()[0];
        const auto shape=item->shape();std::vector<double> direction(4,1),result(4,0);
        const HydroCouple::DifferentialEntry seed{item,HydroCouple::DifferentialRole::StateBefore,HydroCouple::SDK::Differentiation::hostBuffer(direction,shape)};
        const HydroCouple::DifferentialEntry output{item,HydroCouple::DifferentialRole::Output,HydroCouple::SDK::Differentiation::hostBuffer(result,shape)};
        ASSERT_TRUE(d->jvp({&seed,1},{&output,1},&message))<<message;for(double value:result)EXPECT_NEAR(value,1,2e-14);
        component->finish();
    }
    closeLibrary(handle);
}

TEST(FVQualPluginTest, LoadsThroughTheStampedConventionAndRunsItsDefaults)
{
    void *handle = openLibrary(pluginPath());
    ASSERT_NE(handle, nullptr) << libraryError();

    auto abi = reinterpret_cast<HydroCoupleComponentAbiFn>(
        librarySymbol(handle, HYDROCOUPLE_COMPONENT_ABI_SYMBOL));
    ASSERT_NE(abi, nullptr) << libraryError();
    EXPECT_STREQ(abi(), HYDROCOUPLE_COMPONENT_ABI_STAMP);

    auto info = reinterpret_cast<HydroCoupleComponentInfoFn>(
        librarySymbol(handle, HYDROCOUPLE_COMPONENT_INFO_SYMBOL));
    ASSERT_NE(info, nullptr) << libraryError();
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
        ASSERT_EQ(component->inputs().size(), 2u);
        EXPECT_EQ(component->inputs()[0]->id(), "air_temperature");
        ASSERT_GE(component->outputs().size(), 10u);
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
    closeLibrary(handle);
}
