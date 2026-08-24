// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/mesh/sigmaspec.h"

#include <gtest/gtest.h>

using FVQual::Mesh::SigmaSpec;

TEST(SigmaSpec, UniformProducesEvenlySpacedInterfaces)
{
    SigmaSpec spec;
    spec.distribution = SigmaSpec::Distribution::Uniform;
    spec.layerCount   = 4;

    std::vector<double> interfaces;
    std::string         message;

    ASSERT_TRUE(spec.resolve(interfaces, message)) << message;
    ASSERT_EQ(interfaces.size(), 5u);

    EXPECT_DOUBLE_EQ(interfaces.front(), 0.0);
    EXPECT_DOUBLE_EQ(interfaces.back(), 1.0);

    for (size_t k = 0; k + 1 < interfaces.size(); ++k)
        EXPECT_NEAR(interfaces[k + 1] - interfaces[k], 0.25, 1e-12);
}

TEST(SigmaSpec, InterfacesAreStrictlyIncreasing)
{
    for (SigmaSpec::Distribution distribution :
         {SigmaSpec::Distribution::Uniform, SigmaSpec::Distribution::Stretched})
    {
        SigmaSpec spec;
        spec.distribution = distribution;
        spec.layerCount   = 12;
        spec.theta        = 4.0;
        spec.b            = 0.0;

        std::vector<double> interfaces;
        std::string         message;

        ASSERT_TRUE(spec.resolve(interfaces, message)) << message;
        ASSERT_EQ(interfaces.size(), 13u);

        for (size_t k = 0; k + 1 < interfaces.size(); ++k)
            EXPECT_LT(interfaces[k], interfaces[k + 1]) << "interface " << k;
    }
}

TEST(SigmaSpec, SurfaceClusteringPutsThinnerLayersOnTop)
{
    // b = 0 clusters toward the surface, where a stratified reservoir keeps
    // its gradient. The first layer must therefore be thinner than the last.
    SigmaSpec spec;
    spec.distribution = SigmaSpec::Distribution::Stretched;
    spec.layerCount   = 10;
    spec.theta        = 6.0;
    spec.b            = 0.0;

    std::vector<double> interfaces;
    std::string         message;

    ASSERT_TRUE(spec.resolve(interfaces, message)) << message;

    const double firstThickness = interfaces[1] - interfaces[0];
    const double lastThickness =
        interfaces[interfaces.size() - 1] - interfaces[interfaces.size() - 2];

    EXPECT_LT(firstThickness, lastThickness);
}

TEST(SigmaSpec, ExplicitInterfacesArePassedThrough)
{
    SigmaSpec spec;
    spec.distribution       = SigmaSpec::Distribution::Explicit;
    spec.explicitInterfaces = {0.0, 0.1, 0.4, 1.0};

    std::vector<double> interfaces;
    std::string         message;

    ASSERT_TRUE(spec.resolve(interfaces, message)) << message;
    EXPECT_EQ(interfaces, spec.explicitInterfaces);
    EXPECT_EQ(spec.resolvedLayerCount(), 3);
}

TEST(SigmaSpec, MalformedSpecsAreRejected)
{
    std::vector<double> interfaces;
    std::string         message;

    {
        SigmaSpec spec;
        spec.layerCount = 0;
        EXPECT_FALSE(spec.resolve(interfaces, message));
        EXPECT_FALSE(message.empty());
    }

    {
        SigmaSpec spec;
        spec.distribution = SigmaSpec::Distribution::Stretched;
        spec.layerCount   = 5;
        spec.theta        = 0.0; // stretching needs theta > 0
        EXPECT_FALSE(spec.resolve(interfaces, message));
    }

    {
        SigmaSpec spec;
        spec.distribution       = SigmaSpec::Distribution::Explicit;
        spec.explicitInterfaces = {0.0, 0.5, 0.4, 1.0}; // not increasing
        EXPECT_FALSE(spec.resolve(interfaces, message));
    }

    {
        SigmaSpec spec;
        spec.minLayerThickness = -1.0;
        EXPECT_FALSE(spec.resolve(interfaces, message));
    }
}
