// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/mesh/sigmaspec.h"

#include "hydrocouplesdk/tools/sigmagrid.h"

namespace FVQual::Mesh {

namespace Tools = HydroCouple::SDK::Tools;

int SigmaSpec::resolvedLayerCount() const noexcept
{
    if (distribution == Distribution::Explicit)
    {
        return explicitInterfaces.size() >= 2
                   ? static_cast<int>(explicitInterfaces.size()) - 1
                   : 0;
    }

    return layerCount;
}

bool SigmaSpec::resolve(std::vector<double> &interfaces, std::string &message) const
{
    if (minLayerThickness < 0.0)
    {
        message = "minLayerThickness must be non-negative";
        return false;
    }

    switch (distribution)
    {
    case Distribution::Uniform:
        if (layerCount < 1)
        {
            message = "layerCount must be at least 1";
            return false;
        }
        interfaces = Tools::SigmaLayering::uniform(layerCount);
        break;

    case Distribution::Stretched:
        if (layerCount < 1)
        {
            message = "layerCount must be at least 1";
            return false;
        }
        if (theta <= 0.0)
        {
            message = "stretched layering requires theta > 0";
            return false;
        }
        if (b < 0.0 || b > 1.0)
        {
            message = "stretching blend b must lie in [0, 1]";
            return false;
        }
        interfaces = Tools::SigmaLayering::stretched(layerCount, theta, b);
        break;

    case Distribution::Explicit:
        if (explicitInterfaces.size() < 2)
        {
            message = "explicit layering needs at least 2 interfaces";
            return false;
        }
        interfaces = explicitInterfaces;
        break;
    }

    if (!Tools::SigmaLayering::validate(interfaces, message))
        return false;

    message.clear();
    return true;
}

} // namespace FVQual::Mesh
