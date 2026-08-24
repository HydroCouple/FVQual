// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  main.cpp
 * \brief The fvqual command-line runner.
 *
 * Phase 1 exposes the mesh pipeline: read a UGRID mesh, extrude it through a
 * sigma layering, and report the geometry a modeller needs to sanity-check a
 * domain before any physics exists — cell counts, volumes, and the
 * hydrostatic-consistency screening that predicts sigma pressure-gradient
 * trouble. `run` arrives with the hydrodynamic core in Phase 3.
 */

#include "fvqual/io/ugridreader.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/state/statestore.h"
#include "fvqual/version.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <string_view>

namespace {

using FVQual::IO::UgridReader;
using FVQual::Mesh::HorizontalMesh;
using FVQual::Mesh::LayeredMesh;
using FVQual::Mesh::SigmaSpec;

void printUsage()
{
    std::cout
        << "fvqual " << FVQual_FULL_VERSION << "\n"
        << "3D unstructured finite-volume hydrodynamic and water-quality model\n\n"
        << "Usage:\n"
        << "  fvqual mesh <mesh.nc> [options]   Inspect a UGRID mesh and its sigma grid\n"
        << "  fvqual version                    Print the version\n"
        << "  fvqual help                       Print this message\n\n"
        << "Options for 'mesh':\n"
        << "  --layers <n>        Number of sigma layers (default 10)\n"
        << "  --theta <value>     Song-Haidvogel stretching intensity (default 0 = uniform)\n"
        << "  --surface <z>       Uniform initial water-surface elevation, m (default 0)\n"
        << "  --min-thickness <t> Minimum layer thickness, m (default 0.05)\n"
        << "  --haney-limit <r>   Steepness-ratio limit for screening (default 0.2)\n"
        << "  --mesh-name <name>  UGRID mesh topology variable (default: first found)\n";
}

/*!
 * \brief Reads the value following \p index, or reports a missing argument.
 */
bool argumentValue(int argc, char **argv, int index, std::string &value)
{
    if (index + 1 >= argc)
    {
        std::cerr << "error: " << argv[index] << " needs a value\n";
        return false;
    }

    value = argv[index + 1];
    return true;
}

int runMeshCommand(int argc, char **argv)
{
    if (argc < 3)
    {
        std::cerr << "error: 'mesh' needs a NetCDF file path\n";
        return EXIT_FAILURE;
    }

    const std::string path = argv[2];

    SigmaSpec   spec;
    double      surface    = 0.0;
    double      haneyLimit = 0.2;
    std::string meshName;

    spec.layerCount        = 10;
    spec.minLayerThickness = 0.05;

    for (int i = 3; i < argc; ++i)
    {
        const std::string_view flag = argv[i];
        std::string            value;

        if (flag == "--layers")
        {
            if (!argumentValue(argc, argv, i, value))
                return EXIT_FAILURE;
            spec.layerCount = std::stoi(value);
            ++i;
        }
        else if (flag == "--theta")
        {
            if (!argumentValue(argc, argv, i, value))
                return EXIT_FAILURE;
            spec.theta        = std::stod(value);
            spec.distribution = spec.theta > 0.0 ? SigmaSpec::Distribution::Stretched
                                                 : SigmaSpec::Distribution::Uniform;
            ++i;
        }
        else if (flag == "--surface")
        {
            if (!argumentValue(argc, argv, i, value))
                return EXIT_FAILURE;
            surface = std::stod(value);
            ++i;
        }
        else if (flag == "--min-thickness")
        {
            if (!argumentValue(argc, argv, i, value))
                return EXIT_FAILURE;
            spec.minLayerThickness = std::stod(value);
            ++i;
        }
        else if (flag == "--haney-limit")
        {
            if (!argumentValue(argc, argv, i, value))
                return EXIT_FAILURE;
            haneyLimit = std::stod(value);
            ++i;
        }
        else if (flag == "--mesh-name")
        {
            if (!argumentValue(argc, argv, i, value))
                return EXIT_FAILURE;
            meshName = value;
            ++i;
        }
        else
        {
            std::cerr << "error: unknown option '" << flag << "'\n";
            return EXIT_FAILURE;
        }
    }

    if (!UgridReader::available())
    {
        std::cerr << "error: this build has NetCDF support disabled; "
                     "reconfigure with -DFVQual_USE_NETCDF=ON\n";
        return EXIT_FAILURE;
    }

    HydroCouple::SDK::IO::MeshDefinition definition;
    std::string                          message;

    if (!UgridReader::read(path, meshName, definition, message))
    {
        std::cerr << "error: " << message << "\n";
        return EXIT_FAILURE;
    }

    HorizontalMesh horizontal;
    if (!horizontal.build(definition, message))
    {
        std::cerr << "error: " << message << "\n";
        return EXIT_FAILURE;
    }

    if (definition.nodeZ.empty())
    {
        std::cerr << "error: the mesh carries no node elevations; FVQual needs "
                     "bathymetry to build a sigma grid\n";
        return EXIT_FAILURE;
    }

    LayeredMesh layered;
    if (!layered.build(horizontal, spec, surface, message))
    {
        std::cerr << "error: " << message << "\n";
        return EXIT_FAILURE;
    }

    const auto haney = layered.haneyReport(haneyLimit);

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "mesh          " << definition.meshName << "\n"
              << "nodes         " << horizontal.nodeCount() << "\n"
              << "faces         " << horizontal.faceCount() << "\n"
              << "edges         " << horizontal.edgeCount() << "\n"
              << "layers        " << layered.layerCount() << "\n"
              << "cells         " << layered.cellCount() << "\n"
              << "plan area     " << horizontal.totalArea() << " m2\n"
              << "water volume  " << layered.totalVolume() << " m3\n"
              << "min centroid  " << horizontal.minimumCentroidDistance() << " m\n"
              << haney.summary() << "\n";

    return haney.passed() ? EXIT_SUCCESS : EXIT_FAILURE;
}

} // namespace

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        printUsage();
        return EXIT_FAILURE;
    }

    const std::string_view command = argv[1];

    if (command == "help" || command == "--help" || command == "-h")
    {
        printUsage();
        return EXIT_SUCCESS;
    }

    if (command == "version" || command == "--version")
    {
        std::cout << FVQual_FULL_VERSION << "\n";
        return EXIT_SUCCESS;
    }

    if (command == "mesh")
        return runMeshCommand(argc, argv);

    std::cerr << "error: unknown command '" << command << "'\n\n";
    printUsage();
    return EXIT_FAILURE;
}
