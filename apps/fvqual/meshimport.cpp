// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin
#include "hydrocouplesdk/io/meshpackage.h"
#include "fvqual/mesh/horizontalmesh.h"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

int importMeshCommand(int argc, char **argv)
{
    if (argc < 4) {
        std::cerr << "usage: fvqual import-mesh <SWMM-or-SMS-file> <output.mesh.json>"
                     " --crs <projected CRS or LOCAL> --datum <vertical datum>"
                     " --elevation-kind <bed-elevation|terrain-elevation>"
                     " [--xy-units m|ft|us-survey-ft] [--z-units ...] [--depth-units ...]\n";
        return EXIT_FAILURE;
    }
    HydroCouple::SDK::IO::MeshImportOptions options;
    try {
        for (int i = 4; i < argc; i += 2) {
            if (i + 1 >= argc) throw std::invalid_argument("mesh import option needs a value");
            const std::string flag = argv[i], value = argv[i + 1];
            if (flag == "--crs") options.horizontalCrs = value;
            else if (flag == "--datum") options.verticalDatum = value;
            else if (flag == "--elevation-kind") options.elevationKind = value;
            else if (flag == "--xy-units") options.horizontalUnits = value;
            else if (flag == "--z-units") options.verticalUnits = value;
            else if (flag == "--depth-units") options.depthUnits = value;
            else throw std::invalid_argument("unknown mesh import option: " + flag);
        }
        HydroCouple::SDK::IO::MeshPackage package;
        std::string message;
        if (!HydroCouple::SDK::IO::importSurfaceMesh(argv[2], options, package, message))
            throw std::invalid_argument(message);
        FVQual::Mesh::HorizontalMesh horizontal;
        if (!horizontal.build(package.mesh, message)) throw std::invalid_argument(message);
        if (!HydroCouple::SDK::IO::writeMeshPackage(argv[3], package, message))
            throw std::runtime_error(message);
        int64_t triangles = 0, quads = 0;
        for (int64_t f = 0; f < horizontal.faceCount(); ++f) {
            triangles += horizontal.faceNodeCount(f) == 3;
            quads += horizontal.faceNodeCount(f) == 4;
        }
        std::cout << "saved " << argv[3] << "\ntriangles " << triangles
                  << "\nquads " << quads << "\nboundaries " << package.boundaries.size()
                  << "\ncoordinate units m; elevation positive up in " << package.verticalDatum << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}

int inspectMeshCommand(int argc, char **argv)
{
    if (argc < 3) {
        std::cerr << "usage: fvqual inspect-mesh <package.mesh.json> [--report <report.json>]\n";
        return EXIT_FAILURE;
    }
    try {
        std::string reportPath;
        for (int i = 3; i < argc; i += 2) {
            if (std::string(argv[i]) != "--report" || i + 1 >= argc)
                throw std::invalid_argument("inspect-mesh accepts only --report <new report.json>");
            reportPath = argv[i+1];
        }
        HydroCouple::SDK::IO::MeshPackage package;
        std::string message;
        if (!HydroCouple::SDK::IO::readMeshPackage(argv[2], package, message))
            throw std::invalid_argument(message);
        FVQual::Mesh::HorizontalMesh mesh;
        if (!mesh.build(package.mesh, message)) throw std::invalid_argument(message);
        std::map<std::pair<int64_t,int64_t>, int64_t> edges;
        for (int64_t e = 0; e < mesh.edgeCount(); ++e) {
            auto a = mesh.edgeNodes()[2*e], b = mesh.edgeNodes()[2*e+1];
            edges.emplace(std::minmax(a,b), e);
        }
        for (const auto &boundary : package.boundaries) {
            const auto it = edges.find(std::minmax(boundary.node0,boundary.node1));
            if (it == edges.end() || !mesh.isBoundaryEdge(it->second))
                throw std::invalid_argument("declared boundary endpoints must identify an exterior mesh edge");
        }
        const auto r = mesh.geometryReport();
        nlohmann::json report{
            {"schema", "fvqual-mesh-qualification"}, {"version", "1.0"},
            {"package", std::filesystem::absolute(argv[2]).string()},
            {"geometry_valid", true}, {"nodes",mesh.nodeCount()}, {"faces",mesh.faceCount()},
            {"edges",mesh.edgeCount()}, {"triangles",r.triangles}, {"quads",r.quads},
            {"polygons",r.polygons}, {"boundary_edges",r.boundaryEdges},
            {"boundary_loops",r.boundaryLoops}, {"connected_domains",r.connectedDomains},
            {"mapped_boundaries",package.boundaries.size()}, {"area_m2",mesh.totalArea()},
            {"minimum_area_m2",r.minimumArea}, {"maximum_area_m2",r.maximumArea},
            {"minimum_normal_distance_m",r.minimumNormalDistance},
            {"maximum_nonorthogonality_degrees",r.maximumNonOrthogonalityDegrees},
            {"maximum_relative_closure_error",r.maximumRelativeClosureError},
            {"horizontal_crs",package.horizontalCrs}, {"vertical_datum",package.verticalDatum},
            {"elevation_kind",package.elevationKind}, {"coordinate_units","m"},
            {"quad_conversion","none"}, {"provenance",package.provenance},
            {"reservoir_ready",false},
            {"qualification_note","Geometry inspection does not establish bathymetry, wet columns or solver convergence. Boundary loops require comparison with the intended outline; holes are not automatically gaps."}};
        if (!reportPath.empty()) {
            if (std::filesystem::exists(reportPath))
                throw std::invalid_argument("report exists; choose a new path");
            std::ofstream out(reportPath);
            out << report.dump(2) << '\n'; out.close();
            if (!out) throw std::runtime_error("cannot write mesh report");
        }
        std::cout << report.dump(2) << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << "error: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
