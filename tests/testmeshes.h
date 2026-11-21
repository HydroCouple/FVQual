// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  testmeshes.h
 * \brief Small hand-checkable meshes shared by the test suite.
 *
 * Every mesh here has an area, a centroid, and an edge count that can be
 * verified by hand, so a failure points at the code rather than at the
 * fixture.
 */

#ifndef FVQual_TESTS_TESTMESHES_H
#define FVQual_TESTS_TESTMESHES_H

#include "hydrocouplesdk/io/meshdefinition.h"

namespace FVQual::Testing {

/*!
 * \brief One counter-clockwise unit right triangle: (0,0), (1,0), (0,1).
 *        Area 1/2, centroid (1/3, 1/3), 3 boundary edges.
 */
inline HydroCouple::SDK::IO::MeshDefinition unitTriangle()
{
    HydroCouple::SDK::IO::MeshDefinition mesh;
    mesh.meshName        = "mesh2d";
    mesh.nodeX           = {0.0, 1.0, 0.0};
    mesh.nodeY           = {0.0, 0.0, 1.0};
    mesh.faceNodeOffsets = {0, 3};
    mesh.faceNodes       = {0, 1, 2};
    return mesh;
}

/*!
 * \brief The unit square split into two triangles sharing the diagonal 0-2.
 *        Area 1, 5 unique edges of which exactly 1 is interior.
 */
inline HydroCouple::SDK::IO::MeshDefinition twoTriangleSquare()
{
    HydroCouple::SDK::IO::MeshDefinition mesh;
    mesh.meshName        = "mesh2d";
    mesh.nodeX           = {0.0, 1.0, 1.0, 0.0};
    mesh.nodeY           = {0.0, 0.0, 1.0, 1.0};
    mesh.faceNodeOffsets = {0, 3, 6};
    mesh.faceNodes       = {0, 1, 2, 0, 2, 3};
    return mesh;
}

/*!
 * \brief A quad and a triangle sharing one edge — the mixed tri-quad case
 *        that exercises padded UGRID connectivity.
 *
 * Quad  (0,0) (2,0) (2,1) (0,1)  area 2, and
 * tri   (2,0) (3,0) (2,1)        area 1/2, sharing edge 1-2.
 */
inline HydroCouple::SDK::IO::MeshDefinition quadPlusTriangle()
{
    HydroCouple::SDK::IO::MeshDefinition mesh;
    mesh.meshName        = "mesh2d";
    mesh.nodeX           = {0.0, 2.0, 2.0, 0.0, 3.0};
    mesh.nodeY           = {0.0, 0.0, 1.0, 1.0, 0.0};
    mesh.faceNodeOffsets = {0, 4, 7};
    mesh.faceNodes       = {0, 1, 2, 3, 1, 4, 2};
    return mesh;
}

/*!
 * \brief A regular strip of \p columns quads, each 1 x 1, running along x.
 *        Area = columns; interior edges = columns - 1.
 */
inline HydroCouple::SDK::IO::MeshDefinition quadStrip(int columns)
{
    HydroCouple::SDK::IO::MeshDefinition mesh;
    mesh.meshName = "mesh2d";

    for (int i = 0; i <= columns; ++i)
    {
        mesh.nodeX.push_back(static_cast<double>(i));
        mesh.nodeY.push_back(0.0);
        mesh.nodeX.push_back(static_cast<double>(i));
        mesh.nodeY.push_back(1.0);
    }

    mesh.faceNodeOffsets.push_back(0);

    for (int i = 0; i < columns; ++i)
    {
        const int64_t bottomLeft  = 2 * i;
        const int64_t topLeft     = 2 * i + 1;
        const int64_t bottomRight = 2 * (i + 1);
        const int64_t topRight    = 2 * (i + 1) + 1;

        mesh.faceNodes.push_back(bottomLeft);
        mesh.faceNodes.push_back(bottomRight);
        mesh.faceNodes.push_back(topRight);
        mesh.faceNodes.push_back(topLeft);
        mesh.faceNodeOffsets.push_back(static_cast<int64_t>(mesh.faceNodes.size()));
    }

    return mesh;
}

/*!
 * \brief A straight channel of \p columns quads spanning \p length in x and
 *        \p width in y — the domain for the 1-D analytical benchmarks
 *        (seiche, wind setup, dam break).
 *
 * Cell centers sit at (i + 1/2) * length / columns.
 */
inline HydroCouple::SDK::IO::MeshDefinition channelStrip(int columns, double length,
                                                         double width)
{
    HydroCouple::SDK::IO::MeshDefinition mesh;
    mesh.meshName = "mesh2d";

    const double dx = length / static_cast<double>(columns);

    for (int i = 0; i <= columns; ++i)
    {
        const double x = static_cast<double>(i) * dx;

        mesh.nodeX.push_back(x);
        mesh.nodeY.push_back(0.0);
        mesh.nodeX.push_back(x);
        mesh.nodeY.push_back(width);
    }

    mesh.faceNodeOffsets.push_back(0);

    for (int i = 0; i < columns; ++i)
    {
        mesh.faceNodes.push_back(2 * i);
        mesh.faceNodes.push_back(2 * (i + 1));
        mesh.faceNodes.push_back(2 * (i + 1) + 1);
        mesh.faceNodes.push_back(2 * i + 1);
        mesh.faceNodeOffsets.push_back(static_cast<int64_t>(mesh.faceNodes.size()));
    }

    return mesh;
}

/*!
 * \brief A rectangular basin of \p nx by \p ny quads.
 *
 * The channel strip above is one cell wide, which is the right domain for the
 * 1-D analytical benchmarks and the wrong one for anything that circulates. A
 * wind-driven overturning cell in a one-cell-wide chain has no lateral return
 * path, so every return flow is forced through the vertical — which is exactly
 * the suspicion hanging over DeGray's interior warm bias, and exactly what a
 * laterally-averaged schematisation cannot be used to test.
 *
 * Cell centers sit at ((i + 1/2) dx, (j + 1/2) dy).
 */
inline HydroCouple::SDK::IO::MeshDefinition basinGrid(int nx, int ny, double length,
                                                      double width)
{
    HydroCouple::SDK::IO::MeshDefinition mesh;
    mesh.meshName = "mesh2d";

    const double dx = length / static_cast<double>(nx);
    const double dy = width / static_cast<double>(ny);

    const auto node = [&](int i, int j) { return j * (nx + 1) + i; };

    for (int j = 0; j <= ny; ++j)
    {
        for (int i = 0; i <= nx; ++i)
        {
            mesh.nodeX.push_back(static_cast<double>(i) * dx);
            mesh.nodeY.push_back(static_cast<double>(j) * dy);
        }
    }

    mesh.faceNodeOffsets.push_back(0);

    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {
            mesh.faceNodes.push_back(node(i, j));
            mesh.faceNodes.push_back(node(i + 1, j));
            mesh.faceNodes.push_back(node(i + 1, j + 1));
            mesh.faceNodes.push_back(node(i, j + 1));

            mesh.faceNodeOffsets.push_back(static_cast<int64_t>(mesh.faceNodes.size()));
        }
    }

    return mesh;
}

} // namespace FVQual::Testing

#endif // FVQual_TESTS_TESTMESHES_H
