// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/component/sigmaverticalcoordinate.h"

#include "fvqual/mesh/layeredmesh.h"

#include <limits>

namespace FVQual::Component {

using HydroCouple::Spatial::VerticalCoordinateKind;

namespace {

constexpr double kBadElevation = std::numeric_limits<double>::quiet_NaN();

} // namespace

SigmaVerticalCoordinate::SigmaVerticalCoordinate(
    const FVQual::Mesh::LayeredMesh &mesh) noexcept
    : m_mesh(&mesh)
{
}

VerticalCoordinateKind SigmaVerticalCoordinate::kind() const
{
    return VerticalCoordinateKind::Sigma;
}

int64_t SigmaVerticalCoordinate::layerCount() const
{
    return m_mesh ? m_mesh->layerCount() : 0;
}

bool SigmaVerticalCoordinate::isTimeVarying() const
{
    return true;
}

double SigmaVerticalCoordinate::interfaceElevation(int64_t cellIndex,
                                                   int64_t interfaceIndex) const
{
    if (m_mesh == nullptr || cellIndex < 0 || cellIndex >= m_mesh->columnCount() ||
        interfaceIndex < 0 || interfaceIndex >= m_mesh->interfaceCount())
    {
        return kBadElevation;
    }

    // Delegated rather than recomputed from sigma, eta and the bed. The mesh
    // already holds the elevations its own physics integrates over, and a
    // second derivation here would be free to disagree with them — silently,
    // and only for the stretched distributions where the arithmetic differs.
    // What a coupling partner needs is the mesh the model actually solved on.
    return m_mesh->z(cellIndex, static_cast<int>(interfaceIndex));
}

int64_t SigmaVerticalCoordinate::columnCount() const
{
    return m_mesh ? m_mesh->columnCount() : 0;
}

uint64_t SigmaVerticalCoordinate::geometryEpoch() const
{
    return m_epochOffset + (m_mesh ? m_mesh->geometryEpoch() : 0);
}

std::span<const double> SigmaVerticalCoordinate::interfaceElevations() const
{
    return m_mesh ? std::span<const double>(m_mesh->interfaceZ())
                  : std::span<const double>();
}

void SigmaVerticalCoordinate::refresh(const FVQual::Mesh::LayeredMesh &mesh) noexcept
{
    const uint64_t nextEpoch = geometryEpoch() + 1;
    m_mesh = &mesh;
    m_epochOffset = nextEpoch - mesh.geometryEpoch();
}

} // namespace FVQual::Component
