// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/state/statestore.h"

#include <cstring>

namespace FVQual::State {

namespace IO = HydroCouple::SDK::IO;

void StateStore::allocate(const Mesh::LayeredMesh &mesh, const SpeciesRegistry &registry)
{
    m_columnCount    = mesh.columnCount();
    m_cellCount      = mesh.cellCount();
    m_interfaceSlots = mesh.interfaceSlotCount();
    m_layerCount     = mesh.layerCount();
    m_speciesCount   = registry.count();

    const size_t columns    = static_cast<size_t>(m_columnCount);
    const size_t cells      = static_cast<size_t>(m_cellCount);
    const size_t interfaces = static_cast<size_t>(m_interfaceSlots);

    eta.assign(columns, 0.0);
    ubar.assign(columns, 0.0);
    vbar.assign(columns, 0.0);

    u.assign(cells, 0.0);
    v.assign(cells, 0.0);
    density.assign(cells, 0.0);

    w.assign(interfaces, 0.0);
    kz.assign(interfaces, 0.0);
    km.assign(interfaces, 0.0);

    species.assign(static_cast<size_t>(m_speciesCount) * cells, 0.0);

    // Initial free surface comes from the mesh the state was sized against.
    const std::vector<double> &surface = mesh.surfaceElevation();
    if (surface.size() == columns)
        eta.assign(surface.begin(), surface.end());

    for (int s = 0; s < m_speciesCount; ++s)
    {
        const double initial = registry.definition(s).initialValue;

        if (initial == 0.0)
            continue;

        std::span<double> block = speciesBlock(s);
        for (double &value : block)
            value = initial;
    }
}

double StateStore::totalMass(int s, const Mesh::LayeredMesh &mesh) const noexcept
{
    const std::vector<double> &volume = mesh.cellVolume();
    const std::span<const double> block = speciesBlock(s);

    double total = 0.0;

    for (int64_t i = 0; i < m_cellCount; ++i)
        total += block[static_cast<size_t>(i)] * volume[static_cast<size_t>(i)];

    return total;
}

double StateStore::volumeWeightedMean(int s, const Mesh::LayeredMesh &mesh) const noexcept
{
    const double volume = mesh.totalVolume();

    return volume > 0.0 ? totalMass(s, mesh) / volume : 0.0;
}

IO::FieldSlice StateStore::layeredSlice(const std::string &itemId,
                                        std::span<const double> values) const
{
    IO::FieldSlice slice;
    slice.itemId = itemId;
    slice.kind   = HydroCouple::DataKind::Float64;
    slice.shape  = {m_columnCount, static_cast<int64_t>(m_layerCount)};

    slice.data.resize(values.size() * sizeof(double));
    std::memcpy(slice.data.data(), values.data(), values.size() * sizeof(double));

    return slice;
}

IO::FieldSlice StateStore::columnSlice(const std::string &itemId,
                                       std::span<const double> values) const
{
    IO::FieldSlice slice;
    slice.itemId = itemId;
    slice.kind   = HydroCouple::DataKind::Float64;
    slice.shape  = {m_columnCount};

    slice.data.resize(values.size() * sizeof(double));
    std::memcpy(slice.data.data(), values.data(), values.size() * sizeof(double));

    return slice;
}

} // namespace FVQual::State
