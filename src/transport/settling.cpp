// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/transport/settling.h"

#include "fvqual/parallel.h"

#include <algorithm>

namespace FVQual::Transport {

void SettlingSet::initialize(const State::SpeciesRegistry &species)
{
    const size_t count = static_cast<size_t>(species.count());

    m_velocity.assign(count, 0.0);
    m_receiver.assign(count, -1);
    m_delivered.assign(count, 0.0);
    m_lost.assign(count, 0.0);
}

bool SettlingSet::setVelocity(int species, double velocity,
                              const State::SpeciesRegistry &registry,
                              std::string &message)
{
    if (species < 0 || species >= static_cast<int>(m_velocity.size()))
    {
        message = "species index " + std::to_string(species) +
                  " is out of range for the settling set";
        return false;
    }

    if (!(velocity >= 0.0))
    {
        // Upward transport of a particulate is buoyancy, not settling, and the
        // forward substitution below assumes the flux is downward — a negative
        // velocity would make the matrix upper triangular and this solve wrong
        // rather than merely inaccurate.
        message = "settling velocity must be non-negative; buoyant rise is not "
                  "modelled here";
        return false;
    }

    if (!registry.definition(species).isTransported())
    {
        message = "species '" + registry.definition(species).id +
                  "' is benthic and already on the bed; it cannot settle";
        return false;
    }

    m_velocity[static_cast<size_t>(species)] = velocity;

    message.clear();
    return true;
}

bool SettlingSet::setBenthicReceiver(int species, int benthic,
                                     const State::SpeciesRegistry &registry,
                                     std::string &message)
{
    if (species < 0 || species >= static_cast<int>(m_receiver.size()))
    {
        message = "species index " + std::to_string(species) +
                  " is out of range for the settling set";
        return false;
    }

    if (benthic < 0 || benthic >= static_cast<int>(m_receiver.size()))
    {
        message = "benthic species index " + std::to_string(benthic) +
                  " is out of range";
        return false;
    }

    if (registry.definition(benthic).isTransported())
    {
        message = "'" + registry.definition(benthic).id +
                  "' is transported, so material settled into it would be advected "
                  "away again; a benthic receiver must be SpeciesKind::MsxBenthic";
        return false;
    }

    m_receiver[static_cast<size_t>(species)] = benthic;

    message.clear();
    return true;
}

double SettlingSet::velocity(int species) const
{
    if (species < 0 || species >= static_cast<int>(m_velocity.size()))
        return 0.0;

    return m_velocity[static_cast<size_t>(species)];
}

int SettlingSet::benthicReceiver(int species) const
{
    if (species < 0 || species >= static_cast<int>(m_receiver.size()))
        return -1;

    return m_receiver[static_cast<size_t>(species)];
}

bool SettlingSet::any() const noexcept
{
    return std::any_of(m_velocity.begin(), m_velocity.end(),
                       [](double w) { return w > 0.0; });
}

void SettlingSet::apply(const Mesh::LayeredMesh &mesh, State::StateStore &state,
                        double dt)
{
    const int64_t columns = mesh.columnCount();
    const int     layers  = mesh.layerCount();

    const std::vector<double> &thickness = mesh.thickness();
    const std::vector<double> &area      = mesh.horizontal().faceArea();

    for (int s = 0; s < static_cast<int>(m_velocity.size()); ++s)
    {
        const double w = m_velocity[static_cast<size_t>(s)];

        if (!(w > 0.0))
            continue;

        const int    receiver = m_receiver[static_cast<size_t>(s)];
        const double wdt      = w * dt;

        double delivered = 0.0;

        FVQual_PARALLEL_FOR_REDUCE(+ : delivered)
        for (int64_t c = 0; c < columns; ++c)
        {
            if (mesh.hasHypsography())
            {
                // A settling particle crosses the lower interface or strikes
                // the projected banks. Store both in the column's aggregate
                // bed inventory, whose declared areal basis is the plan polygon.
                double incoming = 0.0, bedAmount = 0.0;
                for (int k = 0; k < layers; ++k)
                {
                    const size_t cell = static_cast<size_t>(mesh.cell(c, k));
                    const double volume = mesh.cellVolume()[cell];
                    if (!(volume > 0.0)) continue;
                    const double top = mesh.interfaceArea()[mesh.interfaceSlot(c, k)];
                    const double bottom = mesh.interfaceArea()[mesh.interfaceSlot(c, k+1)];
                    const double updated = (volume * state.speciesAt(s, cell) + incoming)
                                         / (volume + wdt * top);
                    state.speciesAt(s, cell) = updated;
                    incoming = wdt * bottom * updated;
                    bedAmount += wdt * (top-bottom) * updated;
                }
                bedAmount += incoming;
                delivered += bedAmount;
                if (receiver >= 0)
                    state.speciesAt(receiver, mesh.cell(c, layers-1)) += bedAmount / area[static_cast<size_t>(c)];
                continue;
            }

            // Forward substitution from the surface down. Layer 0 is the
            // surface, and nothing settles in from above it, so the first row
            // has no coupling and the sweep needs no starting value.
            double fromAbove = 0.0;

            for (int k = 0; k < layers; ++k)
            {
                const int64_t cell = mesh.cell(c, k);
                const double  h    = thickness[static_cast<size_t>(cell)];

                if (!(h > 0.0))
                    continue;

                const double courant = wdt / h;

                // (1 + wdt/h) C_new = C_old + (wdt/h) C_above_new
                const double updated =
                    (state.speciesAt(s, cell) + courant * fromAbove) / (1.0 + courant);

                state.speciesAt(s, cell) = updated;

                // What this cell passes down, which the next row consumes.
                fromAbove = updated;
            }

            // Out of the bottom cell and onto the sediment. Areal, so no
            // thickness appears: g/m³ x m/s x s = g/m².
            const double perArea = wdt * fromAbove;

            delivered += perArea * area[static_cast<size_t>(c)];

            if (receiver >= 0)
            {
                const int64_t bed = mesh.cell(c, layers - 1);

                state.speciesAt(receiver, bed) += perArea;
            }
        }

        m_delivered[static_cast<size_t>(s)] += delivered;

        if (receiver < 0)
            m_lost[static_cast<size_t>(s)] += delivered;
    }
}

double SettlingSet::deliveredToBed(int species) const
{
    if (species < 0 || species >= static_cast<int>(m_delivered.size()))
        return 0.0;

    return m_delivered[static_cast<size_t>(species)];
}

double SettlingSet::lostToBed(int species) const
{
    if (species < 0 || species >= static_cast<int>(m_lost.size()))
        return 0.0;

    return m_lost[static_cast<size_t>(species)];
}

void SettlingSet::resetTotals()
{
    std::fill(m_delivered.begin(), m_delivered.end(), 0.0);
    std::fill(m_lost.begin(), m_lost.end(), 0.0);
}

} // namespace FVQual::Transport
