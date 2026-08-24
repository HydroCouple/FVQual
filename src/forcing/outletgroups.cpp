// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/forcing/outletgroups.h"

#include <algorithm>
#include <cmath>

namespace FVQual::Forcing {

int OutletGroupSet::add(const OutletGroup &group, const StructureSet &structures,
                        std::string &message)
{
    if (group.ports.empty())
    {
        message = "an outlet group needs at least one port";
        return -1;
    }

    for (int port : group.ports)
    {
        if (port < 0 || port >= structures.count())
        {
            message = "outlet group '" + group.id + "' references a structure that "
                                                    "does not exist";
            return -1;
        }
    }

    // Two ports at the same elevation cannot blend to anything their shared
    // depth does not already hold, but one port listed twice would be given
    // two independent shares of the demand and release double what was asked.
    std::vector<int> sorted = group.ports;
    std::sort(sorted.begin(), sorted.end());

    if (std::adjacent_find(sorted.begin(), sorted.end()) != sorted.end())
    {
        message = "outlet group '" + group.id + "' lists the same port twice";
        return -1;
    }

    if (group.demand < 0.0)
    {
        message = "an outlet group's demand is a magnitude and cannot be negative";
        return -1;
    }

    m_groups.push_back(group);
    m_delivered.push_back(0.0);
    m_achieved.push_back(0.0);

    message.clear();
    return static_cast<int>(m_groups.size()) - 1;
}

bool OutletGroupSet::setDemand(int index, double demand, std::string &message)
{
    if (index < 0 || index >= count())
    {
        message = "outlet group index is out of range";
        return false;
    }

    if (demand < 0.0)
    {
        message = "an outlet group's demand is a magnitude and cannot be negative";
        return false;
    }

    m_groups[static_cast<size_t>(index)].demand = demand;

    message.clear();
    return true;
}

bool OutletGroupSet::setTarget(int index, double temperature, std::string &message)
{
    if (index < 0 || index >= count())
    {
        message = "outlet group index is out of range";
        return false;
    }

    m_groups[static_cast<size_t>(index)].targetTemperature = temperature;

    message.clear();
    return true;
}

bool OutletGroupSet::update(const Mesh::LayeredMesh &mesh, const SourceSet &sources,
                            std::span<const double> density,
                            std::span<const double> temperature, double gravity,
                            double referenceDensity, StructureSet &structures,
                            std::string &message)
{
    for (size_t g = 0; g < m_groups.size(); ++g)
    {
        const OutletGroup &group = m_groups[g];
        const size_t       ports = group.ports.size();

        m_portTemperature.assign(ports, 0.0);
        m_portCapacity.assign(ports, 0.0);
        m_portRelease.assign(ports, 0.0);

        // What each port would deliver, and the most it could pass. Both are
        // probed at the full demand: the withdrawal envelope depends on the
        // discharge, so the temperature a port delivers is only meaningful
        // against a stated flow, and the flow worth asking about is the one
        // that might be sent through it.
        for (size_t p = 0; p < ports; ++p)
        {
            const int        index     = group.ports[p];
            const Structure &structure = structures.structure(index);

            m_portCapacity[p] =
                StructureSet::rate(structure, mesh.z(structure.column, 0));

            m_portTemperature[p] = sources.withdrawalMean(
                mesh, structures.sourceIndex(index), group.demand, density, gravity,
                referenceDensity, temperature);
        }

        if (group.demand <= 0.0)
        {
            for (size_t p = 0; p < ports; ++p)
            {
                if (!structures.setSchedule(group.ports[p], 0.0, message))
                    return false;
            }

            m_delivered[g] = 0.0;
            m_achieved[g]  = 0.0;
            continue;
        }

        if (!group.blend)
        {
            // No temperature objective: share the demand out in proportion to
            // what each port can pass, so the group drains evenly rather than
            // loading one outlet to capacity while another sits idle.
            double capacity = 0.0;

            for (size_t p = 0; p < ports; ++p)
                capacity += m_portCapacity[p];

            for (size_t p = 0; p < ports; ++p)
            {
                m_portRelease[p] = capacity > 0.0
                                       ? group.demand * m_portCapacity[p] / capacity
                                       : 0.0;
            }
        }
        else
        {
            // The tightest pair of ports bracketing the target. Tightest rather
            // than any bracketing pair because a blend of two ports close to
            // the target is far less sensitive to error in either one's
            // delivered temperature than a blend of the coldest and warmest.
            const double target = group.targetTemperature;

            size_t warm = ports;
            size_t cold = ports;

            for (size_t p = 0; p < ports; ++p)
            {
                if (m_portCapacity[p] <= 0.0)
                    continue;

                if (m_portTemperature[p] >= target &&
                    (warm == ports || m_portTemperature[p] < m_portTemperature[warm]))
                {
                    warm = p;
                }

                if (m_portTemperature[p] <= target &&
                    (cold == ports || m_portTemperature[p] > m_portTemperature[cold]))
                {
                    cold = p;
                }
            }

            if (warm < ports && cold < ports && warm != cold)
            {
                const double spread = m_portTemperature[warm] - m_portTemperature[cold];

                // Heat conservation in the blended release.
                const double fraction =
                    spread > 0.0 ? (target - m_portTemperature[cold]) / spread : 0.5;

                const double clamped = std::clamp(fraction, 0.0, 1.0);

                m_portRelease[warm] = clamped * group.demand;
                m_portRelease[cold] = (1.0 - clamped) * group.demand;
            }
            else
            {
                // No pair brackets the target, so the whole reservoir is either
                // warmer or colder than what was asked for. The best available
                // release is the single port closest to it — and the achieved
                // temperature will miss, which is the honest answer.
                size_t closest = ports;
                double best    = 0.0;

                for (size_t p = 0; p < ports; ++p)
                {
                    if (m_portCapacity[p] <= 0.0)
                        continue;

                    const double gap = std::abs(m_portTemperature[p] - target);

                    if (closest == ports || gap < best)
                    {
                        best    = gap;
                        closest = p;
                    }
                }

                if (closest < ports)
                    m_portRelease[closest] = group.demand;
            }
        }

        // Capacity limiting. A port asked for more than its head can pass
        // delivers its capacity; the difference is a shortfall and is reported
        // rather than quietly redistributed, because an operator who cannot
        // meet a demand needs to know that, not to have the model paper over it.
        double delivered = 0.0;
        double weighted  = 0.0;

        for (size_t p = 0; p < ports; ++p)
        {
            m_portRelease[p] = std::min(m_portRelease[p], m_portCapacity[p]);

            if (!structures.setSchedule(group.ports[p], m_portRelease[p], message))
                return false;

            delivered += m_portRelease[p];
            weighted  += m_portRelease[p] * m_portTemperature[p];
        }

        m_delivered[g] = delivered;
        m_achieved[g]  = delivered > 0.0 ? weighted / delivered : 0.0;
    }

    message.clear();
    return true;
}

} // namespace FVQual::Forcing
