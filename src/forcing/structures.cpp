// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/forcing/structures.h"

#include <algorithm>
#include <cmath>

namespace FVQual::Forcing {

namespace {

constexpr double kGravity = 9.80665;

/*!
 * \brief Conventional discharge coefficients, used when none is given.
 */
double defaultCoefficient(StructureType type)
{
    switch (type)
    {
    case StructureType::Weir:    return 1.84; // sharp-crested, SI form
    case StructureType::Orifice: return 0.62;
    case StructureType::Pipe:    return 0.80; // absorbs entrance and friction losses
    case StructureType::Pump:    return 1.0;
    }

    return 1.0;
}

} // namespace

double StructureSet::rate(const Structure &structure, double upstream)
{
    const double coefficient = structure.coefficient > 0.0
                                   ? structure.coefficient
                                   : defaultCoefficient(structure.type);

    const double head = upstream - structure.crestElevation;

    if (structure.type == StructureType::Pump)
        return std::abs(structure.capacity);

    // Nothing to drive the flow: the water is at or below the crest.
    if (head <= 0.0)
        return 0.0;

    const double tail = structure.tailwaterElevation;

    switch (structure.type)
    {
    case StructureType::Weir:
    {
        if (structure.width <= 0.0)
            return 0.0;

        const double free = coefficient * structure.width * std::pow(head, 1.5);

        // Villemonte's submergence correction. A weir drowned by its own
        // tailwater passes less than the free-flow rating, and by the time the
        // downstream level reaches the upstream one it passes nothing.
        const double downstreamHead = tail - structure.crestElevation;

        if (downstreamHead <= 0.0)
            return free;

        if (downstreamHead >= head)
            return 0.0;

        const double ratio = downstreamHead / head;

        return free * std::pow(1.0 - std::pow(ratio, 1.5), 0.385);
    }

    case StructureType::Orifice:
    {
        if (structure.area <= 0.0)
            return 0.0;

        // Free discharge is driven by the head on the orifice centre;
        // once the tailwater submerges it, only the difference across it
        // drives the flow.
        double driving = head;

        if (tail > structure.crestElevation)
            driving = upstream - tail;

        if (driving <= 0.0)
            return 0.0;

        return coefficient * structure.area * std::clamp(structure.opening, 0.0, 1.0) *
               std::sqrt(2.0 * kGravity * driving);
    }

    case StructureType::Pipe:
    {
        if (structure.area <= 0.0)
            return 0.0;

        // A pipe runs full, so the head difference always drives it — there is
        // no free-discharge branch to fall back on.
        const double downstream =
            tail > structure.crestElevation ? tail : structure.crestElevation;

        const double driving = upstream - downstream;

        if (driving <= 0.0)
            return 0.0;

        return coefficient * structure.area * std::clamp(structure.opening, 0.0, 1.0) *
               std::sqrt(2.0 * kGravity * driving);
    }

    case StructureType::Pump:
        break;
    }

    return 0.0;
}

int StructureSet::add(const Structure &structure, SourceSet &sources,
                      std::string &message)
{
    if (structure.column < 0)
    {
        message = "structure column index must be non-negative";
        return -1;
    }

    if (structure.type == StructureType::Weir && structure.width <= 0.0)
    {
        message = "a weir needs a positive crest length";
        return -1;
    }

    if ((structure.type == StructureType::Orifice ||
         structure.type == StructureType::Pipe) &&
        structure.area <= 0.0)
    {
        message = "an orifice or pipe needs a positive flow area";
        return -1;
    }

    if (structure.type == StructureType::Pump && structure.capacity == 0.0)
    {
        message = "a pump needs a non-zero capacity";
        return -1;
    }

    if (structure.controlled && structure.stopElevation > structure.startElevation)
    {
        message = "the stop elevation must not exceed the start elevation, "
                  "or the deadband is inverted and the structure would chatter";
        return -1;
    }

    // The source this structure will drive. Outlet works release, so the
    // placement must be one that describes an outflow.
    Source source;
    source.origin = SourceOrigin::Structure;
    source.column      = structure.column;
    source.discharge   = 0.0;
    source.elevation   = structure.crestElevation;
    source.outletWidth = structure.outletWidth;
    source.temperature = structure.temperature;
    source.concentration = structure.concentration;
    source.id          = structure.id;

    const bool adds = structure.type == StructureType::Pump && structure.capacity > 0.0;

    if (adds)
    {
        // A pumped inflow has a stated composition and no reason to follow an
        // outlet envelope; let it find its own level unless told otherwise.
        source.placement = structure.placement == Placement::SelectiveWithdrawal
                               ? Placement::DensityMatched
                               : structure.placement;
    }
    else
    {
        source.placement = structure.placement;
    }

    const int sourceIndex = sources.add(source, message);

    if (sourceIndex < 0)
        return -1;

    m_structures.push_back(structure);
    m_sourceIndex.push_back(sourceIndex);
    m_discharge.push_back(0.0);
    m_running.push_back(structure.controlled ? char{0} : char{1});

    message.clear();
    return static_cast<int>(m_structures.size()) - 1;
}

bool StructureSet::setOpening(int index, double opening, std::string &message)
{
    if (index < 0 || index >= count())
    {
        message = "structure index is out of range";
        return false;
    }

    if (!(opening >= 0.0 && opening <= 1.0))
    {
        message = "gate opening must lie in [0, 1]";
        return false;
    }

    m_structures[static_cast<size_t>(index)].opening = opening;

    message.clear();
    return true;
}

bool StructureSet::setTailwater(int index, double elevation, std::string &message)
{
    if (index < 0 || index >= count())
    {
        message = "structure index is out of range";
        return false;
    }

    m_structures[static_cast<size_t>(index)].tailwaterElevation = elevation;

    message.clear();
    return true;
}

bool StructureSet::setSchedule(int index, double discharge, std::string &message)
{
    if (index < 0 || index >= count())
    {
        message = "structure index is out of range";
        return false;
    }

    if (discharge < 0.0)
    {
        message = "a scheduled release is a magnitude and cannot be negative; "
                  "direction is a property of the structure, not the schedule";
        return false;
    }

    Structure &structure = m_structures[static_cast<size_t>(index)];

    structure.scheduled          = true;
    structure.scheduledDischarge = discharge;

    message.clear();
    return true;
}

bool StructureSet::clearSchedule(int index, std::string &message)
{
    if (index < 0 || index >= count())
    {
        message = "structure index is out of range";
        return false;
    }

    m_structures[static_cast<size_t>(index)].scheduled = false;

    message.clear();
    return true;
}

bool StructureSet::update(std::span<const double> surfaceElevation, SourceSet &sources,
                          std::string &message)
{
    for (size_t s = 0; s < m_structures.size(); ++s)
    {
        const Structure &structure = m_structures[s];

        if (structure.column >= static_cast<int64_t>(surfaceElevation.size()))
        {
            message = "structure '" + structure.id + "' references a column that does "
                                                     "not exist in this mesh";
            return false;
        }

        const double upstream = surfaceElevation[static_cast<size_t>(structure.column)];

        // Level control with a deadband: start at or above the start level,
        // stop at or below the stop level, and hold state in between. Testing
        // the two thresholds separately rather than against one level is what
        // keeps a pump from cycling every step.
        if (structure.controlled)
        {
            bool &running = reinterpret_cast<bool &>(m_running[s]);

            if (!running && upstream >= structure.startElevation)
                running = true;
            else if (running && upstream <= structure.stopElevation)
                running = false;

            if (!running)
            {
                m_discharge[s] = 0.0;

                if (!sources.setDischarge(m_sourceIndex[s], 0.0, message))
                    return false;

                continue;
            }
        }

        double magnitude = rate(structure, upstream);

        // Under a schedule the rating is the capacity, not the release: an
        // operator opens the gate as far as the demand needs and no further,
        // and cannot open it past what the head will pass.
        if (structure.scheduled)
            magnitude = std::min(magnitude, std::abs(structure.scheduledDischarge));

        // Outlet works release; a pump follows the sign of its capacity.
        double signed_ = -magnitude;

        if (structure.type == StructureType::Pump && structure.capacity > 0.0)
            signed_ = magnitude;

        m_discharge[s] = signed_;

        if (!sources.setDischarge(m_sourceIndex[s], signed_, message))
            return false;
    }

    message.clear();
    return true;
}

double StructureSet::totalDischarge() const noexcept
{
    double total = 0.0;

    for (double value : m_discharge)
        total += value;

    return total;
}

} // namespace FVQual::Forcing
