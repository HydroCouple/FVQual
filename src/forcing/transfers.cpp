// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/forcing/transfers.h"

#include <string>

namespace FVQual::Forcing {

int TransferSet::add(const Transfer &transfer, SourceSet &sources, std::string &message)
{
    if (transfer.fromColumn < 0 || transfer.toColumn < 0)
    {
        message = "transfer column indices must be non-negative";
        return -1;
    }

    if (transfer.fromColumn == transfer.toColumn)
    {
        message = "transfer '" + transfer.id + "' moves water from a column to "
                                               "itself, which does nothing but "
                                               "cycle it through the intake";
        return -1;
    }

    if (transfer.discharge < 0.0)
    {
        message = "a transfer rate is a magnitude; reverse the columns to "
                  "reverse the direction";
        return -1;
    }

    Source intake;
    intake.id          = transfer.id + " intake";
    intake.column      = transfer.fromColumn;
    intake.discharge   = -transfer.discharge;
    intake.placement   = transfer.fromPlacement;
    intake.elevation   = transfer.fromElevation;
    intake.outletWidth = transfer.outletWidth;

    const int intakeIndex = sources.add(intake, message);

    if (intakeIndex < 0)
        return -1;

    Source outfall;
    outfall.id        = transfer.id + " outfall";
    outfall.column    = transfer.toColumn;
    outfall.discharge = transfer.discharge;
    outfall.placement = transfer.toPlacement;
    outfall.elevation = transfer.toElevation;

    const int outfallIndex = sources.add(outfall, message);

    if (outfallIndex < 0)
        return -1;

    m_transfers.push_back(transfer);
    m_intake.push_back(intakeIndex);
    m_outfall.push_back(outfallIndex);
    m_deliveredTemperature.push_back(0.0);

    message.clear();
    return static_cast<int>(m_transfers.size()) - 1;
}

bool TransferSet::setDischarge(int index, double discharge, std::string &message)
{
    if (index < 0 || index >= count())
    {
        message = "transfer index is out of range";
        return false;
    }

    if (discharge < 0.0)
    {
        message = "a transfer rate is a magnitude; reverse the columns to "
                  "reverse the direction";
        return false;
    }

    m_transfers[static_cast<size_t>(index)].discharge = discharge;

    message.clear();
    return true;
}

bool TransferSet::updateDischarges(SourceSet &sources, std::string &message)
{
    for (size_t t = 0; t < m_transfers.size(); ++t)
    {
        const double rate = m_transfers[t].discharge;

        // Equal and opposite by construction. The pair is written together so
        // the two halves cannot drift apart and quietly become a net source.
        if (!sources.setDischarge(m_intake[t], -rate, message))
            return false;

        if (!sources.setDischarge(m_outfall[t], rate, message))
            return false;
    }

    message.clear();
    return true;
}

void TransferSet::carry(const Mesh::LayeredMesh &mesh, int species,
                        int temperatureIndex, std::span<const double> values,
                        SourceSet &sources)
{
    std::string ignored;

    for (size_t t = 0; t < m_transfers.size(); ++t)
    {
        if (m_transfers[t].discharge <= 0.0)
            continue;

        // The mean through the weights the withdrawal is actually removing mass
        // at, so what arrives is exactly what left.
        const double drawn = sources.weightedMean(mesh, m_intake[t], values);

        // The index is this class's own invariant, so the result carries no
        // information a caller could act on.
        (void)sources.setInflowValue(m_outfall[t], species, temperatureIndex, drawn,
                                     ignored);

        if (species == temperatureIndex)
            m_deliveredTemperature[t] = drawn;
    }
}

} // namespace FVQual::Forcing
