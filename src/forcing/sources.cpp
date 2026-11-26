// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/forcing/sources.h"

#include "fvqual/hydro/equationofstate.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace FVQual::Forcing {

namespace {

//! Guard against dividing by a vanishing density contrast.
constexpr double kTiny = 1.0e-10;

/*!
 * \brief Elevation of a layer's centre.
 */
double layerCentre(const Mesh::LayeredMesh &mesh, int64_t column, int layer)
{
    return 0.5 * (mesh.z(column, layer) + mesh.z(column, layer + 1));
}

/*!
 * \brief Bulk buoyancy frequency N (1/s) of a column.
 *
 * N sets the thickness of a withdrawal envelope: strong stratification
 * confines the approaching flow to a thin band, weak stratification lets it
 * draw from most of the column.
 *
 * Deliberately the *column-scale* value rather than the gradient local to the
 * outlet. The local form has a pathology that shows up immediately in a
 * two-layer reservoir: a port set well below a sharp thermocline sits in
 * uniform cold water where the local gradient is zero, so the envelope would
 * be reported as unbounded and the outlet would appear to draw from the
 * surface — the exact opposite of what a deep port does. It is the density
 * contrast the outlet must work against that confines the flow, and that is a
 * property of the column.
 *
 * CE-QUAL-W2 refines this by iterating the zone limits against the gradient
 * inside them; that refinement is worth having later and is not needed for the
 * envelope to behave correctly.
 */
double buoyancyFrequency(const Mesh::LayeredMesh &mesh, std::span<const double> density,
                         int64_t column, double gravity, double referenceDensity)
{
    const int layers = mesh.layerCount();

    if (layers < 2)
        return 0.0;

    const double surfaceDensity = density[static_cast<size_t>(mesh.cell(column, 0))];
    const double bedDensity =
        density[static_cast<size_t>(mesh.cell(column, layers - 1))];

    const double depth = mesh.z(column, 0) - mesh.z(column, layers);

    if (depth <= 0.0)
        return 0.0;

    const double squared =
        gravity / referenceDensity * (bedDensity - surfaceDensity) / depth;

    return squared > 0.0 ? std::sqrt(squared) : 0.0;
}

} // namespace

int SourceSet::add(const Source &source, std::string &message)
{
    if (source.column < 0)
    {
        message = "source column index must be non-negative";
        return -1;
    }

    if (source.placement == Placement::SelectiveWithdrawal && source.discharge > 0.0)
    {
        message = "selective withdrawal describes an outflow; "
                  "use DensityMatched for an inflow that finds its own level";
        return -1;
    }

    if (source.placement == Placement::DensityMatched && source.discharge < 0.0)
    {
        message = "density-matched placement describes an inflow; "
                  "use SelectiveWithdrawal for an outflow";
        return -1;
    }

    if (!std::isfinite(source.discharge))
    {
        message = "source discharge must be finite";
        return -1;
    }

    m_sources.push_back(source);
    message.clear();

    return static_cast<int>(m_sources.size()) - 1;
}

bool SourceSet::setDischarge(int index, double discharge, std::string &message)
{
    if (index < 0 || index >= count())
    {
        message = "source index is out of range";
        return false;
    }

    if (!std::isfinite(discharge))
    {
        message = "source discharge must be finite";
        return false;
    }

    Source &source = m_sources[static_cast<size_t>(index)];

    // A schedule may drive a discharge to zero, but it must not flip a
    // withdrawal into an inflow behind a placement rule that cannot represent
    // one.
    if (source.placement == Placement::SelectiveWithdrawal && discharge > 0.0)
    {
        message = "a selective withdrawal cannot be given a positive discharge";
        return false;
    }

    if (source.placement == Placement::DensityMatched && discharge < 0.0)
    {
        message = "a density-matched inflow cannot be given a negative discharge";
        return false;
    }

    source.discharge = discharge;
    message.clear();

    return true;
}

bool SourceSet::setPlacement(int index, Placement placement, std::string &message)
{
    if (index < 0 || index >= count())
    {
        message = "source index " + std::to_string(index) + " is out of range";
        return false;
    }

    m_sources[static_cast<size_t>(index)].placement = placement;

    message.clear();
    return true;
}

void SourceSet::clear() noexcept
{
    m_sources.clear();
    m_envelope.clear();
    m_placed.clear();
    m_sourceWeights.clear();
}

void SourceSet::inflowConcentration(const Mesh::LayeredMesh &mesh, int species,
                                    int temperatureIndex,
                                    std::vector<double> &out) const
{
    out.assign(static_cast<size_t>(mesh.cellCount()), 0.0);

    if (m_sourceWeights.empty())
        return;

    const int layers = mesh.layerCount();

    // Discharge-weighted, so two inflows meeting in one cell deliver the
    // mixture their flows imply rather than whichever was registered last.
    std::vector<double> weight(static_cast<size_t>(mesh.cellCount()), 0.0);

    for (size_t s = 0; s < m_sources.size(); ++s)
    {
        const Source &source = m_sources[s];

        if (source.discharge <= 0.0 || source.column >= mesh.columnCount())
            continue;

        double incoming = 0.0;

        if (species == temperatureIndex)
        {
            incoming = source.temperature;
        }
        else if (static_cast<size_t>(species) < source.concentration.size())
        {
            incoming = source.concentration[static_cast<size_t>(species)];
        }

        for (int k = 0; k < layers; ++k)
        {
            const double share =
                m_sourceWeights[s * static_cast<size_t>(layers) + static_cast<size_t>(k)];

            if (share <= 0.0)
                continue;

            const size_t cell = static_cast<size_t>(mesh.cell(source.column, k));
            const double flow = source.discharge * share;

            out[cell] += flow * incoming;
            weight[cell] += flow;
        }
    }

    for (size_t i = 0; i < out.size(); ++i)
    {
        if (weight[i] > 0.0)
            out[i] /= weight[i];
    }
}

void SourceSet::scalarLinearization(const Mesh::LayeredMesh &mesh,
    std::vector<double> &withdrawal, std::vector<ScalarSourceEntry> &incoming) const
{
    withdrawal.assign(static_cast<size_t>(mesh.cellCount()),0.0); incoming.clear();
    const size_t layers=static_cast<size_t>(mesh.layerCount());
    for (size_t s=0;s<m_sources.size();++s) {
        const auto &source=m_sources[s];
        if (source.column<0 || source.column>=mesh.columnCount()) continue;
        for (size_t k=0;k<layers;++k) {
            const size_t cell=mesh.cell(source.column,static_cast<int>(k));
            const double rate=source.discharge*m_sourceWeights.at(s*layers+k);
            if (rate>0) incoming.push_back({cell,s,rate});
            else withdrawal[cell]+=rate;
        }
    }
}

void SourceSet::scalarLoad(const Mesh::LayeredMesh &mesh, int species,
                          int temperatureIndex, std::span<const double> ambient,
                          std::vector<double> &out) const
{
    out.assign(static_cast<size_t>(mesh.cellCount()), 0.0);
    const int layers = mesh.layerCount();
    for (size_t s = 0; s < m_sources.size(); ++s) {
        const auto &source = m_sources[s];
        if (source.column < 0 || source.column >= mesh.columnCount()) continue;
        const double incoming = species == temperatureIndex ? source.temperature :
            static_cast<size_t>(species) < source.concentration.size() ? source.concentration[species] : 0;
        for (int k = 0; k < layers; ++k) {
            const size_t cell = static_cast<size_t>(mesh.cell(source.column, k));
            out[cell] += source.discharge * m_sourceWeights[s * layers + k] *
                (source.discharge > 0 ? incoming : ambient[cell]);
        }
    }
}

double SourceSet::netDischarge() const noexcept
{
    double total = 0.0;

    for (const Source &source : m_sources)
        total += source.discharge;

    return total;
}

double SourceSet::withdrawalTimeStep(const Mesh::LayeredMesh &mesh,
                                    std::span<const double> density, double gravity,
                                    double referenceDensity, double courant) const
{
    std::vector<double> withdrawal(static_cast<size_t>(mesh.cellCount()), 0.0), weights;
    for (const auto &source : m_sources) {
        if (source.discharge >= 0 || source.column < 0 || source.column >= mesh.columnCount()) continue;
        double envelope = 0, placed = 0;
        weightsCore(source, mesh, density, gravity, referenceDensity, weights, envelope, placed);
        for (int k = 0; k < mesh.layerCount(); ++k)
            withdrawal[static_cast<size_t>(mesh.cell(source.column, k))] -= source.discharge * weights[k];
    }
    double limit = std::numeric_limits<double>::infinity();
    for (size_t cell = 0; cell < withdrawal.size(); ++cell)
        if (withdrawal[cell] > 0)
            limit = std::min(limit, courant * mesh.cellVolume()[cell] / withdrawal[cell]);
    return limit;
}

void SourceSet::weightsFor(const Source &source, const Mesh::LayeredMesh &mesh,
                           std::span<const double> density, double gravity,
                           double referenceDensity, int index,
                           std::vector<double> &weights)
{
    weightsCore(source, mesh, density, gravity, referenceDensity, weights,
                m_envelope[static_cast<size_t>(index)],
                m_placed[static_cast<size_t>(index)]);
}

void SourceSet::weightsCore(const Source &source, const Mesh::LayeredMesh &mesh,
                            std::span<const double> density, double gravity,
                            double referenceDensity, std::vector<double> &weights,
                            double &envelope, double &placed) const
{
    const int     layers = mesh.layerCount();
    const int64_t column = source.column;

    weights.assign(static_cast<size_t>(layers), 0.0);

    envelope = 0.0;
    placed   = layerCentre(mesh, column, 0);

    switch (source.placement)
    {
    case Placement::Uniform:
    {
        // Proportional to thickness, so a uniform source does not concentrate
        // itself in whichever layers happen to be thin.
        double total = 0.0;

        for (int k = 0; k < layers; ++k)
            total += mesh.thickness()[static_cast<size_t>(mesh.cell(column, k))];

        if (total <= 0.0)
            return;

        for (int k = 0; k < layers; ++k)
        {
            weights[static_cast<size_t>(k)] =
                mesh.thickness()[static_cast<size_t>(mesh.cell(column, k))] / total;
        }

        placed =
            0.5 * (mesh.z(column, 0) + mesh.z(column, layers));
        return;
    }

    case Placement::Surface:
        weights[0] = 1.0;
        placed = layerCentre(mesh, column, 0);
        return;

    case Placement::Bed:
        weights[static_cast<size_t>(layers) - 1] = 1.0;
        placed = layerCentre(mesh, column, layers - 1);
        return;

    case Placement::Elevation:
    {
        int    nearest = 0;
        double best    = 1.0e30;

        for (int k = 0; k < layers; ++k)
        {
            const double gap = std::abs(layerCentre(mesh, column, k) - source.elevation);

            if (gap < best)
            {
                best    = gap;
                nearest = k;
            }
        }

        weights[static_cast<size_t>(nearest)] = 1.0;
        placed = layerCentre(mesh, column, nearest);
        return;
    }

    case Placement::DensityMatched:
    {
        // The inflow descends until the water around it is at least as dense
        // as it is, and settles there. A river colder than the whole reservoir
        // reaches the bed as an underflow; one warmer than the surface stays
        // on top as an overflow. Both fall out of the same scan.
        const double inflowDensity =
            Hydro::EquationOfState::density(source.temperature);

        int settled = layers - 1;

        for (int k = 0; k < layers; ++k)
        {
            if (density[static_cast<size_t>(mesh.cell(column, k))] >= inflowDensity)
            {
                settled = k;
                break;
            }
        }

        weights[static_cast<size_t>(settled)] = 1.0;
        placed = layerCentre(mesh, column, settled);
        return;
    }

    case Placement::SelectiveWithdrawal:
    {
        // Selective withdrawal after Bohan and Grace (1973) and Smith et al.
        // (1987) — the theory CE-QUAL-W2 implements in DOWNSTREAM_WITHDRAWAL,
        // reproduced here at the formulation level (see THIRD_PARTY_LICENSES).
        //
        // The withdrawal zone is found by marching outward from the outlet's
        // own layer. For each candidate layer the buoyancy frequency is taken
        // between *that layer and the outlet layer*, over their separation H,
        //
        //     N(H) = sqrt( g |rho_k - rho_str| / (H rho_str) ) ,
        //
        // and the zone thickness it implies is
        //
        //     point sink:  d = ( c_bi Q / N )^(1/3)
        //     line sink:   d = ( 2 c_bi Q / (W N) )^(1/2) ,
        //
        // the march stopping at the first layer whose separation has caught up
        // with the thickness. That N is *implicit* — it depends on the very
        // separation being tested — and the implicitness is the whole point.
        // Evaluating N once at column scale instead, then keeping these
        // constants, inflates the zone without bound: in a 24 m column it
        // returns a 36 m half-height, which is not a subtle error but a
        // withdrawal that draws from water that is not there. The pairing of a
        // locally evaluated N with these coefficients is not decomposable.
        //
        // c_bi is the boundary-interference coefficient: 2 when the outlet sits
        // in the top or bottom tenth of the water column, where the zone can
        // only develop on one side, and 1 where it develops freely.
        const double magnitude = std::abs(source.discharge);

        if (magnitude <= 0.0)
            return;

        const double surface = mesh.z(column, 0);
        const double bed     = mesh.z(column, layers);
        const double depth   = surface - bed;

        if (depth <= 0.0)
            return;

        // The layer the outlet discharges from.
        int outlet = 0;
        {
            double best = 1.0e30;

            for (int k = 0; k < layers; ++k)
            {
                const double gap =
                    std::abs(layerCentre(mesh, column, k) - source.elevation);

                if (gap < best)
                {
                    best   = gap;
                    outlet = k;
                }
            }
        }

        const size_t outletCell = static_cast<size_t>(mesh.cell(column, outlet));
        const double outletDensity = density[outletCell];

        const double position = (source.elevation - bed) / depth;
        const double interference = (position < 0.1 || position > 0.9) ? 2.0 : 1.0;

        const auto thicknessFor = [&](double separation, double contrast) {
            // N over this separation; a vanishing contrast means nothing
            // confines the flow, so the implied zone is unbounded.
            const double frequency = std::sqrt(gravity * contrast /
                                               (separation * outletDensity + kTiny));

            if (frequency <= kTiny)
                return std::numeric_limits<double>::infinity();

            return source.outletWidth > 0.0
                       ? std::sqrt(2.0 * interference * magnitude /
                                   (source.outletWidth * frequency))
                       : std::cbrt(interference * magnitude / frequency);
        };

        // Upper limit of the zone.
        int top = 0;

        for (int k = outlet - 1; k >= 0; --k)
        {
            const double separation = layerCentre(mesh, column, k) - source.elevation;

            if (separation <= 0.0)
                continue;

            const double contrast =
                std::abs(density[static_cast<size_t>(mesh.cell(column, k))] -
                         outletDensity);

            if (separation >= thicknessFor(separation, contrast))
            {
                top = k;
                break;
            }
        }

        // Lower limit.
        int bottom = layers - 1;

        for (int k = outlet + 1; k < layers; ++k)
        {
            const double separation = source.elevation - layerCentre(mesh, column, k);

            if (separation <= 0.0)
                continue;

            const double contrast =
                std::abs(density[static_cast<size_t>(mesh.cell(column, k))] -
                         outletDensity);

            if (separation >= thicknessFor(separation, contrast))
            {
                bottom = k;
                break;
            }
        }

        envelope = 0.5 * (layerCentre(mesh, column, top) -
                          layerCentre(mesh, column, bottom));
        placed   = source.elevation;

        // Where the fastest approach velocity sits. Bohan and Grace (1973)
        // report that the maximum lies at the port only when the zone develops
        // symmetrically about it; otherwise it shifts, by
        //
        //     y1 = H sin^2( 1.57 z1 / H )                          (eq. 4)
        //
        // with z1 the distance from the port down to the lower limit, H the
        // total zone thickness, and y1 the height of the maximum above that
        // lower limit. Symmetric zones give z1/H = 1/2 and sin^2(pi/4) = 1/2,
        // returning the port itself, so this generalises the earlier behaviour
        // rather than replacing it.
        //
        // Only applied where the zone develops freely. Once it runs into the
        // surface or the bed the relationship no longer holds, and the
        // reference implementation likewise skips it.
        const double above = layerCentre(mesh, column, top) - source.elevation;
        const double below = source.elevation - layerCentre(mesh, column, bottom);

        int    fastest         = outlet;
        double fastestDensity  = outletDensity;

        const bool developsFreely = top > 0 && bottom < layers - 1;

        if (developsFreely && above + below > 0.0)
        {
            const double thickness = above + below;
            const double height =
                thickness * std::pow(std::sin(1.57 * below / thickness), 2.0);

            const double elevation = source.elevation - below + height;

            double best = 1.0e30;

            for (int k = top; k <= bottom; ++k)
            {
                const double gap = std::abs(layerCentre(mesh, column, k) - elevation);

                if (gap < best)
                {
                    best    = gap;
                    fastest = k;
                }
            }

            fastestDensity = density[static_cast<size_t>(mesh.cell(column, fastest))];
        }

        // Velocity profile. Parabolic in the *density* anomaly about the
        // fastest layer rather than in distance from it, normalised separately
        // above and below by the contrast at each limit. In a linearly
        // stratified column the distance and density readings coincide; across
        // a sharp thermocline they do not, and it is the density form that
        // keeps the draw inside the layer the outlet actually sits in.
        const double aboveContrast = std::max(
            std::abs(density[static_cast<size_t>(mesh.cell(column, top))] -
                     fastestDensity),
            kTiny);
        const double belowContrast = std::max(
            std::abs(density[static_cast<size_t>(mesh.cell(column, bottom))] -
                     fastestDensity),
            kTiny);

        double total = 0.0;

        for (int k = top; k <= bottom; ++k)
        {
            const size_t cell = static_cast<size_t>(mesh.cell(column, k));

            const double contrast = k > fastest ? belowContrast : aboveContrast;
            const double ratio    = (density[cell] - fastestDensity) / contrast;

            const double velocity = std::clamp(1.0 - ratio * ratio, 0.0, 1.0);

            if (velocity <= 0.0)
                continue;

            // Flow through a layer is its approach velocity times its area.
            const double weight = velocity * mesh.thickness()[cell];

            weights[static_cast<size_t>(k)] = weight;
            total += weight;
        }

        if (total <= 0.0)
        {
            weights[static_cast<size_t>(outlet)] = 1.0;
            return;
        }

        for (int k = 0; k < layers; ++k)
            weights[static_cast<size_t>(k)] /= total;

        return;
    }
    }
}

double SourceSet::weightedMean(const Mesh::LayeredMesh &mesh, int index,
                               std::span<const double> values) const
{
    const int    layers = mesh.layerCount();
    const size_t offset = static_cast<size_t>(index) * static_cast<size_t>(layers);

    if (offset + static_cast<size_t>(layers) > m_sourceWeights.size())
        return 0.0;

    const int64_t column = m_sources[static_cast<size_t>(index)].column;

    if (column >= mesh.columnCount())
        return 0.0;

    double total    = 0.0;
    double weighted = 0.0;

    for (int k = 0; k < layers; ++k)
    {
        const double share = m_sourceWeights[offset + static_cast<size_t>(k)];

        if (share == 0.0)
            continue;

        weighted += share * values[static_cast<size_t>(mesh.cell(column, k))];
        total    += share;
    }

    return total > 0.0 ? weighted / total : 0.0;
}

bool SourceSet::setInflowValue(int index, int species, int temperatureIndex,
                               double value, std::string &message)
{
    if (index < 0 || index >= count())
    {
        message = "source index is out of range";
        return false;
    }

    Source &source = m_sources[static_cast<size_t>(index)];

    if (species == temperatureIndex)
    {
        source.temperature = value;
    }
    else
    {
        if (static_cast<size_t>(species) >= source.concentration.size())
            source.concentration.resize(static_cast<size_t>(species) + 1, 0.0);

        source.concentration[static_cast<size_t>(species)] = value;
    }

    message.clear();
    return true;
}

double SourceSet::withdrawalMean(const Mesh::LayeredMesh &mesh, int index,
                                 double probeDischarge,
                                 std::span<const double> density, double gravity,
                                 double referenceDensity,
                                 std::span<const double> values) const
{
    const Source &stored = m_sources[static_cast<size_t>(index)];

    if (stored.column >= mesh.columnCount())
        return 0.0;

    // Probed at a stated discharge rather than the source's own, because the
    // question a blending operator asks is what a port *would* deliver if it
    // were opened, and a closed port has no envelope at all: the withdrawal
    // half-height scales as Q^(1/3), so at zero discharge it collapses and the
    // answer would be the temperature of whichever single layer the outlet
    // happens to sit in.
    Source probe     = stored;
    probe.discharge  = -std::abs(probeDischarge);

    std::vector<double> weights;
    double              envelope = 0.0;
    double              placed   = 0.0;

    weightsCore(probe, mesh, density, gravity, referenceDensity, weights, envelope,
                placed);

    const int layers = mesh.layerCount();

    double total     = 0.0;
    double weighted  = 0.0;

    for (int k = 0; k < layers; ++k)
    {
        const double share = weights[static_cast<size_t>(k)];

        if (share == 0.0)
            continue;

        weighted += share * values[static_cast<size_t>(mesh.cell(probe.column, k))];
        total    += share;
    }

    return total > 0.0 ? weighted / total : 0.0;
}

void SourceSet::columnTotals(std::span<double> columnDischarge) const
{
    std::fill(columnDischarge.begin(), columnDischarge.end(), 0.0);

    for (const Source &source : m_sources)
    {
        if (source.column < static_cast<int64_t>(columnDischarge.size()))
            // Boundary water is already carried by the external mode's own
            // boundary flux; reporting it here would move it twice.
            if (!source.externallyAccounted)
                columnDischarge[static_cast<size_t>(source.column)] += source.discharge;
    }
}

void SourceSet::distribute(const Mesh::LayeredMesh &mesh, std::span<const double> density,
                           double gravity, double referenceDensity,
                           std::vector<double> &cellDischarge,
                           std::vector<double> &columnDischarge)
{
    cellDischarge.assign(static_cast<size_t>(mesh.cellCount()), 0.0);
    columnDischarge.assign(static_cast<size_t>(mesh.columnCount()), 0.0);

    m_envelope.assign(m_sources.size(), 0.0);
    m_placed.assign(m_sources.size(), 0.0);

    const int layers = mesh.layerCount();

    m_sourceWeights.assign(m_sources.size() * static_cast<size_t>(layers), 0.0);

    for (size_t s = 0; s < m_sources.size(); ++s)
    {
        const Source &source = m_sources[s];

        if (source.column >= mesh.columnCount())
            continue;

        if (source.discharge == 0.0)
        {
            // Still record where it would have acted, so a diagnostic reading
            // of a shut-off outlet is not misleading.
            weightsFor(source, mesh, density, gravity, referenceDensity,
                       static_cast<int>(s), m_weights);
            continue;
        }

        weightsFor(source, mesh, density, gravity, referenceDensity,
                   static_cast<int>(s), m_weights);

        for (int k = 0; k < layers; ++k)
        {
            m_sourceWeights[s * static_cast<size_t>(layers) + static_cast<size_t>(k)] =
                m_weights[static_cast<size_t>(k)];
        }

        for (int k = 0; k < layers; ++k)
        {
            const double share = m_weights[static_cast<size_t>(k)];

            if (share == 0.0)
                continue;

            cellDischarge[static_cast<size_t>(mesh.cell(source.column, k))] +=
                source.discharge * share;
        }

        // Boundary water is already carried by the external mode's own
        // boundary flux; reporting it here would move it twice.
        if (!source.externallyAccounted)
            columnDischarge[static_cast<size_t>(source.column)] += source.discharge;
    }
}

} // namespace FVQual::Forcing
