// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/transport/euleriantransport.h"

#include "fvqual/parallel.h"

#include <algorithm>
#include <cmath>

namespace FVQual::Transport {

namespace {

/*!
 * \brief Minmod limiter for the vertical slope.
 */
double minmod(double a, double b) noexcept
{
    if (a * b <= 0.0)
        return 0.0;

    return std::abs(a) < std::abs(b) ? a : b;
}

} // namespace

void EulerianTransport::initialize(const Mesh::LayeredMesh &mesh)
{
    const size_t cells = static_cast<size_t>(mesh.cellCount());
    const size_t edgeLayers = static_cast<size_t>(mesh.horizontal().edgeCount()) *
                              static_cast<size_t>(mesh.layerCount());
    const size_t interfaces = static_cast<size_t>(mesh.interfaceSlotCount());

    m_lowHorizontal.assign(edgeLayers, 0.0);
    m_highHorizontal.assign(edgeLayers, 0.0);
    m_antiHorizontal.assign(edgeLayers, 0.0);

    m_lowVertical.assign(interfaces, 0.0);
    m_highVertical.assign(interfaces, 0.0);
    m_antiVertical.assign(interfaces, 0.0);

    m_gradientX.assign(cells, 0.0);
    m_gradientY.assign(cells, 0.0);
    m_gradientZ.assign(cells, 0.0);
    m_limiter.assign(cells, 1.0);

    m_lowValues.assign(cells, 0.0);
    m_risePlus.assign(cells, 1.0);
    m_riseMinus.assign(cells, 1.0);
}

void EulerianTransport::computeLowOrderFluxes(const Mesh::LayeredMesh &mesh,
                                              std::span<const double> edgeLayerFlux,
                                              std::span<const double> verticalFlux,
                                              std::span<const double> values)
{
    const Mesh::HorizontalMesh &horizontal = mesh.horizontal();

    const int     layers  = mesh.layerCount();
    const int64_t columns = mesh.columnCount();

    const std::vector<int64_t> &edgeLeft  = horizontal.edgeLeftFace();
    const std::vector<int64_t> &edgeRight = horizontal.edgeRightFace();
    const std::vector<double>  &area      = horizontal.faceArea();

    FVQual_PARALLEL_FOR
    for (int64_t e = 0; e < horizontal.edgeCount(); ++e)
    {
        const size_t  ei    = static_cast<size_t>(e);
        const size_t  base  = ei * static_cast<size_t>(layers);
        const int64_t left  = edgeLeft[ei];
        const int64_t right = edgeRight[ei];

        if (right == Mesh::kNoFace)
        {
            for (int k = 0; k < layers; ++k)
                m_lowHorizontal[base + static_cast<size_t>(k)] = 0.0;
            continue;
        }

        for (int k = 0; k < layers; ++k)
        {
            const double flux = edgeLayerFlux[base + static_cast<size_t>(k)];
            const size_t donor = static_cast<size_t>(
                flux >= 0.0 ? mesh.cell(left, k) : mesh.cell(right, k));

            m_lowHorizontal[base + static_cast<size_t>(k)] = flux * values[donor];
        }
    }

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < columns; ++c)
    {
        const double columnArea = area[static_cast<size_t>(c)];

        // Surface and bed are material surfaces: no advective flux through
        // either, which the internal mode guarantees by construction.
        m_lowVertical[static_cast<size_t>(mesh.interfaceSlot(c, 0))]      = 0.0;
        m_lowVertical[static_cast<size_t>(mesh.interfaceSlot(c, layers))] = 0.0;

        for (int k = 1; k < layers; ++k)
        {
            const size_t slot = static_cast<size_t>(mesh.interfaceSlot(c, k));
            const double flux = verticalFlux[slot] * columnArea;

            // Positive flux is upward, so it carries the cell below.
            const size_t donor =
                static_cast<size_t>(flux >= 0.0 ? mesh.cell(c, k) : mesh.cell(c, k - 1));

            m_lowVertical[slot] = flux * values[donor];
        }
    }
}

void EulerianTransport::computeGradients(const Mesh::LayeredMesh &mesh,
                                         std::span<const double> values)
{
    const Mesh::HorizontalMesh &horizontal = mesh.horizontal();

    const int     layers  = mesh.layerCount();
    const int64_t columns = mesh.columnCount();

    const std::vector<double>  &centroidX      = horizontal.faceCentroidX();
    const std::vector<double>  &centroidY      = horizontal.faceCentroidY();
    const std::vector<int64_t> &edgeLeft       = horizontal.edgeLeftFace();
    const std::vector<int64_t> &edgeRight      = horizontal.edgeRightFace();
    const std::vector<int64_t> &faceEdgeOffset = horizontal.faceEdgeOffsets();
    const std::vector<int64_t> &faceEdges      = horizontal.faceEdges();

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < columns; ++c)
    {
        const int64_t begin = faceEdgeOffset[static_cast<size_t>(c)];
        const int64_t end   = faceEdgeOffset[static_cast<size_t>(c) + 1];

        for (int k = 0; k < layers; ++k)
        {
            const size_t cell = static_cast<size_t>(mesh.cell(c, k));

            // ── Horizontal: inverse-distance-weighted least squares ───────
            double sumXX = 0.0;
            double sumXY = 0.0;
            double sumYY = 0.0;
            double sumXC = 0.0;
            double sumYC = 0.0;

            double highest = values[cell];
            double lowest  = values[cell];

            for (int64_t index = begin; index < end; ++index)
            {
                const size_t  ei    = static_cast<size_t>(faceEdges[static_cast<size_t>(index)]);
                const bool    owned = edgeLeft[ei] == c;
                const int64_t other = owned ? edgeRight[ei] : edgeLeft[ei];

                if (other == Mesh::kNoFace)
                    continue;

                const size_t neighbour = static_cast<size_t>(mesh.cell(other, k));

                const double dx = centroidX[static_cast<size_t>(other)] -
                                  centroidX[static_cast<size_t>(c)];
                const double dy = centroidY[static_cast<size_t>(other)] -
                                  centroidY[static_cast<size_t>(c)];

                const double squared = dx * dx + dy * dy;

                if (squared <= 0.0)
                    continue;

                const double weight = 1.0 / squared;
                const double dc     = values[neighbour] - values[cell];

                sumXX += weight * dx * dx;
                sumXY += weight * dx * dy;
                sumYY += weight * dy * dy;
                sumXC += weight * dx * dc;
                sumYC += weight * dy * dc;

                highest = std::max(highest, values[neighbour]);
                lowest  = std::min(lowest, values[neighbour]);
            }

            const double determinant = sumXX * sumYY - sumXY * sumXY;

            double gradientX = 0.0;
            double gradientY = 0.0;

            if (std::abs(determinant) > 1.0e-30)
            {
                gradientX = (sumYY * sumXC - sumXY * sumYC) / determinant;
                gradientY = (sumXX * sumYC - sumXY * sumXC) / determinant;
            }
            else if (sumXX > 0.0)
            {
                // Degenerate stencil (a one-dimensional channel): recover the
                // streamwise component rather than giving up entirely.
                gradientX = sumXC / sumXX;
            }

            m_gradientX[cell] = gradientX;
            m_gradientY[cell] = gradientY;

            // ── Barth-Jespersen limiter ───────────────────────────────────
            //
            // Scale the gradient so no face reconstruction leaves the range
            // of this cell and its neighbours. Zalesak will guarantee
            // boundedness afterwards regardless, but limiting here keeps the
            // antidiffusive flux small enough that the limiter rarely has to
            // clip, which is what preserves second-order accuracy in smooth
            // regions.
            double phi = 1.0;

            for (int64_t index = begin; index < end; ++index)
            {
                const size_t  ei    = static_cast<size_t>(faceEdges[static_cast<size_t>(index)]);
                const bool    owned = edgeLeft[ei] == c;
                const int64_t other = owned ? edgeRight[ei] : edgeLeft[ei];

                if (other == Mesh::kNoFace)
                    continue;

                const double dx = horizontal.edgeMidX()[ei] - centroidX[static_cast<size_t>(c)];
                const double dy = horizontal.edgeMidY()[ei] - centroidY[static_cast<size_t>(c)];

                const double delta = gradientX * dx + gradientY * dy;

                if (delta > 1.0e-30)
                    phi = std::min(phi, (highest - values[cell]) / delta);
                else if (delta < -1.0e-30)
                    phi = std::min(phi, (lowest - values[cell]) / delta);
            }

            m_limiter[cell] = std::clamp(phi, 0.0, 1.0);

            // ── Vertical: minmod slope in elevation ───────────────────────
            double slope = 0.0;

            if (k > 0 && k + 1 < layers)
            {
                const size_t above = static_cast<size_t>(mesh.cell(c, k - 1));
                const size_t below = static_cast<size_t>(mesh.cell(c, k + 1));

                const double zAbove  = 0.5 * (mesh.z(c, k - 1) + mesh.z(c, k));
                const double zCentre = 0.5 * (mesh.z(c, k) + mesh.z(c, k + 1));
                const double zBelow  = 0.5 * (mesh.z(c, k + 1) + mesh.z(c, k + 2));

                const double upper = zAbove - zCentre;
                const double lower = zCentre - zBelow;

                if (upper > 0.0 && lower > 0.0)
                {
                    slope = minmod((values[above] - values[cell]) / upper,
                                   (values[cell] - values[below]) / lower);
                }
            }

            m_gradientZ[cell] = slope;
        }
    }
}

void EulerianTransport::computeHighOrderFluxes(const Mesh::LayeredMesh &mesh,
                                               std::span<const double> edgeLayerFlux,
                                               std::span<const double> verticalFlux,
                                               std::span<const double> values)
{
    const Mesh::HorizontalMesh &horizontal = mesh.horizontal();

    const int     layers  = mesh.layerCount();
    const int64_t columns = mesh.columnCount();

    const std::vector<double>  &centroidX = horizontal.faceCentroidX();
    const std::vector<double>  &centroidY = horizontal.faceCentroidY();
    const std::vector<double>  &midX      = horizontal.edgeMidX();
    const std::vector<double>  &midY      = horizontal.edgeMidY();
    const std::vector<int64_t> &edgeLeft  = horizontal.edgeLeftFace();
    const std::vector<int64_t> &edgeRight = horizontal.edgeRightFace();
    const std::vector<double>  &area      = horizontal.faceArea();

    FVQual_PARALLEL_FOR
    for (int64_t e = 0; e < horizontal.edgeCount(); ++e)
    {
        const size_t  ei    = static_cast<size_t>(e);
        const size_t  base  = ei * static_cast<size_t>(layers);
        const int64_t left  = edgeLeft[ei];
        const int64_t right = edgeRight[ei];

        if (right == Mesh::kNoFace)
        {
            for (int k = 0; k < layers; ++k)
                m_highHorizontal[base + static_cast<size_t>(k)] = 0.0;
            continue;
        }

        for (int k = 0; k < layers; ++k)
        {
            const double  flux  = edgeLayerFlux[base + static_cast<size_t>(k)];
            const int64_t donorColumn = flux >= 0.0 ? left : right;
            const size_t  donor = static_cast<size_t>(mesh.cell(donorColumn, k));

            const double dx = midX[ei] - centroidX[static_cast<size_t>(donorColumn)];
            const double dy = midY[ei] - centroidY[static_cast<size_t>(donorColumn)];

            const double reconstructed =
                values[donor] +
                m_limiter[donor] * (m_gradientX[donor] * dx + m_gradientY[donor] * dy);

            m_highHorizontal[base + static_cast<size_t>(k)] = flux * reconstructed;
        }
    }

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < columns; ++c)
    {
        const double columnArea = area[static_cast<size_t>(c)];

        m_highVertical[static_cast<size_t>(mesh.interfaceSlot(c, 0))]      = 0.0;
        m_highVertical[static_cast<size_t>(mesh.interfaceSlot(c, layers))] = 0.0;

        for (int k = 1; k < layers; ++k)
        {
            const size_t slot = static_cast<size_t>(mesh.interfaceSlot(c, k));
            const double flux = verticalFlux[slot] * columnArea;

            const int    donorLayer = flux >= 0.0 ? k : k - 1;
            const size_t donor      = static_cast<size_t>(mesh.cell(c, donorLayer));

            const double interfaceZ = mesh.z(c, k);
            const double donorZ =
                0.5 * (mesh.z(c, donorLayer) + mesh.z(c, donorLayer + 1));

            const double reconstructed =
                values[donor] + m_gradientZ[donor] * (interfaceZ - donorZ);

            m_highVertical[slot] = flux * reconstructed;
        }
    }
}

void EulerianTransport::applyDivergence(const Mesh::LayeredMesh &mesh,
                                        std::span<const double> oldVolume,
                                        std::span<const double> horizontalFlux,
                                        std::span<const double> verticalMassFlux,
                                        std::span<const double> cellSource,
                                        std::span<const double> sourceValue, double dt,
                                        std::span<const double> before,
                                        std::span<double> after) const
{
    const Mesh::HorizontalMesh &horizontal = mesh.horizontal();

    const int     layers  = mesh.layerCount();
    const int64_t columns = mesh.columnCount();

    const std::vector<int64_t> &edgeLeft       = horizontal.edgeLeftFace();
    const std::vector<int64_t> &faceEdgeOffset = horizontal.faceEdgeOffsets();
    const std::vector<int64_t> &faceEdges      = horizontal.faceEdges();

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < columns; ++c)
    {
        const int64_t begin = faceEdgeOffset[static_cast<size_t>(c)];
        const int64_t end   = faceEdgeOffset[static_cast<size_t>(c) + 1];

        for (int k = 0; k < layers; ++k)
        {
            const size_t cell = static_cast<size_t>(mesh.cell(c, k));

            double outflow = 0.0;

            for (int64_t index = begin; index < end; ++index)
            {
                const size_t ei = static_cast<size_t>(faceEdges[static_cast<size_t>(index)]);

                outflow += (edgeLeft[ei] == c ? 1.0 : -1.0) *
                           horizontalFlux[ei * static_cast<size_t>(layers) +
                                          static_cast<size_t>(k)];
            }

            // Leaving through the top, entering through the bottom.
            outflow += verticalMassFlux[static_cast<size_t>(mesh.interfaceSlot(c, k))];
            outflow -= verticalMassFlux[static_cast<size_t>(mesh.interfaceSlot(c, k + 1))];

            double mass = before[cell] * oldVolume[cell] - dt * outflow;

            if (!cellSource.empty())
            {
                const double discharge = cellSource[cell];

                if (discharge != 0.0)
                {
                    // Inflow brings its own composition; outflow leaves at the
                    // ambient value, which keeps the concentration unchanged
                    // while the volume falls.
                    const double arriving =
                        discharge > 0.0 ? sourceValue[cell] : before[cell];

                    mass += dt * discharge * arriving;
                }
            }

            const double volume = mesh.cellVolume()[cell];

            after[cell] = volume > 0.0 ? mass / volume : before[cell];
        }
    }
}

void EulerianTransport::limitAntidiffusiveFluxes(const Mesh::LayeredMesh &mesh,
                                                 std::span<const double> oldValues,
                                                 std::span<const double> lowValues,
                                                 double dt)
{
    const Mesh::HorizontalMesh &horizontal = mesh.horizontal();

    const int     layers  = mesh.layerCount();
    const int64_t columns = mesh.columnCount();

    const std::vector<int64_t> &edgeLeft       = horizontal.edgeLeftFace();
    const std::vector<int64_t> &edgeRight      = horizontal.edgeRightFace();
    const std::vector<int64_t> &faceEdgeOffset = horizontal.faceEdgeOffsets();
    const std::vector<int64_t> &faceEdges      = horizontal.faceEdges();

    // ── Admissible bounds and the antidiffusive budget of every cell ──────
    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < columns; ++c)
    {
        const int64_t begin = faceEdgeOffset[static_cast<size_t>(c)];
        const int64_t end   = faceEdgeOffset[static_cast<size_t>(c) + 1];

        for (int k = 0; k < layers; ++k)
        {
            const size_t cell = static_cast<size_t>(mesh.cell(c, k));

            // Bounds are taken over the old *and* the low-order solution of
            // this cell and its neighbours: Zalesak's choice, and the reason
            // the result cannot exceed what the monotone scheme would allow.
            double highest = std::max(oldValues[cell], lowValues[cell]);
            double lowest  = std::min(oldValues[cell], lowValues[cell]);

            double enteringFlux = 0.0;
            double leavingFlux  = 0.0;

            for (int64_t index = begin; index < end; ++index)
            {
                const size_t  ei    = static_cast<size_t>(faceEdges[static_cast<size_t>(index)]);
                const bool    owned = edgeLeft[ei] == c;
                const int64_t other = owned ? edgeRight[ei] : edgeLeft[ei];

                if (other == Mesh::kNoFace)
                    continue;

                const size_t neighbour = static_cast<size_t>(mesh.cell(other, k));

                highest = std::max({highest, oldValues[neighbour], lowValues[neighbour]});
                lowest  = std::min({lowest, oldValues[neighbour], lowValues[neighbour]});

                const double outward =
                    (owned ? 1.0 : -1.0) *
                    m_antiHorizontal[ei * static_cast<size_t>(layers) +
                                     static_cast<size_t>(k)];

                enteringFlux += std::max(0.0, -outward);
                leavingFlux += std::max(0.0, outward);
            }

            for (int side = 0; side < 2; ++side)
            {
                const int neighbourLayer = side == 0 ? k - 1 : k + 1;

                if (neighbourLayer < 0 || neighbourLayer >= layers)
                    continue;

                const size_t neighbour = static_cast<size_t>(mesh.cell(c, neighbourLayer));

                highest = std::max({highest, oldValues[neighbour], lowValues[neighbour]});
                lowest  = std::min({lowest, oldValues[neighbour], lowValues[neighbour]});

                // Upward flux at interface k leaves this cell; at k+1 it
                // enters from below.
                const double outward =
                    side == 0
                        ? m_antiVertical[static_cast<size_t>(mesh.interfaceSlot(c, k))]
                        : -m_antiVertical[static_cast<size_t>(mesh.interfaceSlot(c, k + 1))];

                enteringFlux += std::max(0.0, -outward);
                leavingFlux += std::max(0.0, outward);
            }

            const double volume = mesh.cellVolume()[cell];
            const double rise   = (highest - lowValues[cell]) * volume / dt;
            const double fall   = (lowValues[cell] - lowest) * volume / dt;

            m_risePlus[cell] =
                enteringFlux > 0.0 ? std::min(1.0, std::max(0.0, rise) / enteringFlux) : 1.0;
            m_riseMinus[cell] =
                leavingFlux > 0.0 ? std::min(1.0, std::max(0.0, fall) / leavingFlux) : 1.0;
        }
    }

    // ── Scale each face by the tighter of its two cells' budgets ──────────
    double admittedSum   = 0.0;
    int64_t admittedFaces = 0;

    FVQual_PARALLEL_FOR_REDUCE(+ : admittedSum, admittedFaces)
    for (int64_t e = 0; e < horizontal.edgeCount(); ++e)
    {
        const size_t  ei    = static_cast<size_t>(e);
        const size_t  base  = ei * static_cast<size_t>(layers);
        const int64_t left  = edgeLeft[ei];
        const int64_t right = edgeRight[ei];

        if (right == Mesh::kNoFace)
            continue;

        for (int k = 0; k < layers; ++k)
        {
            const size_t index = base + static_cast<size_t>(k);
            const double flux  = m_antiHorizontal[index];

            if (flux == 0.0)
                continue;

            const size_t donor    = static_cast<size_t>(mesh.cell(left, k));
            const size_t receiver = static_cast<size_t>(mesh.cell(right, k));

            const double scale = flux > 0.0
                                     ? std::min(m_riseMinus[donor], m_risePlus[receiver])
                                     : std::min(m_riseMinus[receiver], m_risePlus[donor]);

            m_antiHorizontal[index] = flux * scale;

            admittedSum += scale;
            ++admittedFaces;
        }
    }

    for (int64_t c = 0; c < columns; ++c)
    {
        for (int k = 1; k < layers; ++k)
        {
            const size_t slot = static_cast<size_t>(mesh.interfaceSlot(c, k));
            const double flux = m_antiVertical[slot];

            if (flux == 0.0)
                continue;

            // Upward positive: the cell below donates, the cell above receives.
            const size_t below = static_cast<size_t>(mesh.cell(c, k));
            const size_t above = static_cast<size_t>(mesh.cell(c, k - 1));

            const double scale = flux > 0.0
                                     ? std::min(m_riseMinus[below], m_risePlus[above])
                                     : std::min(m_riseMinus[above], m_risePlus[below]);

            m_antiVertical[slot] = flux * scale;

            admittedSum += scale;
            ++admittedFaces;
        }
    }

    m_admittedFraction =
        admittedFaces > 0 ? admittedSum / static_cast<double>(admittedFaces) : 1.0;
}

void EulerianTransport::advect(const Mesh::LayeredMesh &mesh,
                               std::span<const double> oldVolume,
                               std::span<const double> edgeLayerFlux,
                               std::span<const double> verticalFlux,
                               std::span<const double> cellSource,
                               std::span<const double> sourceValue, double dt,
                               AdvectionScheme scheme, std::span<double> values)
{
    if (dt <= 0.0 || mesh.cellCount() <= 0)
        return;

    computeLowOrderFluxes(mesh, edgeLayerFlux, verticalFlux, values);

    if (scheme == AdvectionScheme::Upwind)
    {
        applyDivergence(mesh, oldVolume, m_lowHorizontal, m_lowVertical, cellSource,
                        sourceValue, dt, values, m_lowValues);

        std::copy(m_lowValues.begin(), m_lowValues.end(), values.begin());
        m_admittedFraction = 0.0;
        return;
    }

    // The monotone solution, which both bounds the answer and anchors the
    // correction.
    applyDivergence(mesh, oldVolume, m_lowHorizontal, m_lowVertical, cellSource,
                    sourceValue, dt, values, m_lowValues);

    computeGradients(mesh, values);
    computeHighOrderFluxes(mesh, edgeLayerFlux, verticalFlux, values);

    for (size_t i = 0; i < m_antiHorizontal.size(); ++i)
        m_antiHorizontal[i] = m_highHorizontal[i] - m_lowHorizontal[i];

    for (size_t i = 0; i < m_antiVertical.size(); ++i)
        m_antiVertical[i] = m_highVertical[i] - m_lowVertical[i];

    limitAntidiffusiveFluxes(mesh, values, m_lowValues, dt);

    // Apply the admitted correction to the monotone solution. The volumes
    // here are the *new* ones on both sides, because the correction is a pure
    // redistribution: it moves mass between cells without changing the total.
    const int     layers  = mesh.layerCount();
    const int64_t columns = mesh.columnCount();

    const Mesh::HorizontalMesh &horizontal     = mesh.horizontal();
    const std::vector<int64_t> &edgeLeft       = horizontal.edgeLeftFace();
    const std::vector<int64_t> &faceEdgeOffset = horizontal.faceEdgeOffsets();
    const std::vector<int64_t> &faceEdges      = horizontal.faceEdges();

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < columns; ++c)
    {
        const int64_t begin = faceEdgeOffset[static_cast<size_t>(c)];
        const int64_t end   = faceEdgeOffset[static_cast<size_t>(c) + 1];

        for (int k = 0; k < layers; ++k)
        {
            const size_t cell   = static_cast<size_t>(mesh.cell(c, k));
            const double volume = mesh.cellVolume()[cell];

            if (volume <= 0.0)
            {
                values[cell] = m_lowValues[cell];
                continue;
            }

            double outflow = 0.0;

            for (int64_t index = begin; index < end; ++index)
            {
                const size_t ei = static_cast<size_t>(faceEdges[static_cast<size_t>(index)]);

                outflow += (edgeLeft[ei] == c ? 1.0 : -1.0) *
                           m_antiHorizontal[ei * static_cast<size_t>(layers) +
                                            static_cast<size_t>(k)];
            }

            outflow += m_antiVertical[static_cast<size_t>(mesh.interfaceSlot(c, k))];
            outflow -= m_antiVertical[static_cast<size_t>(mesh.interfaceSlot(c, k + 1))];

            values[cell] = m_lowValues[cell] - dt * outflow / volume;
        }
    }
}

} // namespace FVQual::Transport
