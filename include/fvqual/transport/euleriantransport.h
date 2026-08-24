// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   euleriantransport.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  EulerianTransport — conservative scalar advection with MUSCL
 *         reconstruction and Zalesak flux correction.
 *
 * This is the primary water-quality transport engine. The choice of *strictly
 * conservative Eulerian* over a Lagrangian scheme is deliberate, and the
 * reason is that a reservoir eutrophication or dissolved-oxygen study is a
 * budget: users expect the mass they put in to be accounted for, and a scheme
 * that loses a percent per month cannot answer the question no matter how
 * pretty its fields look.
 *
 * **Why two schemes are computed and blended.** First-order upwind is
 * unconditionally bounded but so diffusive that it destroys the very features
 * a stratified model exists to resolve — it smears a thermocline, and a
 * smeared thermocline weakens the density gradient that drives the flow. That
 * error is not cosmetic: it made the internal seiche period run 19% long and
 * the lock-exchange front 19% slow. A second-order MUSCL reconstruction fixes
 * the diffusion but overshoots at sharp fronts, and an overshoot in a
 * concentration is a negative dissolved oxygen or a negative algal biomass —
 * which a nonlinear kinetics routine will happily turn into nonsense.
 *
 * So both are computed: the upwind flux as a guaranteed-bounded fallback, the
 * MUSCL flux for accuracy, and Zalesak's limiter decides face by face how
 * much of the difference between them can be admitted without creating a new
 * extremum. Where the solution is smooth the limiter admits all of it and the
 * scheme is second order; at a front it admits only what the bounds allow.
 *
 * Boundedness is therefore a property of the algorithm, not of a tolerance:
 * the result is guaranteed to lie within the range of the surrounding
 * upwind solution, so a positive field stays positive.
 *
 * **On the engine interface.** A `TransportEngine`
 * abstraction with Eulerian, ELM, and particle implementations behind it. That
 * interface is deliberately *not* introduced here: extracting an abstraction
 * from a single implementation guesses at the seam, and the second engine is
 * what will reveal where it actually belongs.
 *
 * \threadsafety One instance per thread; the working buffers are owned.
 */

#ifndef FVQual_TRANSPORT_EULERIANTRANSPORT_H
#define FVQual_TRANSPORT_EULERIANTRANSPORT_H

#include "fvqual/fvqual_export.h"
#include "fvqual/mesh/layeredmesh.h"

#include <span>
#include <vector>

namespace FVQual::Transport {

/*!
 * \brief Spatial discretization of the advective flux.
 */
enum class AdvectionScheme
{
    /*!
     * \brief First-order upwind: bounded, cheap, and heavily diffusive.
     *        Retained because it is the fallback the flux correction blends
     *        against, and because it makes the accuracy comparison in the
     *        verification suite an honest one.
     */
    Upwind,

    /*!
     * \brief MUSCL reconstruction limited by Zalesak flux correction.
     */
    MusclFct
};

/*!
 * \brief Conservative scalar advection on the layered mesh.
 */
class FVQual_EXPORT EulerianTransport
{
public:
    /*!
     * \brief Sizes the working buffers for a mesh.
     */
    void initialize(const Mesh::LayeredMesh &mesh);

    /*!
     * \brief Advects one species field over a time step.
     *
     * The caller supplies the fluxes rather than the velocities, because the
     * fluxes are the quantity the mode splitting has already reconciled: they
     * sum, across each face, to the transport the external mode actually
     * delivered. Recomputing them from velocities here would silently break
     * that reconciliation.
     *
     * The volume identity `V_new - V_old = -dt * (net flux)` must hold — it
     * does, by construction of the internal mode — and is what makes a
     * uniform field stay uniform.
     *
     * \param[in]     mesh          Layered mesh, carrying the *new* geometry.
     * \param[in]     oldVolume     Cell volumes before the step (m³).
     * \param[in]     edgeLayerFlux Horizontal volume flux per (edge, layer),
     *                              m³/s, positive from the edge's left face
     *                              toward its right.
     * \param[in]     verticalFlux  Upward volume flux per unit plan area at
     *                              each interface (m/s), indexed by
     *                              `LayeredMesh::interfaceSlot`.
     * \param[in]     cellSource    Volumetric source per cell (m³/s), positive
     *                              adding water; empty when there are none.
     *                              This term belongs *here* rather than as a
     *                              correction applied afterwards: continuity
     *                              already carried the source into the vertical
     *                              flux field, so the divergence has changed
     *                              the concentration whether or not anything
     *                              is done about it. Withdrawing water without
     *                              also withdrawing its constituent load
     *                              concentrates what remains, which looks
     *                              exactly like a source that is not there.
     * \param[in]     sourceValue   Concentration arriving with the inflow, per
     *                              cell. Ignored where \p cellSource is
     *                              negative, since an outflow leaves at the
     *                              ambient value by definition.
     * \param[in]     dt            Time step (s).
     * \param[in]     scheme        Advection scheme.
     * \param[in,out] values        Cell concentrations, updated in place.
     */
    void advect(const Mesh::LayeredMesh &mesh, std::span<const double> oldVolume,
                std::span<const double> edgeLayerFlux,
                std::span<const double> verticalFlux,
                std::span<const double> cellSource, std::span<const double> sourceValue,
                double dt, AdvectionScheme scheme, std::span<double> values);

    /*!
     * \brief Advects one species with no volumetric sources.
     */
    void advect(const Mesh::LayeredMesh &mesh, std::span<const double> oldVolume,
                std::span<const double> edgeLayerFlux,
                std::span<const double> verticalFlux, double dt, AdvectionScheme scheme,
                std::span<double> values)
    {
        advect(mesh, oldVolume, edgeLayerFlux, verticalFlux, {}, {}, dt, scheme, values);
    }

    /*!
     * \brief Fraction of the antidiffusive correction admitted by the limiter
     *        on the last call, averaged over faces: 1 where the solution is
     *        smooth, falling toward 0 at sharp fronts. A useful diagnostic of
     *        how hard the limiter is working.
     */
    [[nodiscard]] double lastAdmittedFraction() const noexcept
    {
        return m_admittedFraction;
    }

private:
    void computeLowOrderFluxes(const Mesh::LayeredMesh &mesh,
                               std::span<const double> edgeLayerFlux,
                               std::span<const double> verticalFlux,
                               std::span<const double> values);

    void computeGradients(const Mesh::LayeredMesh &mesh, std::span<const double> values);

    void computeHighOrderFluxes(const Mesh::LayeredMesh &mesh,
                                std::span<const double> edgeLayerFlux,
                                std::span<const double> verticalFlux,
                                std::span<const double> values);

    void applyDivergence(const Mesh::LayeredMesh &mesh, std::span<const double> oldVolume,
                         std::span<const double> horizontal,
                         std::span<const double> vertical,
                         std::span<const double> cellSource,
                         std::span<const double> sourceValue, double dt,
                         std::span<const double> before, std::span<double> after) const;

    void limitAntidiffusiveFluxes(const Mesh::LayeredMesh &mesh,
                                  std::span<const double> oldValues,
                                  std::span<const double> lowValues, double dt);

    // Mass fluxes: low order, high order, and the limited antidiffusive part.
    std::vector<double> m_lowHorizontal;
    std::vector<double> m_lowVertical;
    std::vector<double> m_highHorizontal;
    std::vector<double> m_highVertical;
    std::vector<double> m_antiHorizontal;
    std::vector<double> m_antiVertical;

    // Reconstruction.
    std::vector<double> m_gradientX;
    std::vector<double> m_gradientY;
    std::vector<double> m_gradientZ;
    std::vector<double> m_limiter;

    // Zalesak accumulators.
    std::vector<double> m_lowValues;
    std::vector<double> m_risePlus;
    std::vector<double> m_riseMinus;

    double m_admittedFraction = 1.0;
};

} // namespace FVQual::Transport

#endif // FVQual_TRANSPORT_EULERIANTRANSPORT_H
