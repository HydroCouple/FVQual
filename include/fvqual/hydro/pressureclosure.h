// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   pressureclosure.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  IPressureClosure — the pressure seam of the internal mode, and its
 *         hydrostatic implementation.
 *
 * **Why this is an interface.** FVQual is hydrostatic: pressure is the
 * weight of the water above, and vertical momentum is diagnostic. That covers
 * the entire CE-QUAL-W2 application space — stratification, wind setup,
 * plunging inflows, selective withdrawal, basin-scale internal seiches. It
 * fails for nonlinear internal solitons, Kelvin-Helmholtz billows, and
 * aeration plumes, where vertical and horizontal scales converge.
 *
 * Adding a non-hydrostatic correction later means a 3-D elliptic solve, which
 * would drag an AMG-preconditioned Krylov stack into a deliberately
 * solver-free design. What it must *not* also mean is rewriting the momentum
 * assembly. So every pressure term reaches the momentum equation through this
 * interface, and the flux assembly never learns where the pressure came from.
 * A future `NonHydrostaticPressureClosure` adds its q-correction behind the
 * same call. That, plus carrying `w` as a full prognostic-capable array, is
 * the whole cost of keeping the option open.
 *
 * **What "baroclinic" must mean here.** The integrand is the density
 * *anomaly* rho - rho0, never the full density. This is a correctness
 * requirement, not a normalization preference, because
 *
 *     integral(z..eta) rho  =  rho0 (eta - z) + integral(z..eta) rho'
 *
 * so integrating the full density silently carries the whole barotropic
 * surface-slope term rho0 d(eta)/dx inside the supposedly baroclinic part.
 * The external mode already computes that term through its own pressure flux,
 * so the duplicate arrives there as a *lagged* body force held fixed across a
 * whole subcycling interval — a restoring force that keeps pushing after the
 * surface has already responded. That is a resonantly forced oscillator, and
 * it drove a free surface that should have tilted by millimetres to several
 * metres. Nothing about it is visible while the surface is flat, which is why
 * every rest-state test passed throughout.
 *
 * **The sigma pressure-gradient problem (risk R1).** Terrain-following layers
 * slope with the bed, so the horizontal pressure gradient becomes the small
 * difference of two large terms, and truncation error appears as currents in
 * water that should be still. The defence is *Green-Gauss with a shared
 * elevation*: every face of a cell is evaluated at that cell's own
 * layer-centre elevation, and each side integrates its own profile down to
 * that same true elevation rather than along the sloping layer. For a
 * horizontally uniform ocean the face pressures are then identical and
 * multiply a closed sum of outward face normals, which vanishes identically —
 * so the residual comes only from the discrete representation of the profile,
 * not from the coordinate itself.
 *
 * Two further mitigations were implemented, measured, and found not to help,
 * recorded here so they are not re-attempted blindly:
 *
 * - *Reference-profile subtraction* (a horizontally uniform rho_ref(z)):
 *   removed. In this Green-Gauss form any function of elevation alone
 *   multiplies the same closed normal sum and cancels identically; measured
 *   effect on the spurious velocity was nil. Note this argument applies only
 *   to a z-dependent profile — subtracting the *constant* rho0 is a different
 *   matter entirely and is mandatory, per the paragraph above.
 * - *Piecewise-linear density reconstruction*: retained behind an option but
 *   measured at roughly 0.3% of the spurious velocity over a steep bed, not
 *   the order of magnitude first assumed. The dominant residual is the
 *   below-bed extrapolation where two columns have different depths, which is
 *   where future work should go.
 *
 * \threadsafety One closure per thread; the working buffers are owned.
 */

#ifndef FVQual_HYDRO_PRESSURECLOSURE_H
#define FVQual_HYDRO_PRESSURECLOSURE_H

#include "fvqual/fvqual_export.h"
#include "fvqual/mesh/layeredmesh.h"

#include <span>
#include <string>
#include <vector>

namespace FVQual::Hydro {

/*!
 * \brief The pressure seam: converts a density field into the horizontal
 *        acceleration the internal momentum equation sees.
 */
class FVQual_EXPORT IPressureClosure
{
public:
    virtual ~IPressureClosure() = default;

    /*!
     * \brief Human-readable closure name for diagnostics.
     */
    [[nodiscard]] virtual std::string name() const = 0;

    /*!
     * \brief Whether this closure carries a non-hydrostatic correction.
     */
    [[nodiscard]] virtual bool isNonHydrostatic() const noexcept = 0;

    /*!
     * \brief Computes the baroclinic horizontal acceleration per cell.
     *
     * The barotropic part (the surface-slope term) is deliberately excluded:
     * it belongs to the external mode of the splitting.
     *
     * \param[in]  mesh             Layered mesh, with current layer geometry.
     * \param[in]  density          Cell densities (kg/m³), one per cell.
     * \param[in]  referenceDensity Boussinesq reference density (kg/m³).
     * \param[in]  gravity          Gravitational acceleration (m/s²).
     * \param[out] accelerationX    Cell-centred acceleration, x (m/s²).
     * \param[out] accelerationY    Cell-centred acceleration, y (m/s²).
     */
    virtual void baroclinicAcceleration(const Mesh::LayeredMesh &mesh,
                                        std::span<const double> density,
                                        double referenceDensity, double gravity,
                                        std::vector<double> &accelerationX,
                                        std::vector<double> &accelerationY) = 0;
};

/*!
 * \brief Hydrostatic pressure closure with sigma-error mitigation.
 */
class FVQual_EXPORT HydrostaticPressureClosure final : public IPressureClosure
{
public:
    /*!
     * \brief Tuning of the sigma-error defences.
     */
    struct Options
    {
        /*!
         * \brief Reconstruct density linearly within each layer, minmod
         *        limited, instead of treating it as piecewise constant.
         *
         * Measured benefit is small — under a percent of the spurious
         * velocity over steep bathymetry, not the order of magnitude
         * originally assumed. Kept because it is nearly free and because it
         * will matter once the below-bed treatment improves.
         */
        bool linearDensityReconstruction = true;
    };

    HydrostaticPressureClosure() = default;
    explicit HydrostaticPressureClosure(const Options &options) : m_options(options) {}

    [[nodiscard]] std::string name() const override;
    [[nodiscard]] bool        isNonHydrostatic() const noexcept override { return false; }

    void baroclinicAcceleration(const Mesh::LayeredMesh &mesh,
                                std::span<const double> density,
                                double referenceDensity, double gravity,
                                std::vector<double> &accelerationX,
                                std::vector<double> &accelerationY) override;

private:
    /*!
     * \brief Precomputes each column's cumulative density integral at its
     *        interfaces, and the limited in-layer density slopes.
     */
    void prepareColumns(const Mesh::LayeredMesh &mesh, std::span<const double> density,
                        double referenceDensity);

    /*!
     * \brief Density integral of a column from \p z up to its free surface
     *        (kg/m²) — the quantity whose horizontal gradient is the
     *        baroclinic force.
     */
    [[nodiscard]] double columnIntegral(const Mesh::LayeredMesh &mesh,
                                        std::span<const double> density,
                                        double referenceDensity, int64_t column,
                                        double z) const;

    Options m_options;

    std::vector<double> m_cumulative; //!< Density *anomaly* integral, per column/interface.
    std::vector<double> m_slope;      //!< Per cell, d(rho)/dz within the layer.
};

} // namespace FVQual::Hydro

#endif // FVQual_HYDRO_PRESSURECLOSURE_H
