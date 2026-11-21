// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   turbulence.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  Two-equation k-epsilon vertical mixing (F5.4).
 *
 * **Why this exists, measured rather than assumed.** Richardson damping was the
 * mixing scheme through Phase 3. Against DeGray's measured casts, matched runs
 * to day 150 at 30 layers, 333 measurements, overall RMSE:
 *
 * | | overall | 5-10 m | 20-30 m | 30-60 m |
 * |---|---|---|---|---|
 * | no wind stress | **2.29 °C** | −3.91 | −0.83 | −0.55 |
 * | Richardson + wind | 3.59 °C | −1.11 | +3.22 | +3.65 |
 * | k-epsilon + wind | **3.42 °C** | −1.23 | +2.95 | +3.39 |
 *
 * Wind stress is what fixes the epilimnion: the 5-10 m band goes from −3.91 to
 * −1.11, and without it the modeled April surface layer is 3 m against a
 * measured 14 m. What wind then costs is the interior, which both schemes
 * over-mix by 3 to 4 °C.
 *
 * The closure is a modest improvement over Richardson damping, better in every
 * band below 10 m and 5% better overall. **An earlier version of this comment
 * claimed the opposite** — that the closure was worse, 8.95 °C against 8.35 —
 * and that comparison was invalid: both runs predated the fix to the vertical
 * Courant limit and were spending their time at the time-step floor. The
 * lesson is the measurement's, not the closure's: a comparison between two
 * schemes is only worth its numbers if the model underneath is stepping
 * correctly.
 *
 * The missing physics is a *memory*. Richardson damping is diagnostic: mixing
 * is a function of the local gradients right now, so a stirred surface layer
 * cannot persist into calm, and turbulence cannot be transported away from
 * where it was generated. A two-equation closure carries turbulent kinetic
 * energy and its dissipation as transported quantities, which is what lets a
 * wind-mixed layer sit on top of an intact thermocline — the configuration
 * every measured DeGray profile shows and neither current setting produces.
 *
 * **Formulation.** Standard k-epsilon, which is the `gen`-family GLS closure at
 * `p = 3, m = 1.5, n = -1`:
 *
 *     Dk/Dt  = d/dz( nu_t/sigma_k  dk/dz )  + P + B - eps
 *     Deps/Dt= d/dz( nu_t/sigma_eps deps/dz) + eps/k (c1 P + c3 B - c2 eps)
 *
 * with shear production `P = nu_t S²`, buoyancy production `B = -K_t N²`, and
 * `nu_t = c_mu k²/eps`. Stability functions are Galperin's quasi-equilibrium
 * form, which is algebraic and bounded — enough for a stratified reservoir and
 * far cheaper than a full second-moment closure.
 *
 * \threadsafety One column at a time; the caller owns any parallelism over
 *               columns, as it does for the existing tridiagonal solves.
 */

#ifndef FVQual_HYDRO_TURBULENCE_H
#define FVQual_HYDRO_TURBULENCE_H

#include "fvqual/fvqual_export.h"

#include <cstdint>
#include <vector>

namespace FVQual::Hydro {

/*!
 * \brief Closure coefficients. The defaults are the standard k-epsilon set.
 */
struct FVQual_EXPORT TurbulenceOptions
{
    double cMu     = 0.09;  //!< Eddy-viscosity coefficient.
    double c1      = 1.44;  //!< Shear-production coefficient in the eps equation.
    double c2      = 1.92;  //!< Dissipation-of-dissipation coefficient.
    double c3Stable = -0.4; //!< Buoyancy coefficient, stable stratification.
    double c3Unstable = 1.0;//!< Buoyancy coefficient, unstable stratification.

    double sigmaK   = 1.0;  //!< Schmidt number for k.
    double sigmaEps = 1.3;  //!< Schmidt number for epsilon.

    //! Von Karman constant, for the wall boundary conditions.
    double vonKarman = 0.41;

    /*!
     * \brief Floors (m²/s, m²/s², m²/s³).
     *
     * A closure needs floors or it divides by a vanishing `eps` in still water.
     * These are the smallest values with physical meaning, not tuning knobs —
     * `minDiffusivity` is molecular, and the k and eps floors are what a lake
     * retains under no forcing at all.
     */
    double minTke          = 1.0e-8;
    double minDissipation  = 1.0e-12;
    double minViscosity    = 1.0e-6;
    double minDiffusivity  = 1.0e-7;

    //! Ceiling, so a bad step cannot produce an absurd diffusivity.
    double maxViscosity = 0.1;

    //! Bed roughness length (m), for the bottom wall condition.
    double bedRoughness = 0.005;

    //! Surface roughness length (m), for the wind-driven wall condition.
    double surfaceRoughness = 0.02;
};

/*!
 * \brief k-epsilon closure over one column.
 *
 * Holds `k` and `eps` at layer interfaces, where the shear and buoyancy
 * gradients naturally live and where the mixing coefficients are consumed.
 */
class FVQual_EXPORT TurbulenceClosure
{
public:
    /*!
     * \brief Sizes the state for a mesh.
     *
     * \param[in] columns Column count.
     * \param[in] layers  Layers per column; interfaces are `layers + 1`.
     */
    void initialize(int64_t columns, int layers);

    /*!
     * \brief Advances one column and writes its mixing coefficients.
     *
     * Inputs are per-interface, index 0 at the surface through `layers` at the
     * bed, matching the mesh's interface ordering.
     *
     * \param[in]  column      Column index.
     * \param[in]  dt          Time step (s).
     * \param[in]  thickness   Layer thicknesses (m), `layers` entries.
     * \param[in]  shearSquared Squared velocity shear at each interface (1/s²).
     * \param[in]  buoyancy    Squared buoyancy frequency N² at each interface
     *                         (1/s²); negative where unstable.
     * \param[in]  frictionSurface Surface friction velocity u* (m/s), from the
     *                         wind stress.
     * \param[in]  frictionBed Bed friction velocity u* (m/s).
     * \param[out] viscosity   Eddy viscosity per interface (m²/s).
     * \param[out] diffusivity Eddy diffusivity per interface (m²/s).
     */
    void advanceColumn(int64_t column, double dt, const double *thickness,
                       const double *shearSquared, const double *buoyancy,
                       double frictionSurface, double frictionBed,
                       double *viscosity, double *diffusivity);

    //! Turbulent kinetic energy per interface (m²/s²), for diagnostics.
    [[nodiscard]] const std::vector<double> &tke() const noexcept { return m_tke; }

    //! Dissipation rate per interface (m²/s³).
    [[nodiscard]] const std::vector<double> &dissipation() const noexcept
    {
        return m_dissipation;
    }

    //! Mutable options.
    [[nodiscard]] TurbulenceOptions &options() noexcept { return m_options; }
    [[nodiscard]] const TurbulenceOptions &options() const noexcept
    {
        return m_options;
    }

private:
    TurbulenceOptions m_options;

    int     m_layers  = 0;
    int64_t m_columns = 0;

    std::vector<double> m_tke;
    std::vector<double> m_dissipation;

    //! Per-column scratch, sized on initialize().
    std::vector<double> m_lower, m_diagonal, m_upper, m_rhs, m_solution, m_sweep;
};

/*!
 * \brief Galperin quasi-equilibrium stability function for momentum.
 *
 * Exposed because it is the part most worth checking independently: it must
 * fall as stratification strengthens and stay positive and bounded, and a sign
 * slip here silently inverts the response to stratification.
 *
 * \param[in] alphaN Dimensionless stratification `N² k²/eps²`.
 */
[[nodiscard]] FVQual_EXPORT double stabilityMomentum(double alphaN) noexcept;

//! Galperin stability function for scalars.
[[nodiscard]] FVQual_EXPORT double stabilityScalar(double alphaN) noexcept;

} // namespace FVQual::Hydro

#endif // FVQual_HYDRO_TURBULENCE_H
