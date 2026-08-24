// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   lightfield.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  Beer-Lambert irradiance down each water column (F6.2a).
 *
 * Photosynthesis is the reason a reservoir water-quality model needs a third
 * spatial dimension at all. Algal growth is limited by light, light is
 * exhausted with depth, and what exhausts it includes the algae themselves —
 * so the light field is both an input to the kinetics and a function of them.
 *
 * **The march.** Attenuation is piecewise constant within a layer, which makes
 * the discrete recursion
 *
 *     I(k+1) = I(k) exp(-gamma_k h_k)
 *
 * *exact* rather than approximate — the exponential's own composition law does
 * the work. That is worth knowing, because it means the interface irradiances
 * can be gated against `I0 exp(-gamma z)` to round-off on a uniform column
 * rather than to a discretization tolerance.
 *
 * **Cell values are layer means, not top-of-layer values.** A finite-volume
 * cell holds an average, and for light the distinction is not cosmetic: over a
 * layer of optical depth `tau = gamma h`, the mean is smaller than the value at
 * the top by `(1 - e^-tau)/tau`, which for a 2 m layer in turbid water is a
 * factor of two. Using the top value would bias growth high everywhere, worst
 * in the surface layer where most production happens. The mean has a closed
 * form here, so the honest choice is also the cheap one.
 *
 * **Surface absorption.** A fraction of incident shortwave is absorbed in the
 * first few millimetres and never participates in photosynthesis. It is
 * subtracted before the march and reported separately by surfaceAbsorbed(),
 * because it is a real heat input that a surface heat budget will want rather
 * than energy this class is entitled to discard.
 *
 * \threadsafety One field per solver; update() is internally parallel over
 *               columns, and each column's march is independent.
 */

#ifndef FVQual_KINETICS_LIGHTFIELD_H
#define FVQual_KINETICS_LIGHTFIELD_H

#include "fvqual/fvqual_export.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/state/speciesregistry.h"
#include "fvqual/state/statestore.h"

#include <string>
#include <vector>

namespace FVQual::Kinetics {

/*!
 * \brief Optical properties and the incident forcing.
 */
struct FVQual_EXPORT LightOptions
{
    /*!
     * \brief Shortwave irradiance reaching the water surface (W/m²).
     *
     * Time-varying forcing: set it per step from meteorology, or leave it
     * constant for a controlled experiment.
     */
    double incidentShortwave = 0.0;

    /*!
     * \brief Fraction absorbed at the surface and unavailable below (-).
     *
     * CE-QUAL-W2's BETA. The default 0.45 is its usual value: roughly the
     * non-photosynthetically-active part of the spectrum, which is absorbed
     * within centimetres.
     */
    double surfaceAbsorption = 0.45;

    /*!
     * \brief Attenuation of clear water, before any constituent (1/m).
     *
     * W2's EXH2O. 0.25 1/m is a typical clear-reservoir value; turbid systems
     * run several times higher.
     */
    double backgroundAttenuation = 0.25;
};

/*!
 * \brief Irradiance per cell and per interface, from a Beer-Lambert march.
 */
class FVQual_EXPORT LightField
{
public:
    /*!
     * \brief Sizes the field for \p mesh and \p species.
     *
     * \param[in] mesh    Layered mesh the field is computed on.
     * \param[in] species Species registry, for the self-shading table.
     */
    void initialize(const Mesh::LayeredMesh &mesh,
                    const State::SpeciesRegistry &species);

    /*!
     * \brief Sets a species' contribution to attenuation (m²/g).
     *
     * Self-shading: the concentration of algae, detritus, or suspended solids
     * multiplies this coefficient and adds to the local attenuation. Zero for
     * every species until set, so a species only shades if it was said to.
     *
     * \param[in]  species Species index.
     * \param[in]  value   Specific attenuation (m²/g), non-negative.
     * \param[out] message Failure description when returning false.
     */
    [[nodiscard]] bool setSpecificAttenuation(int species, double value,
                                              std::string &message);

    //! Specific attenuation of \p species (m²/g).
    [[nodiscard]] double specificAttenuation(int species) const;

    /*!
     * \brief Marches irradiance down every column from the current state.
     *
     * Call after transport and before the reactions, so the shading reflects
     * the constituent field the reactions will see.
     *
     * \param[in] mesh    Layered mesh, for thicknesses.
     * \param[in] state   Current state, for the shading species.
     * \param[in] options Incident forcing and optical properties.
     */
    void update(const Mesh::LayeredMesh &mesh, const State::StateStore &state,
                const LightOptions &options);

    /*!
     * \brief Layer-mean irradiance per cell (W/m²) — what a rate law reads.
     */
    [[nodiscard]] const std::vector<double> &meanIrradiance() const noexcept
    {
        return m_mean;
    }

    /*!
     * \brief Irradiance at each layer interface (W/m²), by interface slot.
     *
     * Slot `interfaceSlot(c, 0)` is the surface, after surface absorption;
     * slot `interfaceSlot(c, layerCount)` is the bed. Exposed because the
     * per-layer absorbed energy is the difference across a cell, which is what
     * a solar heat source needs, and because it is the quantity with an exact
     * closed form to gate against.
     */
    [[nodiscard]] const std::vector<double> &interfaceIrradiance() const noexcept
    {
        return m_interface;
    }

    /*!
     * \brief Irradiance absorbed at the surface of each column (W/m²).
     *
     * The part removed before the march. Reported rather than discarded: it is
     * heat that enters the water, and a heat budget that ignored it would be
     * short by nearly half the incident shortwave.
     */
    [[nodiscard]] const std::vector<double> &surfaceAbsorbed() const noexcept
    {
        return m_surfaceAbsorbed;
    }

    //! Whether initialize() has been called.
    [[nodiscard]] bool ready() const noexcept { return !m_mean.empty(); }

private:
    std::vector<double> m_mean;
    std::vector<double> m_interface;
    std::vector<double> m_surfaceAbsorbed;
    std::vector<double> m_attenuation; //!< Specific attenuation per species.
};

/*!
 * \brief Mean of `exp(-tau s)` over `s` in [0, 1] — the layer-mean factor.
 *
 * `(1 - e^-tau)/tau`, written so it stays accurate as `tau` vanishes. The naive
 * form loses every significant digit there and then divides by zero at
 * `tau = 0`, which is an unlit or vanishingly thin layer — not an exotic case.
 */
[[nodiscard]] FVQual_EXPORT double meanAttenuationFactor(double tau) noexcept;

} // namespace FVQual::Kinetics

#endif // FVQual_KINETICS_LIGHTFIELD_H
