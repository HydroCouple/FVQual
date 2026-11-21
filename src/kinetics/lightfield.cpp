// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/kinetics/lightfield.h"

#include "fvqual/parallel.h"

#include <cmath>

namespace FVQual::Kinetics {

double meanAttenuationFactor(double tau) noexcept
{
    // Below this the series is more accurate than the closed form, which is
    // computing (1 - (1 - tau + ...))/tau out of two nearly equal numbers.
    // expm1 handles the rest of the range without cancellation.
    constexpr double kSmall = 1.0e-8;

    if (std::abs(tau) < kSmall)
        return 1.0 - 0.5 * tau;

    return -std::expm1(-tau) / tau;
}

void LightField::initialize(const Mesh::LayeredMesh &mesh,
                            const State::SpeciesRegistry &species)
{
    m_mean.assign(static_cast<size_t>(mesh.cellCount()), 0.0);
    m_interface.assign(static_cast<size_t>(mesh.interfaceSlotCount()), 0.0);
    m_surfaceAbsorbed.assign(static_cast<size_t>(mesh.columnCount()), 0.0);
    m_attenuation.assign(static_cast<size_t>(species.count()), 0.0);
}

bool LightField::setSpecificAttenuation(int species, double value,
                                        std::string &message)
{
    if (species < 0 || species >= static_cast<int>(m_attenuation.size()))
    {
        message = "species index " + std::to_string(species) +
                  " is out of range for the light field";
        return false;
    }

    if (!(value >= 0.0))
    {
        // A negative coefficient is a constituent that *creates* light on its
        // way down. It would also let the total attenuation go negative and the
        // march grow without bound.
        message = "specific attenuation must be non-negative";
        return false;
    }

    m_attenuation[static_cast<size_t>(species)] = value;

    message.clear();
    return true;
}

double LightField::specificAttenuation(int species) const
{
    if (species < 0 || species >= static_cast<int>(m_attenuation.size()))
        return 0.0;

    return m_attenuation[static_cast<size_t>(species)];
}

void LightField::update(const Mesh::LayeredMesh &mesh, const State::StateStore &state,
                        const LightOptions &options)
{
    const int64_t columns = mesh.columnCount();
    const int     layers  = mesh.layerCount();
    const int     species = state.speciesCount();

    const std::vector<double> &thickness = mesh.thickness();

    FVQual_PARALLEL_FOR
    for (int64_t c = 0; c < columns; ++c)
    {
        // What arrives at this column: one value for the whole mesh, or a
        // per-column field where shading, cloud or a coupled atmosphere vary
        // it. The scalar branch is the original arithmetic untouched.
        const double incident =
            options.incidentPerColumn != nullptr
                ? (*options.incidentPerColumn)[static_cast<size_t>(c)]
                : options.incidentShortwave;

        // What survives the surface, and what the surface kept.
        const double penetrating = incident * (1.0 - options.surfaceAbsorption);

        m_surfaceAbsorbed[static_cast<size_t>(c)] = incident - penetrating;

        double irradiance = penetrating;

        m_interface[static_cast<size_t>(mesh.interfaceSlot(c, 0))] = irradiance;

        for (int k = 0; k < layers; ++k)
        {
            const int64_t cell = mesh.cell(c, k);
            const double  h    = thickness[static_cast<size_t>(cell)];

            // Self-shading. Every species contributes whatever it was given a
            // coefficient for, so which constituents shade is a modeling
            // statement rather than a hard-coded list.
            double gamma = options.backgroundAttenuation;

            for (int s = 0; s < species; ++s)
            {
                const double alpha = m_attenuation[static_cast<size_t>(s)];

                if (alpha != 0.0)
                    gamma += alpha * state.speciesAt(s, cell);
            }

            const double tau = gamma * h;

            m_mean[static_cast<size_t>(cell)] = irradiance * meanAttenuationFactor(tau);

            irradiance *= std::exp(-tau);

            m_interface[static_cast<size_t>(mesh.interfaceSlot(c, k + 1))] = irradiance;
        }
    }
}

} // namespace FVQual::Kinetics
