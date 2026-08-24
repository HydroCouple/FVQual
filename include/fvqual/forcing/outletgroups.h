// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   outletgroups.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  Multi-port outlet operation to a release target (F4.2).
 *
 * A dam with one outlet releases whatever that outlet draws. A dam with ports
 * at several elevations can *choose*, and choosing is the reason selective
 * withdrawal was implemented at all: a deep port draws cold hypolimnetic water,
 * a surface port draws warm epilimnetic water, and blending the two hits a
 * downstream temperature that neither reaches alone. Release temperature is a
 * regulated quantity below many dams, so this is the operational lever the
 * whole structures phase exists to provide.
 *
 * **What is decided here, and what is not.** This class decides *how the demand
 * is split between ports*. It does not decide the demand itself — total release
 * comes from an operating schedule the caller supplies, because it is set by
 * downstream water rights, power generation, and flood control, none of which
 * the reservoir model knows about. Nor does it decide the capacity: each port's
 * rating curve still limits what its head can pass, and a demand beyond the
 * combined capacity is reported as a shortfall rather than met by inventing
 * water (see deliveredRelease()).
 *
 * **The mixing law.** With ports delivering temperatures `T_w` and `T_c` that
 * bracket the target `T*`, conservation of heat in the blended release gives
 *
 *     f_w = (T* - T_c) / (T_w - T_c),
 *
 * the fraction drawn from the warm port. When no pair brackets the target the
 * best achievable release is the single closest port, and the achieved
 * temperature then falls short of the target — correctly, because no
 * combination of outlets can produce water the reservoir does not contain.
 *
 * \threadsafety One set per solver; update() mutates the recorded assignment.
 */

#ifndef FVQual_FORCING_OUTLETGROUPS_H
#define FVQual_FORCING_OUTLETGROUPS_H

#include "fvqual/forcing/structures.h"
#include "fvqual/fvqual_export.h"
#include "fvqual/mesh/layeredmesh.h"

#include <span>
#include <string>
#include <vector>

namespace FVQual::Forcing {

/*!
 * \brief A set of ports operated together to meet one release target.
 */
struct FVQual_EXPORT OutletGroup
{
    std::string id;

    /*!
     * \brief Indices, in the paired StructureSet, of the ports to operate.
     *
     * They need not share a column: a dam with outlet works split across two
     * segments is still operated to one target.
     */
    std::vector<int> ports;

    /*!
     * \brief Total release demanded of the group (m³/s, magnitude).
     */
    double demand = 0.0;

    /*!
     * \brief Target temperature of the blended release (deg C).
     */
    double targetTemperature = 0.0;

    /*!
     * \brief Operate to \c targetTemperature. When false the demand is split
     *        between the ports in proportion to their capacities, which is what
     *        a group with no temperature objective does.
     */
    bool blend = true;
};

/*!
 * \brief The outlet groups acting on a model.
 */
class FVQual_EXPORT OutletGroupSet
{
public:
    /*!
     * \brief Registers a group.
     *
     * \param[in]  group      The group.
     * \param[in]  structures The structure set its ports index into.
     * \param[out] message    Failure description when returning -1.
     * \return Its index, or -1 when the specification is invalid.
     */
    int add(const OutletGroup &group, const StructureSet &structures,
            std::string &message);

    [[nodiscard]] int count() const noexcept
    {
        return static_cast<int>(m_groups.size());
    }

    [[nodiscard]] const OutletGroup &group(int index) const noexcept
    {
        return m_groups[static_cast<size_t>(index)];
    }

    /*!
     * \brief Sets the total release demanded of a group (m³/s).
     */
    [[nodiscard]] bool setDemand(int index, double demand, std::string &message);

    /*!
     * \brief Sets the target release temperature of a group (deg C).
     */
    [[nodiscard]] bool setTarget(int index, double temperature, std::string &message);

    /*!
     * \brief Chooses the port split and writes it into the structures.
     *
     * Called once per internal step, before the structures are rated. Once per
     * *internal* step is the right cadence and the subcycle is not: the choice
     * is driven by the temperature field, which transport advances once a step,
     * so re-deciding it inside the barotropic subcycle would repeat the whole
     * envelope calculation for every port against a temperature field that has
     * not moved. Capacity limiting still happens each substep, where the head
     * does move — StructureSet::update() applies it.
     *
     * \param[in]     mesh             The layered mesh.
     * \param[in]     sources          Sources the structures drive.
     * \param[in]     density          Cell densities (kg/m³).
     * \param[in]     temperature      Cell temperatures (deg C).
     * \param[in]     gravity          Gravitational acceleration (m/s²).
     * \param[in]     referenceDensity Boussinesq reference density (kg/m³).
     * \param[in,out] structures       Structures to schedule.
     * \param[out]    message          Failure description when returning false.
     */
    [[nodiscard]] bool update(const Mesh::LayeredMesh &mesh, const SourceSet &sources,
                              std::span<const double> density,
                              std::span<const double> temperature, double gravity,
                              double referenceDensity, StructureSet &structures,
                              std::string &message);

    /*!
     * \brief Release the group was able to schedule (m³/s).
     *
     * Less than the demand when the ports' combined capacity cannot meet it.
     * Exposed because a model that quietly released less than it was asked to
     * would otherwise look like a mass-balance error somewhere else entirely.
     */
    [[nodiscard]] double deliveredRelease(int index) const noexcept
    {
        return m_delivered[static_cast<size_t>(index)];
    }

    /*!
     * \brief Flow-weighted temperature of the scheduled release (deg C).
     *
     * The target when it is achievable, and the closest the reservoir can come
     * when it is not.
     */
    [[nodiscard]] double achievedTemperature(int index) const noexcept
    {
        return m_achieved[static_cast<size_t>(index)];
    }

private:
    std::vector<OutletGroup> m_groups;
    std::vector<double>      m_delivered;
    std::vector<double>      m_achieved;

    // Scratch, sized to the largest group.
    std::vector<double> m_portTemperature;
    std::vector<double> m_portCapacity;
    std::vector<double> m_portRelease;
};

} // namespace FVQual::Forcing

#endif // FVQual_FORCING_OUTLETGROUPS_H
