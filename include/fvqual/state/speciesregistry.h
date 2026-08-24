// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   speciesregistry.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  SpeciesRegistry — the single naming authority for transported
 *         quantities.
 *
 * Everything the transport engines move is a species, including the two
 * reserved ones: temperature and water age. That is a deliberate choice
 * — heat is not a bolted-on special case with its own advection path,
 * it rides the same kernels as every other scalar, and so does age. One
 * transport path, one set of conservation tests, one place for bugs to hide.
 *
 * The naming convention and species kinds mirror the OpenSWMM reaction
 * engine's registry so a `.rxn` file written for one reads in the other
 * as a transported species.
 *
 * \threadsafety Immutable after setup; concurrent reads are safe.
 */

#ifndef FVQual_STATE_SPECIESREGISTRY_H
#define FVQual_STATE_SPECIESREGISTRY_H

#include "fvqual/fvqual_export.h"

#include <string>
#include <vector>

namespace FVQual::State {

/*!
 * \brief Role a species plays in the model.
 */
enum class SpeciesKind
{
    Pollutant,   //!< Ordinary transported constituent.
    Temperature, //!< Reserved: water temperature (deg C).
    Age,         //!< Reserved: water age (days).
    MsxBulk,     //!< Multispecies bulk species (water column).
    MsxBenthic   //!< Multispecies benthic species (bed-attached, not advected).
};

/*!
 * \brief One species definition.
 */
struct FVQual_EXPORT SpeciesDefinition
{
    std::string id;                       //!< Unique identifier.
    SpeciesKind kind = SpeciesKind::Pollutant;
    std::string units;                    //!< SI units.
    std::string standardName;             //!< CF standard_name, when one exists.
    double      initialValue = 0.0;       //!< Uniform initial condition.

    /*!
     * \brief Whether the transport engines advect this species. Benthic
     *        species live on the bed and are not carried by the flow.
     */
    [[nodiscard]] bool isTransported() const noexcept
    {
        return kind != SpeciesKind::MsxBenthic;
    }
};

/*!
 * \brief Ordered registry of species; the index is the species' identity in
 *        every SoA block.
 */
class FVQual_EXPORT SpeciesRegistry
{
public:
    /*!
     * \brief Reserved id of the temperature species.
     */
    static constexpr const char *kTemperature = "__TEMPERATURE__";

    /*!
     * \brief Reserved id of the water-age species.
     */
    static constexpr const char *kAge = "__AGE__";

    /*!
     * \brief Registers a species.
     * \param[in]  definition Species to add.
     * \param[out] message    Failure description when returning -1.
     * \return The new species index, or -1 when the id is empty or duplicated.
     */
    int add(const SpeciesDefinition &definition, std::string &message);

    /*!
     * \brief Registers temperature with the conventional units and initial
     *        value.
     * \param[in] initialCelsius Uniform initial temperature (deg C).
     */
    int addTemperature(double initialCelsius, std::string &message);

    /*!
     * \brief Registers water age (days), initially zero everywhere.
     */
    int addAge(std::string &message);

    /*!
     * \brief Index of \p id, or -1 when it is not registered.
     */
    [[nodiscard]] int indexOf(const std::string &id) const noexcept;

    [[nodiscard]] int count() const noexcept
    {
        return static_cast<int>(m_definitions.size());
    }

    [[nodiscard]] const SpeciesDefinition &definition(int index) const noexcept
    {
        return m_definitions[static_cast<size_t>(index)];
    }

    [[nodiscard]] const std::vector<SpeciesDefinition> &definitions() const noexcept
    {
        return m_definitions;
    }

    /*!
     * \brief Index of the reserved temperature species, or -1.
     */
    [[nodiscard]] int temperatureIndex() const noexcept { return m_temperatureIndex; }

    /*!
     * \brief Index of the reserved age species, or -1.
     */
    [[nodiscard]] int ageIndex() const noexcept { return m_ageIndex; }

private:
    std::vector<SpeciesDefinition> m_definitions;
    int                            m_temperatureIndex = -1;
    int                            m_ageIndex         = -1;
};

} // namespace FVQual::State

#endif // FVQual_STATE_SPECIESREGISTRY_H
