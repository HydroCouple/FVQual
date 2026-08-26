// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   forcingfield.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  ForcingField and MeteorologyField — forcing that is uniform or
 *         per column, chosen per variable (F3.7).
 *
 * Meteorology was spatially uniform on the argument that a reservoir is small
 * against a weather system. That is true of the weather and false of the water:
 * shading, fetch and sheltering vary along a reservoir, a narrow tree-lined arm
 * does not see the wind the main pool does, and a coupled atmospheric model
 * delivers fields rather than scalars.
 *
 * **The uniform case stays a single value.** A `ForcingField` holds either one
 * double or a vector, so the common configuration costs one number and the
 * arithmetic is identical to what the scalar path did — not merely close, but
 * the same operations on the same values, which is what lets the existing heat
 * budget gates keep their meaning.
 *
 * **Granularity is the surface column, not the cell.** Wind, air temperature
 * and cloud act at the air-water interface, so a per-cell value for them has no
 * physical referent. Where depth resolution does mean something it already
 * exists: `LightField` attenuates per cell with self-shading.
 *
 * **`surfaceHeatFlux` is deliberately untouched.** It remains a pure function
 * of one `Meteorology` and one temperature. Only the supplier changes, and
 * `MeteorologyField::at()` materialises the struct the function already takes.
 *
 * \threadsafety `at()` and the accessors are const and race-free; the setters
 *               are not, and are meant to be called between steps.
 */

#ifndef FVQual_FORCING_FORCINGFIELD_H
#define FVQual_FORCING_FORCINGFIELD_H

#include "fvqual/fvqual_export.h"
#include "fvqual/hydro/heatbudget.h"

#include <cstdint>
#include <string>
#include <vector>

namespace FVQual::Forcing {

/*!
 * \brief One forcing variable, either uniform over the domain or given per
 *        column.
 */
class FVQual_EXPORT ForcingField
{
public:
    ForcingField() = default;

    //! Constructs a uniform field.
    explicit ForcingField(double uniform) noexcept
        : m_uniform(uniform)
    {
    }

    /*!
     * \brief Makes the field uniform, discarding any per-column values.
     */
    void setUniform(double value) noexcept
    {
        m_uniform = value;
        m_values.clear();
    }

    /*!
     * \brief Gives the field a value per column.
     *
     * Refused, with the field left untouched, unless \p values has exactly
     * \p columns entries — a silently recycled or truncated array is a forcing
     * error that looks like a physics error later.
     *
     * \param[in]  values  One value per column.
     * \param[in]  columns Column count the mesh expects.
     * \param[out] message Failure description when returning false.
     */
    [[nodiscard]] bool setPerColumn(std::vector<double> values, int64_t columns,
                                    std::string &message);

    /*!
     * \brief The value at \p column.
     *
     * No bounds check: the field is validated against the mesh's column count
     * when it is set, and every caller loops over that same count. A check here
     * would sit in the innermost loop of the heat budget and never fire.
     */
    [[nodiscard]] double at(int64_t column) const noexcept
    {
        return m_values.empty() ? m_uniform
                                : m_values[static_cast<size_t>(column)];
    }

    //! Whether the field is a single value rather than a per-column array.
    [[nodiscard]] bool isUniform() const noexcept { return m_values.empty(); }

    //! The uniform value; meaningless unless isUniform().
    [[nodiscard]] double uniformValue() const noexcept { return m_uniform; }

    //! The per-column values; empty when uniform.
    [[nodiscard]] const std::vector<double> &values() const noexcept
    {
        return m_values;
    }

private:
    double              m_uniform = 0.0;
    std::vector<double> m_values; //!< Empty means uniform.
};

/*!
 * \brief Meteorological forcing, each variable uniform or per column.
 *
 * The spatially varying variables are those that describe the weather at a
 * place. The rest of `Meteorology` — emissivities, the wind-function
 * coefficients, the psychrometric constant, the drag coefficient — are
 * calibration constants for the *formulation*, not observations of a location,
 * and are carried once in \ref base.
 */
class FVQual_EXPORT MeteorologyField
{
public:
    /*!
     * \brief Non-spatial coefficients, and the template `at()` starts from.
     *
     * Its spatially varying members are ignored; the fields below supply those.
     */
    Hydro::Meteorology base;

    ForcingField airTemperature{20.0};    //!< Dry-bulb air temperature (deg C).
    ForcingField dewPoint{10.0};          //!< Dew point (deg C).
    ForcingField windSpeed{2.0};          //!< Wind speed (m/s).
    ForcingField windDirection{0.0};      //!< Wind direction (radians; see Meteorology).
    ForcingField cloudCover{0.0};         //!< Sky covered, 0 to 1.
    ForcingField shortwave{0.0};          //!< Incident shortwave (W/m²).
    ForcingField pressure{101325.0};      //!< Atmospheric pressure (Pa).
    ForcingField albedo{0.06};            //!< Shortwave reflected at the surface.

    // surfaceAbsorption is deliberately absent. W2's BETA is an optical
    // property of the water and its spectrum, not an observation of the weather
    // at a place, and LightField already owns a copy of it. Making it vary per
    // column without plumbing the same variation into the light field would
    // leave the surface split and the penetrating profile disagreeing about the
    // same number — a silent inconsistency of exactly the kind this module
    // exists to stop. It stays in `base`, where both consumers read one value.
    // Promoting it later is additive.

    /*!
     * \brief Sets every variable uniform from one `Meteorology`, coefficients
     *        included.
     *
     * The migration path for code that used to mutate a single struct, and the
     * configuration most applications want.
     */
    void setUniform(const Hydro::Meteorology &weather);

    /*!
     * \brief Assigning a `Meteorology` makes the whole field uniform.
     *
     * Exists so that code written against the old scalar member — including
     * every existing heat-budget gate — keeps compiling and keeps meaning what
     * it meant. A migration that forces its own tests to be rewritten has
     * discarded the evidence that the behaviour did not change.
     */
    MeteorologyField &operator=(const Hydro::Meteorology &weather)
    {
        setUniform(weather);
        return *this;
    }

    /*!
     * \brief The weather at \p column.
     *
     * Returns by value because `surfaceHeatFlux` takes a `Meteorology`, and
     * assembling one costs less than restructuring a function whose gates are
     * the reason the heat budget is trusted at all.
     */
    [[nodiscard]] Hydro::Meteorology at(int64_t column) const;

    //! Whether every variable is uniform.
    [[nodiscard]] bool isUniform() const noexcept;

    /*!
     * \brief Fraction of variables carrying per-column values, for reporting.
     */
    [[nodiscard]] int perColumnCount() const noexcept;
};

} // namespace FVQual::Forcing

#endif // FVQual_FORCING_FORCINGFIELD_H
