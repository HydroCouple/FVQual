// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   w2timeseries.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  CE-QUAL-W2 time-series and vertical-profile inputs (F8.1).
 *
 * The forcing half of a W2 model: meteorology, inflows, inflow temperatures,
 * outflows, and the initial vertical profile. The format is the same throughout
 * — a title line, a blank, a header naming the columns, then rows whose first
 * column is the Julian day.
 *
 * **Interpolation is linear and clamped.** A met record at three-hour intervals
 * driving a model on a thirty-second step has to be interpolated, and clamping
 * outside the record is what keeps a run that starts slightly before the first
 * sample from extrapolating air temperature into fiction.
 */

#ifndef FVQual_IO_W2TIMESERIES_H
#define FVQual_IO_W2TIMESERIES_H

#include "fvqual/fvqual_export.h"

#include <iosfwd>
#include <string>
#include <vector>

namespace FVQual::IO {

/*!
 * \brief A W2 time-series file: Julian day plus one or more data columns.
 */
class FVQual_EXPORT W2TimeSeries
{
public:
    [[nodiscard]] bool read(const std::string &path, std::string &message);

    /*!
     * \brief Reads a file whose records wrap across several lines.
     *
     * **Why this is needed.** W2's fixed-format files put one value per segment,
     * and a reservoir with more segments than fit a line simply continues onto
     * the next one, with no Julian day to mark it. DeGray's `wsc.npt` carries 31
     * wind-sheltering values per record across four lines, and its header names
     * nine — so neither the line nor the header says how wide a record is. The
     * caller does, from the bathymetry.
     *
     * With \p valuesPerRecord at 0 this is exactly \ref read().
     *
     * \param[in]  path            File to read.
     * \param[in]  valuesPerRecord Values after the Julian day, or 0 to take the
     *                             count from the first data line.
     * \param[out] message         Failure description, with a line number.
     */
    [[nodiscard]] bool read(const std::string &path, int valuesPerRecord,
                            std::string &message);

    //! As above, from an open stream.
    [[nodiscard]] bool parse(std::istream &input, const std::string &origin,
                             int valuesPerRecord, std::string &message);

    [[nodiscard]] bool parse(std::istream &input, const std::string &origin,
                             std::string &message);

    //! Data columns, excluding the Julian-day column.
    [[nodiscard]] int columnCount() const noexcept
    {
        return static_cast<int>(m_names.size());
    }

    [[nodiscard]] int rowCount() const noexcept
    {
        return static_cast<int>(m_time.size());
    }

    //! Column names as the header gave them, excluding `JDAY`.
    [[nodiscard]] const std::vector<std::string> &columnNames() const noexcept
    {
        return m_names;
    }

    //! Julian days, ascending.
    [[nodiscard]] const std::vector<double> &time() const noexcept { return m_time; }

    //! Index of \p name, or -1. Case-insensitive.
    [[nodiscard]] int columnIndex(const std::string &name) const;

    /*!
     * \brief Linearly interpolated value of \p column at \p julianDay.
     *
     * Clamped to the first and last samples rather than extrapolated.
     */
    [[nodiscard]] double valueAt(double julianDay, int column) const;

    //! The record's first and last Julian days.
    [[nodiscard]] double firstDay() const noexcept
    {
        return m_time.empty() ? 0.0 : m_time.front();
    }

    [[nodiscard]] double lastDay() const noexcept
    {
        return m_time.empty() ? 0.0 : m_time.back();
    }

private:
    std::vector<std::string>         m_names;
    std::vector<double>              m_time;
    std::vector<std::vector<double>> m_values; //!< [column][row].
};

/*!
 * \brief A W2 vertical-profile file: named profiles, one value per wetted layer.
 *
 * The number of values is itself a check on the vertical datum. DeGray's
 * profile carries 28 entries, which is exactly the number of layers wetted at
 * its initial surface of 123.80 m given `EBOT = 66.35 m` — an independent
 * confirmation of an elevation the bathymetry file never states.
 */
class FVQual_EXPORT W2VerticalProfile
{
public:
    [[nodiscard]] bool read(const std::string &path, std::string &message);

    [[nodiscard]] bool parse(std::istream &input, const std::string &origin,
                             std::string &message);

    //! Names of the profiles the file carries, in order.
    [[nodiscard]] const std::vector<std::string> &names() const noexcept
    {
        return m_names;
    }

    /*!
     * \brief Values of the named profile, surface first; empty if absent.
     */
    [[nodiscard]] const std::vector<double> &profile(const std::string &name) const;

private:
    std::vector<std::string>         m_names;
    std::vector<std::vector<double>> m_profiles;
};

} // namespace FVQual::IO

#endif // FVQual_IO_W2TIMESERIES_H
