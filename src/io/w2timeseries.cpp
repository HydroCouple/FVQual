// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/io/w2timeseries.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>

namespace FVQual::IO {

namespace {

std::vector<std::string> tokens(const std::string &line)
{
    std::vector<std::string> out;
    std::istringstream       in(line);
    std::string              token;

    while (in >> token)
        out.push_back(token);

    return out;
}

std::string upper(std::string text)
{
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return text;
}

//! Whether every token parses as a number — how a data row is told from a header.
bool allNumeric(const std::vector<std::string> &fields)
{
    if (fields.empty())
        return false;

    for (const std::string &field : fields)
    {
        try
        {
            size_t used = 0;

            (void)std::stod(field, &used);

            if (used != field.size())
                return false;
        }
        catch (const std::exception &)
        {
            return false;
        }
    }

    return true;
}

const std::vector<double> kEmpty;

} // namespace

// ── Time series ───────────────────────────────────────────────────────────

bool W2TimeSeries::read(const std::string &path, std::string &message)
{
    return read(path, 0, message);
}

bool W2TimeSeries::read(const std::string &path, int valuesPerRecord,
                        std::string &message)
{
    std::ifstream input(path);

    if (!input)
    {
        message = "cannot open time-series file '" + path + "'";
        return false;
    }

    return parse(input, path, valuesPerRecord, message);
}

bool W2TimeSeries::parse(std::istream &input, const std::string &origin,
                         std::string &message)
{
    return parse(input, origin, 0, message);
}

bool W2TimeSeries::parse(std::istream &input, const std::string &origin,
                         int valuesPerRecord, std::string &message)
{
    std::vector<std::string>         names;
    std::vector<double>              time;
    std::vector<std::vector<double>> values;

    int  number = 0;
    bool seen   = false;

    std::string line;

    const auto fail = [&](const std::string &why) {
        message = origin + ":" + std::to_string(number) + ": " + why;
        return false;
    };

    while (std::getline(input, line))
    {
        ++number;

        const std::vector<std::string> fields = tokens(line);

        if (fields.empty())
            continue;

        if (!allNumeric(fields))
        {
            // A header. The last one before the data wins, which handles the
            // title line preceding it.
            if (seen)
                continue;

            names.clear();

            for (size_t i = 1; i < fields.size(); ++i)
                names.push_back(fields[i]);

            continue;
        }

        if (fields.size() < 2 && valuesPerRecord <= 0)
            return fail("a data row needs a Julian day and at least one value");

        // A wrapped record: the Julian day is on the first line only, and the
        // rest of its values follow on lines that look like bare numbers.
        std::vector<std::string> record = fields;

        if (valuesPerRecord > 0)
        {
            while (static_cast<int>(record.size()) - 1 < valuesPerRecord)
            {
                std::string more;

                if (!std::getline(input, more))
                    return fail("file ended part way through a record: wanted " +
                                std::to_string(valuesPerRecord) + " values, have " +
                                std::to_string(record.size() - 1));

                ++number;

                const std::vector<std::string> extra = tokens(more);

                if (extra.empty())
                    continue;

                if (!allNumeric(extra))
                    return fail("expected the rest of a record, found text");

                record.insert(record.end(), extra.begin(), extra.end());
            }

            if (static_cast<int>(record.size()) - 1 > valuesPerRecord)
                return fail("record has " + std::to_string(record.size() - 1) +
                            " values but " + std::to_string(valuesPerRecord) +
                            " were expected");
        }

        if (!seen)
        {
            seen = true;

            // A file whose header did not survive gets positional names, so a
            // caller can still address the columns.
            if (names.size() != record.size() - 1)
            {
                names.clear();

                for (size_t i = 1; i < record.size(); ++i)
                    names.push_back("COL" + std::to_string(i));
            }

            values.assign(names.size(), {});
        }

        if (record.size() - 1 != names.size())
        {
            return fail("row has " + std::to_string(record.size() - 1) +
                        " values but the header named " + std::to_string(names.size()));
        }

        const double day = std::stod(record.front());

        if (!time.empty() && day < time.back())
            return fail("Julian days must not go backwards");

        time.push_back(day);

        for (size_t i = 0; i < names.size(); ++i)
            values[i].push_back(std::stod(record[i + 1]));
    }

    if (time.empty())
        return fail("no data rows");

    m_names  = std::move(names);
    m_time   = std::move(time);
    m_values = std::move(values);

    message.clear();
    return true;
}

int W2TimeSeries::columnIndex(const std::string &name) const
{
    const std::string wanted = upper(name);

    for (size_t i = 0; i < m_names.size(); ++i)
    {
        if (upper(m_names[i]) == wanted)
            return static_cast<int>(i);
    }

    return -1;
}

double W2TimeSeries::valueAt(double julianDay, int column) const
{
    if (column < 0 || column >= columnCount() || m_time.empty())
        return 0.0;

    const std::vector<double> &series = m_values[static_cast<size_t>(column)];

    if (julianDay <= m_time.front())
        return series.front();

    if (julianDay >= m_time.back())
        return series.back();

    const auto upperBound = std::upper_bound(m_time.begin(), m_time.end(), julianDay);

    const size_t high = static_cast<size_t>(upperBound - m_time.begin());
    const size_t low  = high - 1;

    const double span = m_time[high] - m_time[low];

    if (span <= 0.0)
        return series[low];

    const double fraction = (julianDay - m_time[low]) / span;

    return series[low] + fraction * (series[high] - series[low]);
}

// ── Vertical profiles ─────────────────────────────────────────────────────

bool W2VerticalProfile::read(const std::string &path, std::string &message)
{
    std::ifstream input(path);

    if (!input)
    {
        message = "cannot open vertical-profile file '" + path + "'";
        return false;
    }

    return parse(input, path, message);
}

bool W2VerticalProfile::parse(std::istream &input, const std::string &origin,
                              std::string &message)
{
    std::vector<std::string>         names;
    std::vector<std::vector<double>> profiles;

    std::string line;
    bool        started = false;

    while (std::getline(input, line))
    {
        const std::vector<std::string> fields = tokens(line);

        if (fields.empty())
            continue;

        if (!allNumeric(fields))
        {
            // A section header: the constituent name, then repeated tags
            // (`T2I`, `C2I`) that carry no information here. The first
            // non-numeric line is the file title, which has no tags after it.
            if (fields.size() < 2)
                continue;

            const std::string tag = upper(fields[1]);

            if (tag != "T2I" && tag != "C2I")
                continue;

            names.push_back(fields.front());
            profiles.emplace_back();

            started = true;
            continue;
        }

        if (!started)
            continue;

        for (const std::string &field : fields)
            profiles.back().push_back(std::stod(field));
    }

    if (names.empty())
    {
        message = origin + ": no profiles found";
        return false;
    }

    m_names    = std::move(names);
    m_profiles = std::move(profiles);

    message.clear();
    return true;
}

const std::vector<double> &W2VerticalProfile::profile(const std::string &name) const
{
    const std::string wanted = upper(name);

    for (size_t i = 0; i < m_names.size(); ++i)
    {
        if (upper(m_names[i]) == wanted)
            return m_profiles[i];
    }

    return kEmpty;
}

} // namespace FVQual::IO
