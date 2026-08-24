// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/kinetics/reactionfile.h"

#include "fvqual/hydro/internalmode.h"
#include "fvqual/kinetics/expression.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>

namespace FVQual::Kinetics {

namespace {

//! Sections the reader understands.
enum class Section
{
    None,
    Species,
    Constants,
    Rates,
    Settling,
    Light
};

std::string upper(std::string text)
{
    std::transform(text.begin(), text.end(), text.begin(),
                   [](unsigned char c) { return static_cast<char>(std::toupper(c)); });
    return text;
}

std::string trim(const std::string &text)
{
    const auto first = text.find_first_not_of(" \t\r\n");

    if (first == std::string::npos)
        return {};

    const auto last = text.find_last_not_of(" \t\r\n");

    return text.substr(first, last - first + 1);
}

//! Strips a trailing comment. `;` follows the SWMM and EPANET convention.
std::string stripComment(const std::string &line)
{
    const auto mark = line.find(';');

    return mark == std::string::npos ? line : line.substr(0, mark);
}

//! Splits on whitespace into tokens, with no field allowed to contain spaces.
std::vector<std::string> split(const std::string &line)
{
    std::vector<std::string> out;
    std::istringstream       in(line);
    std::string              token;

    while (in >> token)
        out.push_back(token);

    return out;
}

//! Splits on whitespace into at most \p limit fields, the last keeping the rest.
std::vector<std::string> fields(const std::string &line, size_t limit)
{
    std::vector<std::string> out;

    size_t i = 0;

    while (i < line.size())
    {
        while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i])))
            ++i;

        if (i >= line.size())
            break;

        if (out.size() + 1 == limit)
        {
            // The final field takes everything remaining, so a rate expression
            // keeps its spaces.
            out.push_back(trim(line.substr(i)));
            break;
        }

        const size_t start = i;

        while (i < line.size() && !std::isspace(static_cast<unsigned char>(line[i])))
            ++i;

        out.push_back(line.substr(start, i - start));
    }

    return out;
}

/*!
 * \brief Evaluates a numeric field, which may be an arithmetic expression.
 *
 * `0.23 / 86400` rather than `2.662037e-6`. The same parser the rate laws use,
 * so there is no second numeric syntax to learn or to get subtly wrong.
 */
bool evaluateNumber(const std::string &text, double &value, std::string &why)
{
    Expression expression;

    if (parse(text, expression) != 0 || !expression.valid)
    {
        why = "'" + text + "' is not a number or arithmetic expression";
        return false;
    }

    // No variables are bound, so any name in it is unresolvable — and would
    // evaluate as zero, which is exactly the silent failure the network's own
    // binding check exists to prevent.
    for (const Token &token : expression.postfix)
    {
        if (token.type == TokenType::VARIABLE)
        {
            why = "'" + text + "' refers to '" + token.var_name +
                  "', but a numeric field cannot use variables";
            return false;
        }
    }

    value = evaluate(expression, nullptr, 0);

    return true;
}

} // namespace

bool ReactionFile::read(const std::string &path, std::string &message)
{
    std::ifstream input(path);

    if (!input)
    {
        message = "cannot open reaction file '" + path + "'";
        return false;
    }

    return parse(input, path, message);
}

bool ReactionFile::parse(std::istream &input, const std::string &origin,
                         std::string &message)
{
    // Parsed into locals and committed only on success, so a file rejected
    // half-way through leaves the object exactly as it was. Clearing the
    // members up front would have left whatever parsed before the error, which
    // is the worst of both: a failed call and a partly-configured model.
    std::vector<Species>     species;
    std::vector<Constant>    constants;
    std::vector<Rate>        rates;
    std::vector<Settling>    settlingEntries;
    std::vector<Attenuation> attenuation;

    Section section = Section::None;

    std::string line;
    int         number = 0;

    // Every failure below reports where it happened. A configuration file that
    // says only "parse error" makes the user bisect their own file by hand.
    const auto fail = [&](const std::string &why) {
        message = origin + ":" + std::to_string(number) + ": " + why;
        return false;
    };

    while (std::getline(input, line))
    {
        ++number;

        const std::string content = trim(stripComment(line));

        if (content.empty())
            continue;

        if (content.front() == '[')
        {
            if (content.back() != ']')
                return fail("section header '" + content + "' is not closed");

            const std::string name = upper(trim(content.substr(1, content.size() - 2)));

            if (name == "SPECIES")
                section = Section::Species;
            else if (name == "CONSTANTS")
                section = Section::Constants;
            else if (name == "RATES")
                section = Section::Rates;
            else if (name == "SETTLING")
                section = Section::Settling;
            else if (name == "LIGHT")
                section = Section::Light;
            else
                return fail("unknown section '" + name + "'");

            continue;
        }

        if (section == Section::None)
            return fail("content before any section header");

        std::string why;

        switch (section)
        {
            case Section::Species:
            {
                const auto parts = split(content);

                if (parts.size() < 3)
                    return fail("a species needs an id, units and initial value");

                if (parts.size() > 4)
                    return fail("a species takes at most four fields");

                Species entry;
                entry.id    = parts[0];
                entry.units = parts[1];

                if (!evaluateNumber(parts[2], entry.initialValue, why))
                    return fail(why);

                if (parts.size() == 4)
                {
                    if (upper(parts[3]) != "BENTHIC")
                        return fail("'" + parts[3] + "' is not a species qualifier; "
                                    "expected BENTHIC");

                    entry.benthic = true;
                }

                const auto clash = std::find_if(
                    species.begin(), species.end(),
                    [&](const Species &s) { return s.id == entry.id; });

                if (clash != species.end())
                    return fail("species '" + entry.id + "' is declared twice");

                species.push_back(std::move(entry));
                break;
            }

            case Section::Constants:
            {
                const auto parts = split(content);

                if (parts.size() != 2)
                    return fail("a constant needs a name and a value, with no "
                                "spaces inside the value");

                Constant constant;
                constant.name = parts[0];

                if (!evaluateNumber(parts[1], constant.value, why))
                    return fail(why);

                constants.push_back(std::move(constant));
                break;
            }

            case Section::Rates:
            {
                const auto parts = fields(content, 3);

                if (parts.size() < 3)
                {
                    return fail("a rate law needs a species, a scope and an "
                                "expression");
                }

                Rate rate;
                rate.id = parts[0];

                // The scope is mandatory, not defaulted. It cannot be optional:
                // an omitted field before a free-form expression is
                // unparseable, because `DO   0 - k * CBOD` is indistinguishable
                // from a scope named `0`. Making it explicit also makes a model
                // file say where each process acts, which is worth the words.
                const std::string scope = upper(parts[1]);

                if (scope == "WATER")
                    rate.scope = RateScope::WaterColumn;
                else if (scope == "SURFACE")
                    rate.scope = RateScope::Surface;
                else if (scope == "BED")
                    rate.scope = RateScope::Bed;
                else
                    return fail("'" + parts[1] +
                                "' is not a scope; expected WATER, SURFACE or BED");

                rate.expression = parts[2];

                if (rate.expression.empty())
                    return fail("rate law for '" + rate.id + "' is empty");

                rates.push_back(std::move(rate));
                break;
            }

            case Section::Settling:
            {
                const auto parts = split(content);

                if (parts.size() < 2 || parts.size() > 3)
                    return fail("a settling entry needs a species and a velocity, "
                                "and optionally a benthic receiver");

                Settling entry;
                entry.id = parts[0];

                if (!evaluateNumber(parts[1], entry.velocity, why))
                    return fail(why);

                if (parts.size() == 3)
                    entry.receiver = parts[2];

                settlingEntries.push_back(std::move(entry));
                break;
            }

            case Section::Light:
            {
                const auto parts = split(content);

                if (parts.size() != 2)
                    return fail("a light entry needs a species and an attenuation");

                Attenuation entry;
                entry.id = parts[0];

                if (!evaluateNumber(parts[1], entry.value, why))
                    return fail(why);

                attenuation.push_back(std::move(entry));
                break;
            }

            case Section::None:
                break;
        }
    }

    m_species   = std::move(species);
    m_constants = std::move(constants);
    m_rates     = std::move(rates);
    m_settling  = std::move(settlingEntries);
    m_light     = std::move(attenuation);

    message.clear();
    return true;
}

bool ReactionFile::registerSpecies(State::SpeciesRegistry &registry,
                                   std::string &message) const
{
    for (const Species &species : m_species)
    {
        State::SpeciesDefinition definition;
        definition.id           = species.id;
        definition.units        = species.units;
        definition.initialValue = species.initialValue;

        if (species.benthic)
            definition.kind = State::SpeciesKind::MsxBenthic;

        if (registry.add(definition, message) < 0)
            return false;
    }

    message.clear();
    return true;
}

bool ReactionFile::configure(Hydro::InternalModeSolver &solver,
                             std::string &message) const
{
    ReactionNetwork &network = solver.reactions();

    for (const Constant &constant : m_constants)
    {
        if (!network.addConstant(constant.name, constant.value, message))
            return false;
    }

    for (const Rate &rate : m_rates)
    {
        if (!network.addRate(rate.id, rate.expression, rate.scope, message))
            return false;
    }

    for (const Settling &entry : m_settling)
    {
        const int species = solver.speciesIndex(entry.id);

        if (species < 0)
        {
            message = "settling refers to species '" + entry.id +
                      "', which is not registered";
            return false;
        }

        if (!solver.settling().setVelocity(species, entry.velocity, solver.species(),
                                           message))
        {
            return false;
        }

        if (entry.receiver.empty())
            continue;

        const int receiver = solver.speciesIndex(entry.receiver);

        if (receiver < 0)
        {
            message = "settling refers to benthic species '" + entry.receiver +
                      "', which is not registered";
            return false;
        }

        if (!solver.settling().setBenthicReceiver(species, receiver, solver.species(),
                                                  message))
        {
            return false;
        }
    }

    for (const Attenuation &entry : m_light)
    {
        const int species = solver.speciesIndex(entry.id);

        if (species < 0)
        {
            message = "light attenuation refers to species '" + entry.id +
                      "', which is not registered";
            return false;
        }

        if (!solver.light().setSpecificAttenuation(species, entry.value, message))
            return false;
    }

    message.clear();
    return true;
}

} // namespace FVQual::Kinetics
