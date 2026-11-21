// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   reactionfile.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  Reading a reaction network from a text file (F6.3).
 *
 * The reaction network is data rather than code so that changing a constituent
 * model does not mean recompiling the engine. That is only half true while the
 * only way to build one is a sequence of C++ calls — this is the other half.
 *
 * **Format.** Sections in brackets, one record per line, `;` to end-of-line for
 * comments. Section names and keywords are case-insensitive; species ids are
 * not, because they have to match the registry exactly.
 *
 * ```
 * [SPECIES]
 * ;  id       units   initial   [BENTHIC]
 *    DO       g/m3    9.0
 *    CBOD     g/m3    5.0
 *    TSS_BED  g/m2    0.0       BENTHIC
 *
 * [CONSTANTS]
 *    kd       0.23/86400            ; per second
 *
 * [RATES]
 * ;  species  scope     expression
 *    DO       WATER     0 - kd * CBOD
 *    DO       SURFACE   2.3e-5 * (9.0 - DO)
 *    DO       BED       0 - 1.16e-5
 *
 * [SETTLING]
 * ;  species  velocity        [receiver]
 *    TSS      1.0/86400       TSS_BED
 *
 * [LIGHT]
 * ;  species  specific attenuation (m2/g)
 *    ALGAE    0.02
 * ```
 *
 * **Everything is SI.** Seconds, meters, grams. That is a deliberate refusal to
 * be clever: a rate law is arbitrary text over constants whose meaning the
 * parser cannot know, so it could convert some numbers and not others, which
 * would be worse than converting none.
 *
 * **But numeric fields are expressions.** `0.23/86400` is how a modeler thinks
 * about a decay coefficient quoted per day, and it stays readable in the file
 * where `2.6620370370e-6` does not. Every numeric field here is evaluated by the
 * same parser the rate laws use, so the arithmetic is visible rather than
 * pre-digested. This is the one place the file can undo the per-day trap without
 * guessing.
 *
 * **Numeric fields carry no spaces.** `0.23/86400`, not `0.23 / 86400`. The rule
 * exists because a record like a settling entry has an optional field *after*
 * the number, and whitespace is the only delimiter — allowing spaces inside the
 * number would make `TSS 1.0 / 86400 TSS_BED` ambiguous, which it duly was until
 * a test caught it. The rate expression in `[RATES]` is the sole exception: it
 * is always the last field on its line, so it may be written out in full.
 *
 * **The scope on a rate law is mandatory.** It could not be optional: an
 * omitted field in front of a free-form expression is unparseable, since
 * `DO   0 - k * CBOD` is indistinguishable from a law scoped to something
 * called `0`. Requiring it also makes a model file state where every process
 * acts, which is worth the extra word.
 *
 * **Two phases, as with the templates**, because species must exist before the
 * state store is allocated and rate laws can only be compiled after the solver
 * exists.
 */

#ifndef FVQual_KINETICS_REACTIONFILE_H
#define FVQual_KINETICS_REACTIONFILE_H

#include "fvqual/fvqual_export.h"
#include "fvqual/kinetics/reactionnetwork.h"
#include "fvqual/state/speciesregistry.h"

#include <iosfwd>
#include <string>
#include <vector>

namespace FVQual::Hydro {
class InternalModeSolver;
}

namespace FVQual::Kinetics {

/*!
 * \brief A reaction network as read from text, before it touches a model.
 *
 * Parsing and applying are separate so that a malformed file is rejected
 * without having half-configured anything.
 */
class FVQual_EXPORT ReactionFile
{
public:
    /*!
     * \brief Reads and parses \p path.
     *
     * \param[out] message Failure description, with a line number, when
     *                     returning false.
     */
    [[nodiscard]] bool read(const std::string &path, std::string &message);

    /*!
     * \brief Parses from an already-open stream.
     *
     * \param[in]  input   Stream to read to end.
     * \param[in]  origin  Name used in error messages; a path, or anything.
     * \param[out] message Failure description, with a line number.
     */
    [[nodiscard]] bool parse(std::istream &input, const std::string &origin,
                             std::string &message);

    //! Registers every species the file declares.
    [[nodiscard]] bool registerSpecies(State::SpeciesRegistry &registry,
                                       std::string &message) const;

    //! Installs constants, rate laws, settling and light attenuation.
    [[nodiscard]] bool configure(Hydro::InternalModeSolver &solver,
                                 std::string &message) const;

    //! Number of species declared.
    [[nodiscard]] int speciesCount() const noexcept
    {
        return static_cast<int>(m_species.size());
    }

    //! Number of rate laws declared, across all scopes.
    [[nodiscard]] int rateCount() const noexcept
    {
        return static_cast<int>(m_rates.size());
    }

    //! Number of constants declared.
    [[nodiscard]] int constantCount() const noexcept
    {
        return static_cast<int>(m_constants.size());
    }

private:
    struct Species
    {
        std::string id;
        std::string units;
        double      initialValue = 0.0;
        bool        benthic      = false;
    };

    struct Constant
    {
        std::string name;
        double      value = 0.0;
    };

    struct Rate
    {
        std::string id;
        std::string expression;
        RateScope   scope = RateScope::WaterColumn;
    };

    struct Settling
    {
        std::string id;
        double      velocity = 0.0;
        std::string receiver;
    };

    struct Attenuation
    {
        std::string id;
        double      value = 0.0;
    };

    std::vector<Species>     m_species;
    std::vector<Constant>    m_constants;
    std::vector<Rate>        m_rates;
    std::vector<Settling>    m_settling;
    std::vector<Attenuation> m_light;
};

} // namespace FVQual::Kinetics

#endif // FVQual_KINETICS_REACTIONFILE_H
