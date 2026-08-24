// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   reactionnetwork.h
 * \author Caleb Buahin
 * \version 0.1.0-alpha.1
 * \brief  Multispecies reaction kinetics, operator-split from transport (F6.1b).
 *
 * A rate law per species, written as an expression over the other species and a
 * set of named constants, integrated per cell over each internal step. The
 * approach is CE-QUAL-W2's and EPANET-MSX's: the *network* is data, not code, so
 * adding a constituent or changing a decay law does not mean recompiling the
 * engine.
 *
 * **Temperature is not a special case.** It is a transported species like any
 * other, so a rate law may reference `T` — or any other species — simply by
 * name. Arrhenius temperature dependence is then something the modeller writes,
 * not something the engine hard-codes.
 *
 * **Splitting.** Transport runs, then reactions run on the transported field:
 * first-order Lie splitting. The error is `O(dt)` in the splitting itself,
 * which is consistent with the rest of the step, and the reaction sub-problem
 * is integrated adaptively so a stiff network does not dictate the hydrodynamic
 * time step.
 *
 * **What this class does *not* guarantee.** Everywhere else in the model,
 * conservation is an identity the engine enforces: transport moves constituent
 * mass without creating any, and an internal transfer cancels term by term. A
 * reaction network is different in kind. Whether mass is conserved depends
 * entirely on whether the rate laws sum to zero across species — `A -> B`
 * conserves, `-k*A` alone does not, and both are legitimate things to ask for.
 * The engine cannot decide which was meant, so it does not try. What it does
 * instead is make the consequence visible: ReactionDiagnostics reports the
 * smallest value any species reached, so a network that is driving
 * concentrations negative says so rather than quietly producing plausible
 * output.
 *
 * \threadsafety One network per solver. integrate() is internally parallel and
 *               gives each thread its own workspace.
 */

#ifndef FVQual_KINETICS_REACTIONNETWORK_H
#define FVQual_KINETICS_REACTIONNETWORK_H

#include "fvqual/fvqual_export.h"
#include "fvqual/kinetics/expression.h"
#include "fvqual/state/speciesregistry.h"
#include "fvqual/state/statestore.h"

#include <string>
#include <vector>

namespace FVQual::Kinetics {

/*!
 * \brief What the last integration did, beyond changing the state.
 */
struct FVQual_EXPORT ReactionDiagnostics
{
    /*!
     * \brief Smallest concentration any reacting species reached (mixed units).
     *
     * Negative means the network drove a species below zero. That is not
     * clamped away: a clamp would restore plausibility while destroying the
     * mass balance, and a modeller who cannot see it happening cannot fix the
     * rate law that caused it.
     */
    double minimumValue = 0.0;

    //! Cells where a species went negative.
    int64_t negativeCells = 0;

    //! Cells where the integrator failed to reach the requested accuracy.
    int64_t unconvergedCells = 0;

    //! Integrator calls made — one per cell with at least one rate law.
    int64_t integrations = 0;
};

/*!
 * \brief A set of rate laws over the registered species.
 */
class FVQual_EXPORT ReactionNetwork
{
public:
    /*!
     * \brief Registers a named constant usable by every rate law.
     *
     * \param[in]  name    Identifier as it appears in expressions.
     * \param[in]  value   Its value.
     * \param[out] message Failure description when returning false.
     */
    [[nodiscard]] bool addConstant(const std::string &name, double value,
                                   std::string &message);

    /*!
     * \brief Registers `dC/dt` for one species.
     *
     * The expression may reference any registered species by id and any
     * constant by name. A species with no rate law is transported and not
     * reacted, which is how a conservative tracer is expressed — by saying
     * nothing about it rather than by writing a zero.
     *
     * \param[in]  speciesId  Species the law applies to.
     * \param[in]  expression Rate law, in concentration units per second.
     * \param[out] message    Failure description when returning false.
     */
    [[nodiscard]] bool addRate(const std::string &speciesId,
                               const std::string &expression, std::string &message);

    /*!
     * \brief Resolves every expression's names against \p registry.
     *
     * Must be called before integrate(). Separating this from addRate() is
     * what allows a rate law to reference a species registered after it, and
     * it is where an unresolvable name is reported — at configuration, rather
     * than as a silent zero on the first step.
     *
     * \param[in]  registry Species registry the ids refer to.
     * \param[out] message  Failure description when returning false.
     */
    [[nodiscard]] bool compile(const State::SpeciesRegistry &registry,
                               std::string &message);

    /*!
     * \brief Advances every reacting species by \p dt, cell by cell.
     *
     * \param[in]     registry  Species registry.
     * \param[in,out] state     State store; reacting species are updated.
     * \param[in]     cellCount Cells to integrate.
     * \param[in]     dt        Interval (s).
     * \param[out]    message   Failure description when returning false.
     */
    [[nodiscard]] bool integrate(const State::SpeciesRegistry &registry,
                                 State::StateStore &state, int64_t cellCount, double dt,
                                 std::string &message);

    //! Number of species carrying a rate law.
    [[nodiscard]] int reactingCount() const noexcept
    {
        return static_cast<int>(m_rates.size());
    }

    //! Whether there is anything to integrate.
    [[nodiscard]] bool empty() const noexcept { return m_rates.empty(); }

    //! Whether compile() has succeeded since the last change.
    [[nodiscard]] bool compiled() const noexcept { return m_compiled; }

    [[nodiscard]] const ReactionDiagnostics &diagnostics() const noexcept
    {
        return m_diagnostics;
    }

private:
    /*!
     * \brief One rate law, after compilation.
     */
    struct Rate
    {
        std::string id;         //!< Species the law applies to.
        std::string text;       //!< The expression as written.
        Expression  expression; //!< Compiled form.
        int         species = -1; //!< Index in the registry, after compile().
    };

    std::vector<Rate>        m_rates;
    std::vector<std::string> m_constantNames;
    std::vector<double>      m_constantValues;

    //! Variable table handed to the evaluator: species, then constants.
    std::vector<std::string> m_variableNames;
    std::vector<const char *> m_variableTable;

    bool                m_compiled = false;
    ReactionDiagnostics m_diagnostics;

    //! Per-thread workspace: variable values, and the integrated subset.
    std::vector<double> m_values;
    std::vector<double> m_y;
};

} // namespace FVQual::Kinetics

#endif // FVQual_KINETICS_REACTIONNETWORK_H
