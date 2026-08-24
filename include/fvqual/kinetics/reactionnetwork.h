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

#include <span>
#include <string>
#include <vector>

namespace FVQual::Kinetics {

/*!
 * \brief Engine-computed per-cell fields a rate law may reference by name.
 *
 * The third category of quantity in an expression, after species and
 * constants: things the model knows about a cell that are not transported.
 * Light is the reason this exists — algal growth cannot be written without it,
 * and it cannot be a species because it is not conserved and not advected.
 *
 * Each span is either empty or exactly `cellCount` long. Empty means "not
 * supplied", and a rate law that references a field which was not supplied is
 * refused rather than reading a zero — a zero light field is a model with no
 * photosynthesis, which runs and looks plausible.
 */
struct FVQual_EXPORT DerivedFields
{
    std::span<const double> light;     //!< Layer-mean irradiance (W/m²).
    std::span<const double> depth;     //!< Depth of the cell centre below the surface (m).
    std::span<const double> thickness; //!< Cell thickness (m).
};

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
 * \brief Where in the column a rate law acts.
 */
enum class RateScope
{
    /*!
     * \brief Everywhere in the water column, per unit volume.
     *
     * Units are concentration per second. The ordinary case: decay,
     * nitrification, algal growth.
     */
    WaterColumn,

    /*!
     * \brief Only in the cell touching the bed, per unit **area**.
     *
     * Units are mass per unit area per second, and the engine divides by the
     * bed cell's thickness to turn that into a concentration rate. Sediment
     * oxygen demand is quoted as a flux — g/m²/d — because that is what it is:
     * a property of the sediment surface, not of the water above it.
     *
     * Doing the division in the engine rather than asking the modeller to write
     * `/ __THICKNESS__` is deliberate. Forgetting it produces a model where the
     * benthic demand silently scales with the layer count, and the value looks
     * plausible either way.
     */
    Bed,

    /*!
     * \brief Only in the cell touching the atmosphere, per unit **area**.
     *
     * The same areal treatment as Bed, at the other end of the column, and for
     * the same reason: gas exchange is a property of the water surface. It is
     * areal without needing to be declared so — reaeration is
     * `K_L (C_sat - C)`, and a transfer velocity in m/s times a concentration
     * in g/m³ is already g/m²/s.
     */
    Surface
};

/*!
 * \brief Where a cell sits in its column, as a bit field.
 *
 * A one-layer column is both surface and bed, which is not a degenerate case
 * to be excluded: it is a shallow pond, where reaeration and sediment demand
 * genuinely act on the same water.
 */
enum CellPosition : unsigned char
{
    CellInterior  = 0,
    CellAtBed     = 1 << 0,
    CellAtSurface = 1 << 1
};

/*!
 * \brief A set of rate laws over the registered species.
 */
class FVQual_EXPORT ReactionNetwork
{
public:
    //! Reserved name for layer-mean irradiance (W/m²).
    static constexpr const char *kLight = "__LIGHT__";

    //! Reserved name for the depth of the cell centre below the surface (m).
    static constexpr const char *kDepth = "__DEPTH__";

    //! Reserved name for cell thickness (m).
    static constexpr const char *kThickness = "__THICKNESS__";

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
     * \brief Registers `dC/dt` for one species over part of the column.
     *
     * A species may carry one law per scope, and where they overlap — a bed cell
     * is also part of the water column — the rates **add** into a single
     * derivative rather than being applied in sequence. That distinction
     * matters: sequential application would introduce a splitting error between
     * two terms of the same ODE, where none is needed.
     *
     * Two laws for the same species *and* scope are still refused. Summing them
     * would be well defined, but writing the same law twice is far more often a
     * mistake than an intent, and the sum is available by writing one
     * expression.
     *
     * \param[in]  speciesId  Species the law applies to.
     * \param[in]  expression Rate law; units depend on \p scope.
     * \param[in]  scope      Where the law acts.
     * \param[out] message    Failure description when returning false.
     */
    [[nodiscard]] bool addRate(const std::string &speciesId,
                               const std::string &expression, RateScope scope,
                               std::string &message);

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
     * \param[in]     fields    Engine-computed per-cell fields. Required only
     *                          for the ones the rate laws actually reference.
     * \param[in]     position  CellPosition flags per cell. Required only when
     *                          a scoped (Bed or Surface) law exists.
     */
    [[nodiscard]] bool integrate(const State::SpeciesRegistry &registry,
                                 State::StateStore &state, int64_t cellCount, double dt,
                                 std::string &message,
                                 const DerivedFields &fields = {},
                                 std::span<const unsigned char> position = {});

    /*!
     * \brief Whether any rate law references the irradiance field.
     *
     * Lets a caller skip computing light for a network that does not use it.
     */
    [[nodiscard]] bool usesLight() const noexcept;

    //! Whether any rate law is scoped, and so needs the cell-position flags.
    [[nodiscard]] bool usesCellPosition() const noexcept;

    //! Number of species carrying at least one rate law.
    [[nodiscard]] int reactingCount() const noexcept
    {
        return static_cast<int>(m_reacting.size());
    }

    //! Number of rate laws, which may exceed reactingCount() across scopes.
    [[nodiscard]] int rateCount() const noexcept
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
        int         slot    = -1; //!< Index into m_reacting, after compile().
        RateScope   scope   = RateScope::WaterColumn;
    };

    [[nodiscard]] bool addRateImpl(const std::string &speciesId,
                                   const std::string &expression, RateScope scope,
                                   std::string &message);

    //! The derived fields, in the order they occupy the variable table.
    enum Derived : int
    {
        DerivedLight = 0,
        DerivedDepth,
        DerivedThickness,
        DerivedCount
    };

    std::vector<Rate>        m_rates;
    std::vector<std::string> m_constantNames;
    std::vector<double>      m_constantValues;

    /*!
     * \brief Registry indices of the species carrying at least one law.
     *
     * The integrated state vector is indexed by position in *this* list, not by
     * rate law, so that a species with both a water-column and a bed law is
     * integrated once with the two rates summed.
     */
    std::vector<int> m_reacting;

    //! Which derived fields the compiled rate laws actually reference.
    bool m_derivedUsed[DerivedCount] = {false, false, false};

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
