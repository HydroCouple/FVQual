// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/kinetics/reactionnetwork.h"

#include "fvqual/kinetics/odesolver.h"
#include "fvqual/parallel.h"

#include <algorithm>
#include <limits>

namespace FVQual::Kinetics {

namespace {

//! Relative accuracy asked of the reaction integrator.
constexpr double kTolerance = 1.0e-8;

} // namespace

bool ReactionNetwork::addConstant(const std::string &name, double value,
                                  std::string &message)
{
    if (name.empty())
    {
        message = "a constant needs a name";
        return false;
    }

    if (std::find(m_constantNames.begin(), m_constantNames.end(), name) !=
        m_constantNames.end())
    {
        message = "constant '" + name + "' is already defined";
        return false;
    }

    m_constantNames.push_back(name);
    m_constantValues.push_back(value);
    m_compiled = false;

    message.clear();
    return true;
}

bool ReactionNetwork::addRate(const std::string &speciesId,
                              const std::string &expression, std::string &message)
{
    if (speciesId.empty())
    {
        message = "a rate law needs a species";
        return false;
    }

    for (const Rate &rate : m_rates)
    {
        if (rate.id == speciesId)
        {
            // Two laws for one species is ambiguous rather than additive: the
            // modeller either meant to sum them into one expression or made a
            // mistake, and guessing which would be worse than refusing.
            message = "species '" + speciesId + "' already has a rate law";
            return false;
        }
    }

    Rate rate;
    rate.id   = speciesId;
    rate.text = expression;

    if (parse(expression, rate.expression) != 0 || !rate.expression.valid)
    {
        message = "rate law for '" + speciesId + "' is not a valid expression: " +
                  expression;
        return false;
    }

    m_rates.push_back(std::move(rate));
    m_compiled = false;

    message.clear();
    return true;
}

bool ReactionNetwork::compile(const State::SpeciesRegistry &registry,
                              std::string &message)
{
    const int species = registry.count();

    // The variable table every expression is bound against: species first, in
    // registry order, then the constants. Species come first so that a bound
    // variable index below `species` is a concentration and can be refreshed
    // from the integrator's state vector without a lookup.
    m_variableNames.clear();
    m_variableNames.reserve(static_cast<size_t>(species) + DerivedCount +
                            m_constantNames.size());

    for (int s = 0; s < species; ++s)
        m_variableNames.push_back(registry.definition(s).id);

    // Then the engine-computed fields, always present in the table whether or
    // not they are used — their positions have to be fixed for the integration
    // loop to fill them without a lookup.
    const char *const derivedNames[DerivedCount] = {kLight, kDepth, kThickness};

    for (const char *name : derivedNames)
    {
        if (std::find(m_variableNames.begin(), m_variableNames.end(), name) !=
            m_variableNames.end())
        {
            message = std::string("species '") + name +
                      "' collides with a reserved derived-field name";
            return false;
        }

        m_variableNames.push_back(name);
    }

    for (const std::string &name : m_constantNames)
    {
        if (std::find(m_variableNames.begin(), m_variableNames.end(), name) !=
            m_variableNames.end())
        {
            message = "constant '" + name + "' collides with a species of the same name";
            return false;
        }

        m_variableNames.push_back(name);
    }

    m_variableTable.clear();
    m_variableTable.reserve(m_variableNames.size());

    for (const std::string &name : m_variableNames)
        m_variableTable.push_back(name.c_str());

    const int variables = static_cast<int>(m_variableTable.size());

    for (bool &used : m_derivedUsed)
        used = false;

    for (Rate &rate : m_rates)
    {
        rate.species = registry.indexOf(rate.id);

        if (rate.species < 0)
        {
            message = "rate law refers to species '" + rate.id +
                      "', which is not registered";
            return false;
        }

        // Every name in the expression must resolve. An unbound name would
        // otherwise evaluate as zero, turning a typo into a rate law that runs
        // and is silently wrong — the failure mode this whole check exists to
        // prevent.
        const int bound = bind_variables(rate.expression, m_variableTable.data(),
                                         variables);

        int used = 0;

        for (const Token &token : rate.expression.postfix)
        {
            if (token.type == TokenType::VARIABLE)
                ++used;
        }

        if (bound != used)
        {
            message = "rate law for '" + rate.id + "' refers to " +
                      std::to_string(used - bound) +
                      " name(s) that are neither a species, a derived field, nor a "
                      "constant: " + rate.text;
            return false;
        }

        // Record which derived fields this law reads, so integrate() can refuse
        // when one of them was not supplied.
        for (const Token &token : rate.expression.postfix)
        {
            if (token.type != TokenType::VARIABLE)
                continue;

            const int offset = token.var_idx - species;

            if (offset >= 0 && offset < DerivedCount)
                m_derivedUsed[offset] = true;
        }
    }

    // Workspace, one slice per thread.
    const size_t threads = static_cast<size_t>(maxThreads());

    m_values.assign(static_cast<size_t>(variables) * threads, 0.0);
    m_y.assign(m_rates.size() * threads, 0.0);

    m_compiled = true;

    message.clear();
    return true;
}

bool ReactionNetwork::usesLight() const noexcept
{
    return m_derivedUsed[DerivedLight];
}

bool ReactionNetwork::integrate(const State::SpeciesRegistry &registry,
                                State::StateStore &state, int64_t cellCount, double dt,
                                std::string &message, const DerivedFields &fields)
{
    if (m_rates.empty())
    {
        m_diagnostics = ReactionDiagnostics{};
        message.clear();
        return true;
    }

    if (!m_compiled)
    {
        message = "the reaction network must be compiled before it is integrated";
        return false;
    }

    if (!(dt > 0.0))
    {
        message = "reaction interval must be positive";
        return false;
    }

    const int species   = registry.count();
    const int variables = static_cast<int>(m_variableTable.size());
    const int reacting  = static_cast<int>(m_rates.size());

    // A referenced field that was not supplied is refused here. Ablating this
    // check does not produce a quiet zero as I first assumed — it segfaults,
    // because the integration loop indexes the span directly. So this is a
    // bounds guard first and a diagnostic second.
    const std::span<const double> derived[DerivedCount] = {fields.light, fields.depth,
                                                           fields.thickness};

    const char *const derivedNames[DerivedCount] = {kLight, kDepth, kThickness};

    for (int d = 0; d < DerivedCount; ++d)
    {
        if (!m_derivedUsed[d])
            continue;

        // One check, two messages. An earlier version tested emptiness
        // separately, which an ablation showed to be dead code: an empty span
        // has size zero and so already fails the length test. The distinction
        // that survives is only in what the caller is told.
        if (static_cast<int64_t>(derived[d].size()) < cellCount)
        {
            message = derived[d].empty()
                          ? std::string("a rate law references ") + derivedNames[d] +
                                " but that field was not supplied"
                          : std::string("the ") + derivedNames[d] + " field covers " +
                                std::to_string(derived[d].size()) + " cells but " +
                                std::to_string(cellCount) + " were asked for";
            return false;
        }
    }

    // The workspace was sized at compile() time from the then-current thread
    // count. A host that raises it afterwards — omp_set_num_threads() from an
    // embedding application or a Python binding — would otherwise have threads
    // indexing past the end of these buffers, which corrupts memory silently
    // rather than failing.
    const size_t threads = static_cast<size_t>(maxThreads());

    if (m_values.size() < static_cast<size_t>(variables) * threads)
        m_values.assign(static_cast<size_t>(variables) * threads, 0.0);

    if (m_y.size() < m_rates.size() * threads)
        m_y.assign(m_rates.size() * threads, 0.0);

    double  minimum    = std::numeric_limits<double>::max();
    int64_t negatives  = 0;
    int64_t unconverged = 0;

    FVQual_PARALLEL_FOR_CLAUSES(reduction(min : minimum)
                                reduction(+ : negatives)
                                reduction(+ : unconverged))
    for (int64_t cell = 0; cell < cellCount; ++cell)
    {
        const size_t slot = static_cast<size_t>(threadIndex());

        double *values = m_values.data() + slot * static_cast<size_t>(variables);
        double *y      = m_y.data() + slot * static_cast<size_t>(reacting);

        // Every species is visible to every rate law, whether or not it reacts:
        // a decay rate that depends on temperature is the normal case, not a
        // special one.
        for (int s = 0; s < species; ++s)
            values[s] = state.speciesAt(s, cell);

        for (int d = 0; d < DerivedCount; ++d)
        {
            values[species + d] =
                m_derivedUsed[d] ? derived[d][static_cast<size_t>(cell)] : 0.0;
        }

        for (int c = 0; c < static_cast<int>(m_constantValues.size()); ++c)
        {
            values[species + DerivedCount + c] =
                m_constantValues[static_cast<size_t>(c)];
        }

        for (int r = 0; r < reacting; ++r)
            y[r] = state.speciesAt(m_rates[static_cast<size_t>(r)].species, cell);

        // The reacting species are refreshed from the integrator's state on
        // every derivative evaluation; the rest stay frozen through the
        // sub-step, which is what makes this a split operator rather than a
        // second transport.
        const auto derivatives = [&](double, const double *current, double *rates) {
            for (int r = 0; r < reacting; ++r)
                values[m_rates[static_cast<size_t>(r)].species] = current[r];

            for (int r = 0; r < reacting; ++r)
                rates[r] = evaluate_fast(m_rates[static_cast<size_t>(r)].expression,
                                         values);
        };

        // The first trial step is the whole interval; the adaptive controller
        // shrinks it as far as a stiff network needs, so a fast reaction costs
        // accuracy in the reaction sub-problem rather than in the whole model's
        // time step.
        const int status =
            Kinetics::integrate(y, reacting, 0.0, dt, kTolerance, dt, derivatives);

        if (status != ODE_OK)
            ++unconverged;

        bool wentNegative = false;

        for (int r = 0; r < reacting; ++r)
        {
            state.speciesAt(m_rates[static_cast<size_t>(r)].species, cell) = y[r];

            minimum = std::min(minimum, y[r]);

            if (y[r] < 0.0)
                wentNegative = true;
        }

        if (wentNegative)
            ++negatives;
    }

    m_diagnostics.minimumValue     = minimum;
    m_diagnostics.negativeCells    = negatives;
    m_diagnostics.unconvergedCells = unconverged;
    m_diagnostics.integrations     = cellCount;

    message.clear();
    return true;
}

} // namespace FVQual::Kinetics
