// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_kinetics_numerics.cpp
 * \brief Verification of the expression evaluator and ODE integrator (F6.1a).
 *
 * These two files were carried over from openswmm.engine, where they are
 * exercised by that project's suite. Ported code is not tested code: the
 * namespace changed, the build changed, and the compiler and flags changed, so
 * what is verified here is *this* copy in *this* build. A port that compiles is
 * evidence of nothing.
 *
 * The gates are chosen so that a wrong answer is impossible to mistake for a
 * right one. Expressions are checked against values computed independently in
 * the test, and the integrator against closed-form solutions of the ODEs it is
 * given — exponential decay, a stiff pair with a known analytical solution, and
 * a conservative reaction where the invariant is exact rather than approximate.
 */

#include "fvqual/kinetics/expression.h"
#include "fvqual/kinetics/odesolver.h"

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

namespace Kinetics = FVQual::Kinetics;

namespace {

/*!
 * \brief Parses \p text or fails the test with the offending expression.
 */
Kinetics::Expression compile(const std::string &text)
{
    Kinetics::Expression expression;

    const int status = Kinetics::parse(text, expression);

    EXPECT_EQ(status, 0) << "failed to parse: " << text;
    EXPECT_TRUE(expression.valid) << text;

    return expression;
}

} // namespace

// ── Expression evaluation ─────────────────────────────────────────────────

TEST(KineticsExpression, ArithmeticMatchesTheSameFormulaInCpp)
{
    // Precedence and associativity, checked against the compiler's own
    // evaluation of the identical formula rather than against a number typed
    // by hand — a typo in the expected value would otherwise read as a parser
    // bug, or worse, hide one.
    struct Case
    {
        const char *text;
        double      expected;
    };

    const Case cases[] = {
        {"2 + 3 * 4", 2 + 3 * 4},
        {"(2 + 3) * 4", (2 + 3) * 4},
        {"2 - 3 - 4", 2 - 3 - 4},              // left-associative
        {"2 ^ 3 ^ 2", std::pow(2.0, 9.0)},     // right-associative
        {"10 / 4", 10.0 / 4.0},                // real division, not integer
        {"1 + 2 * 3 - 4 / 2", 1 + 2 * 3 - 4.0 / 2.0},
    };

    for (const Case &c : cases)
    {
        const Kinetics::Expression expression = compile(c.text);

        EXPECT_NEAR(Kinetics::evaluate(expression, nullptr, 0), c.expected, 1.0e-12)
            << c.text;
    }
}

TEST(KineticsExpression, UnaryMinusBindsTighterThanExponentiation)
{
    // Pinned because it is a real fork in convention, not because it is
    // obviously right. This parser reads "-3 ^ 2" as (-3)^2 = 9, as Excel
    // does. Mathematical notation and Python read it as -(3^2) = -9.
    //
    // The behavior is inherited from the shared openswmm lineage, so changing
    // it here would silently give a different answer to a reaction expression
    // that produces the other one there. It is recorded rather than altered,
    // and anyone writing a rate law is better served by parenthesizing.
    const Kinetics::Expression expression = compile("-3 ^ 2");

    EXPECT_NEAR(Kinetics::evaluate(expression, nullptr, 0), 9.0, 1.0e-12)
        << "unary-minus precedence changed; shared expressions would diverge";

    // Written explicitly, both readings are available and unambiguous.
    EXPECT_NEAR(Kinetics::evaluate(compile("(0 - 3) ^ 2"), nullptr, 0), 9.0, 1.0e-12);
    EXPECT_NEAR(Kinetics::evaluate(compile("0 - 3 ^ 2"), nullptr, 0), -9.0, 1.0e-12);
}

TEST(KineticsExpression, FunctionsMatchTheStandardLibrary)
{
    const double x = 0.7;

    struct Case
    {
        const char *text;
        double      expected;
    };

    const Case cases[] = {
        {"exp(0.7)", std::exp(x)},   {"log(0.7)", std::log(x)},
        {"sqrt(0.7)", std::sqrt(x)}, {"abs(0 - 0.7)", std::abs(-x)},
        {"sin(0.7)", std::sin(x)},   {"cos(0.7)", std::cos(x)},
        {"tan(0.7)", std::tan(x)},   {"tanh(0.7)", std::tanh(x)},
        {"min(0.7, 2)", 0.7},        {"max(0.7, 2)", 2.0},
    };

    for (const Case &c : cases)
    {
        const Kinetics::Expression expression = compile(c.text);

        EXPECT_NEAR(Kinetics::evaluate(expression, nullptr, 0), c.expected, 1.0e-12)
            << c.text;
    }
}

TEST(KineticsExpression, VariablesBindByNameAndEvaluateFromTheArray)
{
    // The path a reaction rate actually takes: names resolved once at
    // configuration, then evaluated per cell against a flat array.
    Kinetics::Expression expression = compile("K * C * exp(-0.05 * (T - 20))");

    const char *const names[] = {"K", "C", "T"};

    // Returns the number of variables bound, not a status code.
    ASSERT_EQ(Kinetics::bind_variables(expression, names, 3), 3);

    const double values[] = {0.25, 4.0, 25.0};

    const double expected = 0.25 * 4.0 * std::exp(-0.05 * (25.0 - 20.0));

    EXPECT_NEAR(Kinetics::evaluate(expression, values, 3), expected, 1.0e-12);

    // The bound fast path must agree with the checked one, or the optimization
    // is silently computing something else.
    EXPECT_NEAR(Kinetics::evaluate_fast(expression, values), expected, 1.0e-12);

    // And it must track the values, not a cached result.
    const double changed[] = {0.25, 8.0, 25.0};

    EXPECT_NEAR(Kinetics::evaluate_fast(expression, changed), 2.0 * expected, 1.0e-12);
}

TEST(KineticsExpression, MalformedExpressionsAreRejected)
{
    // A rate expression is user input. Silently returning zero for a
    // mistyped formula would produce a model that runs and is wrong.
    for (const char *text : {"2 +", "(2 + 3", "2 + 3)", "", "exp(", "* 3"})
    {
        Kinetics::Expression expression;

        const int status = Kinetics::parse(text, expression);

        EXPECT_TRUE(status != 0 || !expression.valid)
            << "accepted malformed expression: '" << text << "'";
    }
}

// ── ODE integration ───────────────────────────────────────────────────────

TEST(KineticsOde, ExponentialDecayMatchesItsClosedForm)
{
    // dC/dt = -k C,  C(t) = C0 exp(-k t). The workhorse of first-order
    // kinetics, and the case where an integrator's error is unambiguous.
    constexpr double k  = 0.35;
    constexpr double c0 = 12.0;
    constexpr double t1 = 6.0;

    double y = c0;

    const auto derivatives = [](double, const double *state, double *rates) {
        rates[0] = -k * state[0];
    };

    const int status = Kinetics::integrate(&y, 1, 0.0, t1, 1.0e-9, 0.1, derivatives);

    ASSERT_EQ(status, Kinetics::ODE_OK);

    const double exact = c0 * std::exp(-k * t1);

    EXPECT_NEAR(y, exact, 1.0e-7 * exact) << "integrated " << y << " against " << exact;
}

TEST(KineticsOde, ACoupledPairMatchesItsAnalyticalSolution)
{
    // A -> B -> C, the classic sequential decay. B rises then falls, and its
    // peak is where a sloppy integrator shows itself:
    //
    //   A(t) = A0 exp(-k1 t)
    //   B(t) = A0 k1 (exp(-k1 t) - exp(-k2 t)) / (k2 - k1)
    constexpr double k1 = 0.8;
    constexpr double k2 = 0.3;
    constexpr double a0 = 5.0;

    for (double t : {0.5, 1.5, 4.0, 10.0})
    {
        double y[3] = {a0, 0.0, 0.0};

        const auto derivatives = [](double, const double *s, double *r) {
            r[0] = -k1 * s[0];
            r[1] = k1 * s[0] - k2 * s[1];
            r[2] = k2 * s[1];
        };

        ASSERT_EQ(Kinetics::integrate(y, 3, 0.0, t, 1.0e-10, 0.05, derivatives),
                  Kinetics::ODE_OK);

        const double exactA = a0 * std::exp(-k1 * t);
        const double exactB = a0 * k1 * (std::exp(-k1 * t) - std::exp(-k2 * t)) / (k2 - k1);

        EXPECT_NEAR(y[0], exactA, 1.0e-7 * a0) << "A at t = " << t;
        EXPECT_NEAR(y[1], exactB, 1.0e-7 * a0) << "B at t = " << t;

        // Nothing leaves the system, so the total is an identity rather than
        // an approximation — and it is the property a reaction network has to
        // keep if constituent budgets are to close.
        EXPECT_NEAR(y[0] + y[1] + y[2], a0, 1.0e-8 * a0) << "mass at t = " << t;
    }
}

TEST(KineticsOde, AStiffPairIsIntegratedWithoutBlowingUp)
{
    // Reaction networks are stiff by nature: fast equilibria alongside slow
    // decay. An adaptive step must shrink through the transient and open up
    // afterwards rather than either exploding or crawling to the step limit.
    constexpr double fast = 500.0;
    constexpr double slow = 0.1;

    double y[2] = {1.0, 0.0};

    const auto derivatives = [](double, const double *s, double *r) {
        r[0] = -fast * s[0];
        r[1] = fast * s[0] - slow * s[1];
    };

    const int status = Kinetics::integrate(y, 2, 0.0, 5.0, 1.0e-8, 1.0e-4, derivatives);

    ASSERT_EQ(status, Kinetics::ODE_OK) << "integrator gave up on a stiff pair";

    // The fast species is long gone. The slow one is *fed* by it, so it is not
    // a bare exponential: the sequential-decay solution
    // y1 = f/(s-f) (exp(-f t) - exp(-s t)) leaves a factor of f/(f-s) = 1.0002
    // that a bare exp(-s t) would miss by more than the tolerance. The
    // integrator had this right and the first version of this expectation did
    // not.
    const double exactSlow =
        fast / (slow - fast) * (std::exp(-fast * 5.0) - std::exp(-slow * 5.0));

    EXPECT_LT(y[0], 1.0e-12) << "fast species did not relax";
    EXPECT_NEAR(y[1], exactSlow, 1.0e-6) << "slow species reached " << y[1];

    EXPECT_GE(y[0], 0.0) << "integrator produced a negative concentration";
    EXPECT_GE(y[1], 0.0) << "integrator produced a negative concentration";
}

TEST(KineticsOde, AnExpressionDrivenRateIntegratesCorrectly)
{
    // The two pieces together, which is how the kinetics engine will use them:
    // a parsed rate expression evaluated inside the derivative callback.
    Kinetics::Expression expression = compile("0 - K * C");

    const char *const names[] = {"K", "C"};
    ASSERT_EQ(Kinetics::bind_variables(expression, names, 2), 2);

    constexpr double k  = 0.42;
    constexpr double c0 = 3.0;
    constexpr double t1 = 4.0;

    double y = c0;

    const auto derivatives = [&expression](double, const double *state, double *rates) {
        const double values[] = {k, state[0]};
        rates[0] = Kinetics::evaluate_fast(expression, values);
    };

    ASSERT_EQ(Kinetics::integrate(&y, 1, 0.0, t1, 1.0e-9, 0.1, derivatives),
              Kinetics::ODE_OK);

    const double exact = c0 * std::exp(-k * t1);

    EXPECT_NEAR(y, exact, 1.0e-7 * exact);
}
