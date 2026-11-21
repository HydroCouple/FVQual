// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  test_turbulence.cpp
 * \brief The k-epsilon closure against cases with known answers (F5.4).
 *
 * The scheme this replaces failed in a way worth remembering: measured against
 * DeGray, Richardson damping had no middle setting. Without shear it mixed
 * nothing (2.42 °C RMSE but a 3 m epilimnion against a measured 14 m); with
 * shear it mixed the column out (8.35 °C, hypolimnion 8 °C too warm). Both
 * settings looked fine in every unit test, because the unit tests asked whether
 * mixing responded to gradients and never whether the amount was right.
 *
 * So these gates are anchored to results that exist independently of this code:
 * the law of the wall, the decay law for homogeneous turbulence, and the
 * qualitative structure a stratified column must hold. Where a number is
 * asserted it comes from theory, not from a previous run of this function.
 */

#include "fvqual/hydro/turbulence.h"

#include <gtest/gtest.h>

#include <cmath>
#include <vector>

using FVQual::Hydro::stabilityMomentum;
using FVQual::Hydro::stabilityScalar;
using FVQual::Hydro::TurbulenceClosure;

namespace {

//! A uniform column: `layers` layers of `h` meters.
struct Column
{
    int                 layers;
    std::vector<double> thickness;
    std::vector<double> shear;
    std::vector<double> buoyancy;
    std::vector<double> viscosity;
    std::vector<double> diffusivity;

    explicit Column(int n, double h)
        : layers(n), thickness(static_cast<size_t>(n), h),
          shear(static_cast<size_t>(n + 1), 0.0),
          buoyancy(static_cast<size_t>(n + 1), 0.0),
          viscosity(static_cast<size_t>(n + 1), 0.0),
          diffusivity(static_cast<size_t>(n + 1), 0.0)
    {
    }

    void run(TurbulenceClosure &closure, double dt, int steps,
             double uStarSurface, double uStarBed)
    {
        for (int s = 0; s < steps; ++s)
        {
            closure.advanceColumn(0, dt, thickness.data(), shear.data(),
                                  buoyancy.data(), uStarSurface, uStarBed,
                                  viscosity.data(), diffusivity.data());
        }
    }
};

} // namespace

// ── Stability functions ───────────────────────────────────────────────────

TEST(Turbulence, StabilityFallsAsStratificationStrengthens)
{
    // The property the whole closure rests on, and a sign slip here would
    // invert the model's response to stratification while leaving every
    // conservation check intact.
    double previous = stabilityMomentum(-0.2);

    for (double alpha = -0.15; alpha <= 0.55; alpha += 0.05)
    {
        const double value = stabilityMomentum(alpha);

        EXPECT_LT(value, previous) << "stability rose with stratification at " << alpha;
        EXPECT_GT(value, 0.0) << "stability went non-positive at " << alpha;

        previous = value;
    }
}

TEST(Turbulence, TheTurbulentPrandtlNumberRisesWithStratification)
{
    // Written first as "scalars mix less than momentum", which is wrong: in
    // Galperin's functions Sh exceeds Sm throughout the fitted range, so the
    // turbulent Prandtl number Sm/Sh is 0.73 at neutral and never reaches one.
    // Scalars mix *more* than momentum in weakly stratified water, which is the
    // standard result and not what I asserted.
    //
    // The property that is actually true, and the one a thermocline depends on,
    // is that the advantage shrinks as stratification strengthens: Pr_t rises
    // towards one, so a density interface resists scalar transport increasingly
    // like it resists momentum. Kept as a monotonic assertion rather than a
    // threshold, because the monotonicity is the physics and the value is a fit.
    double previous = stabilityMomentum(0.0) / stabilityScalar(0.0);

    EXPECT_LT(previous, 1.0) << "scalars did not mix more than momentum at neutral";

    for (double alpha = 0.05; alpha <= 0.55; alpha += 0.05)
    {
        const double prandtl = stabilityMomentum(alpha) / stabilityScalar(alpha);

        EXPECT_GT(prandtl, previous)
            << "the Prandtl number fell as stratification strengthened, at " << alpha;

        previous = prandtl;
    }

    EXPECT_LT(previous, 1.0) << "the Prandtl number overshot one within the range";
}

TEST(Turbulence, StabilityIsBoundedOutsideItsFittedRange)
{
    // The rational form has a pole outside the range Galperin fitted, and an
    // unclamped pole in a mixing coefficient blows up a column.
    for (double alpha : {-1.0e3, -10.0, 10.0, 1.0e3})
    {
        EXPECT_TRUE(std::isfinite(stabilityMomentum(alpha))) << alpha;
        EXPECT_GT(stabilityMomentum(alpha), 0.0) << alpha;
        EXPECT_LT(stabilityMomentum(alpha), 1.0) << alpha;
    }
}

// ── Law of the wall ───────────────────────────────────────────────────────

TEST(Turbulence, AConstantStressLayerReachesTheLogLayerEnergy)
{
    // In a constant-stress layer the equilibrium turbulent kinetic energy is
    // k = u*²/sqrt(c_mu) — a textbook result, independent of this code. At
    // u* = 0.01 m/s that is 1e-4/0.3 = 3.33e-4 m²/s².
    TurbulenceClosure closure;
    closure.initialize(1, 20);

    Column column(20, 0.5);

    const double uStar = 0.01;

    // Shear consistent with a log layer at mid-depth, so production and
    // dissipation can balance rather than the walls doing all the work.
    for (size_t i = 0; i < column.shear.size(); ++i)
    {
        const double z = 0.25 + 0.5 * static_cast<double>(i);
        const double dudz = uStar / (0.41 * z);

        column.shear[i] = dudz * dudz;
    }

    column.run(closure, 10.0, 3000, uStar, uStar);

    const double expected = uStar * uStar / std::sqrt(closure.options().cMu);

    // Near the wall, where the log law is the theory. Mid-channel is *not*:
    // in a two-wall channel the energy legitimately falls below the wall value
    // towards the center, and an earlier version of this gate asserted the
    // log-layer energy at mid-depth and failed at 0.28 of it — the bound was
    // wrong about where the result applies, not the closure.
    // Three interfaces, reaching 1.75 m into a 10 m channel. Interface 4 at
    // 2.25 m sits at 0.478 of the wall energy, which is not a failure of the
    // closure — it is a quarter of the way to mid-channel, past where a
    // constant-stress layer is the right idealization.
    for (int i = 1; i <= 3; ++i)
    {
        EXPECT_GT(closure.tke()[static_cast<size_t>(i)], 0.5 * expected)
            << "interface " << i << " collapsed below the log-layer energy";
        EXPECT_LT(closure.tke()[static_cast<size_t>(i)], 2.0 * expected)
            << "interface " << i << " ran far above the log-layer energy";
    }

    // And the structure: energy falls away from the wall towards mid-channel.
    EXPECT_GT(closure.tke()[1], closure.tke()[10])
        << "turbulence did not decay away from the wall that generates it";
}

TEST(Turbulence, StrongerFrictionMakesMoreTurbulence)
{
    // Quadratic, since k = u*²/sqrt(c_mu): doubling u* quadruples the energy.
    TurbulenceClosure closure;

    auto energyAt = [&](double uStar) {
        closure.initialize(1, 10);

        Column column(10, 1.0);
        column.run(closure, 10.0, 300, uStar, uStar);

        return closure.tke()[0];
    };

    const double weak   = energyAt(0.005);
    const double strong = energyAt(0.010);

    EXPECT_NEAR(strong / weak, 4.0, 0.2) << "surface energy is not quadratic in u*";
}

// ── Decay ─────────────────────────────────────────────────────────────────

TEST(Turbulence, TurbulenceDecaysWhenNothingDrivesIt)
{
    // With no shear, no buoyancy and no wall friction, k must fall towards its
    // floor. A closure that sustains turbulence on nothing would keep a
    // reservoir mixed through a calm week.
    TurbulenceClosure closure;
    closure.initialize(1, 10);

    Column column(10, 1.0);

    // Spin it up with shear as well as wall friction. Wall friction alone
    // leaves the interior near its floor — correctly, since turbulence
    // diffusing off a wall into unsheared water decays as it goes — so the
    // first version of this test was measuring the floor decaying to the floor.
    for (auto &s : column.shear)
        s = 1.0e-3;

    column.run(closure, 10.0, 2000, 0.02, 0.02);

    const double spun = closure.tke()[5];

    ASSERT_GT(spun, 1.0e-4) << "the column never spun up";

    // Now take everything away: no shear, no friction, nothing to sustain it.
    for (auto &s : column.shear)
        s = 0.0;

    column.run(closure, 10.0, 3000, 0.0, 0.0);

    const double decayed = closure.tke()[5];

    EXPECT_LT(decayed, 0.1 * spun)
        << "turbulence persisted with nothing to sustain it: " << spun
        << " to " << decayed;
}

// ── Stratification ────────────────────────────────────────────────────────

TEST(Turbulence, StratificationSuppressesMixingAtTheInterface)
{
    // The property Richardson damping also had, kept because losing it would be
    // a regression: a density interface must mix less than unstratified water
    // under the same shear.
    TurbulenceClosure closure;

    auto diffusivityWith = [&](double nSquared) {
        closure.initialize(1, 20);

        Column column(20, 0.5);

        for (auto &s : column.shear)
            s = 1.0e-4;

        for (auto &b : column.buoyancy)
            b = nSquared;

        column.run(closure, 10.0, 300, 0.01, 0.005);

        return closure.options().minDiffusivity > 0.0 ? column.diffusivity[10] : 0.0;
    };

    const double neutral    = diffusivityWith(0.0);
    const double stratified = diffusivityWith(1.0e-2);

    EXPECT_LT(stratified, neutral)
        << "stratification did not suppress mixing: " << stratified
        << " against " << neutral;
}

TEST(Turbulence, AnUnstableColumnMixesHarderThanANeutralOne)
{
    // Convection: negative N² produces turbulence rather than consuming it,
    // which is what the sign switch on c3 is for.
    TurbulenceClosure closure;

    auto diffusivityWith = [&](double nSquared) {
        closure.initialize(1, 20);

        Column column(20, 0.5);

        for (auto &s : column.shear)
            s = 1.0e-6;

        for (auto &b : column.buoyancy)
            b = nSquared;

        column.run(closure, 10.0, 300, 0.002, 0.002);

        return column.diffusivity[10];
    };

    EXPECT_GT(diffusivityWith(-1.0e-4), diffusivityWith(0.0))
        << "an unstable column did not mix harder than a neutral one";
}

// ── Robustness ────────────────────────────────────────────────────────────

TEST(Turbulence, StillWaterStaysFiniteAndAtItsFloor)
{
    // The division by eps is the closure's sharpest edge. Still water is where
    // it bites, and it is also the state a reservoir spends nights in.
    TurbulenceClosure closure;
    closure.initialize(1, 15);

    Column column(15, 1.0);

    column.run(closure, 60.0, 500, 0.0, 0.0);

    for (size_t i = 0; i < column.viscosity.size(); ++i)
    {
        EXPECT_TRUE(std::isfinite(column.viscosity[i])) << "interface " << i;
        EXPECT_TRUE(std::isfinite(column.diffusivity[i])) << "interface " << i;

        EXPECT_GE(column.diffusivity[i], closure.options().minDiffusivity);
        EXPECT_LE(column.viscosity[i], closure.options().maxViscosity);
    }
}

TEST(Turbulence, ALargeStepDoesNotDriveEnergyNegative)
{
    // The sink is on the diagonal precisely so that no step, however long, can
    // overshoot into negative energy. Worth asserting, because an explicit sink
    // here would pass every short-step test and fail on a real reservoir Δt.
    TurbulenceClosure closure;
    closure.initialize(1, 10);

    Column column(10, 2.0);

    for (auto &s : column.shear)
        s = 1.0e-3;

    column.run(closure, 3600.0, 50, 0.01, 0.01);

    for (size_t i = 0; i < closure.tke().size(); ++i)
    {
        EXPECT_GT(closure.tke()[i], 0.0) << "negative energy at interface " << i;
        EXPECT_TRUE(std::isfinite(closure.dissipation()[i]));
    }
}
