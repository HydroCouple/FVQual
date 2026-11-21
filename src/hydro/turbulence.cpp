// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

#include "fvqual/hydro/turbulence.h"

#include "fvqual/parallel.h"

#include <algorithm>
#include <cmath>

namespace FVQual::Hydro {

namespace {

//! Galperin's bounds on the dimensionless stratification.
constexpr double kAlphaMin = -0.28;
constexpr double kAlphaMax = 0.56;

/*!
 * \brief Thomas algorithm over `n` unknowns.
 *
 * The same forward-elimination/back-substitution the momentum and scalar column
 * solves use; kept local rather than shared because the closure's system is
 * assembled differently enough that a common helper would take more parameters
 * than it saves.
 */
void solveTridiagonal(int n, const double *lower, const double *diagonal,
                      const double *upper, const double *rhs, double *sweep,
                      double *solution)
{
    double beta = diagonal[0];

    solution[0] = rhs[0] / beta;

    for (int i = 1; i < n; ++i)
    {
        sweep[i] = upper[i - 1] / beta;
        beta     = diagonal[i] - lower[i] * sweep[i];
        solution[i] = (rhs[i] - lower[i] * solution[i - 1]) / beta;
    }

    for (int i = n - 2; i >= 0; --i)
        solution[i] -= sweep[i + 1] * solution[i + 1];
}

} // namespace

double stabilityMomentum(double alphaN) noexcept
{
    // Galperin et al. (1988) quasi-equilibrium. Clamped because the rational
    // form has a pole outside the range it was fitted over, and an unclamped
    // pole in a mixing coefficient is a reliable way to blow up a column.
    const double a = std::clamp(alphaN, kAlphaMin, kAlphaMax);

    return 0.0906 / (1.0 + 1.732 * a);
}

double stabilityScalar(double alphaN) noexcept
{
    const double a = std::clamp(alphaN, kAlphaMin, kAlphaMax);

    return 0.1243 / (1.0 + 2.681 * a);
}

void TurbulenceClosure::initialize(int64_t columns, int layers)
{
    m_columns = columns;
    m_layers  = layers;

    const size_t interfaces = static_cast<size_t>(columns) *
                              static_cast<size_t>(layers + 1);

    m_tke.assign(interfaces, m_options.minTke);
    m_dissipation.assign(interfaces, m_options.minDissipation);

    // One slice per thread. advanceColumn() is called from inside a parallel
    // loop over columns, and a single shared workspace is a race that corrupts
    // the solve silently rather than crashing -- every thread would assemble
    // its tridiagonal on top of whatever its neighbors were half-way through.
    //
    // This is the second time this exact bug has been written in this codebase:
    // F7.1b found the same shared-workspace race in the momentum and scalar
    // column solves. It survives unit testing because the gates call the solver
    // from one thread; it only appears in the threaded model.
    const size_t n = static_cast<size_t>(layers + 1) *
                     static_cast<size_t>(std::max(1, maxThreads()));

    m_lower.assign(n, 0.0);
    m_diagonal.assign(n, 0.0);
    m_upper.assign(n, 0.0);
    m_rhs.assign(n, 0.0);
    m_solution.assign(n, 0.0);
    m_sweep.assign(n, 0.0);
}

void TurbulenceClosure::advanceColumn(int64_t column, double dt,
                                      const double *thickness,
                                      const double *shearSquared,
                                      const double *buoyancy,
                                      double frictionSurface, double frictionBed,
                                      double *viscosity, double *diffusivity)
{
    const int    n    = m_layers + 1;
    const size_t base = static_cast<size_t>(column) * static_cast<size_t>(n);
    const auto  &o    = m_options;

    // This thread's slice of the shared workspace.
    const size_t slice = static_cast<size_t>(threadIndex()) * static_cast<size_t>(n);

    double *lower    = m_lower.data() + slice;
    double *diagonal = m_diagonal.data() + slice;
    double *upper    = m_upper.data() + slice;
    double *rhs      = m_rhs.data() + slice;
    double *solution = m_solution.data() + slice;
    double *sweep    = m_sweep.data() + slice;

    // Spacing between interfaces: the thickness of the layer between them.
    auto spacingAt = [&](int i) {
        if (i <= 0)
            return std::max(thickness[0], 1.0e-12);

        if (i >= m_layers)
            return std::max(thickness[m_layers - 1], 1.0e-12);

        return std::max(0.5 * (thickness[i - 1] + thickness[i]), 1.0e-12);
    };

    // ── Production, and the eddy viscosity it is evaluated with ───────────
    //
    // Evaluated from the *previous* step's k and eps, which is what makes the
    // scheme explicit in production and implicit in transport and sink — the
    // standard arrangement, and the one that keeps the sink term on the
    // diagonal where it cannot drive k negative.
    for (int i = 0; i < n; ++i)
    {
        const size_t s = base + static_cast<size_t>(i);

        const double k   = std::max(m_tke[s], o.minTke);
        const double eps = std::max(m_dissipation[s], o.minDissipation);

        const double alphaN = buoyancy[i] * k * k / (eps * eps);

        const double sm = stabilityMomentum(alphaN);
        const double sh = stabilityScalar(alphaN);

        // nu_t = sqrt(2) s_m k^2 / eps in the Galperin scaling; the sqrt(2)
        // folds the k-vs-q² convention difference between the fit and this
        // formulation.
        const double scale = std::sqrt(2.0) * k * k / eps;

        viscosity[i]   = std::clamp(sm * scale, o.minViscosity, o.maxViscosity);
        diffusivity[i] = std::clamp(sh * scale, o.minDiffusivity, o.maxViscosity);
    }

    // ── k equation ────────────────────────────────────────────────────────
    for (int i = 0; i < n; ++i)
    {
        const size_t s = base + static_cast<size_t>(i);

        const double k   = std::max(m_tke[s], o.minTke);
        const double eps = std::max(m_dissipation[s], o.minDissipation);

        const double shear    = viscosity[i] * shearSquared[i];
        const double buoyant  = -diffusivity[i] * buoyancy[i];

        double lowerTerm = 0.0;
        double upperTerm = 0.0;

        const double h = spacingAt(i);

        if (i > 0)
        {
            const double face = 0.5 * (viscosity[i] + viscosity[i - 1]) / o.sigmaK;
            lowerTerm = -dt * face / (h * spacingAt(i - 1));
        }

        if (i + 1 < n)
        {
            const double face = 0.5 * (viscosity[i] + viscosity[i + 1]) / o.sigmaK;
            upperTerm = -dt * face / (h * spacingAt(i + 1));
        }

        lower[i] = lowerTerm;
        upper[i] = upperTerm;

        // Dissipation goes on the diagonal as eps/k, so the sink is implicit in
        // k and cannot overshoot into negative energy however large the step.
        diagonal[i] = 1.0 - lowerTerm - upperTerm + dt * eps / k;

        rhs[i] = k + dt * (shear + buoyant);
    }

    // Wall conditions: k = u*²/sqrt(c_mu) at both boundaries, which is the
    // constant-stress-layer equilibrium and the reason a log layer comes out.
    const double surfaceTke = std::max(frictionSurface * frictionSurface /
                                           std::sqrt(o.cMu),
                                       o.minTke);
    const double bedTke = std::max(frictionBed * frictionBed / std::sqrt(o.cMu),
                                   o.minTke);

    lower[0] = 0.0;
    upper[0] = 0.0;
    diagonal[0] = 1.0;
    rhs[0] = surfaceTke;

    const size_t last = static_cast<size_t>(n - 1);

    lower[last]    = 0.0;
    upper[last]    = 0.0;
    diagonal[last] = 1.0;
    rhs[last]      = bedTke;

    solveTridiagonal(n, lower, diagonal, upper, rhs, sweep, solution);

    for (int i = 0; i < n; ++i)
        m_tke[base + static_cast<size_t>(i)] =
            std::max(solution[i], o.minTke);

    // ── epsilon equation ──────────────────────────────────────────────────
    for (int i = 0; i < n; ++i)
    {
        const size_t s = base + static_cast<size_t>(i);

        const double k   = std::max(m_tke[s], o.minTke);
        const double eps = std::max(m_dissipation[s], o.minDissipation);

        const double shear   = viscosity[i] * shearSquared[i];
        const double buoyant = -diffusivity[i] * buoyancy[i];

        // c3 differs in sign between stable and unstable stratification: in a
        // stable column buoyancy destroys turbulence, in an unstable one it
        // makes it. One coefficient cannot do both.
        const double c3 = buoyancy[i] >= 0.0 ? o.c3Stable : o.c3Unstable;

        double lowerTerm = 0.0;
        double upperTerm = 0.0;

        const double h = spacingAt(i);

        if (i > 0)
        {
            const double face = 0.5 * (viscosity[i] + viscosity[i - 1]) / o.sigmaEps;
            lowerTerm = -dt * face / (h * spacingAt(i - 1));
        }

        if (i + 1 < n)
        {
            const double face = 0.5 * (viscosity[i] + viscosity[i + 1]) / o.sigmaEps;
            upperTerm = -dt * face / (h * spacingAt(i + 1));
        }

        lower[i] = lowerTerm;
        upper[i] = upperTerm;

        diagonal[i] = 1.0 - lowerTerm - upperTerm + dt * o.c2 * eps / k;

        rhs[i] =
            eps + dt * (eps / k) * (o.c1 * shear + c3 * buoyant);
    }

    // Wall dissipation, eps = u*³/(kappa z), evaluated at the distance the wall
    // node actually represents — the roughness length plus half the adjacent
    // layer — rather than at the roughness length alone.
    //
    // Using z0 by itself is a real error and it shows: for z0 = 0.02 m it sets
    // eps 25 times larger than the first interior point should see, which
    // over-dissipates the near-wall energy and puts a dip one node in from the
    // wall. The log-layer gate found it — interface 1 sat at 3.6e-5 while
    // interfaces 2 and 3 were near 1.5e-4, and turbulence does not have a local
    // minimum next to the surface that generates it.
    const double surfaceScale =
        o.vonKarman * (o.surfaceRoughness + 0.5 * std::max(thickness[0], 1.0e-12));

    const double bedScale =
        o.vonKarman *
        (o.bedRoughness + 0.5 * std::max(thickness[m_layers - 1], 1.0e-12));

    const double surfaceEps =
        std::max(frictionSurface * frictionSurface * frictionSurface / surfaceScale,
                 o.minDissipation);

    const double bedEps =
        std::max(frictionBed * frictionBed * frictionBed / bedScale,
                 o.minDissipation);

    lower[0]    = 0.0;
    upper[0]    = 0.0;
    diagonal[0] = 1.0;
    rhs[0]      = surfaceEps;

    lower[last]    = 0.0;
    upper[last]    = 0.0;
    diagonal[last] = 1.0;
    rhs[last]      = bedEps;

    solveTridiagonal(n, lower, diagonal, upper, rhs, sweep, solution);

    for (int i = 0; i < n; ++i)
        m_dissipation[base + static_cast<size_t>(i)] =
            std::max(solution[i], o.minDissipation);

    // ── Coefficients from the updated state ───────────────────────────────
    for (int i = 0; i < n; ++i)
    {
        const size_t s = base + static_cast<size_t>(i);

        const double k   = m_tke[s];
        const double eps = m_dissipation[s];

        const double alphaN = buoyancy[i] * k * k / (eps * eps);

        const double scale = std::sqrt(2.0) * k * k / eps;

        viscosity[i] = std::clamp(stabilityMomentum(alphaN) * scale,
                                  o.minViscosity, o.maxViscosity);
        diffusivity[i] = std::clamp(stabilityScalar(alphaN) * scale,
                                    o.minDiffusivity, o.maxViscosity);
    }
}

} // namespace FVQual::Hydro
