// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  crosssection_agreement.cpp
 * \brief Does HydroCouple's ICrossSection contract survive two different
 *        underlying representations of the same channel? (F8.2f)
 *
 * The interface is built on a claim: that answering *derived* quantities —
 * width, storage area, flow area, wetted perimeter at a stage — lets a
 * station-elevation model (HEC-RAS, SWMM) and a width-elevation model
 * (CE-QUAL-W2) describe the same channel and agree, where an interface built
 * around survey points would have forced W2 to invent stations it does not have.
 *
 * That claim is cheap to assert in a header and worth measuring. This program
 * implements the same trapezoid both ways and compares them against the closed
 * form, at every stage from the invert to the top of bank.
 *
 * It lives in verification/ rather than src/ deliberately: FVQual is a 3-D mesh
 * model and has no use for a cross-section class of its own. This is here to
 * test the interface, not to ship behind it.
 *
 * Build (from the repository root):
 *
 *     c++ -std=c++20 -I../HydroCouple/include \
 *         verification/interfaces/crosssection_agreement.cpp -o crosssection_agreement
 */

#include "hydrocouplespatial.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using HydroCouple::Spatial::CrossSectionKind;
using HydroCouple::Spatial::ICrossSection;

namespace {

// ── The channel under test ────────────────────────────────────────────────
//
// A symmetric trapezoid: invert at elevation 2, bottom width 6, sides rising
// 8 m over 4 m horizontally on each side, so the width grows by exactly 1 m
// per meter of stage and the top of bank is at elevation 10, 14 m wide.

constexpr double kInvert      = 2.0;
constexpr double kTopOfBank   = 10.0;
constexpr double kBottomWidth = 6.0;

//! Closed form, to keep both implementations honest about a third opinion.
double exactWidth(double stage)
{
    // Strictly below, not at: at the invert a flat bottom is wetted with zero
    // depth, so the limiting width is the bottom width. See ICrossSection.
    if (stage < kInvert) return 0.0;
    return kBottomWidth + (std::min(stage, kTopOfBank) - kInvert);
}

double exactArea(double stage)
{
    if (stage < kInvert) return 0.0;
    const double d = std::min(stage, kTopOfBank) - kInvert;
    return kBottomWidth * d + 0.5 * d * d;
}

double exactPerimeter(double stage)
{
    if (stage < kInvert) return 0.0;
    const double d = std::min(stage, kTopOfBank) - kInvert;
    // Each bank retreats 0.5 m horizontally per meter of rise.
    return kBottomWidth + 2.0 * d * std::sqrt(1.0 + 0.25);
}

// ── Representation 1: survey points, as HEC-RAS and SWMM hold them ────────

class StationElevationSection final : public ICrossSection
{
public:
    StationElevationSection(std::vector<double> stations,
                            std::vector<double> elevations,
                            double leftBank, double rightBank)
        : m_station(std::move(stations)), m_elevation(std::move(elevations)),
          m_leftBank(leftBank), m_rightBank(rightBank)
    {
    }

    [[nodiscard]] CrossSectionKind kind() const override
    {
        return CrossSectionKind::StationElevation;
    }

    [[nodiscard]] double invertElevation() const override
    {
        return *std::min_element(m_elevation.begin(), m_elevation.end());
    }

    [[nodiscard]] double topWidth(double stage) const override
    {
        return integrate(stage, Quantity::Width, -1.0e300, 1.0e300);
    }

    [[nodiscard]] double storageArea(double stage) const override
    {
        return integrate(stage, Quantity::Area, -1.0e300, 1.0e300);
    }

    [[nodiscard]] double flowArea(double stage) const override
    {
        return integrate(stage, Quantity::Area, m_leftBank, m_rightBank);
    }

    [[nodiscard]] double wettedPerimeter(double stage) const override
    {
        return integrate(stage, Quantity::Perimeter, m_leftBank, m_rightBank);
    }

    [[nodiscard]] int64_t stationCount() const override
    {
        return static_cast<int64_t>(m_station.size());
    }

    void stations(double *stations, double *elevations) const override
    {
        for (size_t i = 0; i < m_station.size(); ++i)
        {
            stations[i]   = m_station[i];
            elevations[i] = m_elevation[i];
        }
    }

private:
    enum class Quantity
    {
        Width,
        Area,
        Perimeter
    };

    //! Walks the panels between adjacent survey points, clipped to the water
    //! surface and to the conveying limits.
    [[nodiscard]] double integrate(double stage, Quantity quantity, double from,
                                   double to) const
    {
        double total = 0.0;

        for (size_t i = 0; i + 1 < m_station.size(); ++i)
        {
            double x0 = m_station[i], z0 = m_elevation[i];
            double x1 = m_station[i + 1], z1 = m_elevation[i + 1];

            if (z0 > stage && z1 > stage) continue;

            // Clip the dry end of a partially submerged panel.
            if (z0 > stage)
            {
                const double t = (stage - z1) / (z0 - z1);
                x0             = x1 + t * (x0 - x1);
                z0             = stage;
            }
            else if (z1 > stage)
            {
                const double t = (stage - z0) / (z1 - z0);
                x1             = x0 + t * (x1 - x0);
                z1             = stage;
            }

            // Clip to the conveying window.
            const double lo = std::max(std::min(x0, x1), from);
            const double hi = std::min(std::max(x0, x1), to);
            if (hi <= lo)
            {
                continue;
            }

            const double span = x1 - x0;
            const double frac = span != 0.0 ? (hi - lo) / std::abs(span) : 1.0;

            switch (quantity)
            {
            case Quantity::Width:
                total += hi - lo;
                break;
            case Quantity::Area:
            {
                const double d0 = stage - z0, d1 = stage - z1;
                total += 0.5 * (d0 + d1) * (hi - lo);
                break;
            }
            case Quantity::Perimeter:
                total += std::hypot(span, z1 - z0) * frac;
                break;
            }
        }

        return total;
    }

    std::vector<double> m_station;
    std::vector<double> m_elevation;
    double              m_leftBank;
    double              m_rightBank;
};

// ── Representation 2: width against elevation, as CE-QUAL-W2 holds it ─────

class WidthElevationSection final : public ICrossSection
{
public:
    WidthElevationSection(std::vector<double> elevations, std::vector<double> widths)
        : m_elevation(std::move(elevations)), m_width(std::move(widths))
    {
    }

    [[nodiscard]] CrossSectionKind kind() const override
    {
        return CrossSectionKind::WidthElevation;
    }

    [[nodiscard]] double invertElevation() const override { return m_elevation.front(); }

    [[nodiscard]] double topWidth(double stage) const override
    {
        // Strictly below the invert, per the ICrossSection convention: at the
        // invert itself the flat bottom is wetted at zero depth, so the width
        // is the bottom width and not zero.
        if (stage < m_elevation.front()) return 0.0;

        for (size_t i = 0; i + 1 < m_elevation.size(); ++i)
        {
            if (stage <= m_elevation[i + 1])
            {
                const double t = (stage - m_elevation[i]) /
                                 (m_elevation[i + 1] - m_elevation[i]);
                return m_width[i] + t * (m_width[i + 1] - m_width[i]);
            }
        }
        return m_width.back();
    }

    [[nodiscard]] double storageArea(double stage) const override
    {
        // The integral of width with respect to elevation, which is what this
        // representation is naturally good at and the station form is not.
        double area = 0.0;
        for (size_t i = 0; i + 1 < m_elevation.size(); ++i)
        {
            if (stage <= m_elevation[i]) break;
            const double top = std::min(stage, m_elevation[i + 1]);
            const double dz  = top - m_elevation[i];
            if (dz <= 0.0) continue;
            area += 0.5 * (m_width[i] + topWidth(top)) * dz;
        }
        return area;
    }

    //! No ineffective flow area is representable here, and saying so by
    //! returning the storage area is the honest answer for this producer.
    [[nodiscard]] double flowArea(double stage) const override
    {
        return storageArea(stage);
    }

    [[nodiscard]] double wettedPerimeter(double stage) const override
    {
        // Each bank retreats by half the width increment per unit rise, so the
        // sloping length follows even without stations.
        double perimeter = m_width.front();
        for (size_t i = 0; i + 1 < m_elevation.size(); ++i)
        {
            if (stage <= m_elevation[i]) break;
            const double top = std::min(stage, m_elevation[i + 1]);
            const double dz  = top - m_elevation[i];
            if (dz <= 0.0) continue;
            const double dHalfWidth = 0.5 * (topWidth(top) - m_width[i]);
            perimeter += 2.0 * std::hypot(dz, dHalfWidth);
        }
        return perimeter;
    }

    //! The whole point: this producer has no stations, and says so.
    [[nodiscard]] int64_t stationCount() const override { return 0; }

    void stations(double *, double *) const override {}

private:
    std::vector<double> m_elevation;
    std::vector<double> m_width;
};

int    failures = 0;
double worst    = 0.0;

void check(const char *what, double got, double expected, double tolerance)
{
    const double error = std::abs(got - expected);
    worst              = std::max(worst, error);
    if (!(error <= tolerance))
    {
        std::printf("  FAIL %-38s got %12.6f  expected %12.6f\n", what, got, expected);
        ++failures;
    }
}

} // namespace

int main()
{
    // The same trapezoid, described two ways.
    StationElevationSection stationForm({0.0, 4.0, 10.0, 14.0},
                                        {10.0, 2.0, 2.0, 10.0}, 0.0, 14.0);

    WidthElevationSection widthForm({2.0, 4.0, 6.0, 8.0, 10.0},
                                    {6.0, 8.0, 10.0, 12.0, 14.0});

    std::printf("Cross-section agreement across representations (F8.2f)\n");
    std::printf("  station-elevation kind = %d, stations = %lld\n",
                static_cast<int>(stationForm.kind()),
                static_cast<long long>(stationForm.stationCount()));
    std::printf("  width-elevation   kind = %d, stations = %lld\n\n",
                static_cast<int>(widthForm.kind()),
                static_cast<long long>(widthForm.stationCount()));

    check("invert (station form)", stationForm.invertElevation(), kInvert, 1.0e-12);
    check("invert (width form)", widthForm.invertElevation(), kInvert, 1.0e-12);

    for (double stage = kInvert; stage <= kTopOfBank + 1.0e-9; stage += 0.25)
    {
        check("station width vs closed form", stationForm.topWidth(stage),
              exactWidth(stage), 1.0e-9);
        check("width-form width vs closed form", widthForm.topWidth(stage),
              exactWidth(stage), 1.0e-9);

        check("station area vs closed form", stationForm.storageArea(stage),
              exactArea(stage), 1.0e-9);
        check("width-form area vs closed form", widthForm.storageArea(stage),
              exactArea(stage), 1.0e-9);

        check("station perimeter vs closed form", stationForm.wettedPerimeter(stage),
              exactPerimeter(stage), 1.0e-9);
        check("width-form perimeter vs closed form", widthForm.wettedPerimeter(stage),
              exactPerimeter(stage), 1.0e-9);

        // The claim in one line: two representations, one answer.
        check("the two representations agree on area", stationForm.storageArea(stage),
              widthForm.storageArea(stage), 1.0e-9);
    }

    // And the flow/storage split must not be collapsible. A section whose
    // banks sit inside the surveyed extent conveys over less than it stores.
    StationElevationSection withFloodplain({0.0, 4.0, 10.0, 14.0},
                                           {10.0, 2.0, 2.0, 10.0}, 4.0, 10.0);

    const double storage = withFloodplain.storageArea(kTopOfBank);
    const double flow    = withFloodplain.flowArea(kTopOfBank);

    std::printf("Ineffective flow area, at top of bank:\n");
    std::printf("  storage area %.4f m2, flow area %.4f m2, conveying %.1f%%\n\n",
                storage, flow, 100.0 * flow / storage);

    if (!(flow < storage))
    {
        std::printf("  FAIL flow area was not less than storage area\n");
        ++failures;
    }

    std::printf("%s  (%d failures, worst residual %.3e)\n",
                failures == 0 ? "PASS" : "FAILED", failures, worst);
    return failures == 0 ? 0 : 1;
}
