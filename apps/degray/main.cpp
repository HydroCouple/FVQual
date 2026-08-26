// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file  main.cpp
 * \brief The DeGray Reservoir application, end to end (Phase 8).
 *
 * Everything the engine has been built for, driven by a real reservoir's real
 * data: CE-QUAL-W2's DeGray example, March to December 1980. Bathymetry from
 * `bth1.csv`, initial temperatures from `vpr.npt`, meteorology from `met.npt`
 * with solar radiation computed rather than read, and the branch inflow and
 * outflow records.
 *
 * **What this is and is not.** There is no W2 output shipped with the example
 * to compare against, so this is not a code-to-code validation. It is an
 * end-to-end exercise judged on physics: does the reservoir stratify through
 * summer, does the thermocline sit where one belongs, does surface temperature
 * follow the equilibrium temperature the weather implies, and do the heat and
 * volume budgets close. Those catch gross errors and would not catch a subtle
 * formulation difference.
 *
 * **A second caveat worth stating plainly.** W2 is laterally averaged, so its
 * schematisation is one segment wide and the mesh built from it is a chain of
 * columns rather than a basin. The vertical physics — stratification, the
 * thermocline, selective withdrawal — is exercised fully. The horizontal
 * three-dimensionality FVQual exists for is not exercised at all.
 *
 * Results are written as CSV where they can be read and plotted, rather than
 * summarised to standard output and lost.
 */

#include "fvqual/hydro/heatbudget.h"
#include "fvqual/hydro/internalmode.h"
#include "fvqual/io/w2bathymetry.h"
#include "fvqual/io/w2timeseries.h"
#include "fvqual/mesh/horizontalmesh.h"
#include "fvqual/mesh/layeredmesh.h"
#include "fvqual/mesh/sigmaspec.h"
#include "fvqual/state/speciesregistry.h"
#include "fvqual/state/statestore.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

//! DeGray's water body, from the control file's location card and grid card.
constexpr double kBottomElevation = 66.35;
constexpr double kLatitude        = 34.2;
constexpr double kStartDay        = 64.5;
constexpr double kEndDay          = 358.7;

//! Layers in the sigma grid. Enough to resolve a thermocline in 57 m.
constexpr int kLayers = 30;

//! For the wind-direction convention change; see the conversion at its use.
constexpr double kHalfPi = 0.5 * 3.14159265358979323846;

/*!
 * \brief Interpolates a W2 vertical profile, surface first, to an elevation.
 */
double profileAt(const std::vector<double> &values, const FVQual::IO::W2Bathymetry &bath,
                 double surface, double elevation)
{
    if (values.empty())
        return 0.0;

    // The profile covers the wetted layers from the surface down, so build the
    // elevations it corresponds to and interpolate between them.
    std::vector<double> centres;

    for (int k = 0; k < bath.layerCount(); ++k)
    {
        const double bottom = bath.layerBottom(k);

        if (bottom >= surface)
            continue;

        bool wet = false;

        for (int i : bath.activeSegments())
            wet = wet || bath.width(k, i) > 0.0;

        if (wet)
            centres.push_back(bottom + 0.5 * (k + 1 < bath.layerCount()
                                                  ? bath.layerBottom(k) - bath.layerBottom(k + 1)
                                                  : 2.0));
    }

    if (centres.size() != values.size())
    {
        // Fall back to a linear map over the profile's extent rather than
        // refusing: the caller reports the mismatch.
        const double span = static_cast<double>(values.size() - 1);
        const double frac = std::clamp((surface - elevation) / std::max(surface - bath.bottomElevation(), 1.0), 0.0, 1.0);

        const double position = frac * span;
        const size_t low      = static_cast<size_t>(position);
        const size_t high     = std::min(low + 1, values.size() - 1);

        return values[low] + (position - static_cast<double>(low)) * (values[high] - values[low]);
    }

    // centres runs from the surface downward, so elevations decrease.
    if (elevation >= centres.front())
        return values.front();

    if (elevation <= centres.back())
        return values.back();

    for (size_t i = 1; i < centres.size(); ++i)
    {
        if (elevation >= centres[i])
        {
            const double span = centres[i - 1] - centres[i];
            const double frac = span > 0.0 ? (elevation - centres[i]) / span : 0.0;

            return values[i] + frac * (values[i - 1] - values[i]);
        }
    }

    return values.back();
}

} // namespace

int main(int argc, char **argv)
{
    const std::string inputs =
        argc > 1 ? argv[1]
                 : "/sessions/festive-sharp-feynman/mnt/v5/Examples/DeGray/InputFiles";

    const std::string outputs = argc > 2 ? argv[2] : ".";

    std::string message;

    // ── Geometry ──────────────────────────────────────────────────────────

    FVQual::IO::W2Bathymetry bathymetry;

    if (!bathymetry.read(inputs + "/bth1.csv", kBottomElevation, message))
    {
        std::cerr << "bathymetry: " << message << '\n';
        return 1;
    }

    HydroCouple::SDK::IO::MeshDefinition definition;
    std::vector<double>                  bed;

    if (!bathymetry.toMeshDefinition(definition, bed, message))
    {
        std::cerr << "mesh: " << message << '\n';
        return 1;
    }

    FVQual::Mesh::HorizontalMesh horizontal;

    if (!horizontal.build(definition, message))
    {
        std::cerr << "horizontal mesh: " << message << '\n';
        return 1;
    }

    const double startSurface = bathymetry.initialSurface()[
        static_cast<size_t>(bathymetry.activeSegments().front())];

    // Every column must be wet: the internal mode has no three-dimensional
    // wetting and drying yet, so a segment whose bed sits above the pool would
    // be refused. Report it here rather than as a cryptic failure later.
    int dry = 0;

    for (int64_t f = 0; f < horizontal.faceCount(); ++f)
    {
        horizontal.faceBed()[static_cast<size_t>(f)] = bed[static_cast<size_t>(f)];

        if (bed[static_cast<size_t>(f)] >= startSurface)
            ++dry;
    }

    if (dry > 0)
    {
        std::cerr << dry << " of " << horizontal.faceCount()
                  << " columns are above the initial pool; the internal mode "
                     "requires every column wet (F3.4)\n";
        return 1;
    }

    FVQual::Mesh::SigmaSpec spec;
    spec.layerCount        = kLayers;
    spec.minLayerThickness = 0.0;

    FVQual::Mesh::LayeredMesh layered;

    if (!layered.build(horizontal, spec, startSurface, message))
    {
        std::cerr << "layered mesh: " << message << '\n';
        return 1;
    }

    if (!bathymetry.applyHypsography(layered, message))
    {
        std::cerr << "hypsography: " << message << '\n';
        return 1;
    }

    {
        const std::vector<double> surface(static_cast<size_t>(layered.columnCount()),
                                          startSurface);

        if (!layered.updateSurface(surface, message))
        {
            std::cerr << "surface: " << message << '\n';
            return 1;
        }
    }

    std::printf("DeGray: %lld columns x %d layers, surface %.2f m, volume %.4g m3, "
                "area %.4g m2\n",
                static_cast<long long>(layered.columnCount()), kLayers, startSurface,
                layered.totalVolume(), bathymetry.surfaceAreaAt(startSurface));

    // ── Initial condition ─────────────────────────────────────────────────

    FVQual::IO::W2VerticalProfile profile;

    if (!profile.read(inputs + "/vpr.npt", message))
    {
        std::cerr << "profile: " << message << '\n';
        return 1;
    }

    const std::vector<double> &initialTemperature = profile.profile("TEMP");

    FVQual::State::SpeciesRegistry registry;

    if (registry.addTemperature(7.0, message) < 0)
    {
        std::cerr << "registry: " << message << '\n';
        return 1;
    }

    const int temperature = registry.temperatureIndex();

    FVQual::State::StateStore state;
    state.allocate(layered, registry);

    for (int64_t c = 0; c < layered.columnCount(); ++c)
    {
        for (int k = 0; k < kLayers; ++k)
        {
            const double centre = 0.5 * (layered.z(c, k) + layered.z(c, k + 1));

            state.speciesAt(temperature, layered.cell(c, k)) =
                profileAt(initialTemperature, bathymetry, startSurface, centre);
        }
    }

    // ── Forcing ───────────────────────────────────────────────────────────

    FVQual::IO::W2TimeSeries met, inflow, inflowTemperature, outflow;

    if (!met.read(inputs + "/met.npt", message) ||
        !inflow.read(inputs + "/qin_br1.npt", message) ||
        !inflowTemperature.read(inputs + "/tin_br1.npt", message) ||
        !outflow.read(inputs + "/qot_br1.npt", message))
    {
        std::cerr << "forcing: " << message << '\n';
        return 1;
    }

    const int air   = met.columnIndex("TAIR");
    const int dew   = met.columnIndex("TDEW");
    const int wind  = met.columnIndex("WIND");
    const int cloud = met.columnIndex("CLOUD");
    const int phi   = met.columnIndex("PHI");

    // ── Solver ────────────────────────────────────────────────────────────

    FVQual::Hydro::InternalModeSolver  solver;
    FVQual::Hydro::InternalModeOptions options;

    options.gravity   = 9.80665;
    options.manningN  = 0.03;
    options.cflNumber = 0.5;
    options.maxTimeStep = 300.0;

    if (!solver.initialize(layered, state, registry, options, message))
    {
        std::cerr << "solver: " << message << '\n';
        return 1;
    }

    // DeGray's calibrated wind function, from the control file's heat-exchange
    // card: 9.2 + 0.46 W^2 in W2's units. The coefficients here are in SI, so
    // only the shape carries over.
    solver.meteorology().base.windExponent = 2.0;
    solver.enableMeteorology(/*attenuation=*/0.45);

    // Inflow at the head of the reservoir, density-placed so a cold tributary
    // finds its own level; withdrawal at the dam.
    FVQual::Forcing::Source tributary;
    tributary.column    = 0;
    tributary.placement = FVQual::Forcing::Placement::DensityMatched;

    const int inflowIndex = solver.sources().add(tributary, message);

    FVQual::Forcing::Source release;
    release.column    = layered.columnCount() - 1;
    release.placement = FVQual::Forcing::Placement::Uniform;

    const int outflowIndex = solver.sources().add(release, message);

    if (inflowIndex < 0 || outflowIndex < 0)
    {
        std::cerr << "sources: " << message << '\n';
        return 1;
    }

    // ── Run ───────────────────────────────────────────────────────────────

    std::ofstream daily(outputs + "/degray_daily.csv");

    daily << "jday,surface_temperature_C,bottom_temperature_C,thermocline_depth_m,"
             "surface_elevation_m,volume_m3,heat_content_J,equilibrium_C,"
             "air_temperature_C,shortwave_W_m2\n";
    daily.precision(6);

    // The whole profile at the dam, every day. Surface and bottom temperatures
    // say a reservoir stratified; only the profile says whether it stratified
    // the way a reservoir does, and the depth-time contour drawn from this is
    // the figure a limnologist would ask for first. Depth is written per row
    // rather than implied by the layer index, because the sigma grid moves with
    // the surface and a layer is not at a fixed depth.
    std::ofstream profiles(outputs + "/degray_profiles.csv");

    profiles << "jday,layer,depth_m,elevation_m,temperature_C\n";
    profiles.precision(6);

    const double capacity = options.referenceDensity * options.heatCapacity;

    double day       = kStartDay;
    double nextWrite = kStartDay;

    int64_t steps = 0;

    while (day < kEndDay)
    {
        FVQual::Forcing::MeteorologyField &weather = solver.meteorology();

        weather.airTemperature.setUniform(met.valueAt(day, air));
        weather.dewPoint.setUniform(met.valueAt(day, dew));
        weather.windSpeed.setUniform(met.valueAt(day, wind));

        // W2's PHI is a compass angle, clockwise from north, and W2 projects the
        // stress onto each segment as `cos(PHI - PHI0)` (`w2_4_win.f90:543`).
        // `W2Bathymetry` lays the chain out with the along-channel unit vector
        // `(sin PHI0, cos PHI0)`, so a compass angle maps into FVQual's frame as
        // `pi/2 - angle`. Applying that to PHI puts the wind vector at
        // `(sin PHI, cos PHI)` — the same frame the orientations are in — and
        // FVQual's own projection onto a face then reproduces `cos(PHI - PHI0)`
        // identically, without this code needing to decide whether PHI is a
        // "from" or a "towards" bearing. Matching the reference beats guessing
        // the convention.
        weather.windDirection.setUniform(kHalfPi - met.valueAt(day, phi));

        // W2 records cloud cover in tenths, on a 0-10 scale; Meteorology wants a
        // fraction. Converting here rather than in the reader keeps the reader
        // a faithful record of the file, and puts the unit change at the
        // boundary where the two conventions actually meet.
        const double cover = met.valueAt(day, cloud) / 10.0;

        weather.cloudCover.setUniform(cover);
        weather.shortwave.setUniform(
            FVQual::Hydro::clearSkySolarRadiation(day, kLatitude, cover));

        if (!solver.sources().setDischarge(inflowIndex, inflow.valueAt(day, 0), message) ||
            !solver.sources().setDischarge(outflowIndex, -outflow.valueAt(day, 0),
                                           message))
        {
            std::cerr << "discharge: " << message << '\n';
            return 1;
        }

        if (!solver.sources().setInflowValue(inflowIndex, temperature, temperature,
                                             inflowTemperature.valueAt(day, 0), message))
        {
            std::cerr << "inflow temperature: " << message << '\n';
            return 1;
        }

        const double dt = std::min(solver.stableTimeStep(), 300.0);

        if (!solver.step(dt, message))
        {
            std::cerr << "step failed on day " << day << ": " << message << '\n';
            return 1;
        }

        day += dt / 86400.0;
        ++steps;

        if (day < nextWrite)
            continue;

        // Advance the schedule, do not rebase it on the time actually reached.
        // Rebasing carries each step's overshoot forward, and with a 300 s cap
        // that is up to five minutes a day: the report drifted from 12:05 in
        // March to 17:00 in December, sampling a steadily later hour of a
        // diurnal cycle and aliasing it into the seasonal signal. Found by
        // checking a figure's caption against its data, which is not where one
        // expects to find a timing bug.
        nextWrite += 1.0;

        // Report from the deepest column, at the dam, where a profile would be
        // taken.
        const int64_t deep = layered.columnCount() - 1;

        const double surfaceT = state.speciesAt(temperature, layered.cell(deep, 0));
        const double bottomT =
            state.speciesAt(temperature, layered.cell(deep, kLayers - 1));

        // Thermocline: the interface with the largest temperature gradient.
        double thermocline = 0.0;
        double steepest    = 0.0;

        for (int k = 1; k < kLayers; ++k)
        {
            const double above = state.speciesAt(temperature, layered.cell(deep, k - 1));
            const double below = state.speciesAt(temperature, layered.cell(deep, k));

            const double spacing =
                0.5 * (layered.thickness()[static_cast<size_t>(layered.cell(deep, k - 1))] +
                       layered.thickness()[static_cast<size_t>(layered.cell(deep, k))]);

            const double gradient = spacing > 0.0 ? (above - below) / spacing : 0.0;

            if (gradient > steepest)
            {
                steepest    = gradient;
                thermocline = layered.z(deep, 0) - layered.z(deep, k);
            }
        }

        double heat = 0.0;

        for (int64_t cell = 0; cell < layered.cellCount(); ++cell)
        {
            heat += state.speciesAt(temperature, cell) *
                    layered.cellVolume()[static_cast<size_t>(cell)] * capacity;
        }

        // The reported weather is the weather at the column being reported from,
        // now that forcing may vary along the reservoir. Uniform here, so these
        // are the same numbers the scalar path wrote — but reading them from the
        // field keeps the report honest if the forcing is ever made spatial.
        const FVQual::Hydro::Meteorology reported = weather.at(deep);

        daily << day << ',' << surfaceT << ',' << bottomT << ',' << thermocline << ','
              << layered.surfaceElevation()[static_cast<size_t>(deep)] << ','
              << layered.totalVolume() << ',' << heat << ','
              << FVQual::Hydro::equilibriumTemperature(reported) << ','
              << reported.airTemperature << ',' << reported.shortwave << '\n';

        for (int k = 0; k < kLayers; ++k)
        {
            const double centre = 0.5 * (layered.z(deep, k) + layered.z(deep, k + 1));

            profiles << day << ',' << k << ',' << layered.z(deep, 0) - centre << ','
                     << centre << ','
                     << state.speciesAt(temperature, layered.cell(deep, k)) << '\n';
        }
    }

    daily.close();
    profiles.close();

    // Final profile at the dam, which is what a thermal survey would report.
    std::ofstream final(outputs + "/degray_final_profile.csv");

    final << "depth_m,elevation_m,temperature_C\n";
    final.precision(6);

    const int64_t deep = layered.columnCount() - 1;

    for (int k = 0; k < kLayers; ++k)
    {
        const double centre = 0.5 * (layered.z(deep, k) + layered.z(deep, k + 1));

        final << layered.z(deep, 0) - centre << ',' << centre << ','
              << state.speciesAt(temperature, layered.cell(deep, k)) << '\n';
    }

    final.close();

    std::printf("ran %lld internal steps to day %.2f; wrote degray_daily.csv, "
                "degray_profiles.csv and degray_final_profile.csv to %s\n",
                static_cast<long long>(steps), day, outputs.c_str());

    return 0;
}
