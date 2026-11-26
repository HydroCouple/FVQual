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
#include "fvqual/forcing/operationsconfiguration.h"
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
#include <cstdlib>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <stdexcept>
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
constexpr int kDefaultLayers = 30;

/*!
 * \brief Lateral cells per W2 segment, overridable through `FVQUAL_LANES`.
 *
 * W2 is laterally averaged, so the default of 1 reproduces the schematisation
 * exactly. Above 1 the same bathymetry is laid out across the channel — no new
 * survey data, just the width curve's lateral information used as lateral
 * information — which is the only way to ask whether DeGray's interior warm
 * bias is the one-segment grid rather than the physics.
 */
/*!
 * \brief The deepest lane of the dam segment.
 *
 * A cast at the dam is taken where the water is, which in a section laid out
 * across the channel is the middle. Reporting the last column instead would pick
 * the outermost lane — shallower by construction — so widening the mesh would
 * appear to change the profile when it had only changed where the profile was
 * measured.
 */
int64_t deepestDamColumn(const FVQual::Mesh::LayeredMesh &layered, int lanes)
{
    const int64_t first = layered.columnCount() - lanes;

    int64_t best      = first;
    double  deepestSoFar = layered.columnDepth()[static_cast<size_t>(first)];

    for (int64_t c = first + 1; c < layered.columnCount(); ++c)
    {
        const double depth = layered.columnDepth()[static_cast<size_t>(c)];

        if (depth > deepestSoFar)
        {
            deepestSoFar = depth;
            best         = c;
        }
    }

    return best;
}

int lateralCells()
{
    if (const char *value = std::getenv("FVQUAL_LANES"))
    {
        const int lanes = std::atoi(value);

        if (lanes >= 1)
            return lanes;
    }

    return 1;
}

/*!
 * \brief Layer count, overridable through `FVQUAL_LAYERS`.
 *
 * A knob for convergence testing, not for configuration. If an error shrinks as
 * the grid refines it is the discretisation; if it does not, it is the physics
 * or the forcing, and no amount of resolution will help. That distinction is
 * worth an environment variable.
 */
int layerCount()
{
    if (const char *value = std::getenv("FVQUAL_LAYERS"))
    {
        const int layers = std::atoi(value);

        if (layers >= 4)
            return layers;
    }

    return kDefaultLayers;
}

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

    // Optional third argument: stop early. The full year takes long enough that
    // a diagnostic question about spring stratification should not have to wait
    // for December, and a bounded window is also the only way to run this where
    // a wall clock is enforced.
    const double endDay = argc > 3 ? std::atof(argv[3]) : kEndDay;
    if (!std::isfinite(endDay) || endDay <= kStartDay || endDay > kEndDay) {
        std::cerr << "end day must be in (64.5, 358.7]\n";
        return 1;
    }
    std::filesystem::create_directories(outputs);
    double maximumStep = 300.0;
    if (const char *value = std::getenv("FVQUAL_MAX_DT")) {
        char *end = nullptr;
        maximumStep = std::strtod(value, &end);
        if (!end || *end || !std::isfinite(maximumStep) || maximumStep <= 0) {
            std::cerr << "FVQUAL_MAX_DT must be finite and positive\n";
            return 1;
        }
    }

    // Optional fourth argument, "kepsilon", switches vertical mixing from
    // Richardson damping to the two-equation closure (F5.4). A switch rather
    // than a rebuild, because the only way to know which is better on this
    // reservoir is to run both against the same measurements.
    const bool closure = argc > 4 && std::string(argv[4]) == "kepsilon";

    // Bisection: run the same mesh as a prismatic basin, with no storage curve,
    // to separate a hypsography failure from everything else.
    const bool nohyps = argc > 5 && std::string(argv[5]) == "nohyps";

    // Bisection for the residual interior warm bias: keep wind stress, and the
    // circulation it drives, but switch off turbulent vertical mixing. If the
    // interior still warms, the heat is being advected down by the wind-driven
    // cell rather than mixed down by the closure.
    const bool nomix = argc > 5 && std::string(argv[5]) == "nomix";
    const bool scalarNomix = argc > 5 && std::string(argv[5]) == "scalar_nomix";
    const bool viscosityNomix = argc > 5 && std::string(argv[5]) == "viscosity_nomix";

    bool selectiveOutlet = false;
    double outletElevation = 115.0;
    if (const char *value = std::getenv("FVQUAL_OUTLET_ELEVATION")) {
        char *end = nullptr;
        outletElevation = std::strtod(value, &end);
        if (!end || *end || !std::isfinite(outletElevation) || outletElevation <= kBottomElevation) {
            std::cerr << "FVQUAL_OUTLET_ELEVATION must be above the bottom datum\n";
            return 1;
        }
        selectiveOutlet = true;
    }
    const int kLayers = layerCount();

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

    // FVQUAL_LANES is a ceiling; FVQUAL_MIN_CELL, when set, keeps any lateral
    // cell from falling below that width. Uniform refinement makes 19.8 m cells
    // in DeGray's 79 m arms at four lanes and the barotropic step collapses on
    // water that has nothing to do with the question; a minimum width refines
    // the wide reservoir and leaves the arms alone.
    double minCellWidth = 0.0;

    if (const char *value = std::getenv("FVQUAL_MIN_CELL"))
        minCellWidth = std::atof(value);

    const auto plan  = bathymetry.planLateralCells(lateralCells(), minCellWidth);
    const int  lanes     = plan.lanes.back();  //!< Lanes at the dam.
    const int  headLanes = plan.lanes.front(); //!< Lanes at the head, which may differ.

    // Which W2 segment each column came from. Per-segment inputs — wind
    // sheltering, and anything else the control file varies along the
    // reservoir — have to be looked up through this rather than by assuming
    // one column per segment, which stopped being true once lanes existed.
    std::vector<int> segmentOfColumn(static_cast<size_t>(plan.columnCount), 0);

    for (size_t j = 0; j < plan.lanes.size(); ++j)
        for (int l = 0; l < plan.lanes[j]; ++l)
            segmentOfColumn[static_cast<size_t>(plan.offset[j] + l)] =
                static_cast<int>(j);

    if (!bathymetry.toPlannedMeshDefinition(plan, definition, bed, message))
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

    if (!nohyps && !bathymetry.applyPlannedHypsography(layered, plan, message))
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

    std::printf("mixing: %s\n", closure ? "k-epsilon closure" : "Richardson damping");

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

    // Wind sheltering, `wsc.npt`. W2 scales wind SPEED by this before both the
    // surface heat budget and the momentum stress, so a sheltered reservoir
    // neither evaporates nor overturns as though it sat in open water. DeGray's
    // calibrated value is 0.90 everywhere; stress goes as the square of speed,
    // so ignoring it overstates the wind stress by 23%.
    //
    // Optional, because a model without the file is a model with no sheltering.
    FVQual::IO::W2TimeSeries sheltering;
    bool                     hasSheltering = false;
    {
        std::string ignored;
        // The record is as wide as the bathymetry has segments — the file's
        // own header under-declares it, so the count has to come from here.
        hasSheltering = sheltering.read(inputs + "/wsc.npt",
                                        bathymetry.segmentCount(), ignored);
    }

    // Say so out loud. A calibrated input that is silently absent looks exactly
    // like one that is present and equal to one, and the difference here is a
    // 23% error in the wind stress.
    if (hasSheltering)
        std::printf("wind sheltering: wsc.npt, %d segments, %.2f at the dam\n",
                    sheltering.columnCount(),
                    sheltering.valueAt(kStartDay, sheltering.columnCount() - 1));
    else
        std::printf("wind sheltering: none (no wsc.npt) - full open-water stress\n");

    if (!met.read(inputs + "/met.npt", message) ||
        !inflow.read(inputs + "/qin_br1.npt", message) ||
        !inflowTemperature.read(inputs + "/tin_br1.npt", message) ||
        !outflow.read(inputs + "/qot_br1.npt", message))
    {
        std::cerr << "forcing: " << message << '\n';
        return 1;
    }

    for (const auto *series : {&met, &inflow, &inflowTemperature, &outflow}) {
        if (series->firstDay() > kStartDay || series->lastDay() < endDay) {
            std::cerr << "forcing does not cover requested simulation interval\n";
            return 1;
        }
    }
    if (hasSheltering && (sheltering.firstDay() > kStartDay || sheltering.lastDay() < endDay)) {
        std::cerr << "wind sheltering does not cover requested interval\n";
        return 1;
    }
    const int air   = met.columnIndex("TAIR");
    const int dew   = met.columnIndex("TDEW");
    const int wind  = met.columnIndex("WIND");
    const int cloud = met.columnIndex("CLOUD");
    const int phi   = met.columnIndex("PHI");

    if (std::min({air, dew, wind, cloud, phi}) < 0) {
        std::cerr << "meteorology requires TAIR, TDEW, WIND, CLOUD and PHI\n";
        return 1;
    }

    // ── Solver ────────────────────────────────────────────────────────────

    FVQual::Hydro::InternalModeSolver  solver;
    FVQual::Hydro::InternalModeOptions options;
    if (const char *value = std::getenv("FVQUAL_HEAT_AUDIT")) {
        if (std::string(value) != "0" && std::string(value) != "1") {
            std::cerr << "FVQUAL_HEAT_AUDIT must be 0 or 1\n";
            return 1;
        }
        options.auditHeatTransport = std::string(value) == "1";
    }

    options.gravity   = 9.80665;
    options.manningN  = 0.03;
    options.cflNumber = 0.5;
    options.maxTimeStep = maximumStep;
    options.useTurbulenceClosure = closure;
    options.stableDiffusivityScale = scalarNomix ? 0.0 : 1.0;
    options.stableViscosityScale = viscosityNomix ? 0.0 : 1.0;
    options.auditMomentum = true;

    if (nomix)
    {
        // Convective adjustment stays: without it an unstable column is not
        // resolved at all, which is a different experiment.
        options.neutralViscosity      = 0.0;
        options.backgroundDiffusivity = 1.0e-9;
    }
    const bool nowind = argc > 5 && std::string(argv[5]) == "nowind";

    // Diagnostic bisection: run the same geometry and meteorology with the
    // branch inflow and the dam release switched off, to separate a forcing
    // failure from a geometry one.
    const bool nosources = argc > 5 && std::string(argv[5]) == "nosources";


    if (!solver.initialize(layered, state, registry, options, message))
    {
        std::cerr << "solver: " << message << '\n';
        return 1;
    }

    // DeGray's calibrated wind function, from the control file's heat-exchange
    // card: 9.2 + 0.46 W^2 in W2's units. The coefficients here are in SI, so
    // only the shape carries over.
    solver.meteorology().base.windExponent = 2.0;

    // DeGray's control file states WINDH = 10 m, and W2 drives its evaporative
    // wind function from the 2 m wind while the surface stress uses the 10 m
    // one. Applying coefficients calibrated at 2 m to a 10 m wind inflates the
    // quadratic term by about 1.7x, and the model's surface runs 1.9 C cold.
    solver.meteorology().base.windHeight        = 10.0;
    solver.meteorology().base.evaporationHeight = 2.0;
    solver.meteorology().base.windRoughness     = 0.01;
    if (nowind) solver.meteorology().base.dragCoefficient = 0.0;
    solver.enableMeteorology(/*attenuation=*/0.45);

    // Inflow at the head of the reservoir, density-placed so a cold tributary
    // finds its own level; withdrawal at the dam.
    // With more than one lane the head segment is lanes 0..lanes-1 and the dam
    // is the last lanes columns. Attaching the whole discharge to a single lane
    // would drive a jet down one side of the reservoir, so it is split across
    // the lanes in proportion to their surface area — the share of the section
    // each one actually presents to the flow.
    std::vector<int>    inflowIndices;
    std::vector<int>    outflowIndices;
    std::vector<double> inflowShare;
    std::vector<double> outflowShare;

    {
        const std::vector<double> &area = horizontal.faceArea();

        double headArea = 0.0;
        double damArea  = 0.0;

        for (int l = 0; l < headLanes; ++l)
            headArea += area[static_cast<size_t>(l)];

        for (int l = 0; l < lanes; ++l)
            damArea +=
                area[static_cast<size_t>(layered.columnCount() - lanes + l)];

        // Separate loops: with a planned mesh the head and the dam need not have
        // the same number of lanes, and pairing them in one loop would silently
        // attach the release to the wrong columns whenever they differ.
        for (int l = 0; l < headLanes; ++l)
        {
            FVQual::Forcing::Source tributary;
            tributary.column    = l;
            tributary.placement = FVQual::Forcing::Placement::DensityMatched;

            const int index = solver.sources().add(tributary, message);

            if (index < 0)
                break;

            inflowIndices.push_back(index);
            inflowShare.push_back(area[static_cast<size_t>(l)] / headArea);
        }

        for (int l = 0; l < lanes; ++l)
        {
            const int64_t column = layered.columnCount() - lanes + l;

            FVQual::Forcing::Source release;
            release.column    = column;
            release.placement = selectiveOutlet ? FVQual::Forcing::Placement::SelectiveWithdrawal : FVQual::Forcing::Placement::Uniform;
            release.elevation = outletElevation;
            release.withdrawalAngle = 3.14; // DeGray control workbook C168, point sink.

            const int index = solver.sources().add(release, message);

            if (index < 0)
                break;

            outflowIndices.push_back(index);
            outflowShare.push_back(area[static_cast<size_t>(column)] / damArea);
        }
    }

    if (static_cast<int>(inflowIndices.size()) != headLanes ||
        static_cast<int>(outflowIndices.size()) != lanes)
    {
        std::cerr << "sources: " << message << '\n';
        return 1;
    }

    // ── Run ───────────────────────────────────────────────────────────────

    std::ofstream daily(outputs + "/degray_daily.csv");

    daily << "jday,surface_temperature_C,bottom_temperature_C,thermocline_depth_m,"
             "surface_elevation_m,volume_m3,heat_content_J,equilibrium_C,"
             "air_temperature_C,shortwave_W_m2,water_expected_m3,water_independent_residual_m3,source_heat_in_independent_J,source_heat_out_audit_J\n";
    daily.precision(17);

    // The whole profile at the dam, every day. Surface and bottom temperatures
    // say a reservoir stratified; only the profile says whether it stratified
    // the way a reservoir does, and the depth-time contour drawn from this is
    // the figure a limnologist would ask for first. Depth is written per row
    // rather than implied by the layer index, because the sigma grid moves with
    // the surface and a layer is not at a fixed depth.
    // The mesh itself, written once: enough for a viewer to draw the actual
    // cells rather than a schematic. Polygon nodes, because a planned mesh has
    // faces with more than four of them wherever the lane count changes.
    {
        std::ofstream geometry(outputs + "/degray_mesh.csv");

        geometry << "column,node_index,x,y\n";

        const std::vector<int64_t> &offsets = horizontal.faceNodeOffsets();
        const std::vector<int64_t> &nodes   = horizontal.faceNodes();

        for (int64_t c = 0; c < horizontal.faceCount(); ++c)
        {
            for (int64_t i = offsets[static_cast<size_t>(c)];
                 i < offsets[static_cast<size_t>(c) + 1]; ++i)
            {
                const size_t n = static_cast<size_t>(nodes[static_cast<size_t>(i)]);

                geometry << c << ',' << (i - offsets[static_cast<size_t>(c)]) << ','
                         << horizontal.nodeX()[n] << ',' << horizontal.nodeY()[n]
                         << '\n';
            }
        }
    }

    // The whole three-dimensional field, one row per cell per report. This is
    // what a slice-and-animate view needs and what the single-column profile
    // cannot provide: the dam profile is depth against time, not a field.
    std::ofstream field(outputs + "/degray_field.csv");

    field << "jday,column,layer,top_m,bottom_m,temperature_C\n";
    field.precision(17);

    std::ofstream profiles(outputs + "/degray_profiles.csv");

    profiles << "jday,layer,depth_m,elevation_m,temperature_C\n";
    profiles.precision(17);

    const double capacity = options.referenceDensity * options.heatCapacity;

    double day       = kStartDay;
    double nextWrite = kStartDay;

    int64_t steps       = 0;
    bool    warnedFloor = false;

    // Per-internal-step trace of one column, for the window around a failure.
    // Written as CSV rather than printed, so it can be read and plotted after
    // the fact. FVQUAL_TRACE_COLUMN selects the column; the window is the day
    // either side of DeGray's peak wind at 96.75.
    const char *traceColumn = std::getenv("FVQUAL_TRACE_COLUMN");

    const int64_t traced = traceColumn ? std::atol(traceColumn) : -1;

    std::ofstream trace;

    if (traced >= 0)
    {
        trace.open(outputs + "/degray_trace.csv");
        trace << "jday,dt,eta,depth,u,v,ext_dt,ext_requested_dt,ext_steps,"
                 "ext_max_courant,ext_clamped,eta_left,eta_right,wind,"
                 "mesh_eta,mesh_thick_top,mesh_thick_bot,min_thick_anywhere,"
                 "limit_col,limit_layer,limit_speed,limit_vertical\n";
        trace.precision(9);
    }

    std::ofstream budgets(outputs + "/degray_budgets.csv");
    budgets.precision(17);
    bool budgetHeader = false;
    double expectedWater = layered.totalVolume(), heatIn = 0, heatOut = 0;
    std::ofstream heatAudit;
    std::ofstream momentumAudit(outputs + "/degray_momentum.csv");
    momentumAudit.precision(17);
    momentumAudit << "jday,wind_impulse_x_kg_m_s,wind_impulse_y_kg_m_s,bed_impulse_x_kg_m_s,bed_impulse_y_kg_m_s,"
        "vertical_residual_x_kg_m_s,vertical_residual_y_kg_m_s,wind_work_J,bed_dissipation_J,viscosity_dissipation_J,"
        "implicit_dissipation_J,vertical_energy_residual_J,max_relative_impulse_residual,max_relative_energy_residual,"
        "mode_impulse_x_kg_m_s,mode_impulse_y_kg_m_s,mode_energy_J,max_absolute_impulse_residual_kg_m_s,max_absolute_energy_residual_J\n";
    double maximumCellHeatResidual = 0, maximumRelativeCellHeatResidual = 0;
    if (options.auditHeatTransport) {
        heatAudit.open(outputs + "/degray_heat_transport.csv");
        heatAudit.precision(17);
        heatAudit << "jday,column,layer,top_m,bottom_m,volume_m3,heat_J,initial_heat_J,"
                     "horizontal_advection_J,vertical_advection_J,sources_J,"
                     "horizontal_diffusion_J,vertical_diffusion_J,surface_J,other_J,residual_J,"
                     "upward_advection_top_J,upward_advection_bottom_J,"
                     "upward_diffusion_top_J,upward_diffusion_bottom_J,"
                     "u_m_s,v_m_s,upward_relative_top_m_s,kz_top_m2_s,kz_bottom_m2_s,"
                     "baroclinic_x_m_s2,baroclinic_y_m_s2,density_kg_m3,km_top_m2_s,km_bottom_m2_s\n";
    }
    const auto report = [&]() {
        // Report from the deepest column, at the dam, where a profile would be
        // taken.
        // The deepest lane of the dam segment. At one lane this is the dam
        // column and nothing changes; at more, reporting the LAST column would
        // pick an outer lane, which sits on a higher bed and holds a different
        // water mass — so a one-lane and a four-lane run would be compared at
        // different places in the section and the difference read as physics.
        const int64_t deep = deepestDamColumn(layered, lanes);

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
        const FVQual::Hydro::Meteorology reported = solver.meteorology().at(deep);

        daily << day << ',' << surfaceT << ',' << bottomT << ',' << thermocline << ','
              << layered.surfaceElevation()[static_cast<size_t>(deep)] << ','
              << layered.totalVolume() << ',' << heat << ','
              << FVQual::Hydro::equilibriumTemperature(reported) << ','
              << reported.airTemperature << ',' << reported.shortwave << ',' << expectedWater << ','
              << layered.totalVolume() - expectedWater << ',' << heatIn << ',' << heatOut << '\n';

        for (int k = 0; k < kLayers; ++k)
        {
            const double centre = 0.5 * (layered.z(deep, k) + layered.z(deep, k + 1));

            // Every column, not just the dam, so the field can be sliced.
            if (k == 0)
            {
                for (int64_t c = 0; c < layered.columnCount(); ++c)
                    for (int kk = 0; kk < kLayers; ++kk)
                        field << day << ',' << c << ',' << kk << ','
                              << layered.z(c, kk) << ',' << layered.z(c, kk + 1) << ','
                              << state.speciesAt(temperature, layered.cell(c, kk))
                              << '\n';
            }

            profiles << day << ',' << k << ',' << layered.z(deep, 0) - centre << ','
                     << centre << ','
                     << state.speciesAt(temperature, layered.cell(deep, k)) << '\n';
        }
        const auto values = FVQual::Forcing::processBudgetValues(solver, registry);
        if (options.auditHeatTransport) {
            const auto &audit = solver.heatTransportAudit();
            for (int64_t c = 0; c < layered.columnCount(); ++c)
                for (int k = 0; k < kLayers; ++k) {
                    const size_t cell = layered.cell(c,k), top = layered.interfaceSlot(c,k),
                                 bottom = layered.interfaceSlot(c,k+1);
                    const double heat = capacity * layered.cellVolume()[cell] * state.speciesAt(temperature,cell);
                    const double expected = audit.initial[cell] + audit.horizontalAdvection[cell] +
                        audit.verticalAdvection[cell] + audit.sources[cell] + audit.horizontalDiffusion[cell] +
                        audit.verticalDiffusion[cell] + audit.surface[cell] + audit.other[cell];
                    const double residual = heat - expected;
                    maximumCellHeatResidual = std::max(maximumCellHeatResidual,std::abs(residual));
                    maximumRelativeCellHeatResidual = std::max(maximumRelativeCellHeatResidual,
                        std::abs(residual)/std::max({1.0,std::abs(heat),std::abs(audit.initial[cell])}));
                    // Grid-relative upward velocity uses the actual wet area;
                    // verticalFlux() itself is normalized by plan-polygon area.
                    const double area = layered.interfaceArea()[top];
                    const double w = area > 0 ? solver.verticalFlux()[top] * horizontal.faceArea()[c] / area : 0;
                    heatAudit << day << ',' << c << ',' << k << ',' << layered.z(c,k) << ',' << layered.z(c,k+1)
                        << ',' << layered.cellVolume()[cell] << ',' << heat << ',' << audit.initial[cell]
                        << ',' << audit.horizontalAdvection[cell] << ',' << audit.verticalAdvection[cell]
                        << ',' << audit.sources[cell] << ',' << audit.horizontalDiffusion[cell]
                        << ',' << audit.verticalDiffusion[cell] << ',' << audit.surface[cell]
                        << ',' << audit.other[cell] << ',' << residual
                        << ',' << audit.upwardAdvection[top] << ',' << audit.upwardAdvection[bottom]
                        << ',' << audit.upwardDiffusion[top] << ',' << audit.upwardDiffusion[bottom]
                        << ',' << state.u[cell] << ',' << state.v[cell] << ',' << w
                        << ',' << solver.eddyDiffusivity()[top] << ',' << solver.eddyDiffusivity()[bottom]
                        << ',' << solver.baroclinicAccelerationX()[cell] << ',' << solver.baroclinicAccelerationY()[cell]
                        << ',' << state.density[cell]
                        << ',' << solver.eddyViscosity()[top] << ',' << solver.eddyViscosity()[bottom] << '\n';
                }
            heatAudit.flush();
        }
        if (!budgetHeader) {
            budgets << "jday";
            for (const auto &[id, value] : values.items()) budgets << ',' << id;
            budgets << '\n';
            budgetHeader = true;
        }
        budgets << day;
        for (const auto &[id, value] : values.items()) budgets << ',' << value.get<double>();
        budgets << '\n';
        daily.flush(); profiles.flush(); field.flush(); budgets.flush();
        FVQual::Hydro::MomentumColumnAudit momentum;
        for (const auto &a : solver.momentumAudit()) {
            momentum.windImpulseX+=a.windImpulseX; momentum.windImpulseY+=a.windImpulseY;
            momentum.bedImpulseX+=a.bedImpulseX; momentum.bedImpulseY+=a.bedImpulseY;
            momentum.verticalResidualX+=a.verticalResidualX; momentum.verticalResidualY+=a.verticalResidualY;
            momentum.windWork+=a.windWork; momentum.bedDissipation+=a.bedDissipation;
            momentum.viscosityDissipation+=a.viscosityDissipation; momentum.implicitDissipation+=a.implicitDissipation;
            momentum.verticalEnergyResidual+=a.verticalEnergyResidual;
            momentum.maxRelativeImpulseResidual=std::max(momentum.maxRelativeImpulseResidual,a.maxRelativeImpulseResidual);
            momentum.maxRelativeEnergyResidual=std::max(momentum.maxRelativeEnergyResidual,a.maxRelativeEnergyResidual);
            momentum.maxAbsoluteImpulseResidual=std::max(momentum.maxAbsoluteImpulseResidual,a.maxAbsoluteImpulseResidual);
            momentum.maxAbsoluteEnergyResidual=std::max(momentum.maxAbsoluteEnergyResidual,a.maxAbsoluteEnergyResidual);
            momentum.modeImpulseX+=a.modeImpulseX; momentum.modeImpulseY+=a.modeImpulseY; momentum.modeEnergy+=a.modeEnergy;
        }
        momentumAudit << day << ',' << momentum.windImpulseX << ',' << momentum.windImpulseY << ','
            << momentum.bedImpulseX << ',' << momentum.bedImpulseY << ',' << momentum.verticalResidualX << ','
            << momentum.verticalResidualY << ',' << momentum.windWork << ',' << momentum.bedDissipation << ','
            << momentum.viscosityDissipation << ',' << momentum.implicitDissipation << ',' << momentum.verticalEnergyResidual << ','
            << momentum.maxRelativeImpulseResidual << ',' << momentum.maxRelativeEnergyResidual << ','
            << momentum.modeImpulseX << ',' << momentum.modeImpulseY << ',' << momentum.modeEnergy << ','
            << momentum.maxAbsoluteImpulseResidual << ',' << momentum.maxAbsoluteEnergyResidual << '\n';
        momentumAudit.flush();
        std::printf("report day %.5f (%lld steps)\n", day, static_cast<long long>(steps));
        std::fflush(stdout);
    };
    // Baseline weather for the initial report; each integration step refreshes it.
    solver.meteorology().airTemperature.setUniform(met.valueAt(day, air));
    solver.meteorology().dewPoint.setUniform(met.valueAt(day, dew));
    solver.meteorology().shortwave.setUniform(FVQual::Hydro::clearSkySolarRadiation(
        day, kLatitude, met.valueAt(day, cloud) / 10.0));
    if (!daily || !profiles || !field || !budgets || (options.auditHeatTransport && !heatAudit)) {
        std::cerr << "cannot open DeGray output files\n";
        return 1;
    }
    report();
    nextWrite += 1.0;
    while (day < endDay - 1e-12)

    {
        FVQual::Forcing::MeteorologyField &weather = solver.meteorology();

        weather.airTemperature.setUniform(met.valueAt(day, air));
        weather.dewPoint.setUniform(met.valueAt(day, dew));
        // Sheltering is per segment, so it varies across the reservoir; applying
        // it to the speed rather than to the stress keeps the heat budget and
        // the momentum seeing the same wind, which is the failure mode that
        // once had wind evaporating water without moving it.
        if (hasSheltering)
        {
            std::vector<double> sheltered(static_cast<size_t>(layered.columnCount()));

            for (int64_t c = 0; c < layered.columnCount(); ++c)
            {
                const int segment =
                    std::min(segmentOfColumn[static_cast<size_t>(c)],
                             sheltering.columnCount() - 1);

                sheltered[static_cast<size_t>(c)] =
                    met.valueAt(day, wind) * sheltering.valueAt(day, segment);
            }

            std::string local;

            if (!weather.windSpeed.setPerColumn(sheltered, layered.columnCount(), local))
            {
                std::cerr << "wind sheltering: " << local << '\n';
                return 1;
            }
        }
        else
        {
            weather.windSpeed.setUniform(met.valueAt(day, wind));
        }

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

        const double qin  = nosources ? 0.0 : inflow.valueAt(day, 0);
        const double qout = nosources ? 0.0 : outflow.valueAt(day, 0);

        bool sourcesOk = true;

        for (int l = 0; l < headLanes && sourcesOk; ++l)
            sourcesOk = solver.sources().setDischarge(
                inflowIndices[static_cast<size_t>(l)],
                qin * inflowShare[static_cast<size_t>(l)], message);

        for (int l = 0; l < lanes && sourcesOk; ++l)
            sourcesOk = solver.sources().setDischarge(
                outflowIndices[static_cast<size_t>(l)],
                -qout * outflowShare[static_cast<size_t>(l)], message);

        if (!sourcesOk)
        {
            std::cerr << "discharge: " << message << '\n';
            return 1;
        }

        // Every lane, not just the first. Lanes that carried discharge without
        // a temperature would have injected water at the default and quietly
        // cooled the head of the reservoir in proportion to the lane count.
        const double inflowT = inflowTemperature.valueAt(day, 0);

        for (int l = 0; l < headLanes; ++l)
        {
            if (!solver.sources().setInflowValue(inflowIndices[static_cast<size_t>(l)],
                                                 temperature, temperature, inflowT,
                                                 message))
            {
                std::cerr << "inflow temperature: " << message << '\n';
                return 1;
            }
        }

        // Split at forcing knots and reports. W2 forcing is linear in Julian
        // day; sampling at each step's opening preserves that convention and
        // refinement measures its temporal integration error explicitly.
        double nextKnot = endDay;
        for (const auto *series : {&met, &inflow, &inflowTemperature, &outflow, &sheltering}) {
            const auto found = std::upper_bound(series->time().begin(), series->time().end(), day + 1e-10);
            if (found != series->time().end()) nextKnot = std::min(nextKnot, *found);
        }
        const double targetDay = std::min({endDay, nextWrite, nextKnot});
        const double dt = std::min({solver.stableTimeStep(), maximumStep, (targetDay - day) * 86400.0});

        double withdrawalTemperature = 0;
        for (int l = 0; l < lanes; ++l) {
            const int64_t column = layered.columnCount() - lanes + l;
            double weighted = 0, thickness = 0;
            for (int k = 0; k < kLayers; ++k) {
                const size_t cell = static_cast<size_t>(layered.cell(column, k));
                weighted += layered.thickness()[cell] * state.speciesAt(temperature, cell);
                thickness += layered.thickness()[cell];
            }
            if (selectiveOutlet) {
                // Probe the stated point sink on the opening geometry. This
                // uses the shared placement formula, so it is an audit rather
                // than an independent oracle for the selective envelope.
                withdrawalTemperature += outflowShare[l] * solver.sources().withdrawalMean(
                    layered, outflowIndices[l], -qout * outflowShare[l], state.density,
                    options.gravity, options.referenceDensity, state.speciesBlock(temperature));
            } else withdrawalTemperature += outflowShare[l] * weighted / thickness;
        }
        expectedWater += (qin - qout) * dt;
        heatIn += qin * inflowT * capacity * dt;
        heatOut += qout * withdrawalTemperature * capacity * dt;
        if (!solver.step(dt, message))
        {
            std::cerr << "step failed on day " << day << ": " << message << '\n';
            return 1;
        }

        day += dt / 86400.0;
        if (std::abs(day - targetDay) < 1e-10) day = targetDay;
        ++steps;

        // A Courant violation is silent unless something says so. DeGray's
        // narrowest segment tapers to 0.8% of its surface width at the bed, so
        // the storage-limited step for that column falls below the floor and the
        // solver clamps up to it — taking a step the condition forbids.
        if (traced >= 0 && day > 96.0 && day < 97.2)
        {
            const auto  &ext = solver.external();
            const size_t i   = static_cast<size_t>(traced);

            const auto &eta = ext.surfaceElevation();

            trace << day << ',' << dt << ',' << eta[i] << ','
                  << ext.depth()[i] << ',' << ext.velocityX()[i] << ','
                  << ext.velocityY()[i] << ',' << ext.diagnostics().timeStep << ','
                  << ext.diagnostics().requestedTimeStep << ','
                  << ext.diagnostics().steps << ','
                  << ext.diagnostics().maxCourant << ','
                  << ext.diagnostics().clampedSteps << ','
                  << (traced > 0 ? eta[i - 1] : eta[i]) << ','
                  << (traced + 1 < layered.columnCount() ? eta[i + 1] : eta[i]) << ','
                  << weather.windSpeed.at(0) << ',';

            // The mesh's own view of the same column, and the thinnest layer
            // anywhere: the internal step is limited by thickness/|w|, so a
            // vanishing layer in *any* column collapses it globally.
            double thinnest = 1.0e30;

            for (int64_t cc = 0; cc < layered.columnCount(); ++cc)
                for (int k = 0; k < kLayers; ++k)
                    thinnest = std::min(
                        thinnest,
                        layered.thickness()[static_cast<size_t>(layered.cell(cc, k))]);

            trace << layered.surfaceElevation()[i] << ','
                  << layered.thickness()[static_cast<size_t>(layered.cell(traced, 0))]
                  << ','
                  << layered.thickness()[static_cast<size_t>(
                         layered.cell(traced, kLayers - 1))]
                  << ',' << thinnest << ','
                  << solver.diagnostics().limitingColumn << ','
                  << solver.diagnostics().limitingLayer << ','
                  << solver.diagnostics().limitingSpeed << ','
                  << (solver.diagnostics().limitIsVertical ? 1 : 0) << '\n';
        }

        if (solver.external().timeStepFloorWasBinding() && !warnedFloor)
        {
            std::cerr << "warning: the barotropic time-step floor bound at day "
                      << day << "; the Courant condition was violated and the "
                         "results after this point are not trustworthy\n";

            warnedFloor = true;
        }

        if (day < nextWrite - 1e-10 && day < endDay - 1e-10)
            continue;

        // Advance the schedule, do not rebase it on the time actually reached.
        // Rebasing carries each step's overshoot forward, and with a 300 s cap
        // that is up to five minutes a day: the report drifted from 12:05 in
        // March to 17:00 in December, sampling a steadily later hour of a
        // diurnal cycle and aliasing it into the seasonal signal. Found by
        // checking a figure's caption against its data, which is not where one
        // expects to find a timing bug.
        nextWrite += 1.0;

        report();

    }

    daily.close();
    profiles.close();

    // Final profile at the dam, which is what a thermal survey would report.
    std::ofstream final(outputs + "/degray_final_profile.csv");

    final << "depth_m,elevation_m,temperature_C\n";
    final.precision(17);

    const int64_t deep = deepestDamColumn(layered, lanes);

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

    std::ofstream summary(outputs + "/driver.json");
    summary << nlohmann::json({{"start_jday", kStartDay}, {"end_jday", day},
        {"columns", layered.columnCount()}, {"layers", kLayers}, {"dam_column", deep},
        {"outlet_placement", selectiveOutlet ? "selective-point" : "uniform"},
        {"outlet_elevation_m", selectiveOutlet ? nlohmann::json(outletElevation) : nlohmann::json(nullptr)},
        {"withdrawal_heat_audit", selectiveOutlet ? "shared placement probe; not independent envelope" : "independent thickness-weighted uniform mean"},
        {"maximum_step_s", maximumStep}, {"steps", steps}, {"wind_stress", !nowind},
        {"heat_transport_audit", options.auditHeatTransport},
        {"maximum_cell_heat_residual_J", maximumCellHeatResidual},
        {"maximum_relative_cell_heat_residual", maximumRelativeCellHeatResidual},
        {"shear_mixing_reduced", nomix}, {"convective_mixing_retained", true},
        {"stable_scalar_mixing_scale",options.stableDiffusivityScale},
        {"stable_viscosity_scale",options.stableViscosityScale},
        {"momentum_audit_scope","implicit vertical solve; mode correction measured separately; not full mechanical energy closure"},
        {"closure", closure ? "kepsilon" : "richardson"}, {"hypsography", !nohyps},
        {"barotropic_floor_binding", warnedFloor}, {"forcing_clock", "1980 W2 Julian day; linear, opening-step sample"},
        {"geometry_scope", "laterally averaged W2 cross sections; lateral lanes are inferred, not surveyed"}}).dump(2);
    return (daily && profiles && field && budgets && final && summary && momentumAudit &&
            (!options.auditHeatTransport || heatAudit)) ? 0 : 1;
}
