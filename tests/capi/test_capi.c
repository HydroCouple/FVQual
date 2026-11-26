/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Caleb Buahin
 *
 * test_capi.c — the PY0 gate.
 *
 * Deliberately written in C and compiled by a C compiler. The C API's whole
 * purpose is to be an ABI firewall, and a firewall that is only ever exercised
 * from C++ is not a firewall — it is a C++ header with unusual punctuation.
 * If this file stops compiling, the boundary has leaked.
 *
 * The physics check mirrors the C++ barotropic-seiche gate: a closed basin
 * released from a cosine surface must slosh and must conserve mass.
 */

#include "fvqual/capi/fvqual_engine.h"
#include "fvqual/capi/fvqual_mesh.h"
#include "fvqual/capi/fvqual_options.h"
#include "fvqual/capi/fvqual_sigma.h"
#include "fvqual/capi/fvqual_sim.h"
#include "fvqual/capi/fvqual_coupling.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979323846

static int failures = 0;

static void check(int condition, const char *what)
{
    if (condition)
    {
        printf("  ok    %s\n", what);
    }
    else
    {
        printf("  FAIL  %s\n", what);
        ++failures;
    }
}

static void checkCall(FVQ_ErrorCode code, const char *what)
{
    if (code == FVQ_OK)
    {
        printf("  ok    %s\n", what);
    }
    else
    {
        printf("  FAIL  %s: %s (%s)\n", what, fvq_error_message(code), fvq_last_error());
        ++failures;
    }
}

/* ── Error and version surface ─────────────────────────────────────────── */

static void testDiagnosticsSurface(void)
{
    int major = -1, minor = -1, patch = -1;

    printf("error and version surface\n");

    check(fvq_error_message(FVQ_OK) != NULL, "error messages are available");
    check(fvq_error_message(FVQ_ERR_DRY_COLUMN) != NULL, "dry-column message exists");
    check(fvq_version() != NULL && fvq_version()[0] != '\0', "version string is set");

    fvq_version_numbers(&major, &minor, &patch);
    check(major >= 0 && minor >= 0 && patch >= 0, "version numbers are readable");

    /* NULL arguments must be rejected, not crash: a binding will pass NULL. */
    check(fvq_mesh_node_count(NULL) == 0, "null mesh reports zero nodes");
    check(fvq_sim_time(NULL) == 0.0, "null simulation reports zero time");
    fvq_mesh_destroy(NULL);
    fvq_sim_destroy(NULL);
    check(1, "destroying null handles is safe");
}

/* ── Mesh construction ─────────────────────────────────────────────────── */

static void testMesh(void)
{
    FVQ_Mesh *mesh = NULL;
    double    area = 0.0;
    double   *bed  = NULL;
    int64_t   faces = 0;
    int64_t   i;

    printf("mesh construction\n");

    checkCall(fvq_mesh_create_channel(40, 400.0, 20.0, -10.0, &mesh),
              "channel mesh builds");

    faces = fvq_mesh_face_count(mesh);
    check(faces == 40, "face count matches the request");
    check(fvq_mesh_node_count(mesh) == 82, "node count is (columns + 1) * 2");

    area = fvq_mesh_total_area(mesh);
    check(fabs(area - 400.0 * 20.0) < 1e-9, "total area matches length x width");
    check(fabs(fvq_mesh_minimum_spacing(mesh) - 10.0) < 1e-9, "spacing is the cell size");

    bed = (double *)malloc((size_t)faces * sizeof(double));
    check(bed != NULL, "allocation succeeds");

    if (bed != NULL)
    {
        checkCall(fvq_mesh_face_geometry(mesh, NULL, NULL, NULL, bed, faces),
                  "face geometry copies out");
        check(fabs(bed[0] + 10.0) < 1e-12, "bed elevation carried through");

        for (i = 0; i < faces; ++i)
            bed[i] = -20.0;

        checkCall(fvq_mesh_set_bed(mesh, bed, faces), "bed is settable");
        checkCall(fvq_mesh_face_geometry(mesh, NULL, NULL, NULL, bed, faces),
                  "geometry re-reads");
        check(fabs(bed[10] + 20.0) < 1e-12, "new bed took effect");

        /* A wrong-length buffer must be refused rather than overrun. */
        check(fvq_mesh_set_bed(mesh, bed, faces + 1) == FVQ_ERR_INVALID_ARGUMENT,
              "mismatched buffer length is rejected");

        free(bed);
    }

    fvq_mesh_destroy(mesh);
}

/* ── Lifecycle ordering ────────────────────────────────────────────────── */

static void testLifecycle(void)
{
    FVQ_Mesh       *mesh = NULL;
    FVQ_Simulation *sim  = NULL;
    FVQ_SigmaSpec   sigma;
    int             index = -1;

    printf("lifecycle ordering\n");

    fvq_sigma_default(&sigma);
    sigma.layerCount        = 6;
    sigma.minLayerThickness = 0.0;

    checkCall(fvq_sigma_validate(&sigma), "default sigma validates");

    sigma.layerCount = 0;
    check(fvq_sigma_validate(&sigma) == FVQ_ERR_SIGMA, "zero layers is rejected");
    sigma.layerCount = 6;

    checkCall(fvq_mesh_create_channel(10, 100.0, 20.0, -8.0, &mesh), "mesh builds");
    checkCall(fvq_sim_create(mesh, &sigma, NULL, &sim), "simulation is created");

    /* The mesh was copied, so releasing it now must be harmless — this is the
     * lifetime property the C API promises. */
    fvq_mesh_destroy(mesh);
    mesh = NULL;

    check(fvq_sim_step(sim, 1.0) == FVQ_ERR_LIFECYCLE,
          "stepping before initialize is refused");

    checkCall(fvq_sim_add_temperature(sim, 12.0, &index), "temperature registers");
    check(index == 0, "temperature is the first species");
    check(fvq_sim_add_temperature(sim, 15.0, NULL) == FVQ_ERR_SPECIES,
          "duplicate temperature is refused");

    checkCall(fvq_sim_add_species(sim, "DYE", "g/m3", 1.0, &index), "a species registers");
    check(index == 1, "species indices increment");
    check(fvq_sim_species_count(sim) == 2, "species count is right");

    checkCall(fvq_sim_initialize(sim), "initialize succeeds");

    check(fvq_sim_add_species(sim, "LATE", NULL, 0.0, NULL) == FVQ_ERR_LIFECYCLE,
          "registering after initialize is refused");
    check(fvq_sim_initialize(sim) == FVQ_ERR_LIFECYCLE, "double initialize is refused");

    check(fvq_sim_column_count(sim) == 10, "column count is right");
    check(fvq_sim_layer_count(sim) == 6, "layer count is right");
    check(fvq_sim_cell_count(sim) == 60, "cell count is columns x layers");

    fvq_sim_destroy(sim);
}

/* ── A dry column must be refused, not divided by ──────────────────────── */

static void testDryColumnRefused(void)
{
    FVQ_Mesh       *mesh = NULL;
    FVQ_Simulation *sim  = NULL;
    FVQ_SigmaSpec   sigma;
    FVQ_ErrorCode   code = FVQ_OK;

    printf("dry column handling\n");

    fvq_sigma_default(&sigma);
    sigma.layerCount        = 4;
    sigma.minLayerThickness = 0.0;

    /* Bed above the default surface of zero. */
    checkCall(fvq_mesh_create_channel(5, 50.0, 10.0, 2.0, &mesh), "emergent mesh builds");
    checkCall(fvq_sim_create(mesh, &sigma, NULL, &sim), "simulation is created");

    code = fvq_sim_initialize(sim);
    check(code == FVQ_ERR_DRY_COLUMN, "a dry column is refused with its own code");
    check(fvq_last_error()[0] != '\0', "the failure carries a detail message");

    fvq_sim_destroy(sim);
    fvq_mesh_destroy(mesh);
}

/* ── The physics gate: a barotropic seiche ─────────────────────────────── */

static void testSeiche(void)
{
    FVQ_Mesh       *mesh    = NULL;
    FVQ_Simulation *sim     = NULL;
    FVQ_SigmaSpec   sigma;
    FVQ_Options     options;
    FVQ_Diagnostics diagnostics;

    const int    columns   = 80;
    const double length    = 1000.0;
    const double depth     = 10.0;
    const double amplitude = 0.05;
    const double period    = 2.0 * length / sqrt(9.80665 * depth);

    double *surface = NULL;
    double  before  = 0.0;
    double  after   = 0.0;
    double  drift   = 0.0;
    int     i;

    printf("barotropic seiche\n");

    fvq_sigma_default(&sigma);
    sigma.layerCount        = 4;
    sigma.minLayerThickness = 0.0;

    fvq_options_default(&options);
    options.manningN = 0.0;

    checkCall(fvq_mesh_create_channel(columns, length, 20.0, -depth, &mesh),
              "mesh builds");
    checkCall(fvq_sim_create(mesh, &sigma, &options, &sim), "simulation is created");

    surface = (double *)malloc((size_t)columns * sizeof(double));
    check(surface != NULL, "allocation succeeds");

    if (surface != NULL)
    {
        for (i = 0; i < columns; ++i)
        {
            const double x = ((double)i + 0.5) * length / (double)columns;
            surface[i]     = amplitude * cos(PI * x / length);
        }

        checkCall(fvq_sim_set_initial_surface(sim, surface, columns),
                  "initial surface is accepted");
        free(surface);
    }

    checkCall(fvq_sim_initialize(sim), "initialize succeeds");

    before = fvq_sim_total_volume(sim);
    check(before > 0.0, "the basin holds water");

    check(fvq_sim_stable_time_step(sim) > 0.0, "a stable step is reported");

    /* Run a quarter period, where the surface should have flattened and the
     * flow should be at its strongest. */
    checkCall(fvq_sim_advance(sim, 0.25 * period), "a quarter period advances");

    after = fvq_sim_total_volume(sim);
    drift = fabs(after - before) / before;

    printf("        volume drift %.3e over %.0f s\n", drift, fvq_sim_time(sim));
    check(drift < 1e-10, "mass closes below 1e-10");

    checkCall(fvq_sim_diagnostics(sim, &diagnostics), "diagnostics are readable");
    printf("        max |u| = %.4f m/s, continuity residual = %.3e\n",
           diagnostics.maxHorizontalVelocity, diagnostics.surfaceFluxResidual);

    check(diagnostics.maxHorizontalVelocity > 1e-3, "the basin is actually sloshing");
    check(diagnostics.surfaceFluxResidual < 1e-10, "the two modes stay consistent");
    check(fabs(fvq_sim_time(sim) - 0.25 * period) < 1e-9, "the clock lands on target");

    fvq_sim_destroy(sim);
    fvq_mesh_destroy(mesh);
}


static void testCouplingReplay(void)
{
    FVQ_Mesh *mesh = NULL;
    FVQ_Simulation *sim = NULL;
    FVQ_SigmaSpec sigma;
    FVQ_InflowSpecies species = {"tracer", 0.0002};
    FVQ_InflowTransfer transfer = {0, 0.001, 50232.0, FVQ_INFLOW_UNIFORM, 0.0, &species, 1};
    FVQ_TransferReceipt receipt;
    uint64_t token = 0;
    double before = 0.0, mass = 0.0;
    int tracer = -1;
    printf("conservative inflow and complete replay\n");
    fvq_sigma_default(&sigma);
    sigma.layerCount = 2;
    sigma.minLayerThickness = 0.0;
    checkCall(fvq_mesh_create_channel(2, 10.0, 2.0, -5.0, &mesh), "coupling mesh builds");
    checkCall(fvq_sim_create(mesh, &sigma, NULL, &sim), "coupling simulation builds");
    checkCall(fvq_sim_add_temperature(sim, 10.0, NULL), "temperature is registered");
    checkCall(fvq_sim_add_species(sim, "tracer", "kg/m3", 0.0, &tracer), "tracer is registered");
    checkCall(fvq_sim_initialize(sim), "receiver initializes");
    before = fvq_sim_total_volume(sim);
    checkCall(fvq_sim_save_state(sim, &token), "full checkpoint is saved");
    checkCall(fvq_sim_apply_inflow_frame(sim, "OF1/1", 0.0, 0.1, &transfer, 1, &receipt),
              "one-way frame is accepted");
    check(fabs(fvq_sim_total_volume(sim) - before - 0.001) < 1e-10, "accepted water enters reservoir");
    checkCall(fvq_sim_transfer_species_receipt(sim, "OF1/1", "tracer", &mass), "species receipt is readable");
    check(fabs(mass - 0.0002) < 1e-12, "species receipt closes");
    checkCall(fvq_sim_apply_inflow_frame(sim, "OF1/1", 0.0, 0.1, &transfer, 1, &receipt),
              "identical retry succeeds");
    check(receipt.replayed && fabs(fvq_sim_time(sim) - 0.1) < 1e-12, "retry does not advance or debit again");
    transfer.volumeM3 = 0.002;
    check(fvq_sim_apply_inflow_frame(sim, "OF1/1", 0.0, 0.1, &transfer, 1, &receipt) != FVQ_OK,
          "changed retry is refused");
    transfer.volumeM3 = 0.001;
    checkCall(fvq_sim_restore_state(sim, token), "state and receipt ledger restore together");
    check(fabs(fvq_sim_time(sim)) < 1e-12, "restored clock is zero");
    check(fvq_sim_transfer_species_receipt(sim, "OF1/1", "tracer", &mass) != FVQ_OK,
          "future receipt is removed on restore");
    checkCall(fvq_sim_apply_inflow_frame(sim, "OF1/1", 0.0, 0.1, &transfer, 1, &receipt),
              "restored interval can be replayed");
    check(!receipt.replayed, "restored interval is freshly applied");
    checkCall(fvq_sim_release_state(sim, token), "checkpoint is released");
    check(fvq_sim_release_state(sim, token) != FVQ_OK, "double release is refused");
    fvq_sim_destroy(sim);
    fvq_mesh_destroy(mesh);
}

int main(void)
{
    printf("FVQual C API — version %s\n\n", fvq_version());

    testDiagnosticsSurface();
    testMesh();
    testLifecycle();
    testDryColumnRefused();
    testSeiche();
    testCouplingReplay();

    printf("\n%s\n", failures == 0 ? "all checks passed" : "FAILURES PRESENT");
    return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
