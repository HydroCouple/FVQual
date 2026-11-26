/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Caleb Buahin
 * Conservative, synchronized one-way inflow and same-instance replay.
 */
#ifndef FVQUAL_CAPI_COUPLING_H
#define FVQUAL_CAPI_COUPLING_H
#include "fvqual/capi/fvqual_sim.h"
#ifdef __cplusplus
extern "C" {
#endif
typedef enum FVQ_InflowPlacement {
    FVQ_INFLOW_UNIFORM = 0, FVQ_INFLOW_SURFACE, FVQ_INFLOW_BED,
    FVQ_INFLOW_ELEVATION, FVQ_INFLOW_DENSITY_MATCHED
} FVQ_InflowPlacement;
typedef struct FVQ_InflowSpecies {
    const char *id;
    double massKg;
} FVQ_InflowSpecies;
typedef struct FVQ_InflowTransfer {
    int64_t column;
    double volumeM3;
    double heatJ; /* rho*cp*V*T with reference T=0 degC, not an extra heat source. */
    FVQ_InflowPlacement placement;
    double elevationM;
    const FVQ_InflowSpecies *species;
    int64_t speciesCount;
} FVQ_InflowTransfer;
typedef struct FVQ_TransferReceipt {
    double requestedVolumeM3, acceptedVolumeM3;
    double requestedHeatJ, acceptedHeatJ;
    int replayed;
} FVQ_TransferReceipt;

/*! Executes an interval from the current elapsed time to endSeconds.
 * Species must be transported kg/m3 quantities. Declare every incoming species
 * including zero loads. Same ID + identical contents returns the prior receipt
 * without advancing. Changed contents or negative water are refused.
 * Failed numerical intervals roll back. SWMM already owns the donor debit.
 */
FVQual_CAPI FVQ_ErrorCode fvq_sim_apply_inflow_frame(
    FVQ_Simulation *sim, const char *transferId,
    double startSeconds, double endSeconds,
    const FVQ_InflowTransfer *transfers, int64_t count,
    FVQ_TransferReceipt *receipt);
/*! Accepted mass for a named species in an accepted transfer, kg. */
FVQual_CAPI FVQ_ErrorCode fvq_sim_transfer_species_receipt(
    const FVQ_Simulation *sim, const char *transferId,
    const char *speciesId, double *acceptedKg);
/*! In-memory checkpoints belong to one initialized simulation/configuration.
 * They include the receipt ledger and can be restored repeatedly until released.
 * Restore preserves allocation-stable state buffers, but invalidates geometry
 * views. These tokens are not durable restart files.
 */
FVQual_CAPI FVQ_ErrorCode fvq_sim_save_state(FVQ_Simulation *sim, uint64_t *token);
FVQual_CAPI FVQ_ErrorCode fvq_sim_restore_state(FVQ_Simulation *sim, uint64_t token);
FVQual_CAPI FVQ_ErrorCode fvq_sim_release_state(FVQ_Simulation *sim, uint64_t token);
#ifdef __cplusplus
}
#endif
#endif
