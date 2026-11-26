/* SPDX-License-Identifier: MIT */
#ifndef FVQUAL_CAPI_CONFIGURATION_H
#define FVQUAL_CAPI_CONFIGURATION_H
#include "fvqual/capi/fvqual_sim.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Before initialize: one JSON object with meteorology, sources, structures,
 * outlet_groups, open_boundaries, internal_transfers and optional reaction_deck
 * text. SI operation units; species concentrations use registered units.
 * Schedules hold values between knots in elapsed seconds from zero.
 * Boundary node pairs and reservoir columns are zero-based mesh indices.
 * Validates before replacing any configuration; concentrations must refer to
 * existing species or species declared in the optional reaction deck.
 */
FVQual_CAPI FVQ_ErrorCode fvq_sim_configure_json(FVQ_Simulation *sim, const char *json);
/* Named cumulative process amount, inventory basis or conservation residual.
 * Names match component budget_* outputs; heat J, water m3, known species kg.
 * Other species retain declared concentration*m3 (benthic: *m2) units. */
FVQual_CAPI FVQ_ErrorCode fvq_sim_budget_value(const FVQ_Simulation *sim,
    const char *name, double *value);
#ifdef __cplusplus
}
#endif
#endif
