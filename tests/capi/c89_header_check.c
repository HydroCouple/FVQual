/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Caleb Buahin
 *
 * c89_header_check.c — compiles the public C API headers under strict C89.
 *
 * The ABI firewall is a promise about what compilers can consume these
 * headers. A promise nothing checks is a comment, so this translation unit is
 * built with -std=c89 -pedantic-errors and does nothing else: if a C++-ism,
 * a // comment, or a C99 declaration ever creeps into a public header, the
 * build stops here rather than in a downstream user's toolchain.
 */

#include "fvqual/capi/fvqual_engine.h"
#include "fvqual/capi/fvqual_mesh.h"
#include "fvqual/capi/fvqual_options.h"
#include "fvqual/capi/fvqual_sigma.h"
#include "fvqual/capi/fvqual_sim.h"
#include "fvqual/capi/fvqual_coupling.h"
#include "fvqual/capi/fvqual_differentiation.h"
#include "fvqual/capi/fvqual_column.h"
#include "fvqual/capi/fvqual_reactive.h"
#include "fvqual/capi/fvqual_configuration.h"

/* Reference one symbol of each kind so the declarations are actually parsed
   rather than merely tokenized. */
int fvqual_c89_header_check(void);

int fvqual_c89_header_check(void)
{
    FVQ_SigmaSpec sigma;
    FVQ_Options   options;
    FVQ_ErrorCode code;

    fvq_sigma_default(&sigma);
    fvq_options_default(&options);
    code = fvq_sigma_validate(&sigma);

    return (int)code;
}
