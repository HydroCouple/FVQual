// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Caleb Buahin

/*!
 * \file   parallel.h
 * \author Caleb Buahin
 * \brief  Threading pragmas, wrapped so a serial build stays clean (F7.1).
 *
 * Threading is opt-in at configure time (`FVQual_USE_OPENMP`). Wrapping the
 * pragmas keeps a build without OpenMP free of unknown-pragma warnings under
 * `-Wall` and keeps `#ifdef` out of every loop body.
 *
 * **These are not where the parallelism lives.** A loop only carries one of
 * these macros because it was already written to be race-free — every
 * iteration writing to its own element, gathering rather than scattering. That
 * restructuring is the durable part and is what a Kokkos or GPU backend needs;
 * the pragma is the least interesting line in the file. Adding one to a loop
 * that scatters produces a program that is fast and wrong.
 *
 * \internal Implementation detail of the core library.
 */

#ifndef FVQual_PARALLEL_H
#define FVQual_PARALLEL_H

#ifdef _OPENMP

#include <omp.h>

#define FVQual_PRAGMA(x) _Pragma(#x)

//! Parallel over an index range whose iterations touch disjoint memory.
#define FVQual_PARALLEL_FOR FVQual_PRAGMA(omp parallel for schedule(static))

//! As above, with one or more reduction items, e.g. `(+ : total, count)`.
#define FVQual_PARALLEL_FOR_REDUCE(...)                                                  \
    FVQual_PRAGMA(omp parallel for schedule(static) reduction(__VA_ARGS__))

//! Guards a single scattered update. Rare by design; see BarotropicSolver.
#define FVQual_ATOMIC FVQual_PRAGMA(omp atomic)

#else

#define FVQual_PARALLEL_FOR
#define FVQual_PARALLEL_FOR_REDUCE(...)
#define FVQual_ATOMIC

#endif

namespace FVQual {

/*!
 * \brief Index of the calling thread, or 0 in a serial build.
 *
 * Used to hand each thread its own slice of a shared scratch buffer. The
 * per-column vertical solves need `layerCount` doubles of workspace each, and
 * a single shared buffer is a data race the moment the column loop is
 * threaded — one that would corrupt results quietly rather than crash.
 */
[[nodiscard]] inline int threadIndex() noexcept
{
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

/*!
 * \brief Threads this process may use, or 1 in a serial build.
 */
[[nodiscard]] inline int maxThreads() noexcept
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

} // namespace FVQual

#endif // FVQual_PARALLEL_H
