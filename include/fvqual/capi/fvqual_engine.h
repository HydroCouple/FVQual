/* SPDX-License-Identifier: MIT
 * Copyright (c) 2026 Caleb Buahin
 *
 * fvqual_engine.h — error model, version, and export macro for the FVQual C
 * API.
 *
 * WHY A C API EXISTS AT ALL
 *
 * The Python bindings could have been written straight onto the C++ classes
 * with pybind11. This layer is here instead because it buys three things that
 * binding C++ directly cannot:
 *
 *   1. An ABI firewall. C++ headers, templates, and layouts can be refactored
 *      freely behind this boundary without breaking a compiled wheel.
 *   2. A multi-language surface. MATLAB, R, Julia, and Excel all speak C.
 *   3. Instance scoping. Every call takes a handle, so nothing is global and
 *      independent simulations run concurrently in threads.
 *
 * CONVENTIONS
 *
 * - Handles are pointers to incomplete struct types: opaque to callers, but
 *   still type-checked by the compiler, unlike a bare void*.
 * - Every fallible call returns FVQ_ErrorCode; FVQ_OK is zero. Out-parameters
 *   are written only on success.
 * - fvq_error_message() gives static text for a code; fvq_last_error() gives
 *   the specific detail the engine produced, which is usually the useful one.
 * - All quantities are SI: meters, seconds, kilograms, degrees Celsius. There
 *   is no unit-system switch to bind, because the engine has none.
 * - Declarations are C89-compatible, and CI compiles this header with
 *   -std=c89 -pedantic-errors so that stays true.
 */

#ifndef FVQual_CAPI_ENGINE_H
#define FVQual_CAPI_ENGINE_H

#ifdef __cplusplus
extern "C" {
#endif

/* ── Export macro ──────────────────────────────────────────────────────── */

/* FVQual_CAPI_STATIC comes with a static fvqual_core (src/CMakeLists.txt):
 * there is no DLL to import the functions from. */
#if defined(_WIN32) || defined(__CYGWIN__)
#  if defined(FVQual_CAPI_STATIC)
#    define FVQual_CAPI
#  elif defined(FVQual_CAPI_BUILD)
#    define FVQual_CAPI __declspec(dllexport)
#  else
#    define FVQual_CAPI __declspec(dllimport)
#  endif
#elif defined(__GNUC__) || defined(__clang__)
#  define FVQual_CAPI __attribute__((visibility("default")))
#else
#  define FVQual_CAPI
#endif

/* ── Errors ────────────────────────────────────────────────────────────── */

/*!
 * \brief Outcome of a C API call. Zero is success.
 *
 * The codes are coarse on purpose: they exist so a binding can map them onto
 * an exception type, while the human-readable specifics travel in
 * fvq_last_error(). A Python layer raising BadMeshError("face 7 is
 * degenerate") wants exactly this split.
 */
typedef enum FVQ_ErrorCode
{
    FVQ_OK = 0,               /*!< Success. */
    FVQ_ERR_NULL_ARGUMENT,    /*!< A required pointer argument was NULL. */
    FVQ_ERR_INVALID_ARGUMENT, /*!< An argument was out of range or inconsistent. */
    FVQ_ERR_MESH,             /*!< Mesh construction or validation failed. */
    FVQ_ERR_SIGMA,            /*!< The vertical layering specification is invalid. */
    FVQ_ERR_LIFECYCLE,        /*!< Called out of order (e.g. step before initialize). */
    FVQ_ERR_DRY_COLUMN,       /*!< A column is dry; the 3-D mode requires all wet. */
    FVQ_ERR_SPECIES,          /*!< Species registration failed (duplicate, empty id). */
    FVQ_ERR_FILE,             /*!< A file could not be read or written. */
    FVQ_ERR_UNSUPPORTED,      /*!< Built without the feature this call needs. */
    FVQ_ERR_ALLOCATION,       /*!< Out of memory. */
    FVQ_ERR_INTERNAL          /*!< An unexpected failure; please report it. */
} FVQ_ErrorCode;

/*!
 * \brief Static description of an error code.
 * \return A pointer to a string literal; never NULL, never freed.
 */
FVQual_CAPI const char *fvq_error_message(FVQ_ErrorCode code);

/*!
 * \brief The detail message from the most recent failing call on this thread.
 *
 * Thread-local: each thread sees only its own failures, which is what makes
 * concurrent simulations diagnosable. Valid until the next C API call on the
 * same thread; copy it if you need to keep it. Returns an empty string when
 * nothing has failed yet.
 */
FVQual_CAPI const char *fvq_last_error(void);

/*!
 * \brief Clears the thread-local detail message.
 */
FVQual_CAPI void fvq_clear_error(void);

/* ── Version ───────────────────────────────────────────────────────────── */

/*!
 * \brief Full version string, e.g. "0.1.0-alpha.1".
 */
FVQual_CAPI const char *fvq_version(void);

/*!
 * \brief Numeric version components. Any pointer may be NULL.
 */
FVQual_CAPI void fvq_version_numbers(int *major, int *minor, int *patch);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* FVQual_CAPI_ENGINE_H */
