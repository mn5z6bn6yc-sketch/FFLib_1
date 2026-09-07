/*
 * ff_poly_internal.h
 *
 * INTERNAL header for the ff_poly module. Not installed, not part of the
 * public API. Shared by the ff_poly*.c / ff_krylov.c / ff_frobenius.c
 * translation units that resulted from splitting the old monolithic
 * ff_poly.c.
 *
 * Exposes:
 *   - the ff_poly struct layout (so module files can work on coefficients
 *     without accessor overhead),
 *   - the shared zero-test,
 *   - the base-field matrix helpers that Krylov / Frobenius / charpoly all
 *     need (implemented in ff_mat_util.c).
 */
#ifndef FF_POLY_INTERNAL_H
#define FF_POLY_INTERNAL_H

#include "ff_poly.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /* ------------------------------------------------------------------ */
    /* ff_poly struct layout (module-private)                              */
    /* ------------------------------------------------------------------ */
    struct ff_poly
    {
        const ff_prime_t *ctx;
        int deg;      /* degree; 0 for the zero poly (c[0] == 0) */
        uint32_t *c;  /* coefficients, low to high, capacity deg+1 */
    };

    /* True iff p is NULL or all coefficients are zero. */
    bool pkf_poly_is_zero(const ff_poly_t *p);

    /*
     * Modular exponentiation with a 128-bit exponent: base^exp mod `mod`.
     * Needed by factorization code where exponents are (q^k - 1) / 2 and
     * overflow uint64_t. Implemented in ff_poly_arith.c.
     */
    ff_poly_t *pkf_poly_pow_mod_u128(const ff_poly_t *base,
                                     ff_u128_t exp,
                                     const ff_poly_t *mod);

    /* ------------------------------------------------------------------ */
    /* Base-field matrix helpers (ff_mat_util.c)                           */
    /* ------------------------------------------------------------------ */
    bool ffm_size_mul_ok(size_t a, size_t b);

    bool ffm_vector_is_zero_base(const uint32_t *v, size_t n);

    /* y = A * x for square A. */
    bool ffm_vec_mul_base(const ff_mat_t *A, const uint32_t *x, uint32_t *y);

    ff_mat_t *ffm_mat_copy_base(const ff_mat_t *M);

    uint32_t ffm_mat_trace_base(const ff_mat_t *M);

    /* Build an n x ncols matrix whose columns are taken from a
     * column-major buffer (ncols columns of n entries each). */
    ff_mat_t *ffm_matrix_from_columns(const ff_prime_t *ctx,
                                      size_t n,
                                      const uint32_t *colmajor,
                                      size_t ncols);

    /* Return a new basis matrix [B | added columns]; B may be NULL. */
    ff_mat_t *ffm_basis_append_vectors(const ff_mat_t *B,
                                       const ff_prime_t *ctx,
                                       size_t n,
                                       const uint32_t *colmajor_add,
                                       size_t add_cols);

    /* Return a new matrix [G | local columns]; G may be NULL,
     * local_cols may be 0 (plain copy of G). */
    ff_mat_t *ffm_basis_plus_local_matrix(const ff_mat_t *G,
                                          const ff_prime_t *ctx,
                                          size_t n,
                                          const uint32_t *local_colmajor,
                                          size_t local_cols);

    /* True iff v lies in span(columns of B); optionally returns the
     * unique coefficients. B may be NULL (only the zero vector is in
     * the span). */
    bool ffm_span_contains_vector(const ff_mat_t *B,
                                  size_t n,
                                  const uint32_t *v,
                                  uint32_t *coeffs_out);

#ifdef __cplusplus
}
#endif

#endif /* FF_POLY_INTERNAL_H */
