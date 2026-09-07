/*
ff_poly.h
Polynomial, Krylov, and Frobenius layer header.
*/
#ifndef FF_POLY_H
#define FF_POLY_H

#include "fflib.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /* ============================================================ */
    /* ff_poly                                                      */
    /* ============================================================ */
    typedef struct ff_poly ff_poly_t;

    ff_poly_t *ff_poly_new(const ff_prime_t *ctx, int deg);
    void ff_poly_free(ff_poly_t *p);
    int ff_poly_deg(const ff_poly_t *p);
    uint32_t ff_poly_coeff(const ff_poly_t *p, int i);
    void ff_poly_set_coeff(ff_poly_t *p, int i, uint32_t value);
    void ff_poly_zero(ff_poly_t *p);
    void ff_poly_one(ff_poly_t *p);
    ff_poly_t *ff_poly_copy(const ff_poly_t *p);

    ff_poly_t *ff_poly_add_new(const ff_poly_t *a, const ff_poly_t *b);
    ff_poly_t *ff_poly_sub_new(const ff_poly_t *a, const ff_poly_t *b);
    ff_poly_t *ff_poly_mul_new(const ff_poly_t *a, const ff_poly_t *b);
    uint32_t ff_poly_eval(const ff_poly_t *p, uint32_t x_rep);
    const ff_prime_t *ff_poly_ctx(const ff_poly_t *p);
    void ff_poly_trim(ff_poly_t *p);

    /*
     * Characteristic polynomial of a square base-field matrix.
     * Uses Leverrier's algorithm (O(n^4)), so the prime should not divide
     * any integer 1, 2, ..., n. In practice, use p > n.
     */
    ff_poly_t *ff_poly_charpoly_base(const ff_mat_t *A);

    /*
     * Characteristic polynomial of a square base-field matrix using
     * Hessenberg similarity reduction and DP (O(n^3)).
     * Avoids the O(n^4) cost and division requirements of Leverrier's algorithm.
     * Returned polynomial is monic and stored low-to-high:
     *     c[0] + c[1] x + ... + c[n] x^n,  c[n] = 1.
     */
    ff_poly_t *ff_poly_charpoly_hessenberg_base(const ff_mat_t *A);

    /* ============================================================ */
    /* ff_krylov                                                    */
    /* ============================================================ */
    /* Find the first Krylov dependency for v under A. */
    bool ff_krylov_dependency_base(
        const ff_mat_t *A,
        const uint32_t *v,
        ff_poly_t **poly_out);

    /* Build the local Krylov basis together with the annihilating polynomial for v. */
    bool ff_krylov_basis_poly_base(
        const ff_mat_t *A,
        const uint32_t *v,
        ff_mat_t **basis_out,
        ff_poly_t **poly_out);

    /* ============================================================ */
    /* ff_frobenius                                                 */
    /* ============================================================ */
    /* Fill F as the companion matrix of poly. */
    bool ff_frobenius_companion_from_poly_base(
        ff_mat_t *F,
        const ff_poly_t *poly);

    /* Cyclic Frobenius verification. */
    bool ff_frobenius_cyclic_verify_base(
        const ff_mat_t *A,
        const uint32_t *v,
        ff_mat_t **P_out,
        ff_mat_t **F_out);

    /* Provisional greedy multi-block Frobenius decomposition. */
    bool ff_frobenius_greedy_base(
        const ff_mat_t *A,
        ff_mat_t **P_out,
        ff_mat_t **F_out);

    /*
     * Finds all roots of a polynomial over the base field.
     * Returns an array of roots (in active field representation) and sets *count.
     * Caller must free the returned array. Returns NULL if no roots or on error.
     */
    uint32_t *ff_poly_roots_base(const ff_poly_t *f, size_t *count);

    /* Polynomial division with remainder: num = q * den + r, deg r < deg den.
     * Returns r; if q_out is non-NULL, stores the quotient there.
     * Returns NULL if den is NULL or zero. */
    ff_poly_t *ff_poly_divrem_new(const ff_poly_t *num,
                                  const ff_poly_t *den,
                                  ff_poly_t **q_out);

    /* Polynomial modulo: returns (num % den) */
    ff_poly_t *ff_poly_mod_new(const ff_poly_t *num, const ff_poly_t *den);

    /* Polynomial derivative */
    ff_poly_t *ff_poly_deriv(const ff_poly_t *p);

    /* Modular exponentiation: returns (base^exp % mod) */
    ff_poly_t *ff_poly_mod_pow_new(const ff_poly_t *base, uint64_t exp, const ff_poly_t *mod);

    ff_poly_t *ff_poly_gcd(const ff_poly_t *a, const ff_poly_t *b);
    ff_poly_t *ff_poly_div_new(const ff_poly_t *num, const ff_poly_t *den);

    /*
     * Factors a polynomial over F_p using Berlekamp's Algorithm.
     * Best for small primes. Deterministic and matrix-based.
     * Returns an array of irreducible factors. Caller must free the array and polys.
     * Note: Assumes f is square-free.
     */
    ff_poly_t **ff_poly_factor_berlekamp(const ff_poly_t *f, size_t *out_count);
#ifdef __cplusplus
}
#endif

#endif /* FF_POLY_H */
