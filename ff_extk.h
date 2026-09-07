/*
ff_extk.h
Generic Extension Field Towering: F_{p^k} = F_p[x] / <f(x)>
*/
#ifndef FF_EXTK_H
#define FF_EXTK_H

#include "ff_poly.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /* An element of F_{p^k} is simply a polynomial of degree < k */
    typedef ff_poly_t ff_extk_elem_t;

    typedef struct
    {
        const ff_prime_t *base;
        int k;            /* Extension degree */
        ff_poly_t *irred; /* Irreducible polynomial of degree k */
    } ff_extk_t;

    /*
     * Automatically finds an irreducible polynomial of degree k over F_p
     * and constructs the extension field context.
     */
    ff_extk_t *ff_extk_new(const ff_prime_t *base, int k);
    void ff_extk_free(ff_extk_t *ctx);

    const ff_prime_t *ff_extk_base(const ff_extk_t *ctx);
    int ff_extk_deg(const ff_extk_t *ctx);

    /* Element lifecycle */
    ff_extk_elem_t *ff_extk_elem_new(const ff_extk_t *ctx);
    void ff_extk_elem_free(ff_extk_elem_t *e);
    void ff_extk_zero(const ff_extk_t *ctx, ff_extk_elem_t *out);
    void ff_extk_one(const ff_extk_t *ctx, ff_extk_elem_t *out);
    bool ff_extk_is_zero(const ff_extk_elem_t *a);
    bool ff_extk_eq(const ff_extk_elem_t *a, const ff_extk_elem_t *b);

    /* Arithmetic (All outputs are automatically reduced mod irred) */
    void ff_extk_add(const ff_extk_t *ctx, ff_extk_elem_t *out, const ff_extk_elem_t *a, const ff_extk_elem_t *b);
    void ff_extk_sub(const ff_extk_t *ctx, ff_extk_elem_t *out, const ff_extk_elem_t *a, const ff_extk_elem_t *b);
    void ff_extk_neg(const ff_extk_t *ctx, ff_extk_elem_t *out, const ff_extk_elem_t *a); /* <--- ADD THIS */
    void ff_extk_mul(const ff_extk_t *ctx, ff_extk_elem_t *out, const ff_extk_elem_t *a, const ff_extk_elem_t *b);
    bool ff_extk_inv(const ff_extk_t *ctx, ff_extk_elem_t *out, const ff_extk_elem_t *a);
    ff_extk_t *ff_extk_new_with_poly(const ff_prime_t *base, const ff_poly_t *irred);

    /* ============================================================ */
    /* Matrices over F_{p^k}                                        */
    /* ============================================================ */
    typedef struct ff_mat_extk ff_mat_extk_t;

    ff_mat_extk_t *ff_mat_extk_new(const ff_extk_t *ctx, size_t rows, size_t cols);
    void ff_mat_extk_free(ff_mat_extk_t *M);
    ff_poly_t *ff_mat_extk_at(ff_mat_extk_t *M, size_t r, size_t c);
    const ff_poly_t *ff_mat_extk_atc(const ff_mat_extk_t *M, size_t r, size_t c);
    void ff_mat_extk_identity(ff_mat_extk_t *M);
    void ff_mat_extk_mul(ff_mat_extk_t *C, const ff_mat_extk_t *A, const ff_mat_extk_t *B);
    bool ff_mat_extk_inv(ff_mat_extk_t *X, const ff_mat_extk_t *A);
    int ff_mat_extk_nullspace(ff_mat_extk_t **basis_out, const ff_mat_extk_t *M);

    /* ============================================================ */
    /* Generic Jordan Engine for F_{p^k}                            */
    /* ============================================================ */
    typedef struct
    {
        ff_mat_extk_t *J;
        ff_mat_extk_t *P;
        size_t block_count;
        size_t *block_sizes;
        ff_poly_t **block_eigenvalues; /* Array of pointers to ff_poly_t */
    } ff_jordan_extk_result_t;

    bool ff_jordan_extk_with_eigenvalues(
        const ff_extk_t *E,
        const ff_mat_extk_t *A,
        const ff_poly_t *const *eigs, /* Array of pointers to eigenvalues */
        size_t eig_count,
        ff_jordan_extk_result_t *out);

    void ff_jordan_extk_result_clear(ff_jordan_extk_result_t *res);

#ifdef __cplusplus
}
#endif

#endif /* FF_EXTK_H */