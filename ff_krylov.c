/*
 * ff_krylov.c
 *
 * Krylov subspace machinery:
 *   - ff_krylov_matrix_base      : Krylov matrix [v, Av, A^2 v, ...]
 *                                  (implementation MOVED here from
 *                                  fflib.c; declaration stays in fflib.h)
 *   - ff_krylov_dependency_base  : minimal annihilating polynomial of v
 *   - ff_krylov_basis_poly_base  : local Krylov basis + annihilator
 *
 * Split out of the old monolithic ff_poly.c / fflib.c.
 */
#include "ff_poly_internal.h"

#include <stdlib.h>
#include <string.h>

/* -------------------------------------------------------------------- */
/* Krylov matrix (moved from fflib.c)                                    */
/* -------------------------------------------------------------------- */

bool ff_krylov_matrix_base(ff_mat_t *P, const ff_mat_t *A, const uint32_t *v)
{
    if (!P || !A || !v)
        return false;
    if (P->rows != P->cols)
        return false;
    if (A->rows != A->cols)
        return false;
    if (P->rows != A->rows)
        return false;

    size_t n = A->rows;
    const ff_prime_t *ctx = A->ctx;

    for (size_t i = 0; i < n; ++i)
        P->data[i * P->stride + 0] = v[i];

    for (size_t col = 1; col < n; ++col)
    {
        for (size_t i = 0; i < n; ++i)
        {
            uint32_t sum = 0;

            for (size_t k = 0; k < n; ++k)
            {
                uint32_t av = A->data[i * A->stride + k];
                if (av == 0)
                    continue;

                uint32_t bv = P->data[k * P->stride + (col - 1)];
                if (bv == 0)
                    continue;

                sum = ff_add(ctx, sum, ff_mul(ctx, av, bv));
            }

            P->data[i * P->stride + col] = sum;
        }
    }

    return true;
}

/* -------------------------------------------------------------------- */
/* Minimal annihilating (dependency) polynomial                          */
/* -------------------------------------------------------------------- */

bool ff_krylov_dependency_base(const ff_mat_t *A,
                               const uint32_t *v,
                               ff_poly_t **poly_out)
{
    if (!A || !v || !poly_out)
        return false;

    size_t n = ff_mat_rows(A);
    if (ff_mat_cols(A) != n)
        return false;
    if (n == 0 || n > (size_t)INT32_MAX)
        return false;

    const ff_prime_t *ctx = ff_mat_ctx(A);
    if (!ctx)
        return false;

    if (ffm_vector_is_zero_base(v, n))
        return false;

    if (!ffm_size_mul_ok(n + 1, n))
        return false;

    /* Krylov vectors v, Av, ..., A^n v (column-major, n+1 of them). */
    uint32_t *vecs = calloc((n + 1) * n, sizeof(uint32_t));
    if (!vecs)
        return false;

    memcpy(vecs, v, n * sizeof(uint32_t));

    for (size_t k = 1; k <= n; ++k)
    {
        if (!ffm_vec_mul_base(A, &vecs[(k - 1) * n], &vecs[k * n]))
        {
            free(vecs);
            return false;
        }
    }

    /* Smallest d with A^d v in span{v, ..., A^{d-1} v}. */
    bool found = false;
    ff_poly_t *poly = NULL;

    for (size_t d = 1; d <= n; ++d)
    {
        ff_mat_t *M = ffm_matrix_from_columns(ctx, n, vecs, d);
        if (!M)
            continue;

        uint32_t *coeffs = calloc(d, sizeof(uint32_t));
        if (!coeffs)
        {
            ff_mat_free(M);
            continue;
        }

        bool ok = ff_linear_solve_unique_base(M, &vecs[d * n], coeffs);
        ff_mat_free(M);

        if (ok)
        {
            poly = ff_poly_new(ctx, (int)d);
            if (!poly)
            {
                free(coeffs);
                free(vecs);
                return false;
            }

            /* p(x) = x^d - sum coeffs[i] x^i */
            for (size_t i = 0; i < d; ++i)
                poly->c[i] = ff_neg(ctx, coeffs[i]);
            poly->c[d] = ff_one_rep(ctx);

            free(coeffs);
            found = true;
            break;
        }

        free(coeffs);
    }

    free(vecs);

    if (!found)
        return false;

    *poly_out = poly;
    return true;
}

/* -------------------------------------------------------------------- */
/* Local Krylov basis + annihilator                                      */
/* -------------------------------------------------------------------- */

bool ff_krylov_basis_poly_base(const ff_mat_t *A,
                               const uint32_t *v,
                               ff_mat_t **basis_out,
                               ff_poly_t **poly_out)
{
    if (!basis_out)
        return false;

    *basis_out = NULL;

    ff_poly_t *poly = NULL;
    if (!ff_krylov_dependency_base(A, v, &poly))
        return false;

    size_t n = ff_mat_rows(A);
    const ff_prime_t *ctx = ff_mat_ctx(A);
    int d = ff_poly_deg(poly);

    if (d <= 0)
    {
        ff_poly_free(poly);
        return false;
    }

    if (!ffm_size_mul_ok((size_t)d, n))
    {
        ff_poly_free(poly);
        return false;
    }

    uint32_t *vecs = calloc((size_t)d * n, sizeof(uint32_t));
    if (!vecs)
    {
        ff_poly_free(poly);
        return false;
    }

    memcpy(vecs, v, n * sizeof(uint32_t));

    for (int k = 1; k < d; ++k)
    {
        if (!ffm_vec_mul_base(A, &vecs[(size_t)(k - 1) * n],
                              &vecs[(size_t)k * n]))
        {
            free(vecs);
            ff_poly_free(poly);
            return false;
        }
    }

    ff_mat_t *B = ffm_matrix_from_columns(ctx, n, vecs, (size_t)d);
    free(vecs);

    if (!B)
    {
        ff_poly_free(poly);
        return false;
    }

    *basis_out = B;

    if (poly_out)
    {
        *poly_out = poly;
    }
    else
    {
        ff_poly_free(poly);
    }

    return true;
}
