/*
 * ff_frobenius.c
 *
 * Frobenius (rational canonical) normal form:
 *   - ff_companion_from_poly_base        : companion matrix from a raw
 *                                          coefficient buffer (MOVED here
 *                                          from fflib.c; declaration stays
 *                                          in fflib.h)
 *   - ff_frobenius_companion_from_poly_base
 *   - ff_frobenius_cyclic_verify_base    : cyclic case, self-verifying
 *   - ff_frobenius_greedy_base           : multi-block greedy, self-verifying
 *
 * Split out of the old monolithic ff_poly.c / fflib.c. Both verifiers
 * conjugate-check (P^-1 A P == F) before returning true.
 */
#include "ff_poly_internal.h"

#include <stdlib.h>
#include <string.h>

/* -------------------------------------------------------------------- */
/* Companion matrices                                                    */
/* -------------------------------------------------------------------- */

/* Fill F with the companion matrix of the monic polynomial whose first
 * n coefficients (ascending, WITHOUT the leading 1) are coeffs[0..n-1]. */
void ff_companion_from_poly_base(ff_mat_t *F, const uint32_t *coeffs)
{
    if (!F || !coeffs)
        return;
    if (F->rows != F->cols)
        return;

    size_t n = F->rows;
    const ff_prime_t *ctx = F->ctx;

    memset(F->data, 0, F->rows * F->stride * sizeof(uint32_t));

    uint32_t one = ff_one_rep(ctx);

    for (size_t i = 1; i < n; ++i)
        F->data[i * F->stride + (i - 1)] = one;

    for (size_t i = 0; i < n; ++i)
        F->data[i * F->stride + (n - 1)] = ff_neg(ctx, coeffs[i]);
}

bool ff_frobenius_companion_from_poly_base(ff_mat_t *F, const ff_poly_t *poly)
{
    if (!F || !poly)
        return false;

    size_t n = ff_mat_rows(F);
    if (ff_mat_cols(F) != n)
        return false;
    if ((size_t)poly->deg != n)
        return false;

    ff_companion_from_poly_base(F, poly->c);
    return true;
}

/* -------------------------------------------------------------------- */
/* Cyclic Frobenius verification                                         */
/* -------------------------------------------------------------------- */

bool ff_frobenius_cyclic_verify_base(const ff_mat_t *A,
                                     const uint32_t *v,
                                     ff_mat_t **P_out,
                                     ff_mat_t **F_out)
{
    if (!A || !v)
        return false;

    size_t n = ff_mat_rows(A);
    const ff_prime_t *ctx = ff_mat_ctx(A);

    bool ok = false;
    ff_poly_t *poly = NULL;
    ff_mat_t *P = NULL, *F = NULL, *P_inv = NULL, *T1 = NULL, *T2 = NULL;

    if (!ff_krylov_dependency_base(A, v, &poly))
        return false;

    P = ff_mat_new(ctx, n, n);
    F = ff_mat_new(ctx, n, n);
    P_inv = ff_mat_new(ctx, n, n);
    T1 = ff_mat_new(ctx, n, n);
    T2 = ff_mat_new(ctx, n, n);

    if (!P || !F || !P_inv || !T1 || !T2)
        goto cleanup;

    if (!ff_krylov_matrix_base(P, A, v))
        goto cleanup;

    if (!ff_frobenius_companion_from_poly_base(F, poly))
        goto cleanup;

    if (!ff_mat_inv(P_inv, P))
        goto cleanup;

    ff_mat_mul(T1, P_inv, A);
    ff_mat_mul(T2, T1, P);

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            const uint32_t *x = ff_mat_atc(T2, i, j);
            const uint32_t *y = ff_mat_atc(F, i, j);

            if (!x || !y || *x != *y)
                goto cleanup;
        }
    }

    ok = true;

    if (P_out)
    {
        *P_out = P;
        P = NULL;
    }
    if (F_out)
    {
        *F_out = F;
        F = NULL;
    }

cleanup:
    ff_poly_free(poly);
    ff_mat_free(P);
    ff_mat_free(F);
    ff_mat_free(P_inv);
    ff_mat_free(T1);
    ff_mat_free(T2);
    return ok;
}

/* -------------------------------------------------------------------- */
/* Greedy multi-block Frobenius decomposition                            */
/* -------------------------------------------------------------------- */

bool ff_frobenius_greedy_base(const ff_mat_t *A,
                              ff_mat_t **P_out,
                              ff_mat_t **F_out)
{
    if (!A)
        return false;

    size_t n = ff_mat_rows(A);
    if (ff_mat_cols(A) != n || n == 0)
        return false;

    const ff_prime_t *ctx = ff_mat_ctx(A);
    if (!ctx)
        return false;

    bool ok = false;

    ff_mat_t *G = NULL, *P = NULL, *F = NULL;
    ff_mat_t *P_inv = NULL, *T1 = NULL, *T2 = NULL;

    ff_poly_t **blocks = calloc(n, sizeof(ff_poly_t *));
    size_t *offsets = calloc(n, sizeof(size_t));

    size_t block_count = 0;
    size_t dim = 0;

    uint32_t one = ff_one_rep(ctx);

    if (!blocks || !offsets)
        goto cleanup;

    while (dim < n)
    {
        /* Find the first standard basis vector outside span(G). */
        uint32_t *e = calloc(n, sizeof(uint32_t));
        if (!e)
            goto cleanup;

        bool found_start = false;
        for (size_t i = 0; i < n; ++i)
        {
            e[i] = one;
            if (!ffm_span_contains_vector(G, n, e, NULL))
            {
                found_start = true;
                break;
            }
            e[i] = 0;
        }

        if (!found_start)
        {
            free(e);
            goto cleanup;
        }

        if (!ffm_size_mul_ok(n, n))
        {
            free(e);
            goto cleanup;
        }

        uint32_t *local = calloc(n * n, sizeof(uint32_t));
        uint32_t *curr = calloc(n, sizeof(uint32_t));
        uint32_t *next = calloc(n, sizeof(uint32_t));

        if (!local || !curr || !next)
        {
            free(e);
            free(local);
            free(curr);
            free(next);
            goto cleanup;
        }

        memcpy(local, e, n * sizeof(uint32_t));
        memcpy(curr, e, n * sizeof(uint32_t));
        free(e);

        size_t local_dim = 1;
        bool block_done = false;
        bool fail_block = false;

        while (local_dim <= n)
        {
            if (!ffm_vec_mul_base(A, curr, next))
            {
                fail_block = true;
                break;
            }

            /* Is A*curr in the LOCAL cyclic span? */
            ff_mat_t *L = ffm_matrix_from_columns(ctx, n, local, local_dim);
            if (!L)
            {
                fail_block = true;
                break;
            }

            uint32_t *coeffs = calloc(local_dim, sizeof(uint32_t));
            if (!coeffs)
            {
                ff_mat_free(L);
                fail_block = true;
                break;
            }

            bool in_local = ffm_span_contains_vector(L, n, next, coeffs);
            ff_mat_free(L);

            if (in_local)
            {
                /* Local dependency: record the block polynomial. */
                if (dim + local_dim > n)
                {
                    free(coeffs);
                    fail_block = true;
                    break;
                }

                ff_poly_t *poly = ff_poly_new(ctx, (int)local_dim);
                if (!poly)
                {
                    free(coeffs);
                    fail_block = true;
                    break;
                }

                for (size_t i = 0; i < local_dim; ++i)
                    poly->c[i] = ff_neg(ctx, coeffs[i]);
                poly->c[local_dim] = one;

                blocks[block_count] = poly;
                offsets[block_count] = dim;
                block_count++;

                ff_mat_t *newG = ffm_basis_append_vectors(G, ctx, n, local,
                                                          local_dim);
                free(coeffs);

                if (!newG)
                {
                    fail_block = true;
                    break;
                }

                ff_mat_free(G);
                G = newG;

                dim += local_dim;
                block_done = true;
                break;
            }

            free(coeffs);

            /* Not in the local span -- but if it is in span(G | local),
             * the greedy choice fails for this starting vector. */
            ff_mat_t *combined = ffm_basis_plus_local_matrix(G, ctx, n, local,
                                                             local_dim);
            if (!combined)
            {
                fail_block = true;
                break;
            }

            bool in_combined = ffm_span_contains_vector(combined, n, next, NULL);
            ff_mat_free(combined);

            if (in_combined)
            {
                fail_block = true;
                break;
            }

            if (local_dim >= n)
            {
                fail_block = true;
                break;
            }

            memcpy(&local[local_dim * n], next, n * sizeof(uint32_t));
            local_dim++;
            memcpy(curr, next, n * sizeof(uint32_t));
        }

        free(local);
        free(curr);
        free(next);

        if (fail_block || !block_done)
            goto cleanup;
    }

    if (dim != n)
        goto cleanup;

    /* Assemble the block-diagonal Frobenius matrix F. */
    F = ff_mat_new(ctx, n, n);
    if (!F)
        goto cleanup;

    for (size_t b = 0; b < block_count; ++b)
    {
        ff_poly_t *poly = blocks[b];
        size_t d = (size_t)poly->deg;
        size_t off = offsets[b];

        for (size_t i = 1; i < d; ++i)
        {
            uint32_t *dst = ff_mat_at(F, off + i, off + i - 1);
            if (!dst)
                goto cleanup;
            *dst = one;
        }

        for (size_t i = 0; i < d; ++i)
        {
            uint32_t *dst = ff_mat_at(F, off + i, off + d - 1);
            if (!dst)
                goto cleanup;
            *dst = ff_neg(ctx, poly->c[i]);
        }
    }

    P = G;
    G = NULL;

    /* Verify P^-1 A P == F before claiming success. */
    P_inv = ff_mat_new(ctx, n, n);
    T1 = ff_mat_new(ctx, n, n);
    T2 = ff_mat_new(ctx, n, n);

    if (!P_inv || !T1 || !T2)
        goto cleanup;

    if (!ff_mat_inv(P_inv, P))
        goto cleanup;

    ff_mat_mul(T1, P_inv, A);
    ff_mat_mul(T2, T1, P);

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            const uint32_t *x = ff_mat_atc(T2, i, j);
            const uint32_t *y = ff_mat_atc(F, i, j);

            if (!x || !y || *x != *y)
                goto cleanup;
        }
    }

    ok = true;

    if (P_out)
    {
        *P_out = P;
        P = NULL;
    }
    if (F_out)
    {
        *F_out = F;
        F = NULL;
    }

cleanup:
    ff_mat_free(G);
    ff_mat_free(P);
    ff_mat_free(F);
    ff_mat_free(P_inv);
    ff_mat_free(T1);
    ff_mat_free(T2);

    if (blocks)
    {
        for (size_t i = 0; i < block_count; ++i)
            ff_poly_free(blocks[i]);
        free(blocks);
    }

    free(offsets);
    return ok;
}
