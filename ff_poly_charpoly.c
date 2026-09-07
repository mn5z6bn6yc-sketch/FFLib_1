/*
 * ff_poly_charpoly.c
 *
 * Characteristic polynomials of square base-field matrices:
 *   - ff_poly_charpoly_base             : Leverrier / Newton identities,
 *                                         O(n^4), needs p > n
 *   - ff_poly_charpoly_hessenberg_base  : Hessenberg similarity + DP,
 *                                         O(n^3), no division assumptions
 *
 * Split out of the old monolithic ff_poly.c. The Hessenberg reduction
 * here also answers the library-wide convention question: H = P^-1 A P
 * (verified by ff_frobenius_cyclic_verify_base and ff_staircase).
 */
#include "ff_poly_internal.h"

#include <limits.h>
#include <stdlib.h>

/* -------------------------------------------------------------------- */
/* Leverrier (Newton identities), O(n^4)                                 */
/* -------------------------------------------------------------------- */

ff_poly_t *ff_poly_charpoly_base(const ff_mat_t *A)
{
    if (!A)
        return NULL;

    size_t n_sz = ff_mat_rows(A);
    if (ff_mat_cols(A) != n_sz)
        return NULL;
    if (n_sz > (size_t)INT_MAX)
        return NULL;

    int n = (int)n_sz;
    const ff_prime_t *ctx = ff_mat_ctx(A);
    if (!ctx)
        return NULL;

    if (n == 0)
    {
        ff_poly_t *one = ff_poly_new(ctx, 0);
        if (one)
            ff_poly_one(one);
        return one;
    }

    ff_poly_t *poly = ff_poly_new(ctx, n);
    if (!poly)
        return NULL;

    poly->c[n] = ff_one_rep(ctx);

    uint32_t p = ff_prime_modulus(ctx);

    if ((uint64_t)n >= (uint64_t)p)
    {
        /*
         * Leverrier divides by 1, 2, ..., n. If p <= n one of these is
         * zero modulo p. Use ff_poly_charpoly_hessenberg_base instead.
         */
        ff_poly_free(poly);
        return NULL;
    }

    uint32_t *traces = calloc((size_t)n + 1, sizeof(uint32_t));
    ff_mat_t *Power = ffm_mat_copy_base(A);

    if (!traces || !Power)
    {
        free(traces);
        ff_mat_free(Power);
        ff_poly_free(poly);
        return NULL;
    }

    for (int k = 1; k <= n; ++k)
    {
        if (k > 1)
        {
            ff_mat_t *Next = ff_mat_new(ctx, (size_t)n, (size_t)n);
            if (!Next)
            {
                ff_mat_free(Power);
                free(traces);
                ff_poly_free(poly);
                return NULL;
            }

            ff_mat_mul(Next, Power, A);
            ff_mat_free(Power);
            Power = Next;
        }

        traces[k] = ffm_mat_trace_base(Power);

        uint32_t sum = traces[k];
        for (int i = 1; i < k; ++i)
        {
            uint32_t coeff = poly->c[n - i];
            uint32_t tr = traces[k - i];
            sum = ff_add(ctx, sum, ff_mul(ctx, coeff, tr));
        }

        uint32_t k_rep = ff_to_mont(ctx, (uint32_t)k % p);
        uint32_t inv_k = ff_inv(ctx, k_rep);
        if (inv_k == 0)
        {
            ff_mat_free(Power);
            free(traces);
            ff_poly_free(poly);
            return NULL;
        }

        poly->c[n - k] = ff_neg(ctx, ff_mul(ctx, sum, inv_k));
    }

    ff_mat_free(Power);
    free(traces);
    return poly;
}

/* -------------------------------------------------------------------- */
/* Hessenberg reduction + DP charpoly, O(n^3)                            */
/* -------------------------------------------------------------------- */

/*
 * Reduce A to upper Hessenberg form by similarity transformations
 * (H = P^-1 A P; P is not returned here).
 *
 * Uses a subdiagonal pivot search with row/column swaps, so it succeeds
 * whenever any Hessenberg reduction exists. May still produce a REDUCIBLE
 * Hessenberg matrix (zero on the subdiagonal) -- the DP below handles
 * that via the sub_prod == 0 short-circuit.
 */
static ff_mat_t *mat_to_hessenberg_base(const ff_mat_t *A)
{
    if (!A)
        return NULL;

    const ff_prime_t *ctx = ff_mat_ctx(A);
    size_t n = ff_mat_rows(A);
    if (ff_mat_cols(A) != n || n == 0)
        return NULL;

    ff_mat_t *H = ffm_mat_copy_base(A);
    if (!H)
        return NULL;

    if (n <= 2)
        return H;

    for (size_t c = 0; c + 2 < n; c++)
    {
        /* Find a pivot in column c below the subdiagonal. */
        size_t r_pivot = c + 1;
        while (r_pivot < n)
        {
            const uint32_t *pv = ff_mat_atc(H, r_pivot, c);
            if (pv && *pv != 0)
                break;
            r_pivot++;
        }
        if (r_pivot == n)
            continue; /* column already zero below the subdiagonal */

        if (r_pivot != c + 1)
        {
            /* Swap rows c+1 and r_pivot (left multiplication). */
            for (size_t j = 0; j < n; j++)
            {
                uint32_t tmp = *ff_mat_at(H, c + 1, j);
                *ff_mat_at(H, c + 1, j) = *ff_mat_atc(H, r_pivot, j);
                *ff_mat_at(H, r_pivot, j) = tmp;
            }
            /* Swap columns c+1 and r_pivot (right multiplication). */
            for (size_t i = 0; i < n; i++)
            {
                uint32_t tmp = *ff_mat_at(H, i, c + 1);
                *ff_mat_at(H, i, c + 1) = *ff_mat_atc(H, i, r_pivot);
                *ff_mat_at(H, i, r_pivot) = tmp;
            }
        }

        const uint32_t *pivot_ptr = ff_mat_atc(H, c + 1, c);
        if (!pivot_ptr)
        {
            ff_mat_free(H);
            return NULL;
        }

        uint32_t inv_pivot = ff_inv(ctx, *pivot_ptr);
        if (inv_pivot == 0)
        {
            ff_mat_free(H);
            return NULL;
        }

        /* Eliminate entries below the subdiagonal. */
        for (size_t r = c + 2; r < n; r++)
        {
            const uint32_t *hrc = ff_mat_atc(H, r, c);
            if (!hrc || *hrc == 0)
                continue;

            uint32_t factor = ff_mul(ctx, *hrc, inv_pivot);

            /* Row operation: H[r, :] -= factor * H[c+1, :] */
            for (size_t j = c; j < n; j++)
            {
                const uint32_t *h_c1_j = ff_mat_atc(H, c + 1, j);
                uint32_t *h_r_j = ff_mat_at(H, r, j);
                if (!h_c1_j || !h_r_j)
                {
                    ff_mat_free(H);
                    return NULL;
                }
                *h_r_j = ff_sub(ctx, *h_r_j, ff_mul(ctx, factor, *h_c1_j));
            }

            /* Column compensation: H[:, c+1] += factor * H[:, r] */
            for (size_t i = 0; i < n; i++)
            {
                const uint32_t *h_i_r = ff_mat_atc(H, i, r);
                uint32_t *h_i_c1 = ff_mat_at(H, i, c + 1);
                if (!h_i_r || !h_i_c1)
                {
                    ff_mat_free(H);
                    return NULL;
                }
                *h_i_c1 = ff_add(ctx, *h_i_c1, ff_mul(ctx, factor, *h_i_r));
            }
        }
    }

    return H;
}

/*
 * Characteristic polynomial via Hessenberg similarity + dynamic
 * programming in O(n^3). The result is monic, stored low-to-high:
 *     c[0] + c[1] x + ... + c[n] x^n,  c[n] = 1.
 */
ff_poly_t *ff_poly_charpoly_hessenberg_base(const ff_mat_t *A)
{
    if (!A)
        return NULL;

    const ff_prime_t *ctx = ff_mat_ctx(A);
    size_t n = ff_mat_rows(A);
    if (ff_mat_cols(A) != n)
        return NULL;

    if (n == 0)
    {
        ff_poly_t *one = ff_poly_new(ctx, 0);
        if (one)
            ff_poly_one(one);
        return one;
    }

    ff_mat_t *H = mat_to_hessenberg_base(A);
    if (!H)
        return NULL;

    ff_poly_t **P = calloc(n + 1, sizeof(ff_poly_t *));
    if (!P)
    {
        ff_mat_free(H);
        return NULL;
    }

    P[0] = ff_poly_new(ctx, 0);
    if (!P[0])
        goto fail;
    ff_poly_one(P[0]);

    for (size_t i = 1; i <= n; i++)
    {
        /* P_i(x) = (x - H[i-1][i-1]) * P_{i-1}(x) */
        const uint32_t *diag = ff_mat_atc(H, i - 1, i - 1);
        if (!diag)
            goto fail;

        ff_poly_t *x_minus_h = ff_poly_new(ctx, 1);
        if (!x_minus_h)
            goto fail;
        ff_poly_set_coeff(x_minus_h, 0, ff_neg(ctx, *diag));
        ff_poly_set_coeff(x_minus_h, 1, ff_one_rep(ctx));

        P[i] = ff_poly_mul_new(x_minus_h, P[i - 1]);
        ff_poly_free(x_minus_h);
        if (!P[i])
            goto fail;

        /* Subtract the back-tracking subdiagonal products:
         *   sum_k H[k][k-1] ... H[i-1][i-2] * H[k-1][i-1] * P_{k-1}(x) */
        uint32_t sub_prod = ff_one_rep(ctx);
        for (size_t k_idx = 0; k_idx + 1 < i; k_idx++)
        {
            size_t k = i - 1 - k_idx; /* k: i-1 down to 1 */

            const uint32_t *h_k_k1 = ff_mat_atc(H, k, k - 1);
            if (!h_k_k1)
                goto fail;

            sub_prod = ff_mul(ctx, sub_prod, *h_k_k1);
            if (sub_prod == 0)
                break; /* whole tail vanishes */

            const uint32_t *h_k1_i1 = ff_mat_atc(H, k - 1, i - 1);
            if (!h_k1_i1)
                goto fail;

            uint32_t term = ff_mul(ctx, sub_prod, *h_k1_i1);
            if (term == 0)
                continue;

            int deg_prev = ff_poly_deg(P[k - 1]);
            ff_poly_t *term_poly = ff_poly_new(ctx, deg_prev);
            if (!term_poly)
                goto fail;

            for (int d = 0; d <= deg_prev; d++)
            {
                ff_poly_set_coeff(term_poly, d,
                                  ff_mul(ctx, term, ff_poly_coeff(P[k - 1], d)));
            }

            ff_poly_t *temp = ff_poly_sub_new(P[i], term_poly);
            ff_poly_free(term_poly);
            ff_poly_free(P[i]);
            P[i] = temp;
            if (!P[i])
                goto fail;
        }
    }

    ff_poly_t *result = P[n];
    P[n] = NULL;

    for (size_t i = 0; i <= n; i++)
        ff_poly_free(P[i]);
    free(P);
    ff_mat_free(H);
    return result;

fail:
    for (size_t i = 0; i <= n; i++)
        ff_poly_free(P[i]);
    free(P);
    ff_mat_free(H);
    return NULL;
}
