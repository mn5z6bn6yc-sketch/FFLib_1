/*
 * ff_staircase.c
 *
 * Deterministic finite-field staircase / Hessenberg similarity reduction:
 * computes H and P with H = P^{-1} A P using elementary similarity
 * transformations (row operations paired with compensating column
 * operations), accumulating P as it goes.
 *
 * Split out of the old concatenated ff_decomp.c. The result is verified
 * (P^{-1} A P == H) before returning true, so a failed reduction yields
 * false rather than a wrong H.
 *
 * Provisional limitation (kept from the original): only the
 * adjacent-above entry H[r-1][c] is used as a pivot, with no row/column
 * swapping -- the reduction reports failure if that pivot is zero even
 * when another pivot would work. ff_poly_charpoly.c's internal reduction
 * uses subdiagonal pivots WITH swapping; the two should be consolidated
 * eventually.
 */
#include "ff_decomp.h"
#include "ff_poly_internal.h" /* ffm_mat_copy_base */

bool ff_staircase_hessenberg_base(
    const ff_mat_t *A,
    ff_mat_t **P_out,
    ff_mat_t **H_out)
{
    if (!A)
        return false;

    const ff_prime_t *ctx = ff_mat_ctx(A);
    size_t n = ff_mat_rows(A);

    if (!ctx)
        return false;
    if (ff_mat_cols(A) != n)
        return false;
    if (n == 0)
        return false;

    ff_mat_t *H = ffm_mat_copy_base(A);
    ff_mat_t *P = ff_mat_new(ctx, n, n);

    if (!H || !P)
    {
        ff_mat_free(H);
        ff_mat_free(P);
        return false;
    }

    ff_mat_identity(P);

    if (n >= 3)
    {
        for (size_t c = 0; c + 2 < n; ++c)
        {
            for (size_t r = n - 1; r > c + 1; --r)
            {
                const uint32_t *xp = ff_mat_atc(H, r, c);
                if (!xp || *xp == 0)
                    continue; /* already zero below the subdiagonal */

                const uint32_t *pp = ff_mat_atc(H, r - 1, c);
                if (!pp || *pp == 0)
                {
                    /* No adjacent pivot available (see file comment). */
                    ff_mat_free(H);
                    ff_mat_free(P);
                    return false;
                }

                uint32_t inv_pivot = ff_inv(ctx, *pp);
                if (inv_pivot == 0)
                {
                    ff_mat_free(H);
                    ff_mat_free(P);
                    return false;
                }

                uint32_t f = ff_mul(ctx, *xp, inv_pivot);

                /* Row operation: row_r <- row_r - f * row_{r-1}.
                   Clears H[r][c]. */
                for (size_t j = 0; j < n; ++j)
                {
                    uint32_t *hrj = ff_mat_at(H, r, j);
                    const uint32_t *hr1j = ff_mat_atc(H, r - 1, j);

                    if (!hrj || !hr1j)
                    {
                        ff_mat_free(H);
                        ff_mat_free(P);
                        return false;
                    }

                    uint32_t term = ff_mul(ctx, f, *hr1j);
                    *hrj = ff_sub(ctx, *hrj, term);
                }

                /* Column compensation: col_{r-1} <- col_{r-1} + f * col_r.
                   Keeps the transformation a similarity. */
                for (size_t i = 0; i < n; ++i)
                {
                    uint32_t *hc1 = ff_mat_at(H, i, r - 1);
                    const uint32_t *hcr = ff_mat_atc(H, i, r);

                    if (!hc1 || !hcr)
                    {
                        ff_mat_free(H);
                        ff_mat_free(P);
                        return false;
                    }

                    uint32_t term = ff_mul(ctx, f, *hcr);
                    *hc1 = ff_add(ctx, *hc1, term);
                }

                /* Accumulate P <- P * (I + f * E_{r, r-1}):
                   col_{r-1} <- col_{r-1} + f * col_r. */
                for (size_t i = 0; i < n; ++i)
                {
                    uint32_t *pc1 = ff_mat_at(P, i, r - 1);
                    const uint32_t *pcr = ff_mat_atc(P, i, r);

                    if (!pc1 || !pcr)
                    {
                        ff_mat_free(H);
                        ff_mat_free(P);
                        return false;
                    }

                    uint32_t term = ff_mul(ctx, f, *pcr);
                    *pc1 = ff_add(ctx, *pc1, term);
                }
            }
        }
    }

    /* Verification: P^{-1} A P == H. */
    ff_mat_t *P_inv = ff_mat_new(ctx, n, n);
    ff_mat_t *T1 = ff_mat_new(ctx, n, n);
    ff_mat_t *T2 = ff_mat_new(ctx, n, n);

    if (!P_inv || !T1 || !T2)
    {
        ff_mat_free(H);
        ff_mat_free(P);
        ff_mat_free(P_inv);
        ff_mat_free(T1);
        ff_mat_free(T2);
        return false;
    }

    if (!ff_mat_inv(P_inv, P))
    {
        ff_mat_free(H);
        ff_mat_free(P);
        ff_mat_free(P_inv);
        ff_mat_free(T1);
        ff_mat_free(T2);
        return false;
    }

    ff_mat_mul(T1, P_inv, A);
    ff_mat_mul(T2, T1, P);

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            const uint32_t *x = ff_mat_atc(T2, i, j);
            const uint32_t *y = ff_mat_atc(H, i, j);

            if (!x || !y || *x != *y)
            {
                ff_mat_free(H);
                ff_mat_free(P);
                ff_mat_free(P_inv);
                ff_mat_free(T1);
                ff_mat_free(T2);
                return false;
            }
        }
    }

    ff_mat_free(P_inv);
    ff_mat_free(T1);
    ff_mat_free(T2);

    if (P_out)
        *P_out = P;
    else
        ff_mat_free(P);

    if (H_out)
        *H_out = H;
    else
        ff_mat_free(H);

    return true;
}
