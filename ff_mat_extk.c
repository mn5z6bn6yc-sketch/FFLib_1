/*
 * ff_mat_extk.c
 *
 * Dense matrices over F_{p^k}: storage, identity, multiplication,
 * Gauss-Jordan inverse and nullspace basis.
 *
 * Split out of the old monolithic ff_extk.c. The element layout is a
 * flat ff_poly_t** array (rows*cols pointers, each of capacity k-1)
 * with row stride `stride == cols`; the struct itself stays private
 * via ff_extk_internal.h.
 *
 * Ownership contract (mirrors ff_mat_t / ff_mat_ext2_t): *_new
 * allocates, caller frees; nullspace returns its basis as the COLUMNS
 * of a cols x nullity matrix (return value = nullity, 0 = trivial,
 * -1 = error).
 */
#include "ff_extk_internal.h"

#include <stdlib.h>

/* Shared with ff_jordan_extk.c: capacity-preserving in-place copy. */
void pkf_extk_poly_copy_to(ff_poly_t *dst, const ff_poly_t *src)
{
    if (!dst || !src)
        return;

    int dst_deg = ff_poly_deg(dst);
    for (int i = 0; i <= dst_deg; ++i)
        ff_poly_set_coeff(dst, i, 0);

    int d = ff_poly_deg(src);
    for (int i = 0; i <= d; ++i)
        ff_poly_set_coeff(dst, i, ff_poly_coeff(src, i));
}

ff_mat_extk_t *ff_mat_extk_new(const ff_extk_t *ctx, size_t rows, size_t cols)
{
    if (!ctx || rows == 0 || cols == 0)
        return NULL;

    ff_mat_extk_t *M = calloc(1, sizeof(ff_mat_extk_t));
    if (!M)
        return NULL;

    M->ctx = ctx;
    M->rows = rows;
    M->cols = cols;
    M->stride = cols;

    M->data = calloc(rows * cols, sizeof(ff_poly_t *));
    if (!M->data)
    {
        free(M);
        return NULL;
    }

    for (size_t i = 0; i < rows * cols; ++i)
    {
        M->data[i] = ff_poly_new(ctx->base, ctx->k - 1);
        if (!M->data[i])
        {
            /* ff_mat_extk_free tolerates the NULL tail (calloc). */
            ff_mat_extk_free(M);
            return NULL;
        }
    }

    return M;
}

void ff_mat_extk_free(ff_mat_extk_t *M)
{
    if (!M)
        return;

    if (M->data)
    {
        for (size_t i = 0; i < M->rows * M->cols; ++i)
            ff_poly_free(M->data[i]); /* NULL-safe */

        free(M->data);
    }

    free(M);
}

ff_poly_t *ff_mat_extk_at(ff_mat_extk_t *M, size_t r, size_t c)
{
    if (!M || r >= M->rows || c >= M->cols)
        return NULL;
    return M->data[r * M->stride + c];
}

const ff_poly_t *ff_mat_extk_atc(const ff_mat_extk_t *M, size_t r, size_t c)
{
    if (!M || r >= M->rows || c >= M->cols)
        return NULL;
    return M->data[r * M->stride + c];
}

void ff_mat_extk_identity(ff_mat_extk_t *M)
{
    if (!M)
        return;

    for (size_t i = 0; i < M->rows * M->cols; ++i)
    {
        int d = ff_poly_deg(M->data[i]);
        for (int j = 0; j <= d; ++j)
            ff_poly_set_coeff(M->data[i], j, 0);
    }

    size_t n = (M->rows < M->cols) ? M->rows : M->cols;

    for (size_t i = 0; i < n; ++i)
    {
        ff_poly_set_coeff(M->data[i * M->stride + i], 0,
                          ff_one_rep(M->ctx->base));
    }
}

void ff_mat_extk_mul(ff_mat_extk_t *C, const ff_mat_extk_t *A,
                     const ff_mat_extk_t *B)
{
    if (!A || !B || !C || A->cols != B->rows ||
        C->rows != A->rows || C->cols != B->cols)
        return;

    const ff_extk_t *ctx = A->ctx;

    /* Zero C first: C must NOT alias A or B. */
    for (size_t i = 0; i < C->rows * C->cols; ++i)
    {
        int d = ff_poly_deg(C->data[i]);
        for (int j = 0; j <= d; ++j)
            ff_poly_set_coeff(C->data[i], j, 0);
    }

    ff_poly_t *temp = ff_poly_new(ctx->base, ctx->k - 1);
    ff_poly_t *sum = ff_poly_new(ctx->base, ctx->k - 1);

    if (!temp || !sum)
    {
        ff_poly_free(temp);
        ff_poly_free(sum);
        return;
    }

    for (size_t i = 0; i < A->rows; ++i)
    {
        for (size_t k = 0; k < A->cols; ++k)
        {
            const ff_poly_t *av = A->data[i * A->stride + k];
            if (ff_extk_is_zero(av))
                continue;

            for (size_t j = 0; j < B->cols; ++j)
            {
                const ff_poly_t *bv = B->data[k * B->stride + j];
                if (ff_extk_is_zero(bv))
                    continue;

                ff_extk_mul(ctx, temp, av, bv);
                ff_extk_add(ctx, sum, C->data[i * C->stride + j], temp);
                pkf_extk_poly_copy_to(C->data[i * C->stride + j], sum);
            }
        }
    }

    ff_poly_free(temp);
    ff_poly_free(sum);
}

bool ff_mat_extk_inv(ff_mat_extk_t *X, const ff_mat_extk_t *A)
{
    if (!A || !X || A->rows != A->cols ||
        X->rows != A->rows || X->cols != A->cols)
        return false;

    const ff_extk_t *ctx = A->ctx;
    size_t n = A->rows;

    ff_mat_extk_t *aug = ff_mat_extk_new(ctx, n, 2 * n);
    if (!aug)
        return false;

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
            pkf_extk_poly_copy_to(aug->data[i * aug->stride + j],
                                  A->data[i * A->stride + j]);

        ff_poly_set_coeff(aug->data[i * aug->stride + n + i], 0,
                          ff_one_rep(ctx->base));
    }

    ff_poly_t *factor = ff_poly_new(ctx->base, ctx->k - 1);
    ff_poly_t *term = ff_poly_new(ctx->base, ctx->k - 1);
    ff_poly_t *inv_pivot = ff_poly_new(ctx->base, ctx->k - 1);

    if (!factor || !term || !inv_pivot)
    {
        ff_mat_extk_free(aug);
        ff_poly_free(factor);
        ff_poly_free(term);
        ff_poly_free(inv_pivot);
        return false;
    }

    for (size_t col = 0; col < n; ++col)
    {
        size_t pivot = col;
        while (pivot < n && ff_extk_is_zero(aug->data[pivot * aug->stride + col]))
            pivot++;

        if (pivot == n)
        {
            ff_mat_extk_free(aug);
            ff_poly_free(factor);
            ff_poly_free(term);
            ff_poly_free(inv_pivot);
            return false;
        }

        if (pivot != col)
        {
            for (size_t j = 0; j < 2 * n; ++j)
            {
                ff_poly_t *t = aug->data[col * aug->stride + j];
                aug->data[col * aug->stride + j] = aug->data[pivot * aug->stride + j];
                aug->data[pivot * aug->stride + j] = t;
            }
        }

        if (!ff_extk_inv(ctx, inv_pivot, aug->data[col * aug->stride + col]))
        {
            ff_mat_extk_free(aug);
            ff_poly_free(factor);
            ff_poly_free(term);
            ff_poly_free(inv_pivot);
            return false;
        }

        for (size_t j = col; j < 2 * n; ++j)
            ff_extk_mul(ctx, aug->data[col * aug->stride + j],
                        aug->data[col * aug->stride + j], inv_pivot);

        for (size_t r = 0; r < n; ++r)
        {
            if (r == col)
                continue;

            pkf_extk_poly_copy_to(factor, aug->data[r * aug->stride + col]);
            if (ff_extk_is_zero(factor))
                continue;

            for (size_t j = col; j < 2 * n; ++j)
            {
                ff_extk_mul(ctx, term, factor, aug->data[col * aug->stride + j]);
                ff_extk_sub(ctx, aug->data[r * aug->stride + j],
                            aug->data[r * aug->stride + j], term);
            }
        }
    }

    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            pkf_extk_poly_copy_to(X->data[i * X->stride + j],
                                  aug->data[i * aug->stride + n + j]);

    ff_mat_extk_free(aug);
    ff_poly_free(factor);
    ff_poly_free(term);
    ff_poly_free(inv_pivot);
    return true;
}

int ff_mat_extk_nullspace(ff_mat_extk_t **basis_out, const ff_mat_extk_t *M)
{
    if (!basis_out || !M)
        return -1;

    *basis_out = NULL;

    const ff_extk_t *ctx = M->ctx;
    size_t rows = M->rows;
    size_t cols = M->cols;

    ff_mat_extk_t *aug = ff_mat_extk_new(ctx, rows, cols);
    if (!aug)
        return -1;

    for (size_t i = 0; i < rows * cols; ++i)
        pkf_extk_poly_copy_to(aug->data[i], M->data[i]);

    int *pivot_cols = calloc(cols, sizeof(int));
    if (!pivot_cols)
    {
        ff_mat_extk_free(aug);
        return -1;
    }

    for (size_t i = 0; i < cols; ++i)
        pivot_cols[i] = -1;

    ff_poly_t *factor = ff_poly_new(ctx->base, ctx->k - 1);
    ff_poly_t *term = ff_poly_new(ctx->base, ctx->k - 1);
    ff_poly_t *inv_pivot = ff_poly_new(ctx->base, ctx->k - 1);

    if (!factor || !term || !inv_pivot)
        goto fail;

    {
        size_t row = 0;

        for (size_t col = 0; col < cols && row < rows; ++col)
        {
            size_t pivot = row;
            while (pivot < rows && ff_extk_is_zero(aug->data[pivot * cols + col]))
                pivot++;

            if (pivot == rows)
                continue;

            if (pivot != row)
            {
                for (size_t j = 0; j < cols; ++j)
                {
                    ff_poly_t *t = aug->data[row * cols + j];
                    aug->data[row * cols + j] = aug->data[pivot * cols + j];
                    aug->data[pivot * cols + j] = t;
                }
            }

            if (!ff_extk_inv(ctx, inv_pivot, aug->data[row * cols + col]))
                goto fail;

            for (size_t j = col; j < cols; ++j)
                ff_extk_mul(ctx, aug->data[row * cols + j],
                            aug->data[row * cols + j], inv_pivot);

            for (size_t r = 0; r < rows; ++r)
            {
                if (r == row)
                    continue;

                pkf_extk_poly_copy_to(factor, aug->data[r * cols + col]);
                if (ff_extk_is_zero(factor))
                    continue;

                for (size_t j = col; j < cols; ++j)
                {
                    ff_extk_mul(ctx, term, factor, aug->data[row * cols + j]);
                    ff_extk_sub(ctx, aug->data[r * cols + j],
                                aug->data[r * cols + j], term);
                }
            }

            pivot_cols[col] = (int)row++;
        }

        size_t nullity = cols - row;

        if (nullity == 0)
        {
            free(pivot_cols);
            ff_mat_extk_free(aug);
            ff_poly_free(factor);
            ff_poly_free(term);
            ff_poly_free(inv_pivot);
            return 0;
        }

        ff_mat_extk_t *B = ff_mat_extk_new(ctx, cols, nullity);
        if (!B)
            goto fail;

        size_t bidx = 0;

        for (size_t free_col = 0; free_col < cols; ++free_col)
        {
            if (pivot_cols[free_col] != -1)
                continue;

            ff_poly_set_coeff(B->data[free_col * nullity + bidx], 0,
                              ff_one_rep(ctx->base));

            for (size_t pc = 0; pc < cols; ++pc)
            {
                int prow = pivot_cols[pc];
                if (prow == -1)
                    continue;

                ff_extk_neg(ctx, B->data[pc * nullity + bidx],
                            aug->data[(size_t)prow * cols + free_col]);
            }

            bidx++;
        }

        free(pivot_cols);
        ff_mat_extk_free(aug);
        ff_poly_free(factor);
        ff_poly_free(term);
        ff_poly_free(inv_pivot);

        *basis_out = B;
        return (int)nullity;
    }

fail:
    free(pivot_cols);
    ff_mat_extk_free(aug);
    ff_poly_free(factor);
    ff_poly_free(term);
    ff_poly_free(inv_pivot);
    return -1;
}
