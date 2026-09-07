/*
 * ff_mat_util.c
 *
 * Internal base-field matrix helpers shared by the Krylov, Frobenius and
 * characteristic-polynomial translation units. Extracted from the old
 * monolithic ff_poly.c; declarations live in ff_poly_internal.h.
 */
#include "ff_poly_internal.h"

#include <stdlib.h>
#include <string.h>

bool ffm_size_mul_ok(size_t a, size_t b)
{
    return a == 0 || b <= SIZE_MAX / a;
}

bool ffm_vector_is_zero_base(const uint32_t *v, size_t n)
{
    if (!v)
        return false;

    for (size_t i = 0; i < n; ++i)
    {
        if (v[i] != 0)
            return false;
    }

    return true;
}

bool ffm_vec_mul_base(const ff_mat_t *A, const uint32_t *x, uint32_t *y)
{
    if (!A || !x || !y)
        return false;

    if (ff_mat_rows(A) != ff_mat_cols(A))
        return false;

    size_t n = ff_mat_rows(A);
    const ff_prime_t *F = ff_mat_ctx(A);
    if (!F)
        return false;

    uint32_t *tmp = calloc(n, sizeof(uint32_t));
    if (!tmp)
        return false;

    for (size_t i = 0; i < n; ++i)
    {
        uint32_t sum = 0;
        for (size_t j = 0; j < n; ++j)
        {
            const uint32_t *ap = ff_mat_atc(A, i, j);
            if (!ap || *ap == 0)
                continue;
            if (x[j] == 0)
                continue;

            sum = ff_add(F, sum, ff_mul(F, *ap, x[j]));
        }
        tmp[i] = sum;
    }

    memcpy(y, tmp, n * sizeof(uint32_t));
    free(tmp);
    return true;
}

ff_mat_t *ffm_mat_copy_base(const ff_mat_t *M)
{
    if (!M)
        return NULL;

    const ff_prime_t *ctx = ff_mat_ctx(M);
    size_t rows = ff_mat_rows(M);
    size_t cols = ff_mat_cols(M);

    ff_mat_t *B = ff_mat_new(ctx, rows, cols);
    if (!B)
        return NULL;

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            const uint32_t *src = ff_mat_atc(M, i, j);
            uint32_t *dst = ff_mat_at(B, i, j);

            if (!src || !dst)
            {
                ff_mat_free(B);
                return NULL;
            }

            *dst = *src;
        }
    }

    return B;
}

uint32_t ffm_mat_trace_base(const ff_mat_t *M)
{
    if (!M)
        return 0;

    const ff_prime_t *ctx = ff_mat_ctx(M);
    size_t n = ff_mat_rows(M);

    if (ff_mat_cols(M) != n)
        return 0;

    uint32_t sum = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const uint32_t *d = ff_mat_atc(M, i, i);
        if (!d)
            continue;

        sum = ff_add(ctx, sum, *d);
    }

    return sum;
}

ff_mat_t *ffm_matrix_from_columns(const ff_prime_t *ctx,
                                  size_t n,
                                  const uint32_t *colmajor,
                                  size_t ncols)
{
    if (!ctx || !colmajor || ncols == 0 || n == 0)
        return NULL;

    ff_mat_t *M = ff_mat_new(ctx, n, ncols);
    if (!M)
        return NULL;

    for (size_t j = 0; j < ncols; ++j)
    {
        for (size_t i = 0; i < n; ++i)
        {
            uint32_t *dst = ff_mat_at(M, i, j);
            if (!dst)
            {
                ff_mat_free(M);
                return NULL;
            }

            *dst = colmajor[j * n + i];
        }
    }

    return M;
}

ff_mat_t *ffm_basis_append_vectors(const ff_mat_t *B,
                                   const ff_prime_t *ctx,
                                   size_t n,
                                   const uint32_t *colmajor_add,
                                   size_t add_cols)
{
    if (!ctx || n == 0 || add_cols == 0 || !colmajor_add)
        return NULL;

    size_t old_cols = B ? ff_mat_cols(B) : 0;

    if (!ffm_size_mul_ok(add_cols, n))
        return NULL;

    ff_mat_t *N = ff_mat_new(ctx, n, old_cols + add_cols);
    if (!N)
        return NULL;

    for (size_t j = 0; j < old_cols; ++j)
    {
        for (size_t i = 0; i < n; ++i)
        {
            const uint32_t *src = ff_mat_atc(B, i, j);
            uint32_t *dst = ff_mat_at(N, i, j);

            if (!src || !dst)
            {
                ff_mat_free(N);
                return NULL;
            }

            *dst = *src;
        }
    }

    for (size_t j = 0; j < add_cols; ++j)
    {
        for (size_t i = 0; i < n; ++i)
        {
            uint32_t *dst = ff_mat_at(N, i, old_cols + j);
            if (!dst)
            {
                ff_mat_free(N);
                return NULL;
            }

            *dst = colmajor_add[j * n + i];
        }
    }

    return N;
}

ff_mat_t *ffm_basis_plus_local_matrix(const ff_mat_t *G,
                                      const ff_prime_t *ctx,
                                      size_t n,
                                      const uint32_t *local_colmajor,
                                      size_t local_cols)
{
    if (!ctx || n == 0)
        return NULL;

    size_t old_cols = G ? ff_mat_cols(G) : 0;

    if (local_cols == 0)
    {
        if (old_cols == 0)
            return NULL;

        return ffm_mat_copy_base(G);
    }

    if (!local_colmajor)
        return NULL;

    ff_mat_t *C = ff_mat_new(ctx, n, old_cols + local_cols);
    if (!C)
        return NULL;

    for (size_t j = 0; j < old_cols; ++j)
    {
        for (size_t i = 0; i < n; ++i)
        {
            const uint32_t *src = ff_mat_atc(G, i, j);
            uint32_t *dst = ff_mat_at(C, i, j);

            if (!src || !dst)
            {
                ff_mat_free(C);
                return NULL;
            }

            *dst = *src;
        }
    }

    for (size_t j = 0; j < local_cols; ++j)
    {
        for (size_t i = 0; i < n; ++i)
        {
            uint32_t *dst = ff_mat_at(C, i, old_cols + j);
            if (!dst)
            {
                ff_mat_free(C);
                return NULL;
            }

            *dst = local_colmajor[j * n + i];
        }
    }

    return C;
}

bool ffm_span_contains_vector(const ff_mat_t *B,
                              size_t n,
                              const uint32_t *v,
                              uint32_t *coeffs_out)
{
    if (!v)
        return false;

    if (ffm_vector_is_zero_base(v, n))
    {
        if (B && coeffs_out)
        {
            size_t cols = ff_mat_cols(B);
            for (size_t i = 0; i < cols; ++i)
                coeffs_out[i] = 0;
        }
        return true;
    }

    if (!B)
        return false;

    size_t cols = ff_mat_cols(B);
    if (cols == 0)
        return false;

    uint32_t *tmp = NULL;

    if (coeffs_out)
    {
        tmp = coeffs_out;
    }
    else
    {
        tmp = calloc(cols, sizeof(uint32_t));
        if (!tmp)
            return false;
    }

    bool ok = ff_linear_solve_unique_base(B, v, tmp);

    if (!coeffs_out)
        free(tmp);

    return ok;
}
