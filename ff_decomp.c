/*
    ff_decomp.c

    Concatenated provisional decomposition layer:

      - ff_jordan
      - ff_staircase
      - ff_ratrec

    Requires:
      fflib.h / fflib implementation

    Notes:
      - Jordan engine currently requires eigenvalues to be supplied.
      - Staircase engine produces Hessenberg-like similarity form and
        verifies P^{-1} A P = H.
      - Rational reconstruction uses signed/unsigned __int128 and is
        intended for moderate CRT moduli. For very large reconstructions,
        use GMP-backed rational reconstruction.
*/
#include "ff_decomp.h"
#ifndef FF_DECOMP_H
#define FF_DECOMP_H

#include "fflib.h"

#if defined(__GNUC__) || defined(__clang__)
typedef signed __int128 ff_i128_t;
#else
#error "ff_decomp requires signed __int128 support."
#endif

#ifdef __cplusplus
extern "C"
{
#endif

    /* ============================================================ */
    /* ff_jordan                                                     */
    /* ============================================================ */

    typedef struct
    {
        ff_mat_ext2_t *J;
        ff_mat_ext2_t *P;

        size_t block_count;
        size_t *block_sizes;

        ff_ext2_elem_t *block_eigenvalues;
    } ff_jordan_result_t;

    /*
        Jordan decomposition over an extension field.

        The eigenvalues must be supplied. This provisional engine does not
        yet include polynomial factorisation or automatic eigenvalue discovery.

        A must be a square matrix over the extension field.
    */
    bool ff_jordan_ext2_with_eigenvalues(
        const ff_ext2_t *E,
        const ff_mat_ext2_t *A,
        const ff_ext2_elem_t *eigs,
        size_t eig_count,
        ff_jordan_result_t *out);

    /*
        Convenience wrapper:

        Converts a base-field matrix A into an extension-field matrix and then
        calls ff_jordan_ext2_with_eigenvalues().
    */
    bool ff_jordan_base_with_eigenvalues_ext2(
        const ff_ext2_t *E,
        const ff_mat_t *A_base,
        const ff_ext2_elem_t *eigs,
        size_t eig_count,
        ff_jordan_result_t *out);

    void ff_jordan_result_clear(ff_jordan_result_t *res);

    /* ============================================================ */
    /* ff_staircase                                                  */
    /* ============================================================ */

    /*
        Deterministic finite-field staircase / Hessenberg similarity reduction.

        Computes H and P such that:

            H = P^{-1} A P

        using elementary similarity transformations.

        This provisional engine does not use global block reordering.
        It may fail if a required pivot is zero.
    */
    bool ff_staircase_hessenberg_base(
        const ff_mat_t *A,
        ff_mat_t **P_out,
        ff_mat_t **H_out);

    /* ============================================================ */
    /* ff_ratrec                                                     */
    /* ============================================================ */

    /*
        Rational reconstruction.

        Given x mod M, find n/d such that:

            n * d^{-1} == x mod M

        with small numerator and denominator.

        This provisional implementation uses __int128.
    */
    bool ff_ratrec_from_u128(
        ff_u128_t x,
        ff_u128_t M,
        ff_i128_t *num,
        ff_i128_t *den);

    bool ff_ratrec_crt_i128(
        const ff_crt_t *crt,
        ff_i128_t *num,
        ff_i128_t *den);

    bool ff_ratrec_crt_i64(
        const ff_crt_t *crt,
        int64_t *num,
        int64_t *den);

    bool ff_ratrec_array_i128(
        const uint32_t *primes,
        const uint32_t *residues,
        size_t count,
        ff_i128_t *num,
        ff_i128_t *den);

#ifdef __cplusplus
}
#endif

#endif /* FF_DECOMP_H */

/* ============================================================ */
/* Implementation                                               */
/* ============================================================ */

#include <stdlib.h>
#include <string.h>
#include <stdint.h>

/* ============================================================ */
/* Required fflib accessor declarations                         */
/* ============================================================ */

extern const ff_prime_t *ff_mat_ctx(const ff_mat_t *M);
extern size_t ff_mat_rows(const ff_mat_t *M);
extern size_t ff_mat_cols(const ff_mat_t *M);

extern const ff_ext2_t *ff_mat_ext2_ctx(const ff_mat_ext2_t *M);
extern size_t ff_mat_ext2_rows(const ff_mat_ext2_t *M);
extern size_t ff_mat_ext2_cols(const ff_mat_ext2_t *M);

extern ff_u128_t ff_crt_modulus(const ff_crt_t *crt);

/* ============================================================ */
/* Small utility helpers                                         */
/* ============================================================ */

static bool ffd_size_mul_ok(size_t a, size_t b)
{
    return a == 0 || b <= SIZE_MAX / a;
}

static ff_i128_t ffd_i128_abs(ff_i128_t x)
{
    return (x < 0) ? -x : x;
}

static ff_i128_t ffd_i128_gcd(ff_i128_t a, ff_i128_t b)
{
    if (a < 0)
        a = -a;
    if (b < 0)
        b = -b;

    while (b != 0)
    {
        ff_i128_t t = a % b;
        a = b;
        b = t;
    }

    return a;
}

static ff_u128_t ffd_u128_sqrt(ff_u128_t x)
{
    if (x == 0)
        return 0;

    ff_u128_t lo = 0;
    ff_u128_t hi = ((ff_u128_t)1) << 64;

    while (lo < hi)
    {
        ff_u128_t mid = (lo + hi + 1) >> 1;

        if (mid <= x / mid)
        {
            lo = mid;
        }
        else
        {
            hi = mid - 1;
        }
    }

    return lo;
}

/* ============================================================ */
/* ff_jordan local helpers                                       */
/* ============================================================ */

static bool ffd_ext2_vec_is_zero(
    const ff_ext2_elem_t *v,
    size_t n)
{
    if (!v)
        return true;

    for (size_t i = 0; i < n; ++i)
    {
        if (!ff_ext2_is_zero(&v[i]))
            return false;
    }

    return true;
}

static ff_mat_ext2_t *ffd_ext2_mat_copy(const ff_mat_ext2_t *A)
{
    if (!A)
        return NULL;

    const ff_ext2_t *ctx = ff_mat_ext2_ctx(A);
    size_t rows = ff_mat_ext2_rows(A);
    size_t cols = ff_mat_ext2_cols(A);

    ff_mat_ext2_t *B = ff_mat_ext2_new(ctx, rows, cols);
    if (!B)
        return NULL;

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            const ff_ext2_elem_t *src = ff_mat_ext2_atc(A, i, j);
            ff_ext2_elem_t *dst = ff_mat_ext2_at(B, i, j);

            if (!src || !dst)
            {
                ff_mat_ext2_free(B);
                return NULL;
            }

            *dst = *src;
        }
    }

    return B;
}

static ff_mat_ext2_t *ffd_ext2_mat_from_base(
    const ff_ext2_t *E,
    const ff_mat_t *A)
{
    if (!E || !A)
        return NULL;

    size_t rows = ff_mat_rows(A);
    size_t cols = ff_mat_cols(A);

    ff_mat_ext2_t *B = ff_mat_ext2_new(E, rows, cols);
    if (!B)
        return NULL;

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            const uint32_t *src = ff_mat_atc(A, i, j);
            ff_ext2_elem_t *dst = ff_mat_ext2_at(B, i, j);

            if (!src || !dst)
            {
                ff_mat_ext2_free(B);
                return NULL;
            }

            dst->a = *src;
            dst->b = 0;
        }
    }

    return B;
}

static ff_mat_ext2_t *ffd_ext2_mat_sub_scalar_identity(
    const ff_mat_ext2_t *A,
    const ff_ext2_elem_t *lambda)
{
    if (!A || !lambda)
        return NULL;

    const ff_ext2_t *ctx = ff_mat_ext2_ctx(A);
    size_t n = ff_mat_ext2_rows(A);

    if (ff_mat_ext2_cols(A) != n)
        return NULL;

    ff_mat_ext2_t *N = ffd_ext2_mat_copy(A);
    if (!N)
        return NULL;

    for (size_t i = 0; i < n; ++i)
    {
        ff_ext2_elem_t *d = ff_mat_ext2_at(N, i, i);
        if (!d)
        {
            ff_mat_ext2_free(N);
            return NULL;
        }

        ff_ext2_sub(ctx, d, d, lambda);
    }

    return N;
}

static bool ffd_ext2_mat_vec_mul(
    const ff_mat_ext2_t *M,
    const ff_ext2_elem_t *x,
    ff_ext2_elem_t *y)
{
    if (!M || !x || !y)
        return false;

    const ff_ext2_t *ctx = ff_mat_ext2_ctx(M);
    size_t rows = ff_mat_ext2_rows(M);
    size_t cols = ff_mat_ext2_cols(M);

    if (rows != cols)
        return false;

    ff_ext2_elem_t *tmp = calloc(rows, sizeof(ff_ext2_elem_t));
    if (!tmp)
        return false;

    for (size_t i = 0; i < rows; ++i)
    {
        ff_ext2_elem_t sum;
        ff_ext2_zero(ctx, &sum);

        for (size_t j = 0; j < cols; ++j)
        {
            const ff_ext2_elem_t *m = ff_mat_ext2_atc(M, i, j);
            if (!m || ff_ext2_is_zero(m))
                continue;

            if (ff_ext2_is_zero(&x[j]))
                continue;

            ff_ext2_elem_t prod;
            ff_ext2_mul(ctx, &prod, m, &x[j]);

            ff_ext2_add(ctx, &sum, &sum, &prod);
        }

        tmp[i] = sum;
    }

    memcpy(y, tmp, rows * sizeof(ff_ext2_elem_t));
    free(tmp);

    return true;
}

static ff_ext2_elem_t *ffd_ext2_vec_copy_column(
    const ff_mat_ext2_t *M,
    size_t col)
{
    if (!M)
        return NULL;

    size_t rows = ff_mat_ext2_rows(M);

    if (col >= ff_mat_ext2_cols(M))
        return NULL;

    ff_ext2_elem_t *v = calloc(rows, sizeof(ff_ext2_elem_t));
    if (!v)
        return NULL;

    for (size_t i = 0; i < rows; ++i)
    {
        const ff_ext2_elem_t *src = ff_mat_ext2_atc(M, i, col);
        if (!src)
        {
            free(v);
            return NULL;
        }

        v[i] = *src;
    }

    return v;
}

static ff_mat_ext2_t *ffd_ext2_basis_append_vector(
    ff_mat_ext2_t *B,
    const ff_ext2_t *ctx,
    size_t n,
    const ff_ext2_elem_t *v)
{
    if (!ctx || !v || n == 0)
        return NULL;

    size_t old_cols = B ? ff_mat_ext2_cols(B) : 0;

    ff_mat_ext2_t *N = ff_mat_ext2_new(ctx, n, old_cols + 1);
    if (!N)
        return NULL;

    for (size_t j = 0; j < old_cols; ++j)
    {
        for (size_t i = 0; i < n; ++i)
        {
            const ff_ext2_elem_t *src = ff_mat_ext2_atc(B, i, j);
            ff_ext2_elem_t *dst = ff_mat_ext2_at(N, i, j);

            if (!src || !dst)
            {
                ff_mat_ext2_free(N);
                return NULL;
            }

            *dst = *src;
        }
    }

    for (size_t i = 0; i < n; ++i)
    {
        ff_ext2_elem_t *dst = ff_mat_ext2_at(N, i, old_cols);
        if (!dst)
        {
            ff_mat_ext2_free(N);
            return NULL;
        }

        *dst = v[i];
    }

    return N;
}

static bool ffd_ext2_vector_independent(
    const ff_mat_ext2_t *B,
    const ff_ext2_t *ctx,
    size_t n,
    const ff_ext2_elem_t *v)
{
    if (!ctx || !v || n == 0)
        return false;

    if (ffd_ext2_vec_is_zero(v, n))
        return false;

    size_t cols = B ? ff_mat_ext2_cols(B) : 0;

    if (cols == 0)
        return true;

    if (cols + 1 > n)
        return false;

    ff_mat_ext2_t *aug = ff_mat_ext2_new(ctx, n, cols + 1);
    if (!aug)
        return false;

    for (size_t j = 0; j < cols; ++j)
    {
        for (size_t i = 0; i < n; ++i)
        {
            const ff_ext2_elem_t *src = ff_mat_ext2_atc(B, i, j);
            ff_ext2_elem_t *dst = ff_mat_ext2_at(aug, i, j);

            if (!src || !dst)
            {
                ff_mat_ext2_free(aug);
                return false;
            }

            *dst = *src;
        }
    }

    for (size_t i = 0; i < n; ++i)
    {
        ff_ext2_elem_t *dst = ff_mat_ext2_at(aug, i, cols);
        if (!dst)
        {
            ff_mat_ext2_free(aug);
            return false;
        }

        *dst = v[i];
    }

    ff_mat_ext2_t *ns = NULL;
    int nullity = ff_mat_ext2_nullspace(&ns, aug);

    ff_mat_ext2_free(aug);
    ff_mat_ext2_free(ns);

    if (nullity < 0)
        return false;

    return nullity == 0;
}

static ff_mat_ext2_t *ffd_ext2_column_basis(const ff_mat_ext2_t *M)
{
    if (!M)
        return NULL;

    const ff_ext2_t *ctx = ff_mat_ext2_ctx(M);
    size_t rows = ff_mat_ext2_rows(M);
    size_t cols = ff_mat_ext2_cols(M);

    ff_mat_ext2_t *B = NULL;

    for (size_t j = 0; j < cols; ++j)
    {
        ff_ext2_elem_t *v = ffd_ext2_vec_copy_column(M, j);
        if (!v)
        {
            ff_mat_ext2_free(B);
            return NULL;
        }

        if (ffd_ext2_vector_independent(B, ctx, rows, v))
        {
            ff_mat_ext2_t *N = ffd_ext2_basis_append_vector(B, ctx, rows, v);
            ff_mat_ext2_free(B);
            B = N;

            if (!B)
            {
                free(v);
                return NULL;
            }
        }

        free(v);
    }

    return B;
}

static ff_mat_ext2_t *ffd_ext2_matrix_image(
    const ff_mat_ext2_t *A,
    const ff_mat_ext2_t *basis)
{
    if (!A || !basis)
        return NULL;

    const ff_ext2_t *ctx = ff_mat_ext2_ctx(A);
    size_t rows = ff_mat_ext2_rows(A);
    size_t cols = ff_mat_ext2_cols(basis);

    ff_mat_ext2_t *img = ff_mat_ext2_new(ctx, rows, cols);
    if (!img)
        return NULL;

    ff_mat_ext2_mul(img, A, basis);

    return img;
}

static ff_mat_ext2_t *ffd_ext2_basis_append_matrix(
    ff_mat_ext2_t *B,
    const ff_mat_ext2_t *Add)
{
    if (!Add)
        return B;

    const ff_ext2_t *ctx = ff_mat_ext2_ctx(Add);
    size_t rows = ff_mat_ext2_rows(Add);
    size_t cols = ff_mat_ext2_cols(Add);

    for (size_t j = 0; j < cols; ++j)
    {
        ff_ext2_elem_t *v = ffd_ext2_vec_copy_column(Add, j);
        if (!v)
        {
            ff_mat_ext2_free(B);
            return NULL;
        }

        ff_mat_ext2_t *N = ffd_ext2_basis_append_vector(B, ctx, rows, v);

        free(v);
        ff_mat_ext2_free(B);

        B = N;

        if (!B)
            return NULL;
    }

    return B;
}

/* ============================================================ */
/* ff_jordan implementation                                      */
/* ============================================================ */

void ff_jordan_result_clear(ff_jordan_result_t *res)
{
    if (!res)
        return;

    ff_mat_ext2_free(res->J);
    ff_mat_ext2_free(res->P);

    free(res->block_sizes);
    free(res->block_eigenvalues);

    memset(res, 0, sizeof(*res));
}

bool ff_jordan_ext2_with_eigenvalues(
    const ff_ext2_t *E,
    const ff_mat_ext2_t *A,
    const ff_ext2_elem_t *eigs,
    size_t eig_count,
    ff_jordan_result_t *out)
{
    if (!A || !eigs || !out)
        return false;

    memset(out, 0, sizeof(*out));

    const ff_ext2_t *ctx = ff_mat_ext2_ctx(A);

    if (!ctx || ctx != E)
        return false;

    size_t n = ff_mat_ext2_rows(A);

    if (ff_mat_ext2_cols(A) != n)
        return false;
    if (n == 0)
        return false;

    ff_mat_ext2_t *P_basis = NULL;
    ff_mat_ext2_t *J = ff_mat_ext2_new(ctx, n, n);

    size_t *block_sizes = calloc(n, sizeof(size_t));
    ff_ext2_elem_t *block_eigs = calloc(n, sizeof(ff_ext2_elem_t));

    size_t block_count = 0;
    size_t p_col = 0;

    if (!J || !block_sizes || !block_eigs)
    {
        ff_mat_ext2_free(J);
        free(block_sizes);
        free(block_eigs);
        return false;
    }

    ff_ext2_elem_t one;
    ff_ext2_one(ctx, &one);

    for (size_t ei = 0; ei < eig_count; ++ei)
    {
        bool duplicate = false;

        for (size_t ej = 0; ej < ei; ++ej)
        {
            if (ff_ext2_eq(&eigs[ei], &eigs[ej]))
            {
                duplicate = true;
                break;
            }
        }

        if (duplicate)
            continue;

        ff_ext2_elem_t lambda = eigs[ei];

        ff_mat_ext2_t *N = ffd_ext2_mat_sub_scalar_identity(A, &lambda);
        if (!N)
            goto fail;

        ff_mat_ext2_t **K = calloc(n + 2, sizeof(ff_mat_ext2_t *));
        size_t *dims = calloc(n + 2, sizeof(size_t));

        if (!K || !dims)
        {
            free(K);
            free(dims);
            ff_mat_ext2_free(N);
            goto fail;
        }

        ff_mat_ext2_t *Pk = ffd_ext2_mat_copy(N);
        int maxK = 0;

        if (!Pk)
        {
            free(K);
            free(dims);
            ff_mat_ext2_free(N);
            goto fail;
        }

        for (int k = 1; k <= (int)n; ++k)
        {
            ff_mat_ext2_t *basis = NULL;
            int nullity = ff_mat_ext2_nullspace(&basis, Pk);

            if (nullity < 0)
            {
                ff_mat_ext2_free(basis);
                ff_mat_ext2_free(Pk);
                free(K);
                free(dims);
                ff_mat_ext2_free(N);
                goto fail;
            }

            if (k == 1 && nullity == 0)
            {
                ff_mat_ext2_free(basis);
                maxK = 0;
                break;
            }

            if (k > 1 && (size_t)nullity == dims[k - 1])
            {
                ff_mat_ext2_free(basis);
                maxK = k - 1;
                break;
            }

            K[k] = basis;
            dims[k] = (size_t)nullity;
            maxK = k;

            if ((size_t)nullity == n)
            {
                break;
            }

            if (k < (int)n)
            {
                ff_mat_ext2_t *next = ff_mat_ext2_new(ctx, n, n);
                if (!next)
                {
                    ff_mat_ext2_free(Pk);
                    free(K);
                    free(dims);
                    ff_mat_ext2_free(N);
                    goto fail;
                }

                ff_mat_ext2_mul(next, Pk, N);

                ff_mat_ext2_free(Pk);
                Pk = next;
            }
        }

        ff_mat_ext2_free(Pk);

        if (maxK == 0 || dims[maxK] == 0)
        {
            for (int k = 1; k <= maxK; ++k)
            {
                ff_mat_ext2_free(K[k]);
            }

            free(K);
            free(dims);
            ff_mat_ext2_free(N);
            continue;
        }

        for (int k = maxK; k >= 1; --k)
        {
            size_t delta_k = dims[k] - dims[k - 1];
            size_t delta_k1 = (k < maxK) ? (dims[k + 1] - dims[k]) : 0;

            size_t exact = (delta_k > delta_k1) ? (delta_k - delta_k1) : 0;

            if (exact == 0)
                continue;

            ff_mat_ext2_t *S = NULL;

            if (k > 1 && K[k - 1])
            {
                S = ffd_ext2_mat_copy(K[k - 1]);
                if (!S)
                    goto fail_local;
            }

            ff_mat_ext2_t *Kkp1 = (k < maxK) ? K[k + 1] : K[maxK];

            ff_mat_ext2_t *img = ffd_ext2_matrix_image(N, Kkp1);
            ff_mat_ext2_t *img_basis = ffd_ext2_column_basis(img);

            ff_mat_ext2_free(img);

            if (img_basis)
            {
                ff_mat_ext2_t *S2 = ffd_ext2_basis_append_matrix(S, img_basis);
                ff_mat_ext2_free(S);
                ff_mat_ext2_free(img_basis);

                S = S2;

                if (!S)
                    goto fail_local;
            }

            ff_mat_ext2_t *C = NULL;
            ff_mat_ext2_t *S_work = S ? ffd_ext2_mat_copy(S) : NULL;

            size_t cols_Kk = K[k] ? ff_mat_ext2_cols(K[k]) : 0;

            for (size_t col = 0; col < cols_Kk; ++col)
            {
                if (C && ff_mat_ext2_cols(C) >= exact)
                    break;

                ff_ext2_elem_t *v = ffd_ext2_vec_copy_column(K[k], col);
                if (!v)
                {
                    ff_mat_ext2_free(C);
                    ff_mat_ext2_free(S_work);
                    goto fail_local;
                }

                if (ffd_ext2_vector_independent(S_work, ctx, n, v))
                {
                    ff_mat_ext2_t *C2 = ffd_ext2_basis_append_vector(C, ctx, n, v);
                    ff_mat_ext2_t *S3 = ffd_ext2_basis_append_vector(S_work, ctx, n, v);

                    ff_mat_ext2_free(C);
                    ff_mat_ext2_free(S_work);

                    C = C2;
                    S_work = S3;

                    if (!C || !S_work)
                    {
                        free(v);
                        ff_mat_ext2_free(C);
                        ff_mat_ext2_free(S_work);
                        goto fail_local;
                    }
                }

                free(v);
            }

            ff_mat_ext2_free(S_work);

            if (!C || ff_mat_ext2_cols(C) != exact)
            {
                ff_mat_ext2_free(C);
                goto fail_local;
            }

            for (size_t t = 0; t < exact; ++t)
            {
                if (p_col + k > n)
                {
                    ff_mat_ext2_free(C);
                    goto fail_local;
                }

                ff_ext2_elem_t *top = ffd_ext2_vec_copy_column(C, t);
                if (!top)
                {
                    ff_mat_ext2_free(C);
                    goto fail_local;
                }

                ff_ext2_elem_t *chain = calloc((size_t)k * n, sizeof(ff_ext2_elem_t));
                if (!chain)
                {
                    free(top);
                    ff_mat_ext2_free(C);
                    goto fail_local;
                }

                memcpy(&chain[(size_t)(k - 1) * n], top, n * sizeof(ff_ext2_elem_t));

                for (int idx = k - 2; idx >= 0; --idx)
                {
                    if (!ffd_ext2_mat_vec_mul(
                            N,
                            &chain[(size_t)(idx + 1) * n],
                            &chain[(size_t)idx * n]))
                    {
                        free(top);
                        free(chain);
                        ff_mat_ext2_free(C);
                        goto fail_local;
                    }
                }

                for (int s = 0; s < k; ++s)
                {
                    ff_mat_ext2_t *P2 = ffd_ext2_basis_append_vector(
                        P_basis,
                        ctx,
                        n,
                        &chain[(size_t)s * n]);

                    ff_mat_ext2_free(P_basis);
                    P_basis = P2;

                    if (!P_basis)
                    {
                        free(top);
                        free(chain);
                        ff_mat_ext2_free(C);
                        goto fail_local;
                    }
                }

                for (int s = 0; s < k; ++s)
                {
                    ff_ext2_elem_t *diag = ff_mat_ext2_at(J, p_col + (size_t)s, p_col + (size_t)s);
                    if (!diag)
                    {
                        free(top);
                        free(chain);
                        ff_mat_ext2_free(C);
                        goto fail_local;
                    }

                    *diag = lambda;

                    if (s < k - 1)
                    {
                        ff_ext2_elem_t *sup = ff_mat_ext2_at(
                            J,
                            p_col + (size_t)s,
                            p_col + (size_t)s + 1);

                        if (!sup)
                        {
                            free(top);
                            free(chain);
                            ff_mat_ext2_free(C);
                            goto fail_local;
                        }

                        *sup = one;
                    }
                }

                block_sizes[block_count] = (size_t)k;
                block_eigs[block_count] = lambda;
                block_count++;

                p_col += (size_t)k;

                free(top);
                free(chain);
            }

            ff_mat_ext2_free(C);
            ff_mat_ext2_free(S);

            continue;

        fail_local:
            ff_mat_ext2_free(S);

            for (int kk = 1; kk <= maxK; ++kk)
            {
                ff_mat_ext2_free(K[kk]);
            }

            free(K);
            free(dims);
            ff_mat_ext2_free(N);

            goto fail;
        }

        for (int k = 1; k <= maxK; ++k)
        {
            ff_mat_ext2_free(K[k]);
        }

        free(K);
        free(dims);
        ff_mat_ext2_free(N);
    }

    /* Allow partial decomposition: p_col can be <= n */
    if (p_col > n)
    {
        goto fail;
    }

    ff_mat_ext2_t *P_final = P_basis;
    P_basis = NULL;

    /*
        Verification:
            If p_col == n, verify P^{-1} A P == J
            If p_col < n, verify A * P == P * J_sub (partial decomposition)
    */
    bool verify_ok = true;

    if (p_col == n)
    {
        ff_mat_ext2_t *P_inv = ff_mat_ext2_new(ctx, n, n);
        ff_mat_ext2_t *T1 = ff_mat_ext2_new(ctx, n, n);
        ff_mat_ext2_t *T2 = ff_mat_ext2_new(ctx, n, n);
        if (!P_inv || !T1 || !T2)
        {
            verify_ok = false;
        }
        else
        {
            if (!ff_mat_ext2_inv(P_inv, P_final))
            {
                verify_ok = false;
            }
            else
            {
                ff_mat_ext2_mul(T1, P_inv, A);
                ff_mat_ext2_mul(T2, T1, P_final);
                for (size_t i = 0; i < n; ++i)
                {
                    for (size_t j = 0; j < n; ++j)
                    {
                        const ff_ext2_elem_t *x = ff_mat_ext2_atc(T2, i, j);
                        const ff_ext2_elem_t *y = ff_mat_ext2_atc(J, i, j);
                        if (!x || !y || !ff_ext2_eq(x, y))
                        {
                            verify_ok = false;
                            break;
                        }
                    }
                    if (!verify_ok)
                        break;
                }
            }
        }
        ff_mat_ext2_free(P_inv);
        ff_mat_ext2_free(T1);
        ff_mat_ext2_free(T2);
    }
    else
    {
        /* Partial decomposition: verify A * P == P * J_sub */
        ff_mat_ext2_t *AP = ff_mat_ext2_new(ctx, n, p_col);
        ff_mat_ext2_t *J_sub = ff_mat_ext2_new(ctx, p_col, p_col);
        ff_mat_ext2_t *P_J = ff_mat_ext2_new(ctx, n, p_col);

        if (!AP || !J_sub || !P_J)
        {
            verify_ok = false;
        }
        else
        {
            ff_mat_ext2_mul(AP, A, P_final);

            /* Copy the top-left p_col x p_col block of J into J_sub */
            for (size_t i = 0; i < p_col; ++i)
            {
                for (size_t j = 0; j < p_col; ++j)
                {
                    *ff_mat_ext2_at(J_sub, i, j) = *ff_mat_ext2_atc(J, i, j);
                }
            }

            ff_mat_ext2_mul(P_J, P_final, J_sub);

            for (size_t i = 0; i < n; ++i)
            {
                for (size_t j = 0; j < p_col; ++j)
                {
                    const ff_ext2_elem_t *x = ff_mat_ext2_atc(AP, i, j);
                    const ff_ext2_elem_t *y = ff_mat_ext2_atc(P_J, i, j);
                    if (!x || !y || !ff_ext2_eq(x, y))
                    {
                        verify_ok = false;
                        break;
                    }
                }
                if (!verify_ok)
                    break;
            }
        }
        ff_mat_ext2_free(AP);
        ff_mat_ext2_free(J_sub);
        ff_mat_ext2_free(P_J);
    }

    if (!verify_ok)
    {
        ff_mat_ext2_free(P_final);
        ff_mat_ext2_free(J);
        free(block_sizes);
        free(block_eigs);
        return false;
    }

    out->J = J;
    out->P = P_final;
    out->block_count = block_count;
    out->block_sizes = block_sizes;
    out->block_eigenvalues = block_eigs;
    return true;

fail:
    ff_mat_ext2_free(P_basis);
    ff_mat_ext2_free(J);

    free(block_sizes);
    free(block_eigs);

    return false;
}

bool ff_jordan_base_with_eigenvalues_ext2(
    const ff_ext2_t *E,
    const ff_mat_t *A_base,
    const ff_ext2_elem_t *eigs,
    size_t eig_count,
    ff_jordan_result_t *out)
{
    if (!E || !A_base || !eigs || !out)
        return false;

    ff_mat_ext2_t *A_ext = ffd_ext2_mat_from_base(E, A_base);
    if (!A_ext)
        return false;

    bool ok = ff_jordan_ext2_with_eigenvalues(E, A_ext, eigs, eig_count, out);

    ff_mat_ext2_free(A_ext);

    return ok;
}

/* ============================================================ */
/* ff_staircase implementation                                   */
/* ============================================================ */

static ff_mat_t *ffd_base_mat_copy(const ff_mat_t *A)
{
    if (!A)
        return NULL;

    const ff_prime_t *ctx = ff_mat_ctx(A);
    size_t rows = ff_mat_rows(A);
    size_t cols = ff_mat_cols(A);

    ff_mat_t *B = ff_mat_new(ctx, rows, cols);
    if (!B)
        return NULL;

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            const uint32_t *src = ff_mat_atc(A, i, j);
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

    ff_mat_t *H = ffd_base_mat_copy(A);
    ff_mat_t *P = ff_mat_new(ctx, n, n);

    ff_mat_t *P_inv = NULL;
    ff_mat_t *T1 = NULL;
    ff_mat_t *T2 = NULL;

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
                    continue;

                const uint32_t *pp = ff_mat_atc(H, r - 1, c);
                if (!pp || *pp == 0)
                {
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

                /*
                    Row operation:

                        row_r <- row_r - f row_{r-1}
                */
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

                /*
                    Column compensation:

                        col_{r-1} <- col_{r-1} + f col_r
                */
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

                /*
                    Accumulate P:

                        P <- P * (I + f e_{r,r-1})

                    which means:

                        col_{r-1} <- col_{r-1} + f col_r
                */
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

    /*
        Verify:

            P^{-1} A P == H
    */

    P_inv = ff_mat_new(ctx, n, n);
    T1 = ff_mat_new(ctx, n, n);
    T2 = ff_mat_new(ctx, n, n);

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
    {
        *P_out = P;
    }
    else
    {
        ff_mat_free(P);
    }

    if (H_out)
    {
        *H_out = H;
    }
    else
    {
        ff_mat_free(H);
    }

    return true;
}

/* ============================================================ */
/* ff_ratrec implementation                                      */
/* ============================================================ */

bool ff_ratrec_from_u128(
    ff_u128_t x,
    ff_u128_t M,
    ff_i128_t *num,
    ff_i128_t *den)
{
    if (!num || !den)
        return false;

    if (M == 0)
        return false;

    if (x >= M)
        x %= M;

    const ff_u128_t I128_MAX_U = (((ff_u128_t)1) << 127) - 1;

    if (M > I128_MAX_U)
    {
        return false;
    }

    if (x == 0)
    {
        *num = 0;
        *den = 1;
        return true;
    }

    ff_u128_t halfM = M >> 1;
    ff_u128_t limit_u = ffd_u128_sqrt(halfM);

    if (limit_u == 0)
        limit_u = 1;

    ff_i128_t limit = (ff_i128_t)limit_u;

    ff_i128_t r0 = (ff_i128_t)M;
    ff_i128_t r1 = (ff_i128_t)x;

    ff_i128_t s0 = 0;
    ff_i128_t s1 = 1;

    while (r1 > limit)
    {
        ff_i128_t q = r0 / r1;

        ff_i128_t r2 = r0 - q * r1;
        ff_i128_t s2 = s0 - q * s1;

        r0 = r1;
        r1 = r2;

        s0 = s1;
        s1 = s2;
    }

    if (r1 == 0)
    {
        return false;
    }

    ff_i128_t n = r1;
    ff_i128_t d = s1;

    if (d < 0)
    {
        n = -n;
        d = -d;
    }

    if (d <= 0)
        return false;

    if (ffd_i128_abs(n) > limit)
        return false;
    if (d > limit)
        return false;

    if (ffd_i128_gcd(n, d) != 1)
        return false;

    *num = n;
    *den = d;

    return true;
}

bool ff_ratrec_crt_i128(
    const ff_crt_t *crt,
    ff_i128_t *num,
    ff_i128_t *den)
{
    if (!crt)
        return false;

    ff_u128_t x = 0;
    ff_u128_t M = ff_crt_modulus(crt);

    if (!ff_crt_to_u128(&x, crt))
        return false;

    return ff_ratrec_from_u128(x, M, num, den);
}

bool ff_ratrec_crt_i64(
    const ff_crt_t *crt,
    int64_t *num,
    int64_t *den)
{
    if (!num || !den)
        return false;

    ff_i128_t n128 = 0;
    ff_i128_t d128 = 0;

    if (!ff_ratrec_crt_i128(crt, &n128, &d128))
    {
        return false;
    }

    const ff_i128_t I64_MIN = (ff_i128_t)INT64_MIN;
    const ff_i128_t I64_MAX = (ff_i128_t)INT64_MAX;

    if (n128 < I64_MIN || n128 > I64_MAX)
        return false;
    if (d128 < 1 || d128 > I64_MAX)
        return false;

    *num = (int64_t)n128;
    *den = (int64_t)d128;

    return true;
}

bool ff_ratrec_array_i128(
    const uint32_t *primes,
    const uint32_t *residues,
    size_t count,
    ff_i128_t *num,
    ff_i128_t *den)
{
    if (!primes || !residues || count == 0)
        return false;

    ff_crt_t *crt = ff_crt_new();
    if (!crt)
        return false;

    for (size_t i = 0; i < count; ++i)
    {
        if (!ff_crt_add_residue(crt, primes[i], residues[i]))
        {
            ff_crt_free(crt);
            return false;
        }
    }

    bool ok = ff_ratrec_crt_i128(crt, num, den);

    ff_crt_free(crt);

    return ok;
}
