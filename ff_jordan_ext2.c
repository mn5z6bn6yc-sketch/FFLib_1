/*
 * ff_jordan_ext2.c
 *
 * Jordan decomposition over the quadratic extension F_{p^2}:
 * given a square matrix A over F_{p^2} and its eigenvalues (supplied by
 * the caller), builds J (Jordan form) and P (change of basis) with
 *
 *     A * P = P * J          (equivalently J = P^{-1} A P when full)
 *
 * Partial decompositions (fewer eigenvalues supplied than needed) are
 * allowed: then P is n x p_col with p_col < n, J is n x n with the
 * Jordan blocks in its top-left p_col x p_col corner, and the identity
 * A * P == P * J_sub is verified. In the full case P^{-1} A P == J is
 * verified directly. Either way, no unverified result is returned.
 *
 * Split out of the old concatenated ff_decomp.c. Change vs. the
 * original: eig_count == 0 is rejected up front (it used to produce a
 * vacuously-"successful" empty decomposition with P == NULL).
 */
#include "ff_decomp.h"

#include <stdlib.h>
#include <string.h>

/* -------------------------------------------------------------------- */
/* Local helpers                                                         */
/* -------------------------------------------------------------------- */

static bool ffd_ext2_vec_is_zero(const ff_ext2_elem_t *v, size_t n)
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

/* Embed a base-field matrix into F_{p^2} (imaginary part zero).
 * Assumes A's elements are in the base field's active representation. */
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

/* -------------------------------------------------------------------- */
/* Result lifecycle                                                      */
/* -------------------------------------------------------------------- */

/* Only call on results produced by ff_jordan_*_with_eigenvalues or on
 * zero-initialized structs. */
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

/* -------------------------------------------------------------------- */
/* Jordan engine over F_{p^2}                                            */
/* -------------------------------------------------------------------- */

bool ff_jordan_ext2_with_eigenvalues(
    const ff_ext2_t *E,
    const ff_mat_ext2_t *A,
    const ff_ext2_elem_t *eigs,
    size_t eig_count,
    ff_jordan_result_t *out)
{
    if (!A || !eigs || !out)
        return false;
    if (eig_count == 0)
        return false; /* nothing to decompose against */

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
        /* Skip duplicate eigenvalues. */
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

        /* N = A - lambda I */
        ff_mat_ext2_t *N = ffd_ext2_mat_sub_scalar_identity(A, &lambda);
        if (!N)
            goto fail;

        /* K[k] = ker(N^k), dims[k] = dim K[k]. The sequence of dims is
         * nondecreasing and stabilizes at maxK; the number of Jordan
         * blocks of size EXACTLY k is
         *     (dims[k] - dims[k-1]) - (dims[k+1] - dims[k]). */
        ff_mat_ext2_t **K = calloc(n + 2, sizeof(ff_mat_ext2_t *));
        size_t *dims = calloc(n + 2, sizeof(size_t));

        if (!K || !dims)
        {
            free(K);
            free(dims);
            ff_mat_ext2_free(N);
            goto fail;
        }

        ff_mat_ext2_t *Pk = ffd_ext2_mat_copy(N); /* Pk = N^k */
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
                /* lambda is not an eigenvalue at all. */
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
                break;

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
                ff_mat_ext2_free(K[k]);

            free(K);
            free(dims);
            ff_mat_ext2_free(N);
            continue;
        }

        /* Walk block sizes from largest to smallest. */
        for (int k = maxK; k >= 1; --k)
        {
            size_t delta_k = dims[k] - dims[k - 1];
            size_t delta_k1 = (k < maxK) ? (dims[k + 1] - dims[k]) : 0;
            size_t exact = (delta_k > delta_k1) ? (delta_k - delta_k1) : 0;

            if (exact == 0)
                continue;

            /* Build the exclusion space S = K[k-1] + N(K[k+1]) that the
             * chain TOPS must avoid. */
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

            /* Collect `exact` chain tops from K[k] independent of S. */
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

            /* For each top, build the chain downward:
             * chain[k-1] = top, chain[i] = N * chain[i+1],
             * so A * chain[s+1] = chain[s] + lambda * chain[s+1]. */
            for (size_t t = 0; t < exact; ++t)
            {
                if (p_col + (size_t)k > n)
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

                /* Append chain vectors as columns of P (chain[0] first). */
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

                /* J block: lambda on the diagonal, 1 on the superdiagonal. */
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
                ff_mat_ext2_free(K[kk]);

            free(K);
            free(dims);
            ff_mat_ext2_free(N);

            goto fail;
        }

        for (int k = 1; k <= maxK; ++k)
            ff_mat_ext2_free(K[k]);

        free(K);
        free(dims);
        ff_mat_ext2_free(N);
    }

    if (p_col > n)
        goto fail;

    ff_mat_ext2_t *P_final = P_basis;
    P_basis = NULL;

    /* Verification (nothing unverified is ever returned). */
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
        else if (!ff_mat_ext2_inv(P_inv, P_final))
        {
            verify_ok = false;
        }
        else
        {
            ff_mat_ext2_mul(T1, P_inv, A);
            ff_mat_ext2_mul(T2, T1, P_final);

            for (size_t i = 0; i < n && verify_ok; ++i)
                for (size_t j = 0; j < n && verify_ok; ++j)
                {
                    const ff_ext2_elem_t *x = ff_mat_ext2_atc(T2, i, j);
                    const ff_ext2_elem_t *y = ff_mat_ext2_atc(J, i, j);
                    if (!x || !y || !ff_ext2_eq(x, y))
                        verify_ok = false;
                }
        }

        ff_mat_ext2_free(P_inv);
        ff_mat_ext2_free(T1);
        ff_mat_ext2_free(T2);
    }
    else
    {
        /* Partial decomposition: verify A * P == P * J_sub. */
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

            for (size_t i = 0; i < p_col; ++i)
                for (size_t j = 0; j < p_col; ++j)
                {
                    ff_ext2_elem_t *dst = ff_mat_ext2_at(J_sub, i, j);
                    const ff_ext2_elem_t *src = ff_mat_ext2_atc(J, i, j);
                    if (!dst || !src)
                    {
                        verify_ok = false;
                        break;
                    }
                    *dst = *src;
                }

            if (verify_ok)
            {
                ff_mat_ext2_mul(P_J, P_final, J_sub);

                for (size_t i = 0; i < n && verify_ok; ++i)
                    for (size_t j = 0; j < p_col && verify_ok; ++j)
                    {
                        const ff_ext2_elem_t *x = ff_mat_ext2_atc(AP, i, j);
                        const ff_ext2_elem_t *y = ff_mat_ext2_atc(P_J, i, j);
                        if (!x || !y || !ff_ext2_eq(x, y))
                            verify_ok = false;
                    }
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

/* -------------------------------------------------------------------- */
/* Base-field convenience wrapper                                        */
/* -------------------------------------------------------------------- */

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
