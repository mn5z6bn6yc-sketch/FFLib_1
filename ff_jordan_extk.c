/*
 * ff_jordan_extk.c
 *
 * Generic Jordan engine over F_{p^k}: given a square matrix A over
 * F_{p^k} and supplied eigenvalues (as canonical extk elements,
 * deg < k), builds the Jordan form J and basis P with A * P = P * J.
 *
 * Partial decompositions are allowed and verified: when the supplied
 * eigenvalues do not account for the whole space, P is n x p_col with
 * p_col < n and the identity A * P == P * J_sub is checked; the full
 * case checks P^{-1} A P == J. Nothing unverified is returned.
 *
 * Split out of the old monolithic ff_extk.c. Changes vs. the original:
 *   - eig_count == 0 is rejected (previously a vacuous "success" with
 *     P == NULL and a zero J);
 *   - the unused E parameter is now validated (E must be A's context),
 *     mirroring ff_jordan_ext2_with_eigenvalues.
 */
#include "ff_extk_internal.h"

#include <stdlib.h>
#include <string.h>

/* -------------------------------------------------------------------- */
/* Local helpers                                                         */
/* -------------------------------------------------------------------- */

static ff_mat_extk_t *pkf_extk_mat_copy(const ff_mat_extk_t *A)
{
    if (!A)
        return NULL;

    ff_mat_extk_t *B = ff_mat_extk_new(A->ctx, A->rows, A->cols);
    if (!B)
        return NULL;

    for (size_t i = 0; i < A->rows * A->cols; ++i)
        pkf_extk_poly_copy_to(B->data[i], A->data[i]);

    return B;
}

static ff_mat_extk_t *pkf_extk_mat_sub_scalar_identity(
    const ff_mat_extk_t *A,
    const ff_poly_t *lambda)
{
    ff_mat_extk_t *N = pkf_extk_mat_copy(A);
    if (!N)
        return NULL;

    for (size_t i = 0; i < A->rows; ++i)
        ff_extk_sub(A->ctx, N->data[i * N->stride + i],
                    N->data[i * N->stride + i], lambda);

    return N;
}

static bool pkf_extk_mat_vec_mul(const ff_mat_extk_t *M,
                                 const ff_poly_t *const *x,
                                 ff_poly_t **y)
{
    ff_poly_t *temp = ff_poly_new(M->ctx->base, M->ctx->k - 1);
    ff_poly_t *sum = ff_poly_new(M->ctx->base, M->ctx->k - 1);

    if (!temp || !sum)
    {
        ff_poly_free(temp);
        ff_poly_free(sum);
        return false;
    }

    for (size_t i = 0; i < M->rows; ++i)
    {
        int sum_deg = ff_poly_deg(sum);
        for (int j = 0; j <= sum_deg; ++j)
            ff_poly_set_coeff(sum, j, 0);

        for (size_t j = 0; j < M->cols; ++j)
        {
            if (ff_extk_is_zero(M->data[i * M->stride + j]) ||
                ff_extk_is_zero(x[j]))
                continue;

            ff_extk_mul(M->ctx, temp, M->data[i * M->stride + j], x[j]);
            ff_extk_add(M->ctx, sum, sum, temp);
        }

        pkf_extk_poly_copy_to(y[i], sum);
    }

    ff_poly_free(temp);
    ff_poly_free(sum);
    return true;
}

static ff_poly_t **pkf_extk_vec_copy_column(const ff_mat_extk_t *M, size_t col)
{
    ff_poly_t **v = calloc(M->rows, sizeof(ff_poly_t *));
    if (!v)
        return NULL;

    for (size_t i = 0; i < M->rows; ++i)
    {
        v[i] = ff_poly_new(M->ctx->base, M->ctx->k - 1);
        if (!v[i])
        {
            for (size_t j = 0; j < i; ++j)
                ff_poly_free(v[j]);
            free(v);
            return NULL;
        }
        pkf_extk_poly_copy_to(v[i], M->data[i * M->stride + col]);
    }

    return v;
}

static void pkf_extk_vec_free(ff_poly_t **v, size_t n)
{
    if (!v)
        return;
    for (size_t i = 0; i < n; ++i)
        ff_poly_free(v[i]);
    free(v);
}

static ff_mat_extk_t *pkf_extk_basis_append_vector(
    const ff_extk_t *ctx,
    ff_mat_extk_t *B,
    size_t n,
    ff_poly_t **v)
{
    size_t old_cols = B ? B->cols : 0;

    ff_mat_extk_t *N = ff_mat_extk_new(ctx, n, old_cols + 1);
    if (!N)
        return NULL;

    for (size_t j = 0; j < old_cols; ++j)
        for (size_t i = 0; i < n; ++i)
            pkf_extk_poly_copy_to(N->data[i * N->stride + j],
                                  B->data[i * B->stride + j]);

    for (size_t i = 0; i < n; ++i)
        pkf_extk_poly_copy_to(N->data[i * N->stride + old_cols], v[i]);

    return N;
}

static bool pkf_extk_vector_independent(
    const ff_extk_t *ctx,
    ff_mat_extk_t *B,
    size_t n,
    ff_poly_t **v)
{
    size_t cols = B ? B->cols : 0;

    if (cols + 1 > n)
        return false;

    /* Empty basis: independent iff nonzero. */
    if (!B)
    {
        for (size_t i = 0; i < n; ++i)
        {
            if (!ff_extk_is_zero(v[i]))
                return true;
        }
        return false;
    }

    ff_mat_extk_t *aug = ff_mat_extk_new(ctx, n, cols + 1);
    if (!aug)
        return false;

    for (size_t j = 0; j < cols; ++j)
        for (size_t i = 0; i < n; ++i)
            pkf_extk_poly_copy_to(aug->data[i * aug->stride + j],
                                  B->data[i * B->stride + j]);

    for (size_t i = 0; i < n; ++i)
        pkf_extk_poly_copy_to(aug->data[i * aug->stride + cols], v[i]);

    ff_mat_extk_t *ns = NULL;
    int nullity = ff_mat_extk_nullspace(&ns, aug);

    ff_mat_extk_free(aug);
    ff_mat_extk_free(ns);

    return nullity == 0;
}

static ff_mat_extk_t *pkf_extk_column_basis(const ff_mat_extk_t *M)
{
    ff_mat_extk_t *B = NULL;

    for (size_t j = 0; j < M->cols; ++j)
    {
        ff_poly_t **v = pkf_extk_vec_copy_column(M, j);
        if (!v)
        {
            ff_mat_extk_free(B);
            return NULL;
        }

        if (pkf_extk_vector_independent(M->ctx, B, M->rows, v))
        {
            ff_mat_extk_t *N = pkf_extk_basis_append_vector(M->ctx, B, M->rows, v);
            ff_mat_extk_free(B);
            B = N;
        }

        pkf_extk_vec_free(v, M->rows);

        if (!B)
            return NULL;
    }

    return B;
}

static ff_mat_extk_t *pkf_extk_matrix_image(const ff_mat_extk_t *A,
                                            const ff_mat_extk_t *basis)
{
    ff_mat_extk_t *img = ff_mat_extk_new(A->ctx, A->rows, basis->cols);
    if (!img)
        return NULL;

    ff_mat_extk_mul(img, A, basis);
    return img;
}

static ff_mat_extk_t *pkf_extk_basis_append_matrix(ff_mat_extk_t *B,
                                                   const ff_mat_extk_t *Add)
{
    for (size_t j = 0; j < Add->cols; ++j)
    {
        ff_poly_t **v = pkf_extk_vec_copy_column(Add, j);
        if (!v)
        {
            ff_mat_extk_free(B);
            return NULL;
        }

        ff_mat_extk_t *N = pkf_extk_basis_append_vector(Add->ctx, B, Add->rows, v);

        pkf_extk_vec_free(v, Add->rows);
        ff_mat_extk_free(B);
        B = N;

        if (!B)
            return NULL;
    }

    return B;
}

/* -------------------------------------------------------------------- */
/* Result lifecycle                                                      */
/* -------------------------------------------------------------------- */

/* Only call on results produced by ff_jordan_extk_with_eigenvalues or
 * on zero-initialized structs. Frees the eigenvalue polys AND the
 * array, the matrices, and the size array, then zeroes the struct. */
void ff_jordan_extk_result_clear(ff_jordan_extk_result_t *res)
{
    if (!res)
        return;

    ff_mat_extk_free(res->J);
    ff_mat_extk_free(res->P);
    free(res->block_sizes);

    if (res->block_eigenvalues)
    {
        for (size_t i = 0; i < res->block_count; ++i)
            ff_poly_free(res->block_eigenvalues[i]);

        free(res->block_eigenvalues);
    }

    memset(res, 0, sizeof(*res));
}

/* -------------------------------------------------------------------- */
/* Jordan engine over F_{p^k}                                            */
/* -------------------------------------------------------------------- */

bool ff_jordan_extk_with_eigenvalues(
    const ff_extk_t *E,
    const ff_mat_extk_t *A,
    const ff_poly_t *const *eigs,
    size_t eig_count,
    ff_jordan_extk_result_t *out)
{
    if (!E || !A || !eigs || !out)
        return false;
    if (eig_count == 0)
        return false; /* nothing to decompose against */

    memset(out, 0, sizeof(*out));

    const ff_extk_t *ctx = A->ctx;
    if (!ctx || ctx != E)
        return false;

    size_t n = A->rows;
    if (A->cols != n || n == 0)
        return false;

    ff_mat_extk_t *P_basis = NULL;
    ff_mat_extk_t *J = ff_mat_extk_new(ctx, n, n);
    size_t *block_sizes = calloc(n, sizeof(size_t));
    ff_poly_t **block_eigs = calloc(n, sizeof(ff_poly_t *));
    size_t block_count = 0;
    size_t p_col = 0;

    if (!J || !block_sizes || !block_eigs)
    {
        ff_mat_extk_free(J);
        free(block_sizes);
        free(block_eigs);
        return false;
    }

    ff_poly_t *one = ff_poly_new(ctx->base, ctx->k - 1);
    if (!one)
    {
        ff_mat_extk_free(J);
        free(block_sizes);
        free(block_eigs);
        return false;
    }
    ff_poly_set_coeff(one, 0, ff_one_rep(ctx->base));

    for (size_t ei = 0; ei < eig_count; ++ei)
    {
        bool duplicate = false;
        for (size_t ej = 0; ej < ei; ++ej)
        {
            if (ff_extk_eq(eigs[ei], eigs[ej]))
            {
                duplicate = true;
                break;
            }
        }
        if (duplicate)
            continue;

        const ff_poly_t *lambda = eigs[ei];

        /* N = A - lambda I */
        ff_mat_extk_t *N = pkf_extk_mat_sub_scalar_identity(A, lambda);
        if (!N)
            goto fail;

        /* K[k] = ker(N^k); dims stabilizes at maxK. */
        ff_mat_extk_t **K = calloc(n + 2, sizeof(ff_mat_extk_t *));
        size_t *dims = calloc(n + 2, sizeof(size_t));
        ff_mat_extk_t *Pk = pkf_extk_mat_copy(N);
        int maxK = 0;

        if (!K || !dims || !Pk)
        {
            ff_mat_extk_free(Pk);
            free(K);
            free(dims);
            ff_mat_extk_free(N);
            goto fail;
        }

        for (int k = 1; k <= (int)n; ++k)
        {
            ff_mat_extk_t *basis = NULL;
            int nullity = ff_mat_extk_nullspace(&basis, Pk);

            if (k == 1 && nullity == 0)
            {
                ff_mat_extk_free(basis);
                break;
            }

            if (k > 1 && (size_t)nullity == dims[k - 1])
            {
                ff_mat_extk_free(basis);
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
                ff_mat_extk_t *next = ff_mat_extk_new(ctx, n, n);
                if (!next)
                {
                    ff_mat_extk_free(Pk);
                    goto fail_k;
                }
                ff_mat_extk_mul(next, Pk, N);
                ff_mat_extk_free(Pk);
                Pk = next;
            }
        }

        ff_mat_extk_free(Pk);

        for (int k = maxK; k >= 1; --k)
        {
            size_t delta_k = dims[k] - dims[k - 1];
            size_t delta_k1 = (k < maxK) ? (dims[k + 1] - dims[k]) : 0;
            size_t exact = (delta_k > delta_k1) ? (delta_k - delta_k1) : 0;

            if (exact == 0)
                continue;

            /* Exclusion space S = K[k-1] + N(K[k+1]). */
            ff_mat_extk_t *S = (k > 1 && K[k - 1]) ? pkf_extk_mat_copy(K[k - 1]) : NULL;
            ff_mat_extk_t *Kkp1 = (k < maxK) ? K[k + 1] : K[maxK];
            ff_mat_extk_t *img = pkf_extk_matrix_image(N, Kkp1);
            ff_mat_extk_t *img_basis = img ? pkf_extk_column_basis(img) : NULL;
            ff_mat_extk_free(img);

            if (img_basis)
            {
                S = pkf_extk_basis_append_matrix(S, img_basis);
                ff_mat_extk_free(img_basis);
            }

            /* Collect `exact` tops from K[k] independent of S. */
            ff_mat_extk_t *C = NULL;
            ff_mat_extk_t *S_work = S ? pkf_extk_mat_copy(S) : NULL;

            for (size_t col = 0; col < K[k]->cols; ++col)
            {
                if (C && C->cols >= exact)
                    break;

                ff_poly_t **v = pkf_extk_vec_copy_column(K[k], col);
                if (!v)
                {
                    ff_mat_extk_free(C);
                    ff_mat_extk_free(S_work);
                    ff_mat_extk_free(S);
                    goto fail_k;
                }

                if (pkf_extk_vector_independent(ctx, S_work, n, v))
                {
                    C = pkf_extk_basis_append_vector(ctx, C, n, v);
                    S_work = pkf_extk_basis_append_vector(ctx, S_work, n, v);
                }

                pkf_extk_vec_free(v, n);
            }

            ff_mat_extk_free(S_work);
            ff_mat_extk_free(S);

            /* Build chains from the tops: chain[k-1] = top,
             * chain[i] = N * chain[i+1]. */
            for (size_t t = 0; C && t < exact; ++t)
            {
                ff_poly_t **top = pkf_extk_vec_copy_column(C, t);
                if (!top)
                {
                    ff_mat_extk_free(C);
                    goto fail_k;
                }

                ff_poly_t ***chain = calloc((size_t)k, sizeof(ff_poly_t **));
                if (!chain)
                {
                    pkf_extk_vec_free(top, n);
                    ff_mat_extk_free(C);
                    goto fail_k;
                }

                for (int idx = 0; idx < k; ++idx)
                {
                    chain[idx] = calloc(n, sizeof(ff_poly_t *));
                    if (!chain[idx])
                        goto fail_chain_alloc;

                    for (size_t i = 0; i < n; ++i)
                    {
                        chain[idx][i] = ff_poly_new(ctx->base, ctx->k - 1);
                        if (!chain[idx][i])
                            goto fail_chain_alloc;
                    }
                }

                for (size_t i = 0; i < n; ++i)
                    pkf_extk_poly_copy_to(chain[k - 1][i], top[i]);

                for (int idx = k - 2; idx >= 0; --idx)
                    pkf_extk_mat_vec_mul(N,
                                         (const ff_poly_t *const *)chain[idx + 1],
                                         chain[idx]);

                for (int s = 0; s < k; ++s)
                {
                    P_basis = pkf_extk_basis_append_vector(ctx, P_basis, n, chain[s]);
                    pkf_extk_poly_copy_to(J->data[(p_col + (size_t)s) * J->stride + (p_col + (size_t)s)],
                                          lambda);

                    if (s < k - 1)
                        pkf_extk_poly_copy_to(J->data[(p_col + (size_t)s) * J->stride + (p_col + (size_t)s + 1)],
                                              one);
                }

                block_sizes[block_count] = (size_t)k;
                block_eigs[block_count] = ff_poly_copy(lambda);
                block_count++;
                p_col += (size_t)k;

                for (int idx = 0; idx < k; ++idx)
                    pkf_extk_vec_free(chain[idx], n);
                free(chain);
                pkf_extk_vec_free(top, n);
                continue;

            fail_chain_alloc:
                if (chain)
                {
                    for (int idx = 0; idx < k; ++idx)
                        pkf_extk_vec_free(chain[idx], n);
                    free(chain);
                }
                pkf_extk_vec_free(top, n);
                ff_mat_extk_free(C);
                goto fail_k;
            }

            ff_mat_extk_free(C);
        }

    fail_k:
        for (int k = 1; k <= maxK; ++k)
            ff_mat_extk_free(K[k]);
        free(K);
        free(dims);
        ff_mat_extk_free(N);
        goto fail;
    }

    ff_poly_free(one);

    if (p_col > n)
        goto fail_bare;

    /* Verification. */
    bool verify_ok = true;
    ff_mat_extk_t *AP = ff_mat_extk_new(ctx, n, p_col);
    ff_mat_extk_t *J_sub = ff_mat_extk_new(ctx, p_col, p_col);
    ff_mat_extk_t *PJ = ff_mat_extk_new(ctx, n, p_col);

    if (!AP || !J_sub || !PJ)
    {
        verify_ok = false;
    }
    else
    {
        ff_mat_extk_mul(AP, A, P_basis);

        for (size_t i = 0; i < p_col; ++i)
            for (size_t j = 0; j < p_col; ++j)
                pkf_extk_poly_copy_to(J_sub->data[i * J_sub->stride + j],
                                      J->data[i * J->stride + j]);

        ff_mat_extk_mul(PJ, P_basis, J_sub);

        for (size_t i = 0; i < n && verify_ok; ++i)
            for (size_t j = 0; j < p_col && verify_ok; ++j)
            {
                if (!ff_extk_eq(AP->data[i * AP->stride + j],
                                PJ->data[i * PJ->stride + j]))
                    verify_ok = false;
            }
    }

    ff_mat_extk_free(AP);
    ff_mat_extk_free(J_sub);
    ff_mat_extk_free(PJ);

    if (!verify_ok)
        goto fail_bare;

    out->J = J;
    out->P = P_basis;
    out->block_count = block_count;
    out->block_sizes = block_sizes;
    out->block_eigenvalues = block_eigs;
    return true;

fail:
    ff_poly_free(one);

fail_bare:
    ff_mat_extk_free(P_basis);
    ff_mat_extk_free(J);
    free(block_sizes);

    for (size_t i = 0; i < block_count; ++i)
        ff_poly_free(block_eigs[i]);
    free(block_eigs);

    return false;
}
