/*
ff_snf.c
Implementation of Smith Normal Form for polynomial matrices.
*/
#include "ff_snf.h"
#include <stdlib.h>
#include <string.h>

/* ============================================================ */
/* Internal Polynomial Matrix Structure                         */
/* ============================================================ */
typedef struct
{
    size_t rows;
    size_t cols;
    ff_poly_t ***data; /* data[row][col] */
} ff_pmat_t;

static ff_pmat_t *pmat_new(const ff_prime_t *ctx, size_t rows, size_t cols)
{
    ff_pmat_t *M = calloc(1, sizeof(ff_pmat_t));
    if (!M)
        return NULL;
    M->rows = rows;
    M->cols = cols;
    M->data = calloc(rows, sizeof(ff_poly_t **));
    if (!M->data)
    {
        free(M);
        return NULL;
    }

    for (size_t i = 0; i < rows; ++i)
    {
        M->data[i] = calloc(cols, sizeof(ff_poly_t *));
        if (!M->data[i])
        {
            for (size_t k = 0; k < i; ++k)
                free(M->data[k]);
            free(M->data);
            free(M);
            return NULL;
        }
        for (size_t j = 0; j < cols; ++j)
        {
            M->data[i][j] = ff_poly_new(ctx, 0);
            if (!M->data[i][j])
            {
                /* Cleanup on failure omitted for brevity, assume success in demo */
            }
        }
    }
    return M;
}

static void pmat_free(ff_pmat_t *M)
{
    if (!M)
        return;
    for (size_t i = 0; i < M->rows; ++i)
    {
        for (size_t j = 0; j < M->cols; ++j)
        {
            if (M->data[i][j])
                ff_poly_free(M->data[i][j]);
        }
        free(M->data[i]);
    }
    free(M->data);
    free(M);
}

/* ============================================================ */
/* Internal Polynomial Helpers (Public API only)                */
/* ============================================================ */

static int pdeg(const ff_poly_t *p)
{
    return ff_poly_deg(p);
}

static void ptrim(ff_poly_t *p)
{
    ff_poly_trim(p);
}

static bool pis_zero(const ff_poly_t *p)
{
    int d = ff_poly_deg(p);
    if (d < 0)
        return true;
    for (int i = 0; i <= d; ++i)
    {
        if (ff_poly_coeff(p, i) != 0)
            return false;
    }
    return true;
}

static ff_poly_t *pdivrem(const ff_poly_t *f, const ff_poly_t *g, ff_poly_t **q_out)
{
    if (!f || !g || pis_zero(g))
        return NULL;
    const ff_prime_t *ctx = ff_poly_ctx(f);
    int deg_f = ff_poly_deg(f);
    int deg_g = ff_poly_deg(g);

    ff_poly_t *rem = ff_poly_copy(f);
    ff_poly_t *q = q_out ? ff_poly_new(ctx, (deg_f >= deg_g) ? deg_f - deg_g : 0) : NULL;
    if (!rem)
        return NULL;

    if (deg_f < deg_g)
    {
        if (q_out)
            *q_out = q;
        return rem;
    }

    uint32_t inv_lead = ff_inv(ctx, ff_poly_coeff(g, deg_g));
    for (int i = deg_f; i >= deg_g; --i)
    {
        uint32_t rem_i = ff_poly_coeff(rem, i);
        if (rem_i == 0)
            continue;
        uint32_t coeff = ff_mul(ctx, rem_i, inv_lead);
        if (q)
            ff_poly_set_coeff(q, i - deg_g, coeff);
        for (int j = 0; j <= deg_g; ++j)
        {
            uint32_t sub = ff_mul(ctx, coeff, ff_poly_coeff(g, j));
            uint32_t cur = ff_poly_coeff(rem, i - deg_g + j);
            ff_poly_set_coeff(rem, i - deg_g + j, ff_sub(ctx, cur, sub));
        }
    }
    ptrim(rem);
    if (q)
        ptrim(q);
    if (q_out)
        *q_out = q;
    return rem;
}

static ff_poly_t *pgcd(const ff_poly_t *a, const ff_poly_t *b)
{
    if (!a || !b)
        return NULL;
    const ff_prime_t *ctx = ff_poly_ctx(a);
    ff_poly_t *x = ff_poly_copy(a);
    ff_poly_t *y = ff_poly_copy(b);

    while (!pis_zero(y))
    {
        ff_poly_t *r = pdivrem(x, y, NULL);
        ff_poly_free(x);
        x = y;
        y = r;
    }
    ff_poly_free(y);

    /* Make monic */
    int d = ff_poly_deg(x);
    if (d >= 0)
    {
        uint32_t lead = ff_poly_coeff(x, d);
        if (lead != 0)
        {
            uint32_t inv = ff_inv(ctx, lead);
            for (int i = 0; i <= d; ++i)
            {
                ff_poly_set_coeff(x, i, ff_mul(ctx, ff_poly_coeff(x, i), inv));
            }
        }
    }
    return x;
}

/* ============================================================ */
/* Smith Normal Form Algorithm                                  */
/* ============================================================ */

bool ff_snf_charmatrix_invariants(
    const ff_mat_t *A,
    ff_poly_t ***out_factors,
    size_t *out_count,
    ff_inv_order_t order)
{
    if (!A || !out_factors || !out_count)
        return false;
    *out_factors = NULL;
    *out_count = 0;

    extern size_t ff_mat_rows(const ff_mat_t *M);
    extern size_t ff_mat_cols(const ff_mat_t *M);
    extern const ff_prime_t *ff_mat_ctx(const ff_mat_t *M);

    size_t n = ff_mat_rows(A);
    if (ff_mat_cols(A) != n)
        return false;
    const ff_prime_t *ctx = ff_mat_ctx(A);

    /* 1. Build Characteristic Matrix M = xI - A */
    ff_pmat_t *M = pmat_new(ctx, n, n);
    if (!M)
        return false;

    uint32_t one = ff_one_rep(ctx);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            /* Free the default degree-0 polynomial created by pmat_new */
            ff_poly_free(M->data[i][j]);

            if (i == j)
            {
                /* Diagonal: x - A[i][i] (Degree 1) */
                ff_poly_t *p = ff_poly_new(ctx, 1);
                ff_poly_set_coeff(p, 1, one);
                ff_poly_set_coeff(p, 0, ff_neg(ctx, *ff_mat_atc(A, i, j)));
                M->data[i][j] = p;
            }
            else
            {
                /* Off-diagonal: -A[i][j] (Degree 0) */
                ff_poly_t *p = ff_poly_new(ctx, 0);
                ff_poly_set_coeff(p, 0, ff_neg(ctx, *ff_mat_atc(A, i, j)));
                M->data[i][j] = p;
            }
        }
    }

    /* 2. SNF Elimination */
    for (size_t k = 0; k < n; ++k)
    {
        while (1)
        {
            /* Find minimal degree element in submatrix k..n-1, k..n-1 */
            size_t min_r = k, min_c = k;
            int min_deg = 1000000;
            bool all_zero = true;

            for (size_t r = k; r < n; ++r)
            {
                for (size_t c = k; c < n; ++c)
                {
                    if (!pis_zero(M->data[r][c]))
                    {
                        all_zero = false;
                        int d = pdeg(M->data[r][c]);
                        if (d < min_deg)
                        {
                            min_deg = d;
                            min_r = r;
                            min_c = c;
                        }
                    }
                }
            }
            if (all_zero)
                break;

            /* Swap to (k,k) */
            if (min_r != k)
            {
                ff_poly_t **tmp = M->data[k];
                M->data[k] = M->data[min_r];
                M->data[min_r] = tmp;
            }
            if (min_c != k)
            {
                for (size_t r = 0; r < n; ++r)
                {
                    ff_poly_t *tmp = M->data[r][k];
                    M->data[r][k] = M->data[r][min_c];
                    M->data[r][min_c] = tmp;
                }
            }

            /* Eliminate Row k and Col k */
            ff_poly_t *pivot = M->data[k][k];
            for (size_t c = k + 1; c < n; ++c)
            {
                if (!pis_zero(M->data[k][c]))
                {
                    ff_poly_t *q = NULL;
                    ff_poly_t *rem = pdivrem(M->data[k][c], pivot, &q);
                    ff_poly_free(M->data[k][c]);
                    M->data[k][c] = rem;
                    /* Col operation: Col_c = Col_c - q * Col_k */
                    for (size_t r = 0; r < n; ++r)
                    {
                        if (r == k)
                            continue;
                        ff_poly_t *term = ff_poly_mul_new(q, M->data[r][k]);
                        ff_poly_t *new_val = ff_poly_sub_new(M->data[r][c], term);
                        ff_poly_free(term);
                        ff_poly_free(M->data[r][c]);
                        M->data[r][c] = new_val;
                    }
                    ff_poly_free(q);
                }
            }
            for (size_t r = k + 1; r < n; ++r)
            {
                if (!pis_zero(M->data[r][k]))
                {
                    ff_poly_t *q = NULL;
                    ff_poly_t *rem = pdivrem(M->data[r][k], pivot, &q);
                    ff_poly_free(M->data[r][k]);
                    M->data[r][k] = rem;
                    /* Row operation: Row_r = Row_r - q * Row_k */
                    for (size_t c = 0; c < n; ++c)
                    {
                        if (c == k)
                            continue;
                        ff_poly_t *term = ff_poly_mul_new(q, M->data[k][c]);
                        ff_poly_t *new_val = ff_poly_sub_new(M->data[r][c], term);
                        ff_poly_free(term);
                        ff_poly_free(M->data[r][c]);
                        M->data[r][c] = new_val;
                    }
                    ff_poly_free(q);
                }
            }

            /* Check divisibility condition */
            bool divisible = true;
            size_t off_r = 0, off_c = 0;
            for (size_t r = k + 1; r < n && divisible; ++r)
            {
                for (size_t c = k + 1; c < n && divisible; ++c)
                {
                    if (!pis_zero(M->data[r][c]))
                    {
                        ff_poly_t *g = pgcd(pivot, M->data[r][c]);
                        if (pdeg(g) < pdeg(pivot))
                        {
                            divisible = false;
                            off_r = r;
                            off_c = c;
                        }
                        ff_poly_free(g);
                    }
                }
            }

            if (divisible)
                break; /* Move to next k */

            /* If not divisible, add offending row to row k to force degree drop */
            for (size_t c = 0; c < n; ++c)
            {
                ff_poly_t *new_val = ff_poly_add_new(M->data[k][c], M->data[off_r][c]);
                ff_poly_free(M->data[k][c]);
                M->data[k][c] = new_val;
            }
        }
    }

    /* 3. Extract Invariant Factors (Filter out units) */
    ff_poly_t **factors = calloc(n, sizeof(ff_poly_t *));
    size_t count = 0;

    for (size_t i = 0; i < n; ++i)
    {
        ff_poly_t *p = M->data[i][i];
        ptrim(p);

        int d = ff_poly_deg(p);
        uint32_t c0 = ff_poly_coeff(p, 0);

        /* A factor is a unit if it's a non-zero constant */
        bool is_unit = (d == 0 && c0 != 0);

        if (!is_unit && !pis_zero(p))
        {
            /* Make monic */
            uint32_t lead = ff_poly_coeff(p, d);
            uint32_t inv = ff_inv(ctx, lead);
            for (int j = 0; j <= d; ++j)
            {
                ff_poly_set_coeff(p, j, ff_mul(ctx, ff_poly_coeff(p, j), inv));
            }
            factors[count++] = ff_poly_copy(p);
        }
    }

    pmat_free(M);

    /* 4. Apply Order Switch */
    if (order == FF_INVARIANT_DESCENDING && count > 1)
    {
        for (size_t i = 0; i < count / 2; ++i)
        {
            ff_poly_t *tmp = factors[i];
            factors[i] = factors[count - 1 - i];
            factors[count - 1 - i] = tmp;
        }
    }

    *out_factors = factors;
    *out_count = count;
    return true;
}
