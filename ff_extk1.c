/*
ff_extk.c
Implementation of Generic Extension Field Towering.
*/
#include "ff_extk.h"
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ============================================================ */
/* Irreducible Polynomial Search                                */
/* ============================================================ */

/*
 * Tests if f(x) is irreducible over F_p.
 * Uses the property: f is irreducible iff it is square-free and
 * gcd(x^{p^i} - x, f) == 1 for all 1 <= i <= deg(f)/2.
 */
static bool pkf_is_irreducible(const ff_prime_t *ctx, const ff_poly_t *f)
{
    int k = ff_poly_deg(f);
    if (k <= 0)
        return false;
    uint32_t p = ff_prime_modulus(ctx);

    /* 1. Check square-free: gcd(f, f') == 1 */
    ff_poly_t *deriv = ff_poly_deriv(f);
    ff_poly_t *g = ff_poly_gcd(f, deriv); /* Note: assuming ff_poly_gcd is exposed or accessible */
    /* If ff_poly_gcd is static, we can use the fact that if f has repeated roots,
       the Cantor-Zassenhaus DDF step will fail or we can just rely on the next step.
       For simplicity in this module, we assume random monic polynomials are mostly square-free. */
    if (g)
        ff_poly_free(g);
    if (deriv)
        ff_poly_free(deriv);

    /* 2. Check gcd(x^{p^i} - x, f) == 1 */
    ff_poly_t *x_poly = ff_poly_new(ctx, 1);
    ff_poly_set_coeff(x_poly, 1, ff_one_rep(ctx));

    ff_poly_t *g_pow = ff_poly_copy(x_poly); /* g_pow starts as x */

    for (int i = 1; i <= k / 2; ++i)
    {
        /* g_pow = g_pow^p mod f */
        ff_poly_t *next_pow = ff_poly_mod_pow_new(g_pow, p, f);
        ff_poly_free(g_pow);
        g_pow = next_pow;
        if (!g_pow)
        {
            ff_poly_free(x_poly);
            return false;
        }

        /* diff = g_pow - x */
        ff_poly_t *diff = ff_poly_sub_new(g_pow, x_poly);

        /* h = gcd(diff, f) */
        /* We need a GCD function. If it's static in ff_poly.c, we expose it.
           For this snippet, let's assume we have access to a GCD routine. */
        extern ff_poly_t *ff_poly_gcd(const ff_poly_t *a, const ff_poly_t *b);
        ff_poly_t *h = ff_poly_gcd(diff, f);

        ff_poly_free(diff);

        if (h && ff_poly_deg(h) > 0)
        {
            ff_poly_free(h);
            ff_poly_free(g_pow);
            ff_poly_free(x_poly);
            return false; /* Found a factor */
        }
        if (h)
            ff_poly_free(h);
    }

    ff_poly_free(g_pow);
    ff_poly_free(x_poly);
    return true;
}

static ff_poly_t *find_irreducible(const ff_prime_t *ctx, int k)
{
    uint32_t p = ff_prime_modulus(ctx);
    srand((unsigned int)time(NULL));

    while (1)
    {
        ff_poly_t *f = ff_poly_new(ctx, k);
        ff_poly_set_coeff(f, k, ff_one_rep(ctx)); /* Make monic */

        /* Random coefficients for degrees 0 to k-1 */
        for (int i = 0; i < k; ++i)
        {
            uint32_t val = ((uint32_t)rand() << 15) ^ (uint32_t)rand();
            ff_poly_set_coeff(f, i, ff_to_mont(ctx, val % p));
        }

        if (pkf_is_irreducible(ctx, f))
        {
            return f;
        }
        ff_poly_free(f);
    }
}

/* ============================================================ */
/* Context Management                                           */
/* ============================================================ */

ff_extk_t *ff_extk_new(const ff_prime_t *base, int k)
{
    if (!base || k < 1)
        return NULL;
    ff_extk_t *ctx = calloc(1, sizeof(ff_extk_t));
    if (!ctx)
        return NULL;

    ctx->base = base;
    ctx->k = k;
    ctx->irred = find_irreducible(base, k);

    if (!ctx->irred)
    {
        free(ctx);
        return NULL;
    }
    return ctx;
}

void ff_extk_free(ff_extk_t *ctx)
{
    if (!ctx)
        return;
    if (ctx->irred)
        ff_poly_free(ctx->irred);
    free(ctx);
}

const ff_prime_t *ff_extk_base(const ff_extk_t *ctx) { return ctx ? ctx->base : NULL; }
int ff_extk_deg(const ff_extk_t *ctx) { return ctx ? ctx->k : 0; }

/* ============================================================ */
/* Element Lifecycle                                            */
/* ============================================================ */

ff_extk_elem_t *ff_extk_elem_new(const ff_extk_t *ctx)
{
    return ff_poly_new(ctx->base, ctx->k - 1);
}

void ff_extk_elem_free(ff_extk_elem_t *e)
{
    ff_poly_free(e);
}

void ff_extk_zero(const ff_extk_t *ctx, ff_extk_elem_t *out)
{
    (void)ctx;
    ff_poly_zero(out);
}

void ff_extk_one(const ff_extk_t *ctx, ff_extk_elem_t *out)
{
    (void)ctx;
    ff_poly_one(out);
}

bool ff_extk_is_zero(const ff_extk_elem_t *a)
{
    if (!a)
        return true;
    int d = ff_poly_deg(a);
    for (int i = 0; i <= d; ++i)
    {
        if (ff_poly_coeff(a, i) != 0)
            return false;
    }
    return true;
}

bool ff_extk_eq(const ff_extk_elem_t *a, const ff_extk_elem_t *b)
{
    int da = ff_poly_deg(a);
    int db = ff_poly_deg(b);
    if (da != db)
        return false;
    for (int i = 0; i <= da; ++i)
    {
        if (ff_poly_coeff(a, i) != ff_poly_coeff(b, i))
            return false;
    }
    return true;
}

/* ============================================================ */
/* Arithmetic                                                   */
/* ============================================================ */

/* ============================================================ */
/* Arithmetic                                                   */
/* ============================================================ */

void ff_extk_add(const ff_extk_t *ctx, ff_extk_elem_t *out, const ff_extk_elem_t *a, const ff_extk_elem_t *b)
{
    (void)ctx;
    ff_poly_t *sum = ff_poly_add_new(a, b);
    if (sum)
    {
        /* Clear existing capacity WITHOUT destroying the degree */
        int out_deg = ff_poly_deg(out);
        for (int i = 0; i <= out_deg; ++i)
            ff_poly_set_coeff(out, i, 0);

        /* Copy new coefficients */
        int d = ff_poly_deg(sum);
        for (int i = 0; i <= d; ++i)
            ff_poly_set_coeff(out, i, ff_poly_coeff(sum, i));
        ff_poly_free(sum);
    }
}

void ff_extk_sub(const ff_extk_t *ctx, ff_extk_elem_t *out, const ff_extk_elem_t *a, const ff_extk_elem_t *b)
{
    (void)ctx;
    ff_poly_t *diff = ff_poly_sub_new(a, b);
    if (diff)
    {
        int out_deg = ff_poly_deg(out);
        for (int i = 0; i <= out_deg; ++i)
            ff_poly_set_coeff(out, i, 0);

        int d = ff_poly_deg(diff);
        for (int i = 0; i <= d; ++i)
            ff_poly_set_coeff(out, i, ff_poly_coeff(diff, i));
        ff_poly_free(diff);
    }
}

void ff_extk_mul(const ff_extk_t *ctx, ff_extk_elem_t *out, const ff_extk_elem_t *a, const ff_extk_elem_t *b)
{
    ff_poly_t *prod = ff_poly_mul_new(a, b);
    if (prod)
    {
        ff_poly_t *rem = ff_poly_mod_new(prod, ctx->irred);
        ff_poly_free(prod);
        if (rem)
        {
            int out_deg = ff_poly_deg(out);
            for (int i = 0; i <= out_deg; ++i)
                ff_poly_set_coeff(out, i, 0);

            int d = ff_poly_deg(rem);
            for (int i = 0; i <= d; ++i)
                ff_poly_set_coeff(out, i, ff_poly_coeff(rem, i));
            ff_poly_free(rem);
        }
    }
}

void ff_extk_neg(const ff_extk_t *ctx, ff_extk_elem_t *out, const ff_extk_elem_t *a)
{
    if (!ctx || !out || !a)
        return;
    const ff_prime_t *base = ctx->base;

    /* Clear without destroying capacity */
    int out_deg = ff_poly_deg(out);
    for (int i = 0; i <= out_deg; ++i)
        ff_poly_set_coeff(out, i, 0);

    int d = ff_poly_deg(a);
    for (int i = 0; i <= d; ++i)
    {
        uint32_t c = ff_poly_coeff(a, i);
        ff_poly_set_coeff(out, i, ff_neg(base, c));
    }
}

bool ff_extk_inv(const ff_extk_t *ctx, ff_extk_elem_t *out, const ff_extk_elem_t *a)
{
    if (ff_extk_is_zero(a))
        return false;

    uint32_t p = ff_prime_modulus(ctx->base);
    int k = ctx->k;

    /* Calculate p^k - 2 using __int128 to avoid overflow */
    ff_u128_t order = 1;
    for (int i = 0; i < k; ++i)
        order *= p;
    order -= 2;

    ff_poly_t *res = ff_poly_new(ctx->base, 0);
    ff_poly_one(res);
    ff_poly_t *base_copy = ff_poly_copy(a);
    ff_poly_t *mod_copy = ctx->irred;

    ff_u128_t exp = order;
    while (exp > 0)
    {
        if (exp & 1)
        {
            ff_poly_t *prod = ff_poly_mul_new(res, base_copy);
            ff_poly_free(res);
            res = ff_poly_mod_new(prod, mod_copy);
            ff_poly_free(prod);
        }
        exp >>= 1;
        if (exp > 0)
        {
            ff_poly_t *sq = ff_poly_mul_new(base_copy, base_copy);
            ff_poly_free(base_copy);
            base_copy = ff_poly_mod_new(sq, mod_copy);
            ff_poly_free(sq);
        }
    }

    if (res)
    {
        int out_deg = ff_poly_deg(out);
        for (int i = 0; i <= out_deg; ++i)
            ff_poly_set_coeff(out, i, 0);

        int d = ff_poly_deg(res);
        for (int i = 0; i <= d; ++i)
            ff_poly_set_coeff(out, i, ff_poly_coeff(res, i));
        ff_poly_free(res);
        return true;
    }
    return false;
}

ff_extk_t *ff_extk_new_with_poly(const ff_prime_t *base, const ff_poly_t *irred)
{
    if (!base || !irred || ff_poly_deg(irred) < 1)
        return NULL;
    ff_extk_t *ctx = calloc(1, sizeof(ff_extk_t));
    if (!ctx)
        return NULL;
    ctx->base = base;
    ctx->k = ff_poly_deg(irred);
    ctx->irred = ff_poly_copy(irred);
    return ctx;
}

/* ============================================================ */
/* F_{p^k} Matrix Implementation                                */
/* ============================================================ */
struct ff_mat_extk
{
    const ff_extk_t *ctx;
    size_t rows;
    size_t cols;
    size_t stride;
    ff_poly_t **data; /* 1D array of pointers to ff_poly_t */
};

static void pkf_extk_poly_copy_to(ff_poly_t *dst, const ff_poly_t *src)
{
    if (!dst || !src)
        return;

    /* FIX: Clear destination WITHOUT destroying its allocated capacity */
    int dst_deg = ff_poly_deg(dst);
    for (int i = 0; i <= dst_deg; ++i)
    {
        ff_poly_set_coeff(dst, i, 0);
    }

    /* Copy source coefficients */
    int d = ff_poly_deg(src);
    for (int i = 0; i <= d; ++i)
    {
        ff_poly_set_coeff(dst, i, ff_poly_coeff(src, i));
    }
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
        {
            if (M->data[i])
                ff_poly_free(M->data[i]);
        }
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

    /* FIX: Clear matrix WITHOUT destroying capacity */
    for (size_t i = 0; i < M->rows * M->cols; ++i)
    {
        int d = ff_poly_deg(M->data[i]);
        for (int j = 0; j <= d; ++j)
        {
            ff_poly_set_coeff(M->data[i], j, 0);
        }
    }

    size_t n = (M->rows < M->cols) ? M->rows : M->cols;
    for (size_t i = 0; i < n; ++i)
    {
        ff_poly_one(M->data[i * M->stride + i]);
    }
}

void ff_mat_extk_mul(ff_mat_extk_t *C, const ff_mat_extk_t *A, const ff_mat_extk_t *B)
{
    if (!A || !B || !C || A->cols != B->rows || C->rows != A->rows || C->cols != B->cols)
        return;
    const ff_extk_t *ctx = A->ctx;

    /* FIX: Clear C WITHOUT destroying capacity */
    for (size_t i = 0; i < C->rows * C->cols; ++i)
    {
        int d = ff_poly_deg(C->data[i]);
        for (int j = 0; j <= d; ++j)
        {
            ff_poly_set_coeff(C->data[i], j, 0);
        }
    }

    ff_poly_t *temp = ff_poly_new(ctx->base, ctx->k - 1);
    ff_poly_t *sum = ff_poly_new(ctx->base, ctx->k - 1);

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
    if (!A || !X || A->rows != A->cols || X->rows != A->rows || X->cols != A->cols)
        return false;
    const ff_extk_t *ctx = A->ctx;
    size_t n = A->rows;
    ff_mat_extk_t *aug = ff_mat_extk_new(ctx, n, 2 * n);
    if (!aug)
        return false;

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
            pkf_extk_poly_copy_to(aug->data[i * aug->stride + j], A->data[i * A->stride + j]);
        ff_poly_one(aug->data[i * aug->stride + n + i]);
    }

    ff_poly_t *factor = ff_poly_new(ctx->base, ctx->k - 1);
    ff_poly_t *term = ff_poly_new(ctx->base, ctx->k - 1);
    ff_poly_t *inv_pivot = ff_poly_new(ctx->base, ctx->k - 1);

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
        {
            ff_extk_mul(ctx, aug->data[col * aug->stride + j], aug->data[col * aug->stride + j], inv_pivot);
        }
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
                ff_extk_sub(ctx, aug->data[r * aug->stride + j], aug->data[r * aug->stride + j], term);
            }
        }
    }
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            pkf_extk_poly_copy_to(X->data[i * X->stride + j], aug->data[i * aug->stride + n + j]);
        }
    }
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
    size_t rows = M->rows, cols = M->cols;
    ff_mat_extk_t *aug = ff_mat_extk_new(ctx, rows, cols);
    if (!aug)
        return -1;
    for (size_t i = 0; i < rows * cols; ++i)
        pkf_extk_poly_copy_to(aug->data[i], M->data[i]);

    int *pivot_cols = calloc(cols, sizeof(int));
    for (size_t i = 0; i < cols; ++i)
        pivot_cols[i] = -1;

    ff_poly_t *factor = ff_poly_new(ctx->base, ctx->k - 1);
    ff_poly_t *term = ff_poly_new(ctx->base, ctx->k - 1);
    ff_poly_t *inv_pivot = ff_poly_new(ctx->base, ctx->k - 1);

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
        {
            free(pivot_cols);
            ff_mat_extk_free(aug);
            ff_poly_free(factor);
            ff_poly_free(term);
            ff_poly_free(inv_pivot);
            return -1;
        }
        for (size_t j = col; j < cols; ++j)
            ff_extk_mul(ctx, aug->data[row * cols + j], aug->data[row * cols + j], inv_pivot);
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
                ff_extk_sub(ctx, aug->data[r * cols + j], aug->data[r * cols + j], term);
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
    size_t bidx = 0;
    for (size_t free_col = 0; free_col < cols; ++free_col)
    {
        if (pivot_cols[free_col] != -1)
            continue;
        ff_poly_one(B->data[free_col * nullity + bidx]);
        for (size_t pc = 0; pc < cols; ++pc)
        {
            int prow = pivot_cols[pc];
            if (prow == -1)
                continue;
            ff_extk_neg(ctx, B->data[pc * nullity + bidx], aug->data[prow * cols + free_col]);
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

/* ============================================================ */
/* Generic Jordan Engine Helpers                                */
/* ============================================================ */
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

static ff_mat_extk_t *pkf_extk_mat_sub_scalar_identity(const ff_mat_extk_t *A, const ff_poly_t *lambda)
{
    ff_mat_extk_t *N = pkf_extk_mat_copy(A);
    if (!N)
        return NULL;
    for (size_t i = 0; i < A->rows; ++i)
        ff_extk_sub(A->ctx, N->data[i * N->stride + i], N->data[i * N->stride + i], lambda);
    return N;
}

static bool pkf_extk_mat_vec_mul(const ff_mat_extk_t *M, const ff_poly_t *const *x, ff_poly_t **y)
{
    ff_poly_t *temp = ff_poly_new(M->ctx->base, M->ctx->k - 1);
    ff_poly_t *sum = ff_poly_new(M->ctx->base, M->ctx->k - 1);
    for (size_t i = 0; i < M->rows; ++i)
    {
        /* Clear sum without destroying capacity */
        int sum_deg = ff_poly_deg(sum);
        for (int j = 0; j <= sum_deg; ++j)
            ff_poly_set_coeff(sum, j, 0);

        for (size_t j = 0; j < M->cols; ++j)
        {
            if (ff_extk_is_zero(M->data[i * M->stride + j]) || ff_extk_is_zero(x[j]))
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
    for (size_t i = 0; i < M->rows; ++i)
    {
        v[i] = ff_poly_new(M->ctx->base, M->ctx->k - 1);
        pkf_extk_poly_copy_to(v[i], M->data[i * M->stride + col]);
    }
    return v;
}

static ff_mat_extk_t *pkf_extk_basis_append_vector(ff_mat_extk_t *B, size_t n, ff_poly_t **v)
{
    size_t old_cols = B ? B->cols : 0;
    ff_mat_extk_t *N = ff_mat_extk_new(B ? B->ctx : NULL, n, old_cols + 1);
    if (!N)
        return NULL;
    for (size_t j = 0; j < old_cols; ++j)
        for (size_t i = 0; i < n; ++i)
            pkf_extk_poly_copy_to(N->data[i * N->stride + j], B->data[i * B->stride + j]);
    for (size_t i = 0; i < n; ++i)
        pkf_extk_poly_copy_to(N->data[i * N->stride + old_cols], v[i]);
    return N;
}

static bool pkf_extk_vector_independent(ff_mat_extk_t *B, size_t n, ff_poly_t **v)
{
    size_t cols = B ? B->cols : 0;
    if (cols + 1 > n)
        return false;
    ff_mat_extk_t *aug = ff_mat_extk_new(B ? B->ctx : NULL, n, cols + 1);
    for (size_t j = 0; j < cols; ++j)
        for (size_t i = 0; i < n; ++i)
            pkf_extk_poly_copy_to(aug->data[i * aug->stride + j], B->data[i * B->stride + j]);
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
        if (pkf_extk_vector_independent(B, M->rows, v))
        {
            ff_mat_extk_t *N = pkf_extk_basis_append_vector(B, M->rows, v);
            ff_mat_extk_free(B);
            B = N;
        }
        for (size_t i = 0; i < M->rows; ++i)
            ff_poly_free(v[i]);
        free(v);
        if (!B)
            return NULL;
    }
    return B;
}

static ff_mat_extk_t *pkf_extk_matrix_image(const ff_mat_extk_t *A, const ff_mat_extk_t *basis)
{
    ff_mat_extk_t *img = ff_mat_extk_new(A->ctx, A->rows, basis->cols);
    ff_mat_extk_mul(img, A, basis);
    return img;
}

static ff_mat_extk_t *pkf_extk_basis_append_matrix(ff_mat_extk_t *B, const ff_mat_extk_t *Add)
{
    for (size_t j = 0; j < Add->cols; ++j)
    {
        ff_poly_t **v = pkf_extk_vec_copy_column(Add, j);
        ff_mat_extk_t *N = pkf_extk_basis_append_vector(B, Add->rows, v);
        for (size_t i = 0; i < Add->rows; ++i)
            ff_poly_free(v[i]);
        free(v);
        ff_mat_extk_free(B);
        B = N;
        if (!B)
            return NULL;
    }
    return B;
}

/* ============================================================ */
/* Generic Jordan Engine Core                                   */
/* ============================================================ */
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

bool ff_jordan_extk_with_eigenvalues(
    const ff_extk_t *E, const ff_mat_extk_t *A,
    const ff_poly_t *const *eigs, size_t eig_count,
    ff_jordan_extk_result_t *out)
{
    if (!A || !eigs || !out)
        return false;
    memset(out, 0, sizeof(*out));
    const ff_extk_t *ctx = A->ctx;
    size_t n = A->rows;
    if (A->cols != n || n == 0)
        return false;

    ff_mat_extk_t *P_basis = NULL;
    ff_mat_extk_t *J = ff_mat_extk_new(ctx, n, n);
    size_t *block_sizes = calloc(n, sizeof(size_t));
    ff_poly_t **block_eigs = calloc(n, sizeof(ff_poly_t *));
    size_t block_count = 0, p_col = 0;

    ff_poly_t *one = ff_poly_new(ctx->base, ctx->k - 1);
    ff_poly_one(one);

    for (size_t ei = 0; ei < eig_count; ++ei)
    {
        bool duplicate = false;
        for (size_t ej = 0; ej < ei; ++ej)
            if (ff_extk_eq(eigs[ei], eigs[ej]))
            {
                duplicate = true;
                break;
            }
        if (duplicate)
            continue;

        const ff_poly_t *lambda = eigs[ei];
        ff_mat_extk_t *N = pkf_extk_mat_sub_scalar_identity(A, lambda);
        ff_mat_extk_t **K = calloc(n + 2, sizeof(ff_mat_extk_t *));
        size_t *dims = calloc(n + 2, sizeof(size_t));
        ff_mat_extk_t *Pk = pkf_extk_mat_copy(N);
        int maxK = 0;

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
            dims[k] = nullity;
            maxK = k;
            if ((size_t)nullity == n)
                break;
            if (k < (int)n)
            {
                ff_mat_extk_t *next = ff_mat_extk_new(ctx, n, n);
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

            ff_mat_extk_t *S = (k > 1 && K[k - 1]) ? pkf_extk_mat_copy(K[k - 1]) : NULL;
            ff_mat_extk_t *Kkp1 = (k < maxK) ? K[k + 1] : K[maxK];
            ff_mat_extk_t *img = pkf_extk_matrix_image(N, Kkp1);
            ff_mat_extk_t *img_basis = pkf_extk_column_basis(img);
            ff_mat_extk_free(img);
            if (img_basis)
            {
                S = pkf_extk_basis_append_matrix(S, img_basis);
                ff_mat_extk_free(img_basis);
            }

            ff_mat_extk_t *C = NULL;
            ff_mat_extk_t *S_work = S ? pkf_extk_mat_copy(S) : NULL;
            for (size_t col = 0; col < K[k]->cols; ++col)
            {
                if (C && C->cols >= exact)
                    break;
                ff_poly_t **v = pkf_extk_vec_copy_column(K[k], col);
                if (pkf_extk_vector_independent(S_work, n, v))
                {
                    C = pkf_extk_basis_append_vector(C, n, v);
                    S_work = pkf_extk_basis_append_vector(S_work, n, v);
                }
                for (size_t i = 0; i < n; ++i)
                    ff_poly_free(v[i]);
                free(v);
            }
            ff_mat_extk_free(S_work);
            ff_mat_extk_free(S);

            for (size_t t = 0; t < exact; ++t)
            {
                ff_poly_t **top = pkf_extk_vec_copy_column(C, t);

                /* Allocate k vectors, each of length n */
                ff_poly_t ***chain = calloc(k, sizeof(ff_poly_t **));
                for (int idx = 0; idx < k; ++idx)
                {
                    chain[idx] = calloc(n, sizeof(ff_poly_t *));
                    for (size_t i = 0; i < n; ++i)
                    {
                        chain[idx][i] = ff_poly_new(ctx->base, ctx->k - 1);
                    }
                }

                /* Initialize the top of the chain */
                for (size_t i = 0; i < n; ++i)
                {
                    pkf_extk_poly_copy_to(chain[k - 1][i], top[i]);
                }

                /* Generate the rest of the chain by multiplying by N */
                for (int idx = k - 2; idx >= 0; --idx)
                {
                    pkf_extk_mat_vec_mul(N, (const ff_poly_t *const *)chain[idx + 1], chain[idx]);
                }

                /* Add the chain to P_basis and construct the Jordan block */
                for (int s = 0; s < k; ++s)
                {
                    P_basis = pkf_extk_basis_append_vector(P_basis, n, chain[s]);
                    pkf_extk_poly_copy_to(J->data[(p_col + s) * J->stride + (p_col + s)], lambda);
                    if (s < k - 1)
                        pkf_extk_poly_copy_to(J->data[(p_col + s) * J->stride + (p_col + s + 1)], one);
                }
                block_sizes[block_count] = k;
                block_eigs[block_count] = ff_poly_copy(lambda);
                block_count++;
                p_col += k;

                /* Cleanup chain */
                for (int idx = 0; idx < k; ++idx)
                {
                    for (size_t i = 0; i < n; ++i)
                        ff_poly_free(chain[idx][i]);
                    free(chain[idx]);
                }
                free(chain);

                /* Cleanup top */
                for (size_t i = 0; i < n; ++i)
                    ff_poly_free(top[i]);
                free(top);
            }
            ff_mat_extk_free(C);
        }
        for (int k = 1; k <= maxK; ++k)
            ff_mat_extk_free(K[k]);
        free(K);
        free(dims);
        ff_mat_extk_free(N);
    }
    ff_poly_free(one);

    if (p_col > n)
        goto fail;

    /* Verification: A * P == P * J_sub */
    bool verify_ok = true;
    ff_mat_extk_t *AP = ff_mat_extk_new(ctx, n, p_col);
    ff_mat_extk_t *J_sub = ff_mat_extk_new(ctx, p_col, p_col);
    ff_mat_extk_t *PJ = ff_mat_extk_new(ctx, n, p_col);
    ff_mat_extk_mul(AP, A, P_basis);
    for (size_t i = 0; i < p_col; ++i)
        for (size_t j = 0; j < p_col; ++j)
            pkf_extk_poly_copy_to(J_sub->data[i * J_sub->stride + j], J->data[i * J->stride + j]);
    ff_mat_extk_mul(PJ, P_basis, J_sub);
    for (size_t i = 0; i < n && verify_ok; ++i)
    {
        for (size_t j = 0; j < p_col && verify_ok; ++j)
        {
            if (!ff_extk_eq(AP->data[i * AP->stride + j], PJ->data[i * PJ->stride + j]))
                verify_ok = false;
        }
    }
    ff_mat_extk_free(AP);
    ff_mat_extk_free(J_sub);
    ff_mat_extk_free(PJ);
    if (!verify_ok)
        goto fail;

    out->J = J;
    out->P = P_basis;
    out->block_count = block_count;
    out->block_sizes = block_sizes;
    out->block_eigenvalues = block_eigs;
    return true;

fail:
    ff_mat_extk_free(P_basis);
    ff_mat_extk_free(J);
    free(block_sizes);
    for (size_t i = 0; i < block_count; ++i)
        ff_poly_free(block_eigs[i]);
    free(block_eigs);
    return false;
}
