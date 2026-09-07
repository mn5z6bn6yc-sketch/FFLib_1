
/* ============================================================ */
/* Implementation                                               */
/* ============================================================ */
#include "fflib.h"
#include <stdlib.h>
#include <string.h>

#if !defined(__GNUC__) && !defined(__clang__)
#error "Implementation requires unsigned __int128."
#endif

/* ============================================================ */
/* Small helpers                                                */
/* ============================================================ */

static bool size_mul_ok(size_t a, size_t b)
{
    return a == 0 || b <= SIZE_MAX / a;
}

static uint32_t u32_pow_mod(uint32_t a, uint32_t e, uint32_t p)
{
    if (p == 1)
        return 0;

    uint64_t r = 1 % p;
    uint64_t b = a % p;

    while (e)
    {
        if (e & 1u)
        {
            r = (r * b) % p;
        }

        e >>= 1;

        if (e)
        {
            b = (b * b) % p;
        }
    }

    return (uint32_t)r;
}

/* ============================================================ */
/* Prime field implementation                                   */
/* ============================================================ */

struct ff_prime
{
    uint32_t p;
    ff_backend_t backend;

    /* Montgomery constants */
    uint32_t r2;   /* R^2 mod p, R = 2^32 */
    uint32_t ninv; /* -p^{-1} mod 2^32 */
};

static uint32_t inv_mod_2_32(uint32_t p)
{
    uint32_t x = 1;

    /* Newton iteration for inverse modulo 2^32 */
    for (int i = 0; i < 5; ++i)
    {
        x *= 2u - p * x;
    }

    return x;
}

ff_prime_t *ff_prime_new(uint32_t p, ff_backend_t backend)
{
    if (p < 2)
        return NULL;

    if (backend == FF_BACKEND_MONTGOMERY && (p & 1u) == 0)
    {
        return NULL;
    }

    ff_prime_t *ctx = calloc(1, sizeof *ctx);
    if (!ctx)
        return NULL;

    ctx->p = p;
    ctx->backend = backend;

    if (backend == FF_BACKEND_MONTGOMERY)
    {
        uint64_t r = (1ULL << 32) % p;
        ctx->r2 = (uint32_t)((r * r) % p);
        ctx->ninv = 0u - inv_mod_2_32(p);
    }

    return ctx;
}

void ff_prime_free(ff_prime_t *ctx)
{
    free(ctx);
}

uint32_t ff_prime_modulus(const ff_prime_t *ctx)
{
    return ctx ? ctx->p : 0;
}

ff_backend_t ff_prime_backend(const ff_prime_t *ctx)
{
    return ctx ? ctx->backend : FF_BACKEND_NATIVE;
}

uint32_t ff_add(const ff_prime_t *ctx, uint32_t a, uint32_t b)
{
    uint64_t s = (uint64_t)a + (uint64_t)b;
    if (s >= ctx->p)
        s -= ctx->p;
    return (uint32_t)s;
}

uint32_t ff_sub(const ff_prime_t *ctx, uint32_t a, uint32_t b)
{
    if (a >= b)
        return a - b;
    return (uint32_t)((uint64_t)a + ctx->p - b);
}

uint32_t ff_neg(const ff_prime_t *ctx, uint32_t a)
{
    if (a == 0)
        return 0;
    return ctx->p - a;
}

static uint32_t mont_reduce(const ff_prime_t *ctx, uint64_t T)
{
    uint32_t m = (uint32_t)(T * ctx->ninv);

    ff_u128_t sum = (ff_u128_t)T + (ff_u128_t)m * ctx->p;
    uint64_t u = (uint64_t)(sum >> 32);

    if (u >= ctx->p)
        u -= ctx->p;

    return (uint32_t)u;
}

uint32_t ff_mul(const ff_prime_t *ctx, uint32_t a, uint32_t b)
{
    if (ctx->backend == FF_BACKEND_MONTGOMERY)
    {
        return mont_reduce(ctx, (uint64_t)a * b);
    }

    return (uint32_t)((uint64_t)a * b % ctx->p);
}

uint32_t ff_to_mont(const ff_prime_t *ctx, uint32_t x)
{
    x %= ctx->p;

    if (ctx->backend != FF_BACKEND_MONTGOMERY)
    {
        return x;
    }

    return mont_reduce(ctx, (uint64_t)x * ctx->r2);
}

uint32_t ff_from_mont(const ff_prime_t *ctx, uint32_t x_mont)
{
    if (ctx->backend != FF_BACKEND_MONTGOMERY)
    {
        return x_mont % ctx->p;
    }

    return mont_reduce(ctx, (uint64_t)x_mont);
}

uint32_t ff_zero_rep(const ff_prime_t *ctx)
{
    (void)ctx;
    return 0;
}

uint32_t ff_one_rep(const ff_prime_t *ctx)
{
    return ff_to_mont(ctx, 1);
}

uint32_t ff_minus_one_rep(const ff_prime_t *ctx)
{
    return ff_to_mont(ctx, ctx->p - 1);
}

uint32_t ff_inv(const ff_prime_t *ctx, uint32_t a)
{
    if (a == 0)
        return 0;

    uint32_t regular = ff_from_mont(ctx, a);
    if (regular == 0)
        return 0;

    uint32_t inv_regular = u32_pow_mod(regular, ctx->p - 2, ctx->p);

    return ff_to_mont(ctx, inv_regular);
}

uint32_t ff_pow(const ff_prime_t *ctx, uint32_t a, uint32_t e)
{
    uint32_t r = ff_one_rep(ctx);
    uint32_t base = a;

    while (e)
    {
        if (e & 1u)
        {
            r = ff_mul(ctx, r, base);
        }

        e >>= 1;

        if (e)
        {
            base = ff_mul(ctx, base, base);
        }
    }

    return r;
}

/* ============================================================ */
/* Square roots over prime fields                               */
/* ============================================================ */

static int legendre_regular(const ff_prime_t *ctx, uint32_t a)
{
    a %= ctx->p;

    if (a == 0)
        return 0;
    if (ctx->p == 2)
        return (int)a;

    uint32_t x = u32_pow_mod(a, (ctx->p - 1) / 2, ctx->p);

    if (x == 1)
        return 1;
    if (x == ctx->p - 1)
        return -1;

    return 0;
}

static bool sqrt_regular(const ff_prime_t *ctx, uint32_t a, uint32_t *root)
{
    a %= ctx->p;

    if (a == 0)
    {
        *root = 0;
        return true;
    }

    if (ctx->p == 2)
    {
        *root = a;
        return true;
    }

    int leg = legendre_regular(ctx, a);
    if (leg != 1)
        return false;

    uint32_t p = ctx->p;

    if ((p & 3u) == 3u)
    {
        *root = u32_pow_mod(a, (p + 1) / 4, p);
        return true;
    }

    /* Tonelli-Shanks */

    uint32_t q = p - 1;
    int s = 0;

    while ((q & 1u) == 0)
    {
        q >>= 1;
        s++;
    }

    uint32_t z = 2;
    while (z < p && legendre_regular(ctx, z) != -1)
    {
        z++;
    }

    if (z == p)
        return false;

    uint32_t m = (uint32_t)s;
    uint32_t c = u32_pow_mod(z, q, p);
    uint32_t t = u32_pow_mod(a, q, p);
    uint32_t r = u32_pow_mod(a, (q + 1) / 2, p);

    while (t != 1)
    {
        uint32_t tt = t;
        int i = 0;

        for (i = 1; i < (int)m; ++i)
        {
            tt = (uint32_t)((uint64_t)tt * tt % p);
            if (tt == 1)
                break;
        }

        if (i == (int)m)
            return false;

        uint32_t exp = 1u << (m - (uint32_t)i - 1u);
        uint32_t b = u32_pow_mod(c, exp, p);

        m = (uint32_t)i;
        c = (uint32_t)((uint64_t)b * b % p);
        t = (uint32_t)((uint64_t)t * c % p);
        r = (uint32_t)((uint64_t)r * b % p);
    }

    *root = r;
    return true;
}

bool ff_sqrt(const ff_prime_t *ctx, uint32_t a_rep, uint32_t *root_rep)
{
    if (!ctx || !root_rep)
        return false;

    uint32_t a_regular = ff_from_mont(ctx, a_rep);
    uint32_t root_regular = 0;

    if (!sqrt_regular(ctx, a_regular, &root_regular))
    {
        return false;
    }

    *root_rep = ff_to_mont(ctx, root_regular);
    return true;
}

/* ============================================================ */
/* Quadratic extension field implementation                     */
/* ============================================================ */

struct ff_ext2
{
    const ff_prime_t *base;
    uint32_t omega;     /* regular integer omega */
    uint32_t omega_rep; /* omega in active base representation */
};

ff_ext2_t *ff_ext2_new(const ff_prime_t *base, uint32_t omega)
{
    if (!base)
        return NULL;

    uint32_t p = ff_prime_modulus(base);

    if (p < 3)
        return NULL;

    omega %= p;

    if (omega == 0)
        return NULL;

    /*
        Require x^2 - omega irreducible over F_p.
        For odd p, omega must be a quadratic non-residue.
    */
    if (legendre_regular(base, omega) != -1)
    {
        return NULL;
    }

    ff_ext2_t *ctx = calloc(1, sizeof *ctx);
    if (!ctx)
        return NULL;

    ctx->base = base;
    ctx->omega = omega;
    ctx->omega_rep = ff_to_mont(base, omega);

    return ctx;
}

void ff_ext2_free(ff_ext2_t *ctx)
{
    free(ctx);
}

const ff_prime_t *ff_ext2_base(const ff_ext2_t *ctx)
{
    return ctx ? ctx->base : NULL;
}

void ff_ext2_zero(const ff_ext2_t *ctx, ff_ext2_elem_t *out)
{
    (void)ctx;
    out->a = 0;
    out->b = 0;
}

void ff_ext2_one(const ff_ext2_t *ctx, ff_ext2_elem_t *out)
{
    out->a = ff_one_rep(ctx->base);
    out->b = 0;
}

void ff_ext2_add(
    const ff_ext2_t *ctx,
    ff_ext2_elem_t *out,
    const ff_ext2_elem_t *a,
    const ff_ext2_elem_t *b)
{
    out->a = ff_add(ctx->base, a->a, b->a);
    out->b = ff_add(ctx->base, a->b, b->b);
}

void ff_ext2_sub(
    const ff_ext2_t *ctx,
    ff_ext2_elem_t *out,
    const ff_ext2_elem_t *a,
    const ff_ext2_elem_t *b)
{
    out->a = ff_sub(ctx->base, a->a, b->a);
    out->b = ff_sub(ctx->base, a->b, b->b);
}

void ff_ext2_neg(
    const ff_ext2_t *ctx,
    ff_ext2_elem_t *out,
    const ff_ext2_elem_t *a)
{
    out->a = ff_neg(ctx->base, a->a);
    out->b = ff_neg(ctx->base, a->b);
}

void ff_ext2_mul(
    const ff_ext2_t *ctx,
    ff_ext2_elem_t *out,
    const ff_ext2_elem_t *a,
    const ff_ext2_elem_t *b)
{
    const ff_prime_t *F = ctx->base;

    uint32_t ac = ff_mul(F, a->a, b->a);
    uint32_t bd = ff_mul(F, a->b, b->b);

    uint32_t ad = ff_mul(F, a->a, b->b);
    uint32_t bc = ff_mul(F, a->b, b->a);

    uint32_t real = ff_add(F, ac, ff_mul(F, bd, ctx->omega_rep));
    uint32_t imag = ff_add(F, ad, bc);

    out->a = real;
    out->b = imag;
}

bool ff_ext2_inv(
    const ff_ext2_t *ctx,
    ff_ext2_elem_t *out,
    const ff_ext2_elem_t *a)
{
    if (ff_ext2_is_zero(a))
        return false;

    const ff_prime_t *F = ctx->base;

    /*
        For w^2 = omega:

          (a + b w)^{-1} = (a - b w) / (a^2 - omega b^2)
    */

    uint32_t a2 = ff_mul(F, a->a, a->a);
    uint32_t b2 = ff_mul(F, a->b, a->b);
    uint32_t omega_b2 = ff_mul(F, b2, ctx->omega_rep);

    uint32_t norm = ff_sub(F, a2, omega_b2);

    if (norm == 0)
        return false;

    uint32_t inv_norm = ff_inv(F, norm);

    out->a = ff_mul(F, a->a, inv_norm);
    out->b = ff_neg(F, ff_mul(F, a->b, inv_norm));

    return true;
}

bool ff_ext2_is_zero(const ff_ext2_elem_t *a)
{
    return a->a == 0 && a->b == 0;
}

bool ff_ext2_eq(const ff_ext2_elem_t *a, const ff_ext2_elem_t *b)
{
    return a->a == b->a && a->b == b->b;
}

/* ============================================================ */
/* Base-field matrices                                          */
/* ============================================================ */

ff_mat_t *ff_mat_new(const ff_prime_t *ctx, size_t rows, size_t cols)
{
    if (!ctx || rows == 0 || cols == 0)
        return NULL;

    if (!size_mul_ok(rows, cols))
        return NULL;

    ff_mat_t *M = calloc(1, sizeof *M);
    if (!M)
        return NULL;

    M->ctx = ctx;
    M->rows = rows;
    M->cols = cols;
    M->stride = cols;

    M->data = calloc(rows * cols, sizeof(uint32_t));
    if (!M->data)
    {
        free(M);
        return NULL;
    }

    return M;
}

void ff_mat_free(ff_mat_t *M)
{
    if (!M)
        return;
    free(M->data);
    free(M);
}

uint32_t *ff_mat_at(ff_mat_t *M, size_t r, size_t c)
{
    if (!M || r >= M->rows || c >= M->cols)
        return NULL;
    return &M->data[r * M->stride + c];
}

const uint32_t *ff_mat_atc(const ff_mat_t *M, size_t r, size_t c)
{
    if (!M || r >= M->rows || c >= M->cols)
        return NULL;
    return &M->data[r * M->stride + c];
}

void ff_mat_identity(ff_mat_t *M)
{
    if (!M)
        return;

    memset(M->data, 0, M->rows * M->stride * sizeof(uint32_t));

    uint32_t one = ff_one_rep(M->ctx);

    size_t n = (M->rows < M->cols) ? M->rows : M->cols;

    for (size_t i = 0; i < n; ++i)
    {
        M->data[i * M->stride + i] = one;
    }
}

void ff_mat_mul(
    ff_mat_t *C,
    const ff_mat_t *A,
    const ff_mat_t *B)
{
    if (!A || !B || !C)
        return;
    if (A->cols != B->rows)
        return;
    if (C->rows != A->rows || C->cols != B->cols)
        return;

    const ff_prime_t *ctx = A->ctx;

    memset(C->data, 0, C->rows * C->stride * sizeof(uint32_t));

    for (size_t i = 0; i < A->rows; ++i)
    {
        for (size_t k = 0; k < A->cols; ++k)
        {
            uint32_t av = A->data[i * A->stride + k];
            if (av == 0)
                continue;

            for (size_t j = 0; j < B->cols; ++j)
            {
                uint32_t bv = B->data[k * B->stride + j];
                if (bv == 0)
                    continue;

                uint32_t prod = ff_mul(ctx, av, bv);
                uint32_t *cv = &C->data[i * C->stride + j];

                *cv = ff_add(ctx, *cv, prod);
            }
        }
    }
}

bool ff_mat_inv(ff_mat_t *X, const ff_mat_t *A)
{
    if (!A || !X)
        return false;
    if (A->rows != A->cols)
        return false;
    if (X->rows != A->rows || X->cols != A->cols)
        return false;

    const ff_prime_t *ctx = A->ctx;
    size_t n = A->rows;
    size_t n2 = 2 * n;

    uint32_t *aug = calloc(n * n2, sizeof(uint32_t));
    if (!aug)
        return false;

    uint32_t one = ff_one_rep(ctx);

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            aug[i * n2 + j] = A->data[i * A->stride + j];
        }
        aug[i * n2 + n + i] = one;
    }

    for (size_t col = 0; col < n; ++col)
    {
        size_t pivot = col;

        while (pivot < n && aug[pivot * n2 + col] == 0)
        {
            pivot++;
        }

        if (pivot == n)
        {
            free(aug);
            return false;
        }

        if (pivot != col)
        {
            for (size_t j = 0; j < n2; ++j)
            {
                uint32_t t = aug[col * n2 + j];
                aug[col * n2 + j] = aug[pivot * n2 + j];
                aug[pivot * n2 + j] = t;
            }
        }

        uint32_t inv_pivot = ff_inv(ctx, aug[col * n2 + col]);
        if (inv_pivot == 0)
        {
            free(aug);
            return false;
        }

        for (size_t j = col; j < n2; ++j)
        {
            aug[col * n2 + j] = ff_mul(ctx, aug[col * n2 + j], inv_pivot);
        }

        for (size_t r = 0; r < n; ++r)
        {
            if (r == col)
                continue;

            uint32_t factor = aug[r * n2 + col];
            if (factor == 0)
                continue;

            for (size_t j = col; j < n2; ++j)
            {
                uint32_t term = ff_mul(ctx, factor, aug[col * n2 + j]);
                aug[r * n2 + j] = ff_sub(ctx, aug[r * n2 + j], term);
            }
        }
    }

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            X->data[i * X->stride + j] = aug[i * n2 + n + j];
        }
    }

    free(aug);
    return true;
}

int ff_mat_nullspace(ff_mat_t **basis_out, const ff_mat_t *M)
{
    if (!basis_out || !M)
        return -1;

    *basis_out = NULL;

    const ff_prime_t *ctx = M->ctx;
    size_t rows = M->rows;
    size_t cols = M->cols;

    uint32_t *aug = calloc(rows * cols, sizeof(uint32_t));
    int *pivot_cols = calloc(cols, sizeof(int));

    if (!aug || !pivot_cols)
    {
        free(aug);
        free(pivot_cols);
        return -1;
    }

    for (size_t i = 0; i < cols; ++i)
        pivot_cols[i] = -1;

    memcpy(aug, M->data, rows * cols * sizeof(uint32_t));

    size_t row = 0;

    for (size_t col = 0; col < cols && row < rows; ++col)
    {
        size_t pivot = row;

        while (pivot < rows && aug[pivot * cols + col] == 0)
        {
            pivot++;
        }

        if (pivot == rows)
            continue;

        if (pivot != row)
        {
            for (size_t j = 0; j < cols; ++j)
            {
                uint32_t t = aug[row * cols + j];
                aug[row * cols + j] = aug[pivot * cols + j];
                aug[pivot * cols + j] = t;
            }
        }

        uint32_t inv_pivot = ff_inv(ctx, aug[row * cols + col]);
        if (inv_pivot == 0)
        {
            free(aug);
            free(pivot_cols);
            return -1;
        }

        for (size_t j = col; j < cols; ++j)
        {
            aug[row * cols + j] = ff_mul(ctx, aug[row * cols + j], inv_pivot);
        }

        for (size_t r = 0; r < rows; ++r)
        {
            if (r == row)
                continue;

            uint32_t factor = aug[r * cols + col];
            if (factor == 0)
                continue;

            for (size_t j = col; j < cols; ++j)
            {
                uint32_t term = ff_mul(ctx, factor, aug[row * cols + j]);
                aug[r * cols + j] = ff_sub(ctx, aug[r * cols + j], term);
            }
        }

        pivot_cols[col] = (int)row;
        row++;
    }

    size_t rank = row;
    size_t nullity = cols - rank;

    if (nullity == 0)
    {
        free(aug);
        free(pivot_cols);
        return 0;
    }

    ff_mat_t *B = ff_mat_new(ctx, cols, nullity);
    if (!B)
    {
        free(aug);
        free(pivot_cols);
        return -1;
    }

    uint32_t one = ff_one_rep(ctx);
    size_t bidx = 0;

    for (size_t free_col = 0; free_col < cols; ++free_col)
    {
        if (pivot_cols[free_col] != -1)
            continue;

        B->data[free_col * nullity + bidx] = one;

        for (size_t pc = 0; pc < cols; ++pc)
        {
            int prow = pivot_cols[pc];
            if (prow == -1)
                continue;

            uint32_t val = aug[(size_t)prow * cols + free_col];
            B->data[pc * nullity + bidx] = ff_sub(ctx, 0, val);
        }

        bidx++;
    }

    free(aug);
    free(pivot_cols);

    *basis_out = B;
    return (int)nullity;
}

bool ff_linear_solve_unique_base(
    const ff_mat_t *M,
    const uint32_t *rhs,
    uint32_t *x)
{
    if (!M || !rhs || !x)
        return false;

    const ff_prime_t *ctx = M->ctx;
    size_t rows = M->rows;
    size_t cols = M->cols;
    size_t stride = cols + 1;

    uint32_t *aug = calloc(rows * stride, sizeof(uint32_t));
    int *pivot_row = calloc(cols, sizeof(int));

    if (!aug || !pivot_row)
    {
        free(aug);
        free(pivot_row);
        return false;
    }

    for (size_t i = 0; i < cols; ++i)
        pivot_row[i] = -1;

    for (size_t i = 0; i < rows; ++i)
    {
        for (size_t j = 0; j < cols; ++j)
        {
            aug[i * stride + j] = M->data[i * M->stride + j];
        }
        aug[i * stride + cols] = rhs[i];
    }

    size_t row = 0;

    for (size_t col = 0; col < cols && row < rows; ++col)
    {
        size_t pivot = row;

        while (pivot < rows && aug[pivot * stride + col] == 0)
        {
            pivot++;
        }

        if (pivot == rows)
            continue;

        if (pivot != row)
        {
            for (size_t j = 0; j < stride; ++j)
            {
                uint32_t t = aug[row * stride + j];
                aug[row * stride + j] = aug[pivot * stride + j];
                aug[pivot * stride + j] = t;
            }
        }

        uint32_t inv_pivot = ff_inv(ctx, aug[row * stride + col]);
        if (inv_pivot == 0)
        {
            free(aug);
            free(pivot_row);
            return false;
        }

        for (size_t j = col; j < stride; ++j)
        {
            aug[row * stride + j] = ff_mul(ctx, aug[row * stride + j], inv_pivot);
        }

        for (size_t r = 0; r < rows; ++r)
        {
            if (r == row)
                continue;

            uint32_t factor = aug[r * stride + col];
            if (factor == 0)
                continue;

            for (size_t j = col; j < stride; ++j)
            {
                uint32_t term = ff_mul(ctx, factor, aug[row * stride + j]);
                aug[r * stride + j] = ff_sub(ctx, aug[r * stride + j], term);
            }
        }

        pivot_row[col] = (int)row;
        row++;
    }

    /* Check consistency. */
    for (size_t r = row; r < rows; ++r)
    {
        if (aug[r * stride + cols] != 0)
        {
            free(aug);
            free(pivot_row);
            return false;
        }
    }

    /* Require uniqueness. */
    if (row < cols)
    {
        free(aug);
        free(pivot_row);
        return false;
    }

    for (size_t col = 0; col < cols; ++col)
    {
        int prow = pivot_row[col];

        if (prow == -1)
        {
            free(aug);
            free(pivot_row);
            return false;
        }

        x[col] = aug[(size_t)prow * stride + cols];
    }

    free(aug);
    free(pivot_row);
    return true;
}

/* ============================================================ */
/* Extension-field matrices                                     */
/* ============================================================ */

ff_mat_ext2_t *ff_mat_ext2_new(const ff_ext2_t *ctx, size_t rows, size_t cols)
{
    if (!ctx || rows == 0 || cols == 0)
        return NULL;

    if (!size_mul_ok(rows, cols))
        return NULL;

    ff_mat_ext2_t *M = calloc(1, sizeof *M);
    if (!M)
        return NULL;

    M->ctx = ctx;
    M->rows = rows;
    M->cols = cols;
    M->stride = cols;

    M->data = calloc(rows * cols, sizeof(ff_ext2_elem_t));
    if (!M->data)
    {
        free(M);
        return NULL;
    }

    return M;
}

void ff_mat_ext2_free(ff_mat_ext2_t *M)
{
    if (!M)
        return;
    free(M->data);
    free(M);
}

ff_ext2_elem_t *ff_mat_ext2_at(ff_mat_ext2_t *M, size_t r, size_t c)
{
    if (!M || r >= M->rows || c >= M->cols)
        return NULL;
    return &M->data[r * M->stride + c];
}

const ff_ext2_elem_t *ff_mat_ext2_atc(const ff_mat_ext2_t *M, size_t r, size_t c)
{
    if (!M || r >= M->rows || c >= M->cols)
        return NULL;
    return &M->data[r * M->stride + c];
}

void ff_mat_ext2_identity(ff_mat_ext2_t *M)
{
    if (!M)
        return;

    memset(M->data, 0, M->rows * M->stride * sizeof(ff_ext2_elem_t));

    ff_ext2_elem_t one;
    ff_ext2_one(M->ctx, &one);

    size_t n = (M->rows < M->cols) ? M->rows : M->cols;

    for (size_t i = 0; i < n; ++i)
    {
        M->data[i * M->stride + i] = one;
    }
}

void ff_mat_ext2_mul(
    ff_mat_ext2_t *C,
    const ff_mat_ext2_t *A,
    const ff_mat_ext2_t *B)
{
    if (!A || !B || !C)
        return;
    if (A->cols != B->rows)
        return;
    if (C->rows != A->rows || C->cols != B->cols)
        return;

    const ff_ext2_t *ctx = A->ctx;

    memset(C->data, 0, C->rows * C->stride * sizeof(ff_ext2_elem_t));

    for (size_t i = 0; i < A->rows; ++i)
    {
        for (size_t k = 0; k < A->cols; ++k)
        {
            ff_ext2_elem_t av = A->data[i * A->stride + k];

            if (ff_ext2_is_zero(&av))
                continue;

            for (size_t j = 0; j < B->cols; ++j)
            {
                ff_ext2_elem_t bv = B->data[k * B->stride + j];

                if (ff_ext2_is_zero(&bv))
                    continue;

                ff_ext2_elem_t prod;
                ff_ext2_mul(ctx, &prod, &av, &bv);

                ff_ext2_elem_t *cv = &C->data[i * C->stride + j];
                ff_ext2_add(ctx, cv, cv, &prod);
            }
        }
    }
}

bool ff_mat_ext2_inv(ff_mat_ext2_t *X, const ff_mat_ext2_t *A)
{
    if (!A || !X)
        return false;
    if (A->rows != A->cols)
        return false;
    if (X->rows != A->rows || X->cols != A->cols)
        return false;

    const ff_ext2_t *ctx = A->ctx;
    size_t n = A->rows;
    size_t n2 = 2 * n;

    ff_ext2_elem_t *aug = calloc(n * n2, sizeof(ff_ext2_elem_t));
    if (!aug)
        return false;

    ff_ext2_elem_t one;
    ff_ext2_one(ctx, &one);

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            aug[i * n2 + j] = A->data[i * A->stride + j];
        }
        aug[i * n2 + n + i] = one;
    }

    for (size_t col = 0; col < n; ++col)
    {
        size_t pivot = col;

        while (pivot < n && ff_ext2_is_zero(&aug[pivot * n2 + col]))
        {
            pivot++;
        }

        if (pivot == n)
        {
            free(aug);
            return false;
        }

        if (pivot != col)
        {
            for (size_t j = 0; j < n2; ++j)
            {
                ff_ext2_elem_t t = aug[col * n2 + j];
                aug[col * n2 + j] = aug[pivot * n2 + j];
                aug[pivot * n2 + j] = t;
            }
        }

        ff_ext2_elem_t inv_pivot;

        if (!ff_ext2_inv(ctx, &inv_pivot, &aug[col * n2 + col]))
        {
            free(aug);
            return false;
        }

        for (size_t j = col; j < n2; ++j)
        {
            ff_ext2_mul(ctx, &aug[col * n2 + j], &aug[col * n2 + j], &inv_pivot);
        }

        for (size_t r = 0; r < n; ++r)
        {
            if (r == col)
                continue;

            ff_ext2_elem_t factor = aug[r * n2 + col];

            if (ff_ext2_is_zero(&factor))
                continue;

            for (size_t j = col; j < n2; ++j)
            {
                ff_ext2_elem_t term;
                ff_ext2_mul(ctx, &term, &factor, &aug[col * n2 + j]);

                ff_ext2_sub(ctx, &aug[r * n2 + j], &aug[r * n2 + j], &term);
            }
        }
    }

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            X->data[i * X->stride + j] = aug[i * n2 + n + j];
        }
    }

    free(aug);
    return true;
}

int ff_mat_ext2_nullspace(ff_mat_ext2_t **basis_out, const ff_mat_ext2_t *M)
{
    if (!basis_out || !M)
        return -1;

    *basis_out = NULL;

    const ff_ext2_t *ctx = M->ctx;
    size_t rows = M->rows;
    size_t cols = M->cols;

    ff_ext2_elem_t *aug = calloc(rows * cols, sizeof(ff_ext2_elem_t));
    int *pivot_cols = calloc(cols, sizeof(int));

    if (!aug || !pivot_cols)
    {
        free(aug);
        free(pivot_cols);
        return -1;
    }

    for (size_t i = 0; i < cols; ++i)
        pivot_cols[i] = -1;

    memcpy(aug, M->data, rows * cols * sizeof(ff_ext2_elem_t));

    size_t row = 0;

    for (size_t col = 0; col < cols && row < rows; ++col)
    {
        size_t pivot = row;

        while (pivot < rows && ff_ext2_is_zero(&aug[pivot * cols + col]))
        {
            pivot++;
        }

        if (pivot == rows)
            continue;

        if (pivot != row)
        {
            for (size_t j = 0; j < cols; ++j)
            {
                ff_ext2_elem_t t = aug[row * cols + j];
                aug[row * cols + j] = aug[pivot * cols + j];
                aug[pivot * cols + j] = t;
            }
        }

        ff_ext2_elem_t inv_pivot;

        if (!ff_ext2_inv(ctx, &inv_pivot, &aug[row * cols + col]))
        {
            free(aug);
            free(pivot_cols);
            return -1;
        }

        for (size_t j = col; j < cols; ++j)
        {
            ff_ext2_mul(ctx, &aug[row * cols + j], &aug[row * cols + j], &inv_pivot);
        }

        for (size_t r = 0; r < rows; ++r)
        {
            if (r == row)
                continue;

            ff_ext2_elem_t factor = aug[r * cols + col];

            if (ff_ext2_is_zero(&factor))
                continue;

            for (size_t j = col; j < cols; ++j)
            {
                ff_ext2_elem_t term;
                ff_ext2_mul(ctx, &term, &factor, &aug[row * cols + j]);

                ff_ext2_sub(ctx, &aug[r * cols + j], &aug[r * cols + j], &term);
            }
        }

        pivot_cols[col] = (int)row;
        row++;
    }

    size_t rank = row;
    size_t nullity = cols - rank;

    if (nullity == 0)
    {
        free(aug);
        free(pivot_cols);
        return 0;
    }

    ff_mat_ext2_t *B = ff_mat_ext2_new(ctx, cols, nullity);
    if (!B)
    {
        free(aug);
        free(pivot_cols);
        return -1;
    }

    ff_ext2_elem_t one;
    ff_ext2_one(ctx, &one);

    ff_ext2_elem_t zero;
    ff_ext2_zero(ctx, &zero);

    size_t bidx = 0;

    for (size_t free_col = 0; free_col < cols; ++free_col)
    {
        if (pivot_cols[free_col] != -1)
            continue;

        B->data[free_col * nullity + bidx] = one;

        for (size_t pc = 0; pc < cols; ++pc)
        {
            int prow = pivot_cols[pc];
            if (prow == -1)
                continue;

            ff_ext2_elem_t val = aug[(size_t)prow * cols + free_col];

            ff_ext2_neg(ctx, &B->data[pc * nullity + bidx], &val);
        }

        bidx++;
    }

    free(aug);
    free(pivot_cols);

    *basis_out = B;
    return (int)nullity;
}

/* ============================================================ */
/* Frobenius / Krylov helpers                                   */
/* ============================================================ */

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
    {
        F->data[i * F->stride + (i - 1)] = one;
    }

    for (size_t i = 0; i < n; ++i)
    {
        F->data[i * F->stride + (n - 1)] = ff_neg(ctx, coeffs[i]);
    }
}

bool ff_krylov_matrix_base(
    ff_mat_t *P,
    const ff_mat_t *A,
    const uint32_t *v)
{
    if (!P || !A || !v)
        return false;
    if (P->rows != P->cols)
        return false;
    if (A->rows != A->cols)
        return false;
    if (P->rows != A->rows)
        return false;

    size_t n = A->rows;
    const ff_prime_t *ctx = A->ctx;

    for (size_t i = 0; i < n; ++i)
    {
        P->data[i * P->stride + 0] = v[i];
    }

    for (size_t col = 1; col < n; ++col)
    {
        for (size_t i = 0; i < n; ++i)
        {
            uint32_t sum = 0;

            for (size_t k = 0; k < n; ++k)
            {
                uint32_t av = A->data[i * A->stride + k];
                if (av == 0)
                    continue;

                uint32_t bv = P->data[k * P->stride + (col - 1)];
                if (bv == 0)
                    continue;

                uint32_t prod = ff_mul(ctx, av, bv);
                sum = ff_add(ctx, sum, prod);
            }

            P->data[i * P->stride + col] = sum;
        }
    }

    return true;
}

/* ============================================================ */
/* CRT accumulator                                              */
/* ============================================================ */

struct ff_crt
{
    ff_u128_t val;
    ff_u128_t mod;
};

ff_crt_t *ff_crt_new(void)
{
    ff_crt_t *crt = calloc(1, sizeof *crt);
    if (!crt)
        return NULL;

    crt->val = 0;
    crt->mod = 1;

    return crt;
}

void ff_crt_free(ff_crt_t *crt)
{
    free(crt);
}

void ff_crt_reset(ff_crt_t *crt)
{
    if (!crt)
        return;

    crt->val = 0;
    crt->mod = 1;
}

bool ff_crt_add_residue(
    ff_crt_t *crt,
    uint32_t prime,
    uint32_t residue)
{
    if (!crt || prime < 2)
        return false;

    residue %= prime;

    if ((uint32_t)(crt->mod % prime) == 0)
    {
        return false;
    }

    uint32_t val_mod = (uint32_t)(crt->val % prime);

    uint64_t rhs64;

    if (residue >= val_mod)
    {
        rhs64 = (uint64_t)residue - val_mod;
    }
    else
    {
        rhs64 = (uint64_t)residue + prime - val_mod;
    }

    uint32_t rhs = (uint32_t)(rhs64 % prime);

    uint32_t m_mod = (uint32_t)(crt->mod % prime);
    uint32_t inv = u32_pow_mod(m_mod, prime - 2, prime);

    if (inv == 0)
        return false;

    uint32_t k = (uint32_t)((uint64_t)rhs * inv % prime);

    ff_u128_t max = (ff_u128_t)(-1);

    if (k != 0 && crt->mod > max / k)
    {
        return false;
    }

    ff_u128_t add = crt->mod * k;

    if (crt->val > max - add)
    {
        return false;
    }

    if (crt->mod > max / prime)
    {
        return false;
    }

    crt->val += add;
    crt->mod *= prime;

    return true;
}

bool ff_crt_to_u128(ff_u128_t *out, const ff_crt_t *crt)
{
    if (!out || !crt)
        return false;

    *out = crt->val;
    return true;
}

/* ============================================================ */
/* Opaque Struct Accessors                                      */
/* ============================================================ */

/* Base field matrix accessors */
const ff_prime_t *ff_mat_ctx(const ff_mat_t *M)
{
    return M ? M->ctx : NULL;
}

size_t ff_mat_rows(const ff_mat_t *M)
{
    return M ? M->rows : 0;
}

size_t ff_mat_cols(const ff_mat_t *M)
{
    return M ? M->cols : 0;
}

/* Extension field matrix accessors */
const ff_ext2_t *ff_mat_ext2_ctx(const ff_mat_ext2_t *M)
{
    return M ? M->ctx : NULL;
}

size_t ff_mat_ext2_rows(const ff_mat_ext2_t *M)
{
    return M ? M->rows : 0;
}

size_t ff_mat_ext2_cols(const ff_mat_ext2_t *M)
{
    return M ? M->cols : 0;
}

/* CRT accumulator accessor */
ff_u128_t ff_crt_modulus(const ff_crt_t *crt)
{
    return crt ? crt->mod : 0;
}
