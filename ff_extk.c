/*
 * ff_extk.c
 *
 * Field core of the generic extension tower F_{p^k} = F_p[x] / <irred>:
 * irreducible-polynomial search, context management, element lifecycle
 * and element arithmetic.
 *
 * Split out of the old monolithic ff_extk.c; the matrices live in
 * ff_mat_extk.c and the Jordan engine in ff_jordan_extk.c (both share
 * ff_extk_internal.h).
 *
 * Canonical-form invariant: every element is a polynomial of degree
 * < k, reduced mod `irred`. ff_extk_zero/one and all arithmetic
 * preserve the element's allocated CAPACITY (they never reallocate);
 * outputs must already be allocated by ff_extk_elem_new.
 */
#include "ff_extk.h"

#include <stdlib.h>

/* -------------------------------------------------------------------- */
/* Irreducible polynomial search                                         */
/* -------------------------------------------------------------------- */

/* f is irreducible over F_p iff:
 *   - for i = 1 .. k/2: gcd(f, x^{p^i} - x) == 1   (no factor of degree <= k/2)
 *   - f divides x^{p^k} - x                        (implied once deg f == k
 *     and no small factor exists)
 * Only the gcd chain is checked here, which is sufficient for monic f
 * of degree k. */
static bool pkf_is_irreducible(const ff_prime_t *ctx, const ff_poly_t *f)
{
    int k = ff_poly_deg(f);
    if (k <= 0)
        return false;

    uint32_t p = ff_prime_modulus(ctx);

    ff_poly_t *x_poly = ff_poly_new(ctx, 1);
    if (!x_poly)
        return false;
    ff_poly_set_coeff(x_poly, 1, ff_one_rep(ctx));

    ff_poly_t *g_pow = ff_poly_copy(x_poly); /* x^{p^i}, kept reduced mod f */
    if (!g_pow)
    {
        ff_poly_free(x_poly);
        return false;
    }

    for (int i = 1; i <= k / 2; ++i)
    {
        ff_poly_t *next_pow = ff_poly_mod_pow_new(g_pow, p, f);
        ff_poly_free(g_pow);
        g_pow = next_pow;

        if (!g_pow)
        {
            ff_poly_free(x_poly);
            return false;
        }

        ff_poly_t *diff = ff_poly_sub_new(g_pow, x_poly);
        if (!diff)
        {
            ff_poly_free(g_pow);
            ff_poly_free(x_poly);
            return false;
        }

        ff_poly_t *h = ff_poly_gcd(diff, f);
        ff_poly_free(diff);

        if (h && ff_poly_deg(h) > 0)
        {
            ff_poly_free(h);
            ff_poly_free(g_pow);
            ff_poly_free(x_poly);
            return false;
        }

        ff_poly_free(h);
    }

    ff_poly_free(g_pow);
    ff_poly_free(x_poly);
    return true;
}

/* Random monic polynomial of degree k until one is irreducible.
 * Expected tries ~ k; loops forever on allocation failure or a broken
 * RNG (rand() is unseeded by this library -- seed it in the caller). */
static ff_poly_t *find_irreducible(const ff_prime_t *ctx, int k)
{
    uint32_t p = ff_prime_modulus(ctx);

    for (;;)
    {
        ff_poly_t *f = ff_poly_new(ctx, k);
        if (!f)
            return NULL;

        ff_poly_set_coeff(f, k, ff_one_rep(ctx));

        for (int i = 0; i < k; ++i)
        {
            uint32_t val = ((uint32_t)rand() << 15) ^ (uint32_t)rand();
            ff_poly_set_coeff(f, i, ff_to_mont(ctx, val % p));
        }

        if (pkf_is_irreducible(ctx, f))
            return f;

        ff_poly_free(f);
    }
}

/* -------------------------------------------------------------------- */
/* Context management                                                    */
/* -------------------------------------------------------------------- */

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

    ff_poly_free(ctx->irred);
    free(ctx);
}

const ff_prime_t *ff_extk_base(const ff_extk_t *ctx)
{
    return ctx ? ctx->base : NULL;
}

int ff_extk_deg(const ff_extk_t *ctx)
{
    return ctx ? ctx->k : 0;
}

/* -------------------------------------------------------------------- */
/* Element lifecycle                                                     */
/* -------------------------------------------------------------------- */

ff_extk_elem_t *ff_extk_elem_new(const ff_extk_t *ctx)
{
    if (!ctx)
        return NULL;
    return ff_poly_new(ctx->base, ctx->k - 1);
}

void ff_extk_elem_free(ff_extk_elem_t *e)
{
    ff_poly_free(e);
}

/* Zero `out` in place; capacity is preserved. */
void ff_extk_zero(const ff_extk_t *ctx, ff_extk_elem_t *out)
{
    (void)ctx;
    if (!out)
        return;

    int d = ff_poly_deg(out);
    for (int i = 0; i <= d; ++i)
        ff_poly_set_coeff(out, i, 0);
}

/* Set `out` to 1 in place; capacity is preserved. */
void ff_extk_one(const ff_extk_t *ctx, ff_extk_elem_t *out)
{
    (void)ctx;
    if (!out)
        return;

    int d = ff_poly_deg(out);
    for (int i = 0; i <= d; ++i)
        ff_poly_set_coeff(out, i, 0);

    ff_poly_set_coeff(out, 0, ff_one_rep(ff_poly_ctx(out)));
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

/* Equality assumes canonical form (true for everything this module
 * produces). */
bool ff_extk_eq(const ff_extk_elem_t *a, const ff_extk_elem_t *b)
{
    if (!a || !b)
        return a == b;

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

/* -------------------------------------------------------------------- */
/* Element arithmetic (all outputs reduced mod irred)                    */
/* -------------------------------------------------------------------- */

/* Copy the coefficients of src into out in place. */
static void extk_assign(ff_poly_t *out, const ff_poly_t *src)
{
    int out_deg = ff_poly_deg(out);
    for (int i = 0; i <= out_deg; ++i)
        ff_poly_set_coeff(out, i, 0);

    int d = ff_poly_deg(src);
    for (int i = 0; i <= d; ++i)
        ff_poly_set_coeff(out, i, ff_poly_coeff(src, i));
}

void ff_extk_add(const ff_extk_t *ctx, ff_extk_elem_t *out,
                 const ff_extk_elem_t *a, const ff_extk_elem_t *b)
{
    (void)ctx;
    if (!out)
        return;

    ff_poly_t *sum = ff_poly_add_new(a, b);
    if (sum)
    {
        extk_assign(out, sum);
        ff_poly_free(sum);
    }
}

void ff_extk_sub(const ff_extk_t *ctx, ff_extk_elem_t *out,
                 const ff_extk_elem_t *a, const ff_extk_elem_t *b)
{
    (void)ctx;
    if (!out)
        return;

    ff_poly_t *diff = ff_poly_sub_new(a, b);
    if (diff)
    {
        extk_assign(out, diff);
        ff_poly_free(diff);
    }
}

void ff_extk_mul(const ff_extk_t *ctx, ff_extk_elem_t *out,
                 const ff_extk_elem_t *a, const ff_extk_elem_t *b)
{
    if (!ctx || !out)
        return;

    ff_poly_t *prod = ff_poly_mul_new(a, b);
    if (!prod)
        return;

    ff_poly_t *rem = ff_poly_mod_new(prod, ctx->irred);
    ff_poly_free(prod);

    if (rem)
    {
        extk_assign(out, rem);
        ff_poly_free(rem);
    }
}

void ff_extk_neg(const ff_extk_t *ctx, ff_extk_elem_t *out,
                 const ff_extk_elem_t *a)
{
    if (!ctx || !out || !a)
        return;

    const ff_prime_t *base = ctx->base;

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

/* a^{-1} = a^{p^k - 2} mod irred (Fermat in the multiplicative group).
 * Slow compared to an extended-Euclid inverse; see review notes. */
bool ff_extk_inv(const ff_extk_t *ctx, ff_extk_elem_t *out,
                 const ff_extk_elem_t *a)
{
    if (!ctx || !out)
        return false;
    if (ff_extk_is_zero(a))
        return false;

    uint32_t p = ff_prime_modulus(ctx->base);
    int k = ctx->k;

    /* order of F_{p^k}^* is p^k - 1; exponent is p^k - 2. */
    ff_u128_t exp = 1;
    for (int i = 0; i < k; ++i)
        exp *= p;
    exp -= 2;

    ff_poly_t *res = ff_poly_new(ctx->base, 0);
    if (!res)
        return false;
    ff_poly_one(res);

    ff_poly_t *base_copy = ff_poly_copy(a);
    if (!base_copy)
    {
        ff_poly_free(res);
        return false;
    }

    while (exp > 0)
    {
        if (exp & 1)
        {
            ff_poly_t *prod = ff_poly_mul_new(res, base_copy);
            ff_poly_free(res);
            res = prod ? ff_poly_mod_new(prod, ctx->irred) : NULL;
            ff_poly_free(prod);

            if (!res)
            {
                ff_poly_free(base_copy);
                return false;
            }
        }

        exp >>= 1;

        if (exp > 0)
        {
            ff_poly_t *sq = ff_poly_mul_new(base_copy, base_copy);
            ff_poly_free(base_copy);
            base_copy = sq ? ff_poly_mod_new(sq, ctx->irred) : NULL;
            ff_poly_free(sq);

            if (!base_copy)
            {
                ff_poly_free(res);
                return false;
            }
        }
    }

    extk_assign(out, res);
    ff_poly_free(res);
    ff_poly_free(base_copy);
    return true;
}

/* -------------------------------------------------------------------- */
/* Context from a caller-supplied irreducible                            */
/* -------------------------------------------------------------------- */

/* Copies `irred`; the caller keeps ownership. Irreducibility is NOT
 * re-checked -- a composite modulus silently breaks inverses. */
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

    if (!ctx->irred)
    {
        free(ctx);
        return NULL;
    }

    return ctx;
}
