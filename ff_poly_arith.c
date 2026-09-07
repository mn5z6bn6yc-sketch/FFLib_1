/*
 * ff_poly_arith.c
 *
 * Polynomial division, gcd, modular exponentiation and derivative.
 * Split out of the old monolithic ff_poly.c.
 *
 * This file is the SINGLE home of modular exponentiation. The old tree
 * had five near-identical implementations (pkf_poly_mod_pow,
 * pkf_poly_mod_pow_128, ff_poly_pow_mod, ff_poly_x_pow_mod and the public
 * wrapper); everything now routes through pkf_poly_pow_mod_u128.
 */
#include "ff_poly_internal.h"

#include <stdlib.h>

/* -------------------------------------------------------------------- */
/* Division                                                              */
/* -------------------------------------------------------------------- */

/*
 * Polynomial long division: f = q * g + r, deg r < deg g.
 * Returns the remainder r; if q_out is non-NULL, stores the quotient.
 * Returns NULL if g is NULL or zero. Both results are monic-free (raw)
 * and trimmed.
 *
 * NEW PUBLIC ENTRY POINT (was the static pkf_poly_divrem): ff_snf.c and
 * other callers no longer need private copies of this routine.
 */
ff_poly_t *ff_poly_divrem_new(const ff_poly_t *f,
                              const ff_poly_t *g,
                              ff_poly_t **q_out)
{
    if (q_out)
        *q_out = NULL;

    if (!f || !g || pkf_poly_is_zero(g))
        return NULL;

    const ff_prime_t *ctx = f->ctx;
    int deg_f = f->deg;
    int deg_g = g->deg;

    ff_poly_t *rem = ff_poly_copy(f);
    if (!rem)
        return NULL;

    ff_poly_t *q = NULL;
    if (q_out)
    {
        q = ff_poly_new(ctx, (deg_f >= deg_g) ? deg_f - deg_g : 0);
        if (!q)
        {
            ff_poly_free(rem);
            return NULL;
        }
    }

    if (deg_f < deg_g)
    {
        if (q_out)
            *q_out = q;
        return rem;
    }

    uint32_t inv_lead = ff_inv(ctx, g->c[deg_g]);
    if (inv_lead == 0)
    {
        ff_poly_free(rem);
        ff_poly_free(q);
        return NULL;
    }

    for (int i = deg_f; i >= deg_g; --i)
    {
        if (rem->c[i] == 0)
            continue;

        uint32_t coeff = ff_mul(ctx, rem->c[i], inv_lead);
        if (q)
            q->c[i - deg_g] = coeff;

        for (int j = 0; j <= deg_g; ++j)
        {
            uint32_t sub = ff_mul(ctx, coeff, g->c[j]);
            rem->c[i - deg_g + j] = ff_sub(ctx, rem->c[i - deg_g + j], sub);
        }
    }

    ff_poly_trim(rem);
    if (q)
        ff_poly_trim(q);

    if (q_out)
        *q_out = q;

    return rem;
}

/* Quotient only: f / g (exact or truncated). */
ff_poly_t *ff_poly_div_new(const ff_poly_t *a, const ff_poly_t *b)
{
    ff_poly_t *q = NULL;
    ff_poly_t *r = ff_poly_divrem_new(a, b, &q);
    ff_poly_free(r);
    return q;
}

/* Remainder only: f mod g. */
ff_poly_t *ff_poly_mod_new(const ff_poly_t *num, const ff_poly_t *den)
{
    return ff_poly_divrem_new(num, den, NULL);
}

/* -------------------------------------------------------------------- */
/* GCD                                                                   */
/* -------------------------------------------------------------------- */

/* Euclidean GCD, normalized to be monic. gcd(0, b) = monic(b). */
ff_poly_t *ff_poly_gcd(const ff_poly_t *a, const ff_poly_t *b)
{
    if (!a || !b)
        return NULL;

    const ff_prime_t *ctx = a->ctx;

    ff_poly_t *x = ff_poly_copy(a);
    ff_poly_t *y = ff_poly_copy(b);
    if (!x || !y)
    {
        ff_poly_free(x);
        ff_poly_free(y);
        return NULL;
    }

    while (!pkf_poly_is_zero(y))
    {
        ff_poly_t *r = ff_poly_divrem_new(x, y, NULL);
        ff_poly_free(x);
        x = y;
        y = r;
    }
    ff_poly_free(y);

    /* Make monic. */
    if (x->deg >= 0 && x->c[x->deg] != 0)
    {
        uint32_t inv = ff_inv(ctx, x->c[x->deg]);
        for (int i = 0; i <= x->deg; ++i)
            x->c[i] = ff_mul(ctx, x->c[i], inv);
    }

    return x;
}

/* -------------------------------------------------------------------- */
/* Modular exponentiation (single implementation)                        */
/* -------------------------------------------------------------------- */

ff_poly_t *pkf_poly_pow_mod_u128(const ff_poly_t *base,
                                 ff_u128_t exp,
                                 const ff_poly_t *mod)
{
    if (!base || !mod)
        return NULL;

    if (ff_poly_deg(mod) <= 0)
        return NULL; /* modulus must have degree >= 1 */

    const ff_prime_t *ctx = ff_poly_ctx(base);

    ff_poly_t *res = ff_poly_new(ctx, 0);
    if (!res)
        return NULL;
    ff_poly_one(res);

    ff_poly_t *b = ff_poly_mod_new(base, mod);
    if (!b)
    {
        ff_poly_free(res);
        return NULL;
    }

    while (exp > 0)
    {
        if (exp & 1)
        {
            ff_poly_t *prod = ff_poly_mul_new(res, b);
            ff_poly_free(res);
            res = prod ? ff_poly_mod_new(prod, mod) : NULL;
            ff_poly_free(prod);

            if (!res)
            {
                ff_poly_free(b);
                return NULL;
            }
        }

        exp >>= 1;

        if (exp > 0)
        {
            ff_poly_t *sq = ff_poly_mul_new(b, b);
            ff_poly_free(b);
            b = sq ? ff_poly_mod_new(sq, mod) : NULL;
            ff_poly_free(sq);

            if (!b)
            {
                ff_poly_free(res);
                return NULL;
            }
        }
    }

    ff_poly_free(b);
    return res;
}

/* Public wrapper (exponent widened from uint32_t to uint64_t). */
ff_poly_t *ff_poly_mod_pow_new(const ff_poly_t *base,
                               uint64_t exp,
                               const ff_poly_t *mod)
{
    return pkf_poly_pow_mod_u128(base, (ff_u128_t)exp, mod);
}

/* -------------------------------------------------------------------- */
/* Derivative                                                            */
/* -------------------------------------------------------------------- */

ff_poly_t *ff_poly_deriv(const ff_poly_t *p)
{
    if (!p)
        return NULL;

    if (p->deg <= 0)
    {
        ff_poly_t *z = ff_poly_new(p->ctx, 0);
        if (z)
            ff_poly_zero(z);
        return z;
    }

    const ff_prime_t *ctx = p->ctx;
    uint32_t mod = ff_prime_modulus(ctx);

    ff_poly_t *d = ff_poly_new(ctx, p->deg - 1);
    if (!d)
        return NULL;

    for (int i = 1; i <= p->deg; ++i)
    {
        uint32_t coeff = p->c[i];
        uint32_t i_rep = ff_to_mont(ctx, (uint32_t)i % mod);
        d->c[i - 1] = ff_mul(ctx, coeff, i_rep);
    }

    ff_poly_trim(d);
    return d;
}
