// ff_poly_factor_equal_degree.c
#include <stdlib.h>
#include "ff_poly_factor_equal_degree.h"
#include "poly_random_impl.h"
#include "poly_pow_mod_impl.h"

bool ff_poly_split_equal_degree(
    const ff_poly_t *g_in,
    int k,
    ff_poly_t ***out_factors,
    size_t *out_count)
{
    if (!g_in || !out_factors || !out_count || k <= 0)
        return false;

    const ff_prime_t *ctx = ff_poly_ctx(g_in);
    if (!ctx)
        return false;

    int deg_g = ff_poly_deg(g_in);
    if (deg_g <= 0 || deg_g % k != 0)
        return false;

    uint32_t q = ff_prime_modulus(ctx);

    ff_poly_t *g = ff_poly_copy(g_in);
    if (!g)
        return false;

    size_t cap = 4;
    size_t count = 0;
    ff_poly_t **factors = calloc(cap, sizeof(ff_poly_t *));
    if (!factors)
    {
        ff_poly_free(g);
        return false;
    }

    while (ff_poly_deg(g) > k)
    {
        int retries = 0;
        ff_poly_t *a = NULL;

        while (retries < 64)
        {
            ff_poly_t *h = poly_random_impl(ctx, q, ff_poly_deg(g));
            if (!h)
            {
                ff_poly_free(g);
                free(factors);
                return false;
            }

            uint64_t e = 1;
            for (int i = 0; i < k; ++i)
                e *= (uint64_t)q;
            e = (e - 1) / 2;

            ff_poly_t *pow = poly_pow_mod_impl(h, e, g);
            ff_poly_free(h);
            if (!pow)
            {
                ff_poly_free(g);
                free(factors);
                return false;
            }

            ff_poly_t *pow_minus_one = ff_poly_copy(pow);
            ff_poly_free(pow);
            if (!pow_minus_one)
            {
                ff_poly_free(g);
                free(factors);
                return false;
            }

            uint32_t one = ff_one_rep(ctx);
            uint32_t c0 = ff_poly_coeff(pow_minus_one, 0);
            ff_poly_set_coeff(pow_minus_one, 0, ff_sub(ctx, c0, one));
            ff_poly_trim(pow_minus_one);

            a = ff_poly_gcd(pow_minus_one, g);
            ff_poly_free(pow_minus_one);

            if (!a)
            {
                ff_poly_free(g);
                free(factors);
                return false;
            }

            int deg_a = ff_poly_deg(a);

            if (deg_a > 0 && deg_a < ff_poly_deg(g))
                break;

            ff_poly_free(a);
            a = NULL;
            retries++;
        }

        if (!a)
        {
            ff_poly_free(g);
            free(factors);
            return false;
        }

        ff_poly_t *b = ff_poly_div_new(g, a);
        if (!b)
        {
            ff_poly_free(a);
            ff_poly_free(g);
            free(factors);
            return false;
        }

        if (count == cap)
        {
            cap *= 2;
            ff_poly_t **tmp = realloc(factors, cap * sizeof(ff_poly_t *));
            if (!tmp)
            {
                ff_poly_free(a);
                ff_poly_free(b);
                ff_poly_free(g);
                free(factors);
                return false;
            }
            factors = tmp;
        }

        factors[count++] = a;
        ff_poly_free(g);
        g = b;
    }

    if (count == cap)
    {
        cap *= 2;
        ff_poly_t **tmp = realloc(factors, cap * sizeof(ff_poly_t *));
        if (!tmp)
        {
            ff_poly_free(g);
            free(factors);
            return false;
        }
        factors = tmp;
    }

    factors[count++] = g;

    *out_factors = factors;
    *out_count = count;
    return true;
}
