// poly_pow_mod_impl.c
#include "poly_pow_mod_impl.h"

ff_poly_t *poly_pow_mod_impl(const ff_poly_t *h, uint64_t e, const ff_poly_t *mod)
{
    const ff_prime_t *ctx = ff_poly_ctx(mod);
    if (!ctx)
        return NULL;

    ff_poly_t *res = ff_poly_new(ctx, 0);
    if (!res)
        return NULL;
    ff_poly_one(res);

    ff_poly_t *base = ff_poly_copy(h);
    if (!base)
    {
        ff_poly_free(res);
        return NULL;
    }

    while (e > 0)
    {
        if (e & 1)
        {
            ff_poly_t *tmp = ff_poly_mul_new(res, base);
            ff_poly_free(res);
            res = ff_poly_mod_new(tmp, mod);
            ff_poly_free(tmp);
            if (!res)
            {
                ff_poly_free(base);
                return NULL;
            }
        }

        ff_poly_t *tmp2 = ff_poly_mul_new(base, base);
        ff_poly_free(base);
        base = ff_poly_mod_new(tmp2, mod);
        ff_poly_free(tmp2);
        if (!base)
        {
            ff_poly_free(res);
            return NULL;
        }

        e >>= 1;
    }

    ff_poly_free(base);
    return res;
}
