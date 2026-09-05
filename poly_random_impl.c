// poly_random_impl.c
#include <stdlib.h>
#include "poly_random_impl.h"

ff_poly_t *poly_random_impl(const ff_prime_t *ctx, uint32_t q, int deg)
{
    ff_poly_t *r = ff_poly_new(ctx, deg - 1);
    if (!r)
        return NULL;

    for (int i = 0; i < deg; ++i)
    {
        uint32_t coeff = (uint32_t)(rand() % q);
        ff_poly_set_coeff(r, i, coeff);
    }
    ff_poly_trim(r);
    return r;
}