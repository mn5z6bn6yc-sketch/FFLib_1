// poly_random_impl.h
#ifndef POLY_RANDOM_IMPL_H
#define POLY_RANDOM_IMPL_H

#include "ff_poly.h"
#include "fflib.h"

ff_poly_t *poly_random_impl(const ff_prime_t *ctx, uint32_t q, int deg);

#endif