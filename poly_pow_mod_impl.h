// poly_pow_mod_impl.h
#ifndef POLY_POW_MOD_IMPL_H
#define POLY_POW_MOD_IMPL_H

#include "ff_poly.h"

ff_poly_t *poly_pow_mod_impl(const ff_poly_t *h,
                             uint64_t e,
                             const ff_poly_t *mod);

#endif