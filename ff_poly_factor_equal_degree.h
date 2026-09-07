// ff_poly_factor_equal_degree.h
#ifndef FF_POLY_FACTOR_EQUAL_DEGREE_H
#define FF_POLY_FACTOR_EQUAL_DEGREE_H

#include <stdbool.h>
#include <stddef.h>
#include "ff_poly.h"

bool ff_poly_split_equal_degree(
    const ff_poly_t *g_in,
    int k,
    ff_poly_t ***out_factors,
    size_t *out_count);

#endif
