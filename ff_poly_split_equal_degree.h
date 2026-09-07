/*
 * ff_poly_split_equal_degree.h
 *
 * Equal-Degree Factorization (Cantor-Zassenhaus) over F_p.
 *
 * Renamed from ff_poly_factor_equal_degree.h so the file name matches the
 * symbol it declares (ff_poly_split_equal_degree, implemented in
 * ff_poly_split_equal_degree.c). Implementation now lives in
 * ff_poly_factor.c.
 */
#ifndef FF_POLY_SPLIT_EQUAL_DEGREE_H
#define FF_POLY_SPLIT_EQUAL_DEGREE_H

#include <stdbool.h>
#include <stddef.h>
#include "ff_poly.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /*
     * Splits g into its individual irreducible factors, all of degree k.
     *
     * Preconditions (NOT checked; violating them yields wrong factors):
     *   - g is square-free,
     *   - EVERY irreducible factor of g has degree exactly k
     *     (hence k must divide deg(g) -- that part IS checked).
     *
     * For odd p this is the classical Cantor-Zassenhaus split using
     * h = a^((p^k - 1)/2) mod g. Requires p to be odd.
     *
     * On success: *out_factors = heap array of *out_count factor
     * pointers; caller frees each factor and the array (see the shared
     * ff_poly_factor_list_free helper once it exists).
     * On failure: returns false and sets *out_factors = NULL,
     * *out_count = 0.
     */
    bool ff_poly_split_equal_degree(
        const ff_poly_t *g,
        int k,
        ff_poly_t ***out_factors,
        size_t *out_count);

#ifdef __cplusplus
}
#endif

#endif /* FF_POLY_SPLIT_EQUAL_DEGREE_H */
