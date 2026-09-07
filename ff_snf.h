/*
ff_snf.h
Smith Normal Form for polynomial matrices over F_p[x].
Specifically computes the invariant factors of the characteristic matrix (xI - A).
*/
#ifndef FF_SNF_H
#define FF_SNF_H

#include "fflib.h"
#include "ff_poly.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /* Order switch for invariant factors */
    typedef enum
    {
        FF_INVARIANT_DESCENDING = 0, /* [d_k, ..., d_1] (Minimal poly first) */
        FF_INVARIANT_ASCENDING = 1   /* [d_1, ..., d_k] (Standard divisibility chain) */
    } ff_inv_order_t;

    /*
     * Computes the Smith Normal Form of the polynomial matrix (xI - A).
     * Extracts the non-unit invariant factors d_1(x) | d_2(x) | ... | d_k(x).
     *
     * out_factors: Array of ff_poly_t pointers. Caller must free each poly and the array.
     * out_count:   Number of invariant factors found.
     * order:       FF_INVARIANT_ASCENDING or FF_INVARIANT_DESCENDING.
     *
     * Returns true on success, false on memory/algorithm failure.
     */
    bool ff_snf_charmatrix_invariants(
        const ff_mat_t *A,
        ff_poly_t ***out_factors,
        size_t *out_count,
        ff_inv_order_t order);

#ifdef __cplusplus
}
#endif

#endif /* FF_SNF_H */
