/*
 * ff_extk_internal.h
 *
 * INTERNAL header for the F_{p^k} extension module. Not installed, not
 * part of the public API. Shared by the three translation units that
 * resulted from splitting the old monolithic ff_extk.c:
 *
 *   ff_extk.c        field core (irreducible search, elements, arithmetic)
 *   ff_mat_extk.c    matrices over F_{p^k}
 *   ff_jordan_extk.c generic Jordan engine over F_{p^k}
 *
 * The matrix struct must stay private (callers go through the opaque
 * ff_mat_extk_t and its accessors), but the module's own files need the
 * layout.
 */
#ifndef FF_EXTK_INTERNAL_H
#define FF_EXTK_INTERNAL_H

#include "ff_extk.h"

#ifdef __cplusplus
extern "C"
{
#endif

    struct ff_mat_extk
    {
        const ff_extk_t *ctx;
        size_t rows;
        size_t cols;
        size_t stride;      /* entries in data per row (== cols) */
        ff_poly_t **data;   /* rows*cols element pointers, each of
                               capacity ctx->k - 1 */
    };

    /*
     * Copy src into dst IN PLACE, preserving dst's allocated capacity.
     * dst's coefficients are zeroed up to its current degree first.
     * Used by every extk routine that must not disturb capacity.
     * Implemented in ff_mat_extk.c.
     */
    void pkf_extk_poly_copy_to(ff_poly_t *dst, const ff_poly_t *src);

#ifdef __cplusplus
}
#endif

#endif /* FF_EXTK_INTERNAL_H */
