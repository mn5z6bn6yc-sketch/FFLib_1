/*
ff_bigint.h
Arbitrary-precision integer backend (GMP wrapper).
Provides Big CRT and Rational Reconstruction for lifting
finite-field computations to exact rational numbers (Q).
*/
#ifndef FF_BIGINT_H
#define FF_BIGINT_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#ifdef __cplusplus
extern "C"
{
#endif

    /* Opaque type: The rest of the library never sees mpz_t */
    typedef struct ff_bigint ff_bigint_t;

    /* Lifecycle */
    ff_bigint_t *ff_bigint_new(void);
    void ff_bigint_free(ff_bigint_t *n);

    /* Setters */
    void ff_bigint_set_u32(ff_bigint_t *n, uint32_t val);
    void ff_bigint_set_i64(ff_bigint_t *n, int64_t val);
    void ff_bigint_set_str(ff_bigint_t *n, const char *str, int base);

    /* Getters */
    void ff_bigint_get_str(char *str, size_t *len, int base, const ff_bigint_t *n);
    int64_t ff_bigint_get_i64(const ff_bigint_t *n);
    bool ff_bigint_is_zero(const ff_bigint_t *n);

    /* ============================================================ */
    /* Big CRT Accumulator (Over Z)                                 */
    /* ============================================================ */
    typedef struct ff_crt_gmp ff_crt_gmp_t;

    ff_crt_gmp_t *ff_crt_gmp_new(void);
    void ff_crt_gmp_free(ff_crt_gmp_t *crt);
    void ff_crt_gmp_reset(ff_crt_gmp_t *crt);

    /* Add a residue from a 32-bit prime field computation */
    bool ff_crt_gmp_add_residue(ff_crt_gmp_t *crt, uint32_t prime, uint32_t residue);

    /* Add a residue from an arbitrary big-int prime */
    bool ff_crt_gmp_add_residue_big(ff_crt_gmp_t *crt, const ff_bigint_t *prime, const ff_bigint_t *residue);

    /* Extract the combined result */
    void ff_crt_gmp_get_modulus(ff_crt_gmp_t *crt, ff_bigint_t *out);
    void ff_crt_gmp_get_value(ff_crt_gmp_t *crt, ff_bigint_t *out);

    /* ============================================================ */
    /* Big Rational Reconstruction                                  */
    /* ============================================================ */
    /*
     * Given x mod M, find n/d such that n * d^{-1} == x mod M
     * with |n|, d <= sqrt(M/2).
     */
    bool ff_ratrec_gmp(const ff_bigint_t *x, const ff_bigint_t *M,
                       ff_bigint_t *num, ff_bigint_t *den);

    /* Convenience: Reconstruct directly from a Big CRT accumulator */
    bool ff_ratrec_gmp_crt(ff_crt_gmp_t *crt, ff_bigint_t *num, ff_bigint_t *den);

#ifdef __cplusplus
}
#endif

#endif /* FF_BIGINT_H */