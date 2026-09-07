/*
 * ff_ratrec.c
 *
 * Rational reconstruction over the integers, fixed-width variant:
 * given x mod M, recover n/d with n * d^{-1} == x (mod M) and
 * |n|, d <= sqrt(M/2), which is unique when it exists.
 *
 * Split out of the old concatenated ff_decomp.c. Uses the half-GCD /
 * continued-fraction stopping rule on the extended Euclidean sequence
 * (r0, r1, ...) / (s0, s1, ...), which maintains s_i * x == r_i (mod M).
 *
 * For moduli beyond 2^127 (or results outside the i128 range) use the
 * GMP-backed ff_ratrec_gmp in ff_bigint.c instead.
 *
 * Note: depends on ff_crt_modulus(), which must be declared in fflib.h
 * (see fflib.patch -- the accessor existed in fflib.c but was missing
 * from the header).
 */
#include "ff_decomp.h"

#include <stdlib.h>
#include <stdint.h>

/* -------------------------------------------------------------------- */
/* Small 128-bit helpers                                                 */
/* -------------------------------------------------------------------- */

static ff_i128_t ffd_i128_abs(ff_i128_t x)
{
    return (x < 0) ? -x : x;
}

static ff_i128_t ffd_i128_gcd(ff_i128_t a, ff_i128_t b)
{
    if (a < 0)
        a = -a;
    if (b < 0)
        b = -b;

    while (b != 0)
    {
        ff_i128_t t = a % b;
        a = b;
        b = t;
    }

    return a;
}

/* Floor(sqrt(x)) via binary search; the `mid <= x / mid` form avoids
 * squaring, so no overflow for any ff_u128_t x. */
static ff_u128_t ffd_u128_sqrt(ff_u128_t x)
{
    if (x == 0)
        return 0;

    ff_u128_t lo = 0;
    ff_u128_t hi = ((ff_u128_t)1) << 64;

    while (lo < hi)
    {
        ff_u128_t mid = (lo + hi + 1) >> 1;

        if (mid <= x / mid)
            lo = mid;
        else
            hi = mid - 1;
    }

    return lo;
}

/* -------------------------------------------------------------------- */
/* Core reconstruction                                                   */
/* -------------------------------------------------------------------- */

bool ff_ratrec_from_u128(
    ff_u128_t x,
    ff_u128_t M,
    ff_i128_t *num,
    ff_i128_t *den)
{
    if (!num || !den)
        return false;
    if (M == 0)
        return false;

    if (x >= M)
        x %= M;

    /* Everything downstream is signed 128-bit. */
    const ff_u128_t I128_MAX_U = (((ff_u128_t)1) << 127) - 1;
    if (M > I128_MAX_U)
        return false;

    if (x == 0)
    {
        *num = 0;
        *den = 1;
        return true;
    }

    ff_u128_t limit_u = ffd_u128_sqrt(M >> 1);
    if (limit_u == 0)
        limit_u = 1;

    ff_i128_t limit = (ff_i128_t)limit_u;

    /* Extended Euclid on (M, x), stopping once the remainder drops to
     * the rational-reconstruction bound. Invariant: s_i * x == r_i. */
    ff_i128_t r0 = (ff_i128_t)M;
    ff_i128_t r1 = (ff_i128_t)x;

    ff_i128_t s0 = 0;
    ff_i128_t s1 = 1;

    while (r1 > limit)
    {
        ff_i128_t q = r0 / r1;

        ff_i128_t r2 = r0 - q * r1;
        ff_i128_t s2 = s0 - q * s1;

        r0 = r1;
        r1 = r2;

        s0 = s1;
        s1 = s2;
    }

    if (r1 == 0)
        return false;

    ff_i128_t n = r1;
    ff_i128_t d = s1;

    /* Normalize the sign onto the numerator. */
    if (d < 0)
    {
        n = -n;
        d = -d;
    }

    if (d <= 0)
        return false;

    if (ffd_i128_abs(n) > limit)
        return false;
    if (d > limit)
        return false;

    /* Uniqueness requires the pair to be coprime. */
    if (ffd_i128_gcd(n, d) != 1)
        return false;

    *num = n;
    *den = d;
    return true;
}

/* -------------------------------------------------------------------- */
/* Convenience wrappers                                                  */
/* -------------------------------------------------------------------- */

bool ff_ratrec_crt_i128(
    const ff_crt_t *crt,
    ff_i128_t *num,
    ff_i128_t *den)
{
    if (!crt)
        return false;

    ff_u128_t x = 0;
    ff_u128_t M = ff_crt_modulus(crt);

    if (!ff_crt_to_u128(&x, crt))
        return false;

    return ff_ratrec_from_u128(x, M, num, den);
}

bool ff_ratrec_crt_i64(
    const ff_crt_t *crt,
    int64_t *num,
    int64_t *den)
{
    if (!num || !den)
        return false;

    ff_i128_t n128 = 0;
    ff_i128_t d128 = 0;

    if (!ff_ratrec_crt_i128(crt, &n128, &d128))
        return false;

    const ff_i128_t I64_MIN = (ff_i128_t)INT64_MIN;
    const ff_i128_t I64_MAX = (ff_i128_t)INT64_MAX;

    if (n128 < I64_MIN || n128 > I64_MAX)
        return false;
    if (d128 < 1 || d128 > I64_MAX)
        return false;

    *num = (int64_t)n128;
    *den = (int64_t)d128;
    return true;
}

bool ff_ratrec_array_i128(
    const uint32_t *primes,
    const uint32_t *residues,
    size_t count,
    ff_i128_t *num,
    ff_i128_t *den)
{
    if (!primes || !residues || count == 0)
        return false;

    ff_crt_t *crt = ff_crt_new();
    if (!crt)
        return false;

    for (size_t i = 0; i < count; ++i)
    {
        if (!ff_crt_add_residue(crt, primes[i], residues[i]))
        {
            ff_crt_free(crt);
            return false;
        }
    }

    bool ok = ff_ratrec_crt_i128(crt, num, den);

    ff_crt_free(crt);
    return ok;
}
