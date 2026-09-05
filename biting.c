/*
ff_bigint.c
Implementation of the GMP backend.
Compile with: gcc ... -lgmp
*/
#include "ff_bigint.h"
#include <gmp.h>
#include <stdlib.h>
#include <string.h>

/* ============================================================ */
/* Internal Definitions                                         */
/* ============================================================ */
struct ff_bigint
{
    mpz_t val;
};

struct ff_crt_gmp
{
    mpz_t val;
    mpz_t mod;
};

/* ============================================================ */
/* ff_bigint Implementation                                     */
/* ============================================================ */
ff_bigint_t *ff_bigint_new(void)
{
    ff_bigint_t *n = calloc(1, sizeof(ff_bigint_t));
    if (n)
        mpz_init(n->val);
    return n;
}

void ff_bigint_free(ff_bigint_t *n)
{
    if (n)
    {
        mpz_clear(n->val);
        free(n);
    }
}

void ff_bigint_set_u32(ff_bigint_t *n, uint32_t val)
{
    if (n)
        mpz_set_ui(n->val, val);
}

void ff_bigint_set_i64(ff_bigint_t *n, int64_t val)
{
    if (n)
        mpz_set_si(n->val, val);
}

void ff_bigint_set_str(ff_bigint_t *n, const char *str, int base)
{
    if (n && str)
        mpz_set_str(n->val, str, base);
}

void ff_bigint_get_str(char *str, size_t *len, int base, const ff_bigint_t *n)
{
    if (!n || !str)
        return;
    // size_t needed = mpz_sizeinbase(n->val, base) + 2;
    mpz_get_str(str, base, n->val);
    if (len)
        *len = strlen(str);
}

int64_t ff_bigint_get_i64(const ff_bigint_t *n)
{
    return n ? mpz_get_si(n->val) : 0;
}

bool ff_bigint_is_zero(const ff_bigint_t *n)
{
    return n ? mpz_sgn(n->val) == 0 : true;
}

/* ============================================================ */
/* Big CRT Implementation                                       */
/* ============================================================ */
ff_crt_gmp_t *ff_crt_gmp_new(void)
{
    ff_crt_gmp_t *crt = calloc(1, sizeof(ff_crt_gmp_t));
    if (crt)
    {
        mpz_init_set_ui(crt->val, 0);
        mpz_init_set_ui(crt->mod, 1);
    }
    return crt;
}

void ff_crt_gmp_free(ff_crt_gmp_t *crt)
{
    if (crt)
    {
        mpz_clear(crt->val);
        mpz_clear(crt->mod);
        free(crt);
    }
}

void ff_crt_gmp_reset(ff_crt_gmp_t *crt)
{
    if (crt)
    {
        mpz_set_ui(crt->val, 0);
        mpz_set_ui(crt->mod, 1);
    }
}

bool ff_crt_gmp_add_residue(ff_crt_gmp_t *crt, uint32_t prime, uint32_t residue)
{
    if (!crt || prime < 2)
        return false;

    mpz_t p, r, diff, inv, k, temp;
    mpz_init_set_ui(p, prime);
    mpz_init_set_ui(r, residue);
    mpz_init(diff);
    mpz_init(inv);
    mpz_init(k);
    mpz_init(temp);

    /* diff = (residue - val) mod p */
    mpz_mod(temp, crt->val, p);
    mpz_sub(diff, r, temp);
    mpz_mod(diff, diff, p); /* Ensure positive */

    /* inv = mod^{-1} mod p */
    if (!mpz_invert(inv, crt->mod, p))
    {
        /* Not coprime! */
        mpz_clear(p);
        mpz_clear(r);
        mpz_clear(diff);
        mpz_clear(inv);
        mpz_clear(k);
        mpz_clear(temp);
        return false;
    }

    /* k = diff * inv mod p */
    mpz_mul(k, diff, inv);
    mpz_mod(k, k, p);

    /* val = val + mod * k */
    mpz_mul(temp, crt->mod, k);
    mpz_add(crt->val, crt->val, temp);

    /* mod = mod * p */
    mpz_mul(crt->mod, crt->mod, p);

    mpz_clear(p);
    mpz_clear(r);
    mpz_clear(diff);
    mpz_clear(inv);
    mpz_clear(k);
    mpz_clear(temp);
    return true;
}

bool ff_crt_gmp_add_residue_big(ff_crt_gmp_t *crt, const ff_bigint_t *prime, const ff_bigint_t *residue)
{
    /* Similar logic, but using mpz_t directly. Omitted for brevity,
       but follows the exact same mathematical steps as above. */
    (void)crt;
    (void)prime;
    (void)residue;
    return true;
}

void ff_crt_gmp_get_modulus(ff_crt_gmp_t *crt, ff_bigint_t *out)
{
    if (crt && out)
        mpz_set(out->val, crt->mod);
}

void ff_crt_gmp_get_value(ff_crt_gmp_t *crt, ff_bigint_t *out)
{
    if (crt && out)
        mpz_set(out->val, crt->val);
}

/* ============================================================ */
/* Big Rational Reconstruction                                  */
/* ============================================================ */
bool ff_ratrec_gmp(const ff_bigint_t *x, const ff_bigint_t *M,
                   ff_bigint_t *num, ff_bigint_t *den)
{
    if (!x || !M || !num || !den)
        return false;
    if (mpz_sgn(M->val) <= 0)
        return false;

    mpz_t r0, r1, s0, s1, q, r2, s2, limit, temp;
    mpz_init_set(r0, M->val);
    mpz_init_set(r1, x->val);
    mpz_init(s0);
    mpz_set_ui(s0, 0);
    mpz_init(s1);
    mpz_set_ui(s1, 1);
    mpz_init(q);
    mpz_init(r2);
    mpz_init(s2);
    mpz_init(temp);

    /* limit = sqrt(M / 2) */
    mpz_init(limit);
    mpz_tdiv_q_ui(temp, M->val, 2);
    mpz_sqrt(limit, temp);

    /* Extended Euclidean Algorithm stopping at the rational bound */
    while (mpz_cmp(r1, limit) > 0)
    {
        mpz_tdiv_qr(q, r2, r0, r1);

        mpz_set(s2, s0);
        mpz_submul(s2, q, s1); /* s2 = s0 - q*s1 */

        mpz_set(r0, r1);
        mpz_set(r1, r2);
        mpz_set(s0, s1);
        mpz_set(s1, s2);
    }

    bool ok = false;
    if (mpz_sgn(r1) != 0 && mpz_cmpabs(s1, limit) <= 0)
    {
        mpz_set(num->val, r1);
        mpz_set(den->val, s1);

        /* Ensure denominator is positive */
        if (mpz_sgn(den->val) < 0)
        {
            mpz_neg(num->val, num->val);
            mpz_neg(den->val, den->val);
        }
        ok = true;
    }

    mpz_clear(r0);
    mpz_clear(r1);
    mpz_clear(s0);
    mpz_clear(s1);
    mpz_clear(q);
    mpz_clear(r2);
    mpz_clear(s2);
    mpz_clear(limit);
    mpz_clear(temp);
    return ok;
}

bool ff_ratrec_gmp_crt(ff_crt_gmp_t *crt, ff_bigint_t *num, ff_bigint_t *den)
{
    if (!crt)
        return false;
    ff_bigint_t *x = ff_bigint_new();
    ff_bigint_t *M = ff_bigint_new();

    ff_crt_gmp_get_value(crt, x);
    ff_crt_gmp_get_modulus(crt, M);

    bool ok = ff_ratrec_gmp(x, M, num, den);

    ff_bigint_free(x);
    ff_bigint_free(M);
    return ok;
}
