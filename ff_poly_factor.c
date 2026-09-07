/*
 * ff_poly_factor.c
 *
 * Factorization over F_p:
 *   - ff_poly_split_equal_degree : Cantor-Zassenhaus EDF (the single,
 *     correct recursive implementation; replaces THREE old copies:
 *     the iterative ff_poly_edf_base and ff_poly_split_equal_degree,
 *     and the static pkf_edf_cantor_zassenhaus)
 *   - ff_poly_roots_base         : roots via DDF (gcd with x^p - x) + EDF
 *   - ff_poly_factor_berlekamp   : deterministic factorization
 *
 * Split out of the old monolithic ff_poly.c.
 *
 * FIXES applied while moving:
 *   - The old iterative EDF stored gcd(h^e - 1, g) as a final factor even
 *     when it was a product of several degree-k irreducibles, returning
 *     REDUCIBLE factors for inputs with >= 3 factors. The recursive form
 *     splits both halves until degree k.
 *   - Exponents (q^k - 1)/2 are now computed in ff_u128 (the old uint64_t
 *     arithmetic overflowed for large q, k >= 3).
 *   - Berlekamp's matrix Q is built with COLUMNS = coefficients of
 *     x^(q j). The old live code filled ROWS, i.e. built Q^T: the nullity
 *     count survived (rank is transpose-invariant) but the splitting
 *     phase used left-nullspace vectors as if they were Frobenius fixed
 *     points, which is unsound.
 *   - All failure paths free already-collected factors (the old code
 *     leaked them).
 */
#include "ff_poly_internal.h"
#include "ff_poly_split_equal_degree.h"
#include "poly_random_impl.h"

#include <stdlib.h>

/* -------------------------------------------------------------------- */
/* Shared helpers                                                        */
/* -------------------------------------------------------------------- */

/* Append poly to a growable factor list; takes ownership of poly. */
static bool factor_list_push(ff_poly_t ***list, size_t *count, size_t *cap,
                             ff_poly_t *poly)
{
    if (*count == *cap)
    {
        size_t ncap = (*cap == 0) ? 4 : (*cap * 2);
        ff_poly_t **tmp = realloc(*list, ncap * sizeof *tmp);
        if (!tmp)
            return false;
        *list = tmp;
        *cap = ncap;
    }

    (*list)[(*count)++] = poly;
    return true;
}

static void factor_list_free(ff_poly_t **list, size_t count)
{
    if (!list)
        return;
    for (size_t i = 0; i < count; ++i)
        ff_poly_free(list[i]);
    free(list);
}

/* -------------------------------------------------------------------- */
/* Equal-degree factorization (Cantor-Zassenhaus)                        */
/* -------------------------------------------------------------------- */

/*
 * Precondition: f is square-free and EVERY irreducible factor of f has
 * degree exactly d. Violating this yields wrong factors, not an error.
 *
 * Classical C-Z requires odd p (exponent (p^d - 1) / 2). For p == 2 use
 * the trace-based Ben-Or variant instead (pending).
 */
static bool edf_recurse(const ff_poly_t *f, int d,
                        ff_poly_t ***factors, size_t *count, size_t *cap)
{
    if (f->deg == d)
    {
        /* Already irreducible (by the equal-degree precondition). */
        ff_poly_t *copy = ff_poly_copy(f);
        if (!copy)
            return false;
        if (!factor_list_push(factors, count, cap, copy))
        {
            ff_poly_free(copy);
            return false;
        }
        return true;
    }

    const ff_prime_t *ctx = f->ctx;
    uint32_t p = ff_prime_modulus(ctx);

    /* e = (p^d - 1) / 2, in 128 bits to avoid overflow. */
    ff_u128_t pd = 1;
    for (int i = 0; i < d; ++i)
        pd *= (ff_u128_t)p;
    ff_u128_t e = (pd - 1) / 2;

    for (int attempt = 0; attempt < 256; ++attempt)
    {
        /* Random g of degree < deg(f). poly_random_impl(ctx, p, deg)
         * allocates capacity deg-1 and fills coefficients 0..deg-1. */
        ff_poly_t *g = poly_random_impl(ctx, p, f->deg);
        if (!g)
            return false;

        if (pkf_poly_is_zero(g))
        {
            ff_poly_free(g);
            continue;
        }

        /* h = g^e mod f; for each irreducible factor of f, h is
         * congruent to 0, 1 or -1 modulo that factor. */
        ff_poly_t *h = pkf_poly_pow_mod_u128(g, e, f);
        ff_poly_free(g);
        if (!h)
            continue;

        /* h - 1 */
        uint32_t c0 = ff_poly_coeff(h, 0);
        ff_poly_set_coeff(h, 0, ff_sub(ctx, c0, ff_one_rep(ctx)));
        ff_poly_trim(h);

        ff_poly_t *gcd1 = ff_poly_gcd(f, h);
        ff_poly_free(h);
        if (!gcd1)
            continue;

        int dg = ff_poly_deg(gcd1);

        if (dg > 0 && dg < f->deg)
        {
            /* Non-trivial split: f = gcd1 * rem; split BOTH halves. */
            ff_poly_t *rem = ff_poly_div_new(f, gcd1);
            if (!rem)
            {
                ff_poly_free(gcd1);
                continue;
            }

            bool ok = edf_recurse(gcd1, d, factors, count, cap) &&
                      edf_recurse(rem, d, factors, count, cap);

            ff_poly_free(gcd1);
            ff_poly_free(rem);
            return ok;
        }

        ff_poly_free(gcd1);
    }

    return false; /* exhausted retries */
}

bool ff_poly_split_equal_degree(const ff_poly_t *g,
                                int k,
                                ff_poly_t ***out_factors,
                                size_t *out_count)
{
    if (!g || !out_factors || !out_count || k <= 0)
        return false;

    *out_factors = NULL;
    *out_count = 0;

    int n = ff_poly_deg(g);
    if (n <= 0 || n % k != 0)
        return false;

    ff_poly_t **factors = NULL;
    size_t count = 0;
    size_t cap = 0;

    if (!edf_recurse(g, k, &factors, &count, &cap))
    {
        factor_list_free(factors, count);
        return false;
    }

    *out_factors = factors;
    *out_count = count;
    return true;
}

/* -------------------------------------------------------------------- */
/* Root finding                                                          */
/* -------------------------------------------------------------------- */

uint32_t *ff_poly_roots_base(const ff_poly_t *f, size_t *count)
{
    if (count)
        *count = 0;
    if (!f || !count || pkf_poly_is_zero(f))
        return NULL;

    const ff_prime_t *ctx = f->ctx;
    uint32_t p = ff_prime_modulus(ctx);

    if (f->deg == 0)
        return NULL;

    /*
     * Distinct-degree step for degree 1 only:
     *   f_1 = gcd(f, x^p - x)
     * collects exactly the linear factors (with multiplicity 1, since
     * x^p - x is square-free).
     */
    ff_poly_t *x_poly = ff_poly_new(ctx, 1);
    if (!x_poly)
        return NULL;
    ff_poly_set_coeff(x_poly, 1, ff_one_rep(ctx));

    ff_poly_t *xp = ff_poly_mod_pow_new(x_poly, p, f);
    if (!xp)
    {
        ff_poly_free(x_poly);
        return NULL;
    }

    ff_poly_t *xp_minus_x = ff_poly_sub_new(xp, x_poly);
    ff_poly_free(xp);
    ff_poly_free(x_poly);
    if (!xp_minus_x)
        return NULL;

    ff_poly_t *f1 = ff_poly_gcd(f, xp_minus_x);
    ff_poly_free(xp_minus_x);

    if (!f1 || f1->deg == 0)
    {
        ff_poly_free(f1);
        return NULL; /* no roots in the base field */
    }

    /* Split the square-free product of linears with the public EDF. */
    ff_poly_t **lin = NULL;
    size_t nlin = 0;

    if (!ff_poly_split_equal_degree(f1, 1, &lin, &nlin))
    {
        ff_poly_free(f1);
        return NULL;
    }
    ff_poly_free(f1);

    uint32_t *roots = calloc(nlin, sizeof(uint32_t));
    if (!roots)
    {
        factor_list_free(lin, nlin);
        return NULL;
    }

    for (size_t i = 0; i < nlin; ++i)
    {
        /* Monic linear factor x - r: root is -c[0]. */
        roots[i] = ff_neg(ctx, ff_poly_coeff(lin[i], 0));
    }

    factor_list_free(lin, nlin);
    *count = nlin;
    return roots;
}

/* -------------------------------------------------------------------- */
/* Berlekamp's algorithm (deterministic, small primes)                   */
/* -------------------------------------------------------------------- */

/* Build B = Q - I where column j of Q holds the coefficients of
 * x^(q j) mod f. (The transpose-invariant nullity of B counts the
 * irreducible factors; its nullspace is spanned by Frobenius fixed
 * points, which is what the splitting phase relies on.) */
static ff_mat_t *berlekamp_B(const ff_poly_t *f)
{
    const ff_prime_t *ctx = ff_poly_ctx(f);
    int n = ff_poly_deg(f);

    if (!ctx || n <= 0)
        return NULL;

    uint32_t q = ff_prime_modulus(ctx);

    /* R = x^q mod f */
    ff_poly_t *x_poly = ff_poly_new(ctx, 1);
    if (!x_poly)
        return NULL;
    ff_poly_set_coeff(x_poly, 1, ff_one_rep(ctx));

    ff_poly_t *R = ff_poly_mod_pow_new(x_poly, q, f);
    ff_poly_free(x_poly);
    if (!R)
        return NULL;

    ff_mat_t *Qm = ff_mat_new(ctx, (size_t)n, (size_t)n);
    ff_poly_t *curr = ff_poly_new(ctx, 0);
    if (!Qm || !curr)
    {
        ff_mat_free(Qm);
        ff_poly_free(curr);
        ff_poly_free(R);
        return NULL;
    }
    ff_poly_one(curr); /* i = 0: x^0 = 1 */

    for (int i = 0; i < n; ++i)
    {
        /* COLUMN i = coefficients of x^(q i) mod f */
        for (int j = 0; j < n; ++j)
        {
            uint32_t *dst = ff_mat_at(Qm, (size_t)j, (size_t)i);
            if (!dst)
            {
                ff_mat_free(Qm);
                ff_poly_free(curr);
                ff_poly_free(R);
                return NULL;
            }
            *dst = ff_poly_coeff(curr, j);
        }

        if (i < n - 1)
        {
            ff_poly_t *prod = ff_poly_mul_new(curr, R);
            ff_poly_free(curr);
            curr = prod ? ff_poly_mod_new(prod, f) : NULL;
            ff_poly_free(prod);

            if (!curr)
            {
                ff_mat_free(Qm);
                ff_poly_free(R);
                return NULL;
            }
        }
    }

    ff_poly_free(curr);
    ff_poly_free(R);

    /* B = Q - I */
    uint32_t one = ff_one_rep(ctx);
    for (int i = 0; i < n; ++i)
    {
        uint32_t *e = ff_mat_at(Qm, (size_t)i, (size_t)i);
        if (!e)
        {
            ff_mat_free(Qm);
            return NULL;
        }
        *e = ff_sub(ctx, *e, one);
    }

    return Qm;
}

/* Reconstruct a polynomial from a column of the nullspace basis. */
static ff_poly_t *poly_from_mat_col(const ff_mat_t *M, size_t col,
                                    const ff_prime_t *ctx)
{
    size_t rows = ff_mat_rows(M);

    ff_poly_t *p = ff_poly_new(ctx, (int)rows - 1);
    if (!p)
        return NULL;

    for (size_t i = 0; i < rows; ++i)
    {
        const uint32_t *src = ff_mat_atc(M, i, col);
        if (!src)
        {
            ff_poly_free(p);
            return NULL;
        }
        ff_poly_set_coeff(p, (int)i, *src);
    }

    ff_poly_trim(p);
    return p;
}

/* Core square-free Berlekamp engine. */
static bool berlekamp_squarefree(const ff_poly_t *f,
                                 ff_poly_t ***out_factors, size_t *out_count,
                                 size_t *cap)
{
    const ff_prime_t *ctx = f->ctx;
    uint32_t q = ff_prime_modulus(ctx);
    int n = ff_poly_deg(f);

    if (n <= 1)
        return false;

    ff_mat_t *B = berlekamp_B(f);
    if (!B)
        return false;

    /* Nullspace of B = Frobenius fixed points. */
    ff_mat_t *null_basis = NULL;
    int nullity = ff_mat_nullspace(&null_basis, B);
    ff_mat_free(B);

    if (nullity < 0)
    {
        ff_mat_free(null_basis);
        return false;
    }

    if (nullity <= 1)
    {
        /* Only the constants: f is irreducible. */
        ff_mat_free(null_basis);

        ff_poly_t *copy = ff_poly_copy(f);
        if (!copy)
            return false;
        if (!factor_list_push(out_factors, out_count, cap, copy))
        {
            ff_poly_free(copy);
            return false;
        }
        return true;
    }

    /* Splitting phase: queue of polynomials still to split.
     * Basis column 0 is the constant 1, so start at j = 1. */
    ff_poly_t **to_process = NULL;
    size_t proc_count = 0;
    size_t proc_cap = 0;

    ff_poly_t *f_copy = ff_poly_copy(f);
    if (!f_copy || !factor_list_push(&to_process, &proc_count, &proc_cap, f_copy))
    {
        ff_poly_free(f_copy);
        factor_list_free(to_process, 0);
        ff_mat_free(null_basis);
        return false;
    }

    bool ok = true;

    while (proc_count > 0 && ok)
    {
        ff_poly_t *g = to_process[--proc_count];

        if (ff_poly_deg(g) <= 1)
        {
            if (!factor_list_push(out_factors, out_count, cap, g))
            {
                ff_poly_free(g);
                ok = false;
            }
            continue;
        }

        bool split_found = false;

        for (int j = 1; j < nullity && !split_found; ++j)
        {
            ff_poly_t *vj = poly_from_mat_col(null_basis, (size_t)j, ctx);
            if (!vj)
                continue;

            for (uint32_t s = 0; s < q && !split_found; ++s)
            {
                /* vj - s */
                uint32_t s_rep = ff_to_mont(ctx, s);
                uint32_t vc0 = ff_poly_coeff(vj, 0);
                ff_poly_set_coeff(vj, 0, ff_sub(ctx, vc0, s_rep));
                ff_poly_trim(vj);

                ff_poly_t *factor = ff_poly_gcd(g, vj);
                if (!factor)
                    continue;

                int deg_f = ff_poly_deg(factor);

                if (deg_f > 0 && deg_f < ff_poly_deg(g))
                {
                    ff_poly_t *cofactor = ff_poly_div_new(g, factor);
                    if (!cofactor)
                    {
                        ff_poly_free(factor);
                        ok = false;
                        break;
                    }

                    if (!factor_list_push(&to_process, &proc_count, &proc_cap, factor) ||
                        !factor_list_push(&to_process, &proc_count, &proc_cap, cofactor))
                    {
                        ff_poly_free(factor);
                        ff_poly_free(cofactor);
                        ok = false;
                        break;
                    }

                    ff_poly_free(g);
                    split_found = true;
                }
                else
                {
                    ff_poly_free(factor);
                }
            }

            ff_poly_free(vj);
        }

        if (!ok)
            break;

        if (!split_found)
        {
            /* No basis vector splits g: g is irreducible. */
            if (!factor_list_push(out_factors, out_count, cap, g))
            {
                ff_poly_free(g);
                ok = false;
            }
        }
    }

    /* Drain anything still queued on failure. */
    if (!ok)
    {
        for (size_t i = 0; i < proc_count; ++i)
            ff_poly_free(to_process[i]);
    }

    free(to_process);
    ff_mat_free(null_basis);
    return ok;
}

/* Public wrapper with square-free pre-check. Requires f square-free and
 * of degree >= 2; otherwise returns NULL with *out_count = 0. */
ff_poly_t **ff_poly_factor_berlekamp(const ff_poly_t *f, size_t *out_count)
{
    if (out_count)
        *out_count = 0;
    if (!f || !out_count)
        return NULL;

    /* Berlekamp needs a square-free input: check gcd(f, f'). */
    ff_poly_t *deriv = ff_poly_deriv(f);
    ff_poly_t *g = deriv ? ff_poly_gcd(f, deriv) : NULL;
    bool squarefree = g && ff_poly_deg(g) == 0;
    ff_poly_free(g);
    ff_poly_free(deriv);

    if (!squarefree)
        return NULL;

    ff_poly_t **factors = NULL;
    size_t count = 0;
    size_t cap = 0;

    if (!berlekamp_squarefree(f, &factors, &count, &cap))
    {
        factor_list_free(factors, count);
        return NULL;
    }

    *out_count = count;
    return factors;
}
