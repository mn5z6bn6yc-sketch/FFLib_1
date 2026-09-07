/*
    demo.c — Clean, modular demonstration of the finite‑field library.

    This version breaks the demo into logical sections:
      1. Base field setup
      2. Rectangular matrix + nullspace
      3. Characteristic polynomial
      4. Krylov + Frobenius
      5. Eigenvalue discovery + lifting
      6. Smith Normal Form
      7. Extension‑tower auto‑Jordan wrapper
      8. Big‑integer CRT + rational reconstruction
      9. Berlekamp factorisation

    Each section is isolated into its own function.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#include "fflib.h"
#include "ff_decomp.h"
#include "ff_poly.h"
#include "ff_snf.h"
#include "ff_extk.h"
#include "ff_poly_factor_equal_degree.h"
#include "ff_bigint.h"

#define PRIME 1000000007u

/* Forward declarations */

static void print_poly(const char *name, const ff_poly_t *p, const ff_prime_t *ctx);

static void auto_jordan_tower_wrapper(const ff_mat_t *A, const ff_prime_t *Fp);
static void process_irreducible_factor(const ff_poly_t *g, int k,
                                       const ff_mat_t *A, const ff_prime_t *Fp);
static void demo_circulant_cyclotomic(const ff_prime_t *Fp);

static void print_signed(const char *name, const ff_poly_t *p, const ff_prime_t *ctx)
{
    uint32_t mod = ff_prime_modulus(ctx);
    printf("%s (deg %d):  ", name, ff_poly_deg(p));
    for (int i = ff_poly_deg(p); i >= 0; --i)
    {
        uint32_t c = ff_from_mont(ctx, ff_poly_coeff(p, i));
        int32_t signed_c = (c > mod / 2) ? (int32_t)c - (int32_t)mod : (int32_t)c;

        if (i == ff_poly_deg(p))
        {
            if (signed_c == 1 && i > 0)
                printf("x^%d", i);
            else if (signed_c == -1 && i > 0)
                printf("-x^%d", i);
            else
                printf("%d*x^%d", signed_c, i);
        }
        else if (i == 0)
        {
            printf(" + %d", signed_c);
        }
        else
        {
            if (signed_c == 1)
                printf(" + x^%d", i);
            else if (signed_c == -1)
                printf(" - x^%d", i);
            else if (signed_c > 0)
                printf(" + %d*x^%d", signed_c, i);
            else
                printf(" - %d*x^%d", -signed_c, i);
        }
    }
    printf("\n\n");
}

/* -------------------------------------------------------------
[12] Rational Jordan Form (Base-Field Block Decomposition)
------------------------------------------------------------- */

/* Helper: Fill a k x k Companion Matrix block for q(x) */
static void fill_companion_block(ff_mat_t *R, size_t start_row, size_t start_col, const ff_poly_t *q, const ff_prime_t *Fp)
{
    int k = ff_poly_deg(q);
    /* Subdiagonal 1s */
    for (int i = 1; i < k; ++i)
    {
        *ff_mat_at(R, start_row + i, start_col + i - 1) = ff_one_rep(Fp);
    }
    /* Last column: -c_i */
    for (int i = 0; i < k; ++i)
    {
        uint32_t c = ff_poly_coeff(q, i);
        *ff_mat_at(R, start_row + i, start_col + k - 1) = ff_neg(Fp, c);
    }
}

/* Helper: Fill an (m*k) x (m*k) Rational Jordan Block for q(x)^m */
static void fill_rational_jordan_block(ff_mat_t *R, size_t start, const ff_poly_t *q, int m, const ff_prime_t *Fp)
{
    int k = ff_poly_deg(q);
    for (int i = 0; i < m; ++i)
    {
        /* 1. Place Companion Matrix C(q) on the diagonal */
        fill_companion_block(R, start + i * k, start + i * k, q, Fp);

        /* 2. Place Identity Matrix I_k on the superdiagonal (if m > 1) */
        if (i < m - 1)
        {
            for (int j = 0; j < k; ++j)
            {
                *ff_mat_at(R, start + i * k + j, start + (i + 1) * k + j) = ff_one_rep(Fp);
            }
        }
    }
}

static void demo_rational_jordan(const ff_prime_t *Fp)
{
    printf("--- [12] Rational Jordan Form (Base-Field Block Decomposition) ---\n");

    /* Re-use the 5x5 matrix from demo_tower */
    ff_mat_t *A = ff_mat_new(Fp, 5, 5);
    *ff_mat_at(A, 0, 0) = 0;
    *ff_mat_at(A, 0, 1) = ff_neg(Fp, ff_one_rep(Fp));
    *ff_mat_at(A, 1, 0) = ff_one_rep(Fp);
    *ff_mat_at(A, 1, 1) = 0;
    *ff_mat_at(A, 2, 4) = ff_neg(Fp, ff_one_rep(Fp));
    *ff_mat_at(A, 3, 2) = ff_one_rep(Fp);
    *ff_mat_at(A, 3, 4) = ff_neg(Fp, ff_one_rep(Fp));
    *ff_mat_at(A, 4, 3) = ff_one_rep(Fp);

    /* 1. Compute Characteristic Polynomial */
    ff_poly_t *charpoly = ff_poly_charpoly_hessenberg_base(A);
    int n = ff_poly_deg(charpoly);
    printf("Computed Characteristic Polynomial of degree %d.\n", n);

    /* 2. The Square-Free Trick: Strip multiplicities to find DISTINCT irreducible factors */
    ff_poly_t *deriv = ff_poly_deriv(charpoly);
    ff_poly_t *g = ff_poly_gcd(charpoly, deriv);
    ff_poly_t *sq_free = ff_poly_div_new(charpoly, g);
    ff_poly_free(deriv);
    ff_poly_free(g);

    /* 3. Run DDF/EDF on the square-free part to collect distinct irreducible factors */
    /* 3. Run DDF/EDF on the square-free part to collect distinct irreducible factors */
    typedef struct
    {
        ff_poly_t *q;
        int m;
    } RationalBlock;

    ff_poly_t *rem = ff_poly_copy(sq_free);
    ff_poly_t *h = ff_poly_new(Fp, 1);
    ff_poly_set_coeff(h, 1, ff_one_rep(Fp));
    uint32_t p = ff_prime_modulus(Fp);

    size_t cap = 4, count = 0;
    RationalBlock *blocks = calloc(cap, sizeof(RationalBlock));

    for (int k = 1; k <= n / 2; ++k)
    {
        /* FIX: Early exit if rem is fully factored (degree 0) */
        if (ff_poly_deg(rem) <= 0)
            break;

        ff_poly_t *h_p = ff_poly_mod_pow_new(h, p, rem);
        ff_poly_free(h);
        h = h_p;

        ff_poly_t *x_poly = ff_poly_new(Fp, 1);
        ff_poly_set_coeff(x_poly, 1, ff_one_rep(Fp));
        ff_poly_t *diff = ff_poly_sub_new(h, x_poly);
        ff_poly_free(x_poly);

        ff_poly_t *g_k = ff_poly_gcd(diff, rem);
        ff_poly_free(diff);

        if (g_k && ff_poly_deg(g_k) > 0)
        {
            /* FIX: Divide rem by g_k BEFORE potentially freeing g_k */
            ff_poly_t *quot = ff_poly_div_new(rem, g_k);
            ff_poly_free(rem);
            rem = quot;

            if (ff_poly_deg(g_k) == k)
            {
                if (count == cap)
                {
                    cap *= 2;
                    blocks = realloc(blocks, cap * sizeof(RationalBlock));
                }
                blocks[count].q = g_k;
                blocks[count].m = 0;
                count++;
            }
            else
            {
                ff_poly_t **edf = NULL;
                size_t edf_count = 0;
                if (ff_poly_split_equal_degree(g_k, k, &edf, &edf_count))
                {
                    for (size_t i = 0; i < edf_count; ++i)
                    {
                        if (count == cap)
                        {
                            cap *= 2;
                            blocks = realloc(blocks, cap * sizeof(RationalBlock));
                        }
                        blocks[count].q = edf[i];
                        blocks[count].m = 0;
                        count++;
                    }
                    free(edf);
                }
                ff_poly_free(g_k); /* Now safe to free */
            }
        }
        else
        {
            if (g_k)
                ff_poly_free(g_k);
        }
    }

    /* 4. Count multiplicities in the ORIGINAL charpoly */
    ff_poly_t *curr = ff_poly_copy(charpoly);
    for (size_t i = 0; i < count; ++i)
    {
        int m = 0;
        while (ff_poly_deg(curr) >= ff_poly_deg(blocks[i].q))
        {
            ff_poly_t *rem_check = ff_poly_mod_new(curr, blocks[i].q);
            /* If remainder is exactly 0, it divides perfectly */
            bool is_zero = (ff_poly_deg(rem_check) < 0) || (ff_poly_deg(rem_check) == 0 && ff_poly_coeff(rem_check, 0) == 0);
            ff_poly_free(rem_check);

            if (is_zero)
            {
                ff_poly_t *quot = ff_poly_div_new(curr, blocks[i].q);
                ff_poly_free(curr);
                curr = quot;
                m++;
            }
            else
            {
                break;
            }
        }
        blocks[i].m = m;
    }
    ff_poly_free(curr);
    ff_poly_free(charpoly);

    /* 5. Build the Rational Jordan Matrix R */
    ff_mat_t *R = ff_mat_new(Fp, n, n);
    size_t offset = 0;

    printf("\nConstructing Rational Jordan Form over F_p (No field extensions!):\n");
    for (size_t i = 0; i < count; ++i)
    {
        int k = ff_poly_deg(blocks[i].q);
        int m = blocks[i].m;
        int block_size = k * m;

        printf("  -> Irreducible factor of degree %d, multiplicity %d. Building %dx%d block.\n", k, m, block_size, block_size);
        fill_rational_jordan_block(R, offset, blocks[i].q, m, Fp);
        offset += block_size;

        ff_poly_free(blocks[i].q);
    }
    free(blocks);

    printf("\nRational Jordan Matrix R (%zux%zu):\n", n, n);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            uint32_t val = ff_from_mont(Fp, *ff_mat_atc(R, i, j));
            /* Print signed residues for readability */
            int32_t signed_val = (val > PRIME / 2) ? (int32_t)val - (int32_t)PRIME : (int32_t)val;
            printf("%4d ", signed_val);
        }
        printf("\n");
    }
    printf("\n");

    ff_mat_free(R);
    ff_mat_free(A);
}

/* -------------------------------------------------------------
   Helper: process a single irreducible factor of degree k
   (build F_{p^k}, report success; hook point for Jordan engine)
------------------------------------------------------------- */

void process_irreducible_factor(const ff_poly_t *irred, int k, const ff_mat_t *A, const ff_prime_t *Fp)
{
    printf("         -> Irreducible factor of degree %d:\n", k);
    ff_extk_t *Fpk = ff_extk_new_with_poly(Fp, irred);
    /* the irreducible factor */
    print_poly("            Factor", irred, Fp);

    printf("            Constructing F_{p^%d}...\n", k);
    if (!Fpk)
    {
        printf("         -> FAILED to construct extension field.\n");
        return;
    }
    else
    {
        printf("         -> Irreducible factor of degree %d. Constructing F_{p^%d}...\n", k, k);

        if (!Fpk)
            return;

        /* 1. Lift Base Matrix A to F_{p^k} */
        size_t N = ff_mat_rows(A);
        ff_mat_extk_t *A_extk = ff_mat_extk_new(Fpk, N, N);
        for (size_t i = 0; i < N; ++i)
        {
            for (size_t j = 0; j < N; ++j)
            {
                uint32_t val = ff_from_mont(Fp, *ff_mat_atc(A, i, j));
                ff_poly_t *elem = ff_mat_extk_at(A_extk, i, j);
                ff_poly_set_coeff(elem, 0, ff_to_mont(Fp, val));
            }
        }

        /* 2. The Eigenvalue is the formal variable 'x' */
        ff_poly_t *formal_x = ff_poly_new(Fp, k - 1);
        ff_poly_set_coeff(formal_x, 1, ff_one_rep(Fp));

        const ff_poly_t *eigs[2];
        size_t eig_count = 1;
        eigs[0] = formal_x;

        ff_poly_t *conjugate = NULL;

        /* 3. If k=2, compute the Galois conjugate using Vieta's formulas: beta = -b - alpha */
        if (k == 2)
        {
            uint32_t b = ff_poly_coeff(irred, 1);
            uint32_t neg_b = ff_neg(Fp, b);

            conjugate = ff_poly_new(Fp, 1);
            ff_poly_set_coeff(conjugate, 0, neg_b);
            ff_poly_set_coeff(conjugate, 1, ff_neg(Fp, ff_one_rep(Fp))); /* -x */

            eigs[1] = conjugate;
            eig_count = 2;
            printf("         -> Computed Galois conjugate root for quadratic.\n");
        }

        /* 4. Run the Generic Jordan Engine */
        printf("         -> Running Generic Jordan Engine over F_{p^%d}...\n", k);
        ff_jordan_extk_result_t j_res;
        if (ff_jordan_extk_with_eigenvalues(Fpk, A_extk, eigs, eig_count, &j_res))
        {
            printf("         -> SUCCESS! Found %zu Jordan blocks.\n", j_res.block_count);
            for (size_t i = 0; i < j_res.block_count; ++i)
            {
                printf("            Block %zu: size %zu, eigenvalue = ", i, j_res.block_sizes[i]);

                /* Print the eigenvalue as a polynomial in x */
                /*int deg = ff_poly_deg(j_res.block_eigenvalues[i]);
                for (int d = deg; d >= 0; --d)
                {
                    uint32_t c = ff_from_mont(Fp, ff_poly_coeff(j_res.block_eigenvalues[i], d));
                    if (d == deg)
                        printf("%u*x^%d", c, d);
                    else if (d == 0)
                        printf(" + %u", c);
                    else
                        printf(" + %u*x^%d", c, d);
                }*/
                int deg = ff_poly_deg(j_res.block_eigenvalues[i]);
                for (int d = deg; d >= 0; --d)
                {
                    uint32_t c = ff_from_mont(Fp, ff_poly_coeff(j_res.block_eigenvalues[i], d));
                    int32_t signed_c = (c > PRIME / 2) ? (int32_t)c - (int32_t)PRIME : (int32_t)c;

                    if (d == deg)
                    {
                        if (signed_c == 1 && d > 0)
                            printf("x^%d", d);
                        else if (signed_c == -1 && d > 0)
                            printf("-x^%d", d);
                        else
                            printf("%d*x^%d", signed_c, d);
                    }
                    else if (d == 0)
                    {
                        printf(" + %d", signed_c);
                    }
                    else
                    {
                        if (signed_c == 1)
                            printf(" + x^%d", d);
                        else if (signed_c == -1)
                            printf(" - x^%d", d);
                        else if (signed_c > 0)
                            printf(" + %d*x^%d", signed_c, d);
                        else
                            printf(" - %d*x^%d", -signed_c, d);
                    }
                }
                printf("\n");
            }
            ff_jordan_extk_result_clear(&j_res);
        }
        else
        {
            printf("         -> FAILED! Jordan engine could not decompose.\n");
        }

        if (conjugate)
            ff_poly_free(conjugate);
        ff_poly_free(formal_x);
        ff_mat_extk_free(A_extk);
        ff_extk_free(Fpk);
    }
}

void process_irreducible_factor1(const ff_poly_t *irred, int k, const ff_mat_t *A, const ff_prime_t *Fp)
{
    printf("         -> Irreducible factor of degree %d. Constructing F_{p^%d}...\n", k, k);

    ff_extk_t *Fpk = ff_extk_new_with_poly(Fp, irred);
    if (!Fpk)
    {
        printf("         -> FAILED to construct extension field F_{p^%d}.\n", k);
        return;
    }

    /* For a demo, we don't need the full Jordan engine on the big matrix.
       We just need to show that the polynomial splits in F_{p^k}. */
    ff_extk_elem_t *root = ff_extk_elem_new(Fpk);
    /* root is already zero-initialized by ff_poly_new.
       DO NOT call ff_extk_zero, as it sets deg=0 and breaks ff_poly_set_coeff! */
    if (k > 1)
    {
        ff_poly_set_coeff(root, 1, ff_one_rep(Fp));
    }

    ff_extk_elem_t *result = ff_extk_elem_new(Fpk);
    /* result is already zero-initialized. */

    int deg = ff_poly_deg(irred);
    ff_extk_elem_t *temp = ff_extk_elem_new(Fpk);

    for (int i = deg; i >= 0; --i)
    {
        ff_extk_mul(Fpk, result, result, root);

        /* Clear temp safely without destroying its degree capacity */
        for (int j = 0; j <= ff_poly_deg(temp); ++j)
        {
            ff_poly_set_coeff(temp, j, 0);
        }
        ff_poly_set_coeff(temp, 0, ff_poly_coeff(irred, i));

        ff_extk_add(Fpk, result, result, temp);
    }
    ff_extk_elem_free(temp);
    if (ff_extk_is_zero(result))
    {
        printf("         -> VERIFIED: The polynomial splits in F_{p^%d}!\n", k);
        printf("            The matrix is diagonalizable over this field.\n");
    }
    else
    {
        printf("         -> ERROR: Polynomial did not split.\n");
    }

    ff_extk_elem_free(root);
    ff_extk_elem_free(result);
    ff_extk_free(Fpk);
}

/* -------------------------------------------------------------
   [8] Ultimate auto‑discovery & extension towering wrapper
------------------------------------------------------------- */
static void auto_jordan_tower_wrapper(const ff_mat_t *A, const ff_prime_t *Fp)
{
    printf("\n--- [8] Ultimate Auto-Discovery & Extension Towering Wrapper ---\n");

    /* 1. Characteristic polynomial */
    ff_poly_t *charpoly = ff_poly_charpoly_hessenberg_base(A);
    if (!charpoly)
    {
        printf("Failed to compute characteristic polynomial.\n");
        return;
    }

    int n = ff_poly_deg(charpoly);
    printf("Computed Characteristic Polynomial of degree %d.\n", n);
    print_poly("  Char Poly", charpoly, Fp);
    /* 2. Distinct-Degree Factorization (DDF) */
    ff_poly_t *rem = ff_poly_copy(charpoly);
    ff_poly_t *h = ff_poly_new(Fp, 1);
    ff_poly_set_coeff(h, 1, ff_one_rep(Fp)); /* h = x */

    uint32_t p = ff_prime_modulus(Fp);

    for (int k = 1; k <= n / 2; ++k)
    {
        /* h = h^p mod rem */
        ff_poly_t *h_p = ff_poly_mod_pow_new(h, p, rem);
        ff_poly_free(h);
        h = h_p;

        /* diff = h - x */
        ff_poly_t *x_poly = ff_poly_new(Fp, 1);
        ff_poly_set_coeff(x_poly, 1, ff_one_rep(Fp));
        ff_poly_t *diff = ff_poly_sub_new(h, x_poly);
        ff_poly_free(x_poly);

        /* g = gcd(diff, rem) */
        ff_poly_t *g = ff_poly_gcd(diff, rem);
        ff_poly_free(diff);

        if (!g)
        {
            printf("GCD failed.\n");
            break;
        }

        int deg_g = ff_poly_deg(g);

        if (deg_g > 0)
        {
            printf("\n  [Degree %d] Found factor(s) of degree %d.\n", k, k);

            if (deg_g == k)
            {
                /* Single irreducible factor */
                if (k == 1)
                {
                    uint32_t root = ff_neg(Fp, ff_poly_coeff(g, 0));
                    printf("    -> Root in F_p: %u\n", ff_from_mont(Fp, root));
                }
                else
                {
                    process_irreducible_factor(g, k, A, Fp);
                }
            }
            else
            {
                /* Equal-degree factorization */
                printf("    -> Multiple factors of degree %d found (deg=%d). Running Cantor-Zassenhaus EDF...\n",
                       k, deg_g);

                ff_poly_t **edf_factors = NULL;
                size_t edf_count = 0;

                if (ff_poly_split_equal_degree(g, k, &edf_factors, &edf_count))
                {
                    printf("       Successfully split into %zu irreducible factors!\n", edf_count);

                    for (size_t idx = 0; idx < edf_count; ++idx)
                    {
                        printf("       Processing irreducible factor %zu/%zu...\n", idx + 1, edf_count);
                        process_irreducible_factor(edf_factors[idx], k, A, Fp);
                        ff_poly_free(edf_factors[idx]);
                    }
                    free(edf_factors);
                }
                else
                {
                    printf("       EDF failed to split the factors.\n");
                }
            }

            /* rem = rem / g */
            ff_poly_t *quot = ff_poly_div_new(rem, g);
            ff_poly_free(rem);
            rem = quot;
        }

        ff_poly_free(g);
    }

    /* Remaining factor */
    if (ff_poly_deg(rem) > 0)
    {
        int k = ff_poly_deg(rem);
        printf("\n  [Degree %d] Found remaining irreducible factor of degree %d.\n", k, k);

        ff_extk_t *Fpk = ff_extk_new_with_poly(Fp, rem);
        if (Fpk)
        {
            printf("    -> SUCCESS! Constructed F_{p^%d} for the remaining factor.\n", k);
            ff_extk_free(Fpk);
        }
        else
        {
            printf("    -> FAILED to construct extension field.\n");
        }
    }

    ff_poly_free(h);
    ff_poly_free(rem);
    ff_poly_free(charpoly);
    printf("\n");
}

/* -------------------------------------------------------------
   Utility printing helpers
------------------------------------------------------------- */

static void print_base_mat(const char *name, const ff_mat_t *M, const ff_prime_t *ctx)
{
    printf("%s (%zux%zu):\n", name, ff_mat_rows(M), ff_mat_cols(M));
    for (size_t i = 0; i < ff_mat_rows(M); ++i)
    {
        for (size_t j = 0; j < ff_mat_cols(M); ++j)
        {
            uint32_t val = ff_from_mont(ctx, *ff_mat_atc(M, i, j));
            printf("%10u ", val);
        }
        printf("\n");
    }
    printf("\n");
}

static void print_poly(const char *name, const ff_poly_t *p, const ff_prime_t *ctx)
{
    printf("%s (deg %d): ", name, ff_poly_deg(p));
    for (int i = ff_poly_deg(p); i >= 0; --i)
    {
        uint32_t c = ff_from_mont(ctx, ff_poly_coeff(p, i));
        if (i == ff_poly_deg(p))
            printf("%u*x^%d", c, i);
        else if (i == 0)
            printf(" + %u", c);
        else
            printf(" + %u*x^%d", c, i);
    }
    printf("\n\n");
}

/* -------------------------------------------------------------
   [2] Rectangular matrix + nullspace
------------------------------------------------------------- */
static void demo_nullspace(const ff_prime_t *Fp)
{
    printf("--- [2] Rectangular Matrix & Nullspace ---\n");

    ff_mat_t *M = ff_mat_new(Fp, 8, 12);
    uint32_t one = ff_one_rep(Fp);
    uint32_t two = ff_add(Fp, one, one);
    uint32_t three = ff_add(Fp, two, one);

    for (size_t i = 0; i < 8; ++i)
        for (size_t j = 0; j < 12; ++j)
            *ff_mat_at(M, i, j) =
                (j == i) ? one : (j == i + 1) ? two
                             : (j == i + 2)   ? three
                                              : 0;

    print_base_mat("Rectangular Matrix M", M, Fp);

    ff_mat_t *N = NULL;
    int nullity = ff_mat_nullspace(&N, M);
    printf("Nullity of M: %d\n", nullity);

    if (N)
    {
        print_base_mat("Nullspace Basis", N, Fp);
        ff_mat_free(N);
    }

    ff_mat_free(M);
    printf("\n");
}

/* -------------------------------------------------------------
   [3] Characteristic polynomial
------------------------------------------------------------- */
static void demo_charpoly(const ff_prime_t *Fp)
{
    printf("--- [3] Square Matrix & Characteristic Polynomial ---\n");

    size_t n = 10;
    ff_mat_t *A = ff_mat_new(Fp, n, n);

    srand(42);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            *ff_mat_at(A, i, j) = ff_to_mont(Fp, rand() % PRIME);

    clock_t t1 = clock();
    ff_poly_t *hess = ff_poly_charpoly_hessenberg_base(A);
    clock_t t2 = clock();

    print_poly("Char Poly (Hessenberg)", hess, Fp);
    printf("Time: %.6f sec\n\n", (double)(t2 - t1) / CLOCKS_PER_SEC);

    ff_poly_free(hess);
    ff_mat_free(A);
}

/* -------------------------------------------------------------
   [4] Krylov + Frobenius
------------------------------------------------------------- */
static void demo_frobenius(const ff_prime_t *Fp)
{
    printf("--- [4] Krylov & Frobenius Decomposition ---\n");

    size_t n = 10;
    ff_mat_t *A = ff_mat_new(Fp, n, n);

    srand(42);
    for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
            *ff_mat_at(A, i, j) = ff_to_mont(Fp, rand() % PRIME);

    uint32_t *v = calloc(n, sizeof(uint32_t));
    v[0] = ff_one_rep(Fp);

    ff_poly_t *poly = NULL;
    ff_krylov_dependency_base(A, v, &poly);
    print_poly("Krylov Annihilating Poly", poly, Fp);

    ff_mat_t *P = NULL, *F = NULL;
    if (ff_frobenius_greedy_base(A, &P, &F))
    {
        printf("Greedy Frobenius decomposition successful.\n");
        print_base_mat("P", P, Fp);
        print_base_mat("F", F, Fp);
        ff_mat_free(P);
        ff_mat_free(F);
    }

    ff_poly_free(poly);
    free(v);
    ff_mat_free(A);
    printf("\n");
}

/* -------------------------------------------------------------
   [5] Eigenvalue discovery + lifting
------------------------------------------------------------- */
static void demo_eigen_lifting(const ff_prime_t *Fp)
{
    printf("--- [5] Automatic Eigenvalue Discovery & Lifting ---\n");

    ff_mat_t *A = ff_mat_new(Fp, 2, 2);
    *ff_mat_at(A, 0, 0) = 0;
    *ff_mat_at(A, 0, 1) = ff_neg(Fp, ff_one_rep(Fp));
    *ff_mat_at(A, 1, 0) = ff_one_rep(Fp);
    *ff_mat_at(A, 1, 1) = 0;

    ff_poly_t *charpoly = ff_poly_charpoly_hessenberg_base(A);
    print_poly("Char Poly", charpoly, Fp);

    size_t root_count = 0;
    uint32_t *roots = ff_poly_roots_base(charpoly, &root_count);

    if (!roots || root_count == 0)
    {
        printf("No roots in F_p. Lifting to F_{p^2}...\n");

        uint32_t b = ff_poly_coeff(charpoly, 1);
        uint32_t c = ff_poly_coeff(charpoly, 0);

        uint32_t two = ff_to_mont(Fp, 2);
        uint32_t four = ff_to_mont(Fp, 4);

        uint32_t Delta = ff_sub(Fp, ff_mul(Fp, b, b), ff_mul(Fp, four, c));
        printf("Discriminant Delta = %u\n", ff_from_mont(Fp, Delta));

        ff_ext2_t *Fp2 = ff_ext2_new(Fp, ff_from_mont(Fp, Delta));

        uint32_t inv2 = ff_inv(Fp, two);
        uint32_t real = ff_mul(Fp, ff_neg(Fp, b), inv2);
        uint32_t imag = ff_mul(Fp, ff_one_rep(Fp), inv2);

        ff_ext2_elem_t eigs[2] = {
            {real, imag},
            {real, ff_neg(Fp, imag)}};

        ff_jordan_result_t res;
        if (ff_jordan_base_with_eigenvalues_ext2(Fp2, A, eigs, 2, &res))
        {
            printf("Jordan decomposition successful.\n");
            for (size_t i = 0; i < res.block_count; ++i)
            {
                printf("  Block %zu: size %zu, eigenvalue = %u + %u*w\n",
                       i,
                       res.block_sizes[i],
                       ff_from_mont(Fp, res.block_eigenvalues[i].a),
                       ff_from_mont(Fp, res.block_eigenvalues[i].b));
            }
            ff_jordan_result_clear(&res);
        }

        ff_ext2_free(Fp2);
    }

    free(roots);
    ff_poly_free(charpoly);
    ff_mat_free(A);
    printf("\n");
}

/* -------------------------------------------------------------
   [6] Smith Normal Form
------------------------------------------------------------- */
static void demo_snf(const ff_prime_t *Fp)
{
    printf("--- [6] Smith Normal Form ---\n");

    ff_mat_t *A = ff_mat_new(Fp, 4, 4);
    srand(99);
    for (size_t i = 0; i < 4; ++i)
        for (size_t j = 0; j < 4; ++j)
            *ff_mat_at(A, i, j) = ff_to_mont(Fp, rand() % PRIME);

    ff_poly_t **inv = NULL;
    size_t count = 0;

    if (ff_snf_charmatrix_invariants(A, &inv, &count, FF_INVARIANT_DESCENDING))
    {
        printf("Found %zu invariant factors:\n", count);
        for (size_t i = 0; i < count; ++i)
        {
            print_poly("Invariant", inv[i], Fp);
            ff_poly_free(inv[i]);
        }
        free(inv);
    }

    ff_mat_free(A);
    printf("\n");
}

/* -------------------------------------------------------------
   [8] Extension‑tower auto‑Jordan wrapper demo
------------------------------------------------------------- */
static void demo_tower(const ff_prime_t *Fp)
{
    printf("--- [8] Auto‑Discovery & Extension Towering ---\n");

    ff_mat_t *A = ff_mat_new(Fp, 5, 5);

    /* 2×2 rotation block */
    *ff_mat_at(A, 0, 0) = 0;
    *ff_mat_at(A, 0, 1) = ff_neg(Fp, ff_one_rep(Fp));
    *ff_mat_at(A, 1, 0) = ff_one_rep(Fp);
    *ff_mat_at(A, 1, 1) = 0;

    /* 3×3 companion block for x^3 + x + 1 */
    *ff_mat_at(A, 2, 4) = ff_neg(Fp, ff_one_rep(Fp));
    *ff_mat_at(A, 3, 2) = ff_one_rep(Fp);
    *ff_mat_at(A, 3, 4) = ff_neg(Fp, ff_one_rep(Fp));
    *ff_mat_at(A, 4, 3) = ff_one_rep(Fp);

    auto_jordan_tower_wrapper(A, Fp);
    ff_mat_free(A);
}

/* -------------------------------------------------------------
   [9] Big‑integer CRT + rational reconstruction
------------------------------------------------------------- */
static void demo_bigint(void)
{
    printf("--- [9] Big Integer CRT & Rational Reconstruction ---\n");

    uint32_t primes[] = {1000000007, 1000000009, 999999937};
    uint32_t residues[3];

    for (int i = 0; i < 3; ++i)
    {
        ff_prime_t *Fp = ff_prime_new(primes[i], FF_BACKEND_MONTGOMERY);
        uint32_t num = ff_to_mont(Fp, 123456789 % primes[i]);
        uint32_t den = ff_to_mont(Fp, 987654321 % primes[i]);
        residues[i] = ff_from_mont(Fp, ff_mul(Fp, num, ff_inv(Fp, den)));
        ff_prime_free(Fp);
    }

    ff_crt_gmp_t *crt = ff_crt_gmp_new();
    for (int i = 0; i < 3; ++i)
        ff_crt_gmp_add_residue(crt, primes[i], residues[i]);

    ff_bigint_t *num = ff_bigint_new();
    ff_bigint_t *den = ff_bigint_new();

    if (ff_ratrec_gmp_crt(crt, num, den))
    {
        char ns[256], ds[256];
        ff_bigint_get_str(ns, NULL, 10, num);
        ff_bigint_get_str(ds, NULL, 10, den);
        printf("Exact Rational Result: %s / %s\n", ns, ds);
    }

    ff_bigint_free(num);
    ff_bigint_free(den);
    ff_crt_gmp_free(crt);
    printf("\n");
}

/* -------------------------------------------------------------
   [10] Berlekamp factorisation
------------------------------------------------------------- */
static void demo_berlekamp(void)
{
    printf("--- [10] Berlekamp Factorisation ---\n");

    ff_prime_t *Fp7 = ff_prime_new(7, FF_BACKEND_NATIVE);

    ff_poly_t *f = ff_poly_new(Fp7, 3);
    ff_poly_set_coeff(f, 3, ff_one_rep(Fp7));
    ff_poly_set_coeff(f, 1, ff_to_mont(Fp7, 2));
    ff_poly_set_coeff(f, 0, ff_to_mont(Fp7, 4));

    size_t count = 0;
    ff_poly_t **fac = ff_poly_factor_berlekamp(f, &count);

    printf("Found %zu factors:\n", count);
    for (size_t i = 0; i < count; ++i)
    {
        print_poly("Factor", fac[i], Fp7);
        ff_poly_free(fac[i]);
    }

    free(fac);
    ff_poly_free(f);
    ff_prime_free(Fp7);
    printf("\n");
}

/* -------------------------------------------------------------
[11] Circulant Matrix & Cyclotomic Splitting
------------------------------------------------------------- */
static void demo_circulant_cyclotomic(const ff_prime_t *Fp)
{
    printf("--- [11] Circulant Matrix & Cyclotomic Splitting ---\n");
    size_t n = 9;
    ff_mat_t *C = ff_mat_new(Fp, n, n);

    /* Build a 7x7 Circulant Matrix with first row: 1, 2, 3, 4, 5, 6, 7 */
    uint32_t row[9];
    for (int i = 0; i < n; ++i)
    {
        row[i] = ff_to_mont(Fp, i + 1);
    }

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            /* C[i][j] = row[(j - i) mod n] */
            int idx = ((int)j - (int)i) % (int)n;
            if (idx < 0)
                idx += n;
            *ff_mat_at(C, i, j) = row[idx];
        }
    }

    printf("Constructed 7x7 Circulant Matrix (first row: 1..7).\n");
    printf("Eigenvalues are p(omega^j) where p(x) = 1 + 2x + ... + 7x^6 and omega^7 = 1.\n");
    printf("Since p = 10^9+7 ≡ -1 (mod 7), the 7th roots of unity require F_{p^2}.\n\n");

    /* Run the ultimate auto-discovery wrapper */
    auto_jordan_tower_wrapper(C, Fp);

    ff_mat_free(C);
    printf("\n");
}

/* -------------------------------------------------------------
   Main driver
------------------------------------------------------------- */
int main(void)
{
    printf("======================================================\n");
    printf("=== Finite Field Library Demonstration             ===\n");
    printf("======================================================\n\n");

    ff_prime_t *Fp = ff_prime_new(PRIME, FF_BACKEND_MONTGOMERY);

    demo_nullspace(Fp);
    demo_charpoly(Fp);
    demo_frobenius(Fp);
    demo_eigen_lifting(Fp);
    demo_snf(Fp);
    demo_tower(Fp);
    demo_rational_jordan(Fp);
    demo_bigint();
    demo_berlekamp();

    demo_circulant_cyclotomic(Fp);

    ff_prime_free(Fp);

    printf("======================================================\n");
    printf("=== Demo Complete                                  ===\n");
    printf("======================================================\n");
    return 0;
}