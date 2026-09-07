#ifndef FFLIB_H
#define FFLIB_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#if defined(__GNUC__) || defined(__clang__)
typedef unsigned __int128 ff_u128_t;
#else
typedef unsigned __int128 ff_u128_t; /* VS will need a fallback */
#endif

#ifdef __cplusplus
extern "C"
{
#endif

   /* ============================================================
      Base prime field
      ============================================================ */

   typedef enum
   {
      FF_BACKEND_NATIVE = 0,
      FF_BACKEND_MONTGOMERY = 1
   } ff_backend_t;

   typedef struct ff_prime ff_prime_t;

   /* Prime field API */
   ff_prime_t *ff_prime_new(uint32_t p, ff_backend_t backend);
   void ff_prime_free(ff_prime_t *ctx);
   uint32_t ff_prime_modulus(const ff_prime_t *ctx);
   ff_backend_t ff_prime_backend(const ff_prime_t *ctx);

   uint32_t ff_add(const ff_prime_t *ctx, uint32_t a, uint32_t b);
   uint32_t ff_sub(const ff_prime_t *ctx, uint32_t a, uint32_t b);
   uint32_t ff_mul(const ff_prime_t *ctx, uint32_t a, uint32_t b);
   uint32_t ff_neg(const ff_prime_t *ctx, uint32_t a);
   uint32_t ff_inv(const ff_prime_t *ctx, uint32_t a);
   uint32_t ff_pow(const ff_prime_t *ctx, uint32_t a, uint32_t e);

   uint32_t ff_to_mont(const ff_prime_t *ctx, uint32_t x);
   uint32_t ff_from_mont(const ff_prime_t *ctx, uint32_t x_mont);

   uint32_t ff_zero_rep(const ff_prime_t *ctx);
   uint32_t ff_one_rep(const ff_prime_t *ctx);
   uint32_t ff_minus_one_rep(const ff_prime_t *ctx);

   bool ff_sqrt(const ff_prime_t *ctx, uint32_t a_rep, uint32_t *root_rep);

   /* ============================================================
      Base-field matrices
      ============================================================ */

   typedef struct ff_mat
   {
      const ff_prime_t *ctx;
      size_t rows;
      size_t cols;
      size_t stride;
      uint32_t *data;
   } ff_mat_t;

   ff_mat_t *ff_mat_new(const ff_prime_t *ctx, size_t rows, size_t cols);
   void ff_mat_free(ff_mat_t *M);

   uint32_t *ff_mat_at(ff_mat_t *M, size_t r, size_t c);
   const uint32_t *ff_mat_atc(const ff_mat_t *M, size_t r, size_t c);

   const ff_prime_t *ff_mat_ctx(const ff_mat_t *M);
   size_t ff_mat_rows(const ff_mat_t *M);
   size_t ff_mat_cols(const ff_mat_t *M);

   void ff_mat_identity(ff_mat_t *M);
   void ff_mat_mul(ff_mat_t *C, const ff_mat_t *A, const ff_mat_t *B);
   bool ff_mat_inv(ff_mat_t *X, const ff_mat_t *A);
   int ff_mat_nullspace(ff_mat_t **basis_out, const ff_mat_t *M);
   bool ff_linear_solve_unique_base(const ff_mat_t *M, const uint32_t *rhs, uint32_t *x);

   /* ============================================================
      Quadratic extension field F_{p^2}
      ============================================================ */

   typedef struct
   {
      uint32_t a;
      uint32_t b;
   } ff_ext2_elem_t;

   typedef struct ff_ext2 ff_ext2_t;

   ff_ext2_t *ff_ext2_new(const ff_prime_t *base, uint32_t omega);
   void ff_ext2_free(ff_ext2_t *ctx);
   const ff_prime_t *ff_ext2_base(const ff_ext2_t *ctx);

   void ff_ext2_zero(const ff_ext2_t *ctx, ff_ext2_elem_t *out);
   void ff_ext2_one(const ff_ext2_t *ctx, ff_ext2_elem_t *out);
   void ff_ext2_add(const ff_ext2_t *ctx, ff_ext2_elem_t *out, const ff_ext2_elem_t *a, const ff_ext2_elem_t *b);
   void ff_ext2_sub(const ff_ext2_t *ctx, ff_ext2_elem_t *out, const ff_ext2_elem_t *a, const ff_ext2_elem_t *b);
   void ff_ext2_mul(const ff_ext2_t *ctx, ff_ext2_elem_t *out, const ff_ext2_elem_t *a, const ff_ext2_elem_t *b);
   void ff_ext2_neg(const ff_ext2_t *ctx, ff_ext2_elem_t *out, const ff_ext2_elem_t *a);
   bool ff_ext2_inv(const ff_ext2_t *ctx, ff_ext2_elem_t *out, const ff_ext2_elem_t *a);
   bool ff_ext2_is_zero(const ff_ext2_elem_t *a);
   bool ff_ext2_eq(const ff_ext2_elem_t *a, const ff_ext2_elem_t *b);

   /* ============================================================
      Extension-field matrices
      ============================================================ */

   typedef struct ff_mat_ext2
   {
      const ff_ext2_t *ctx;
      size_t rows;
      size_t cols;
      size_t stride;
      ff_ext2_elem_t *data;
   } ff_mat_ext2_t;

   ff_mat_ext2_t *ff_mat_ext2_new(const ff_ext2_t *ctx, size_t rows, size_t cols);
   void ff_mat_ext2_free(ff_mat_ext2_t *M);

   ff_ext2_elem_t *ff_mat_ext2_at(ff_mat_ext2_t *M, size_t r, size_t c);
   const ff_ext2_elem_t *ff_mat_ext2_atc(const ff_mat_ext2_t *M, size_t r, size_t c);

   void ff_mat_ext2_identity(ff_mat_ext2_t *M);
   void ff_mat_ext2_mul(ff_mat_ext2_t *C, const ff_mat_ext2_t *A, const ff_mat_ext2_t *B);
   bool ff_mat_ext2_inv(ff_mat_ext2_t *X, const ff_mat_ext2_t *A);
   int ff_mat_ext2_nullspace(ff_mat_ext2_t **basis_out, const ff_mat_ext2_t *M);

   /* ============================================================
      Frobenius / Krylov
      ============================================================ */

   void ff_companion_from_poly_base(ff_mat_t *F, const uint32_t *coeffs);
   bool ff_krylov_matrix_base(ff_mat_t *P, const ff_mat_t *A, const uint32_t *v);

   /* ============================================================
      CRT accumulator
      ============================================================ */

   typedef struct ff_crt ff_crt_t;

   ff_crt_t *ff_crt_new(void);
   void ff_crt_free(ff_crt_t *crt);
   void ff_crt_reset(ff_crt_t *crt);
   bool ff_crt_add_residue(ff_crt_t *crt, uint32_t prime, uint32_t residue);
   bool ff_crt_to_u128(ff_u128_t *out, const ff_crt_t *crt);

#ifdef __cplusplus
}
#endif

#endif /* FFLIB_H */
