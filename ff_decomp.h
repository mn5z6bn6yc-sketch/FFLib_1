/*
ff_decomp.h
Decomposition layer header:
  - ff_jordan
  - ff_staircase
  - ff_ratrec
*/

#ifndef FF_DECOMP_H
#define FF_DECOMP_H

#include "fflib.h"

#if defined(__GNUC__) || defined(__clang__)
typedef signed __int128 ff_i128_t;
#else
#error "ff_decomp requires signed __int128 support."
#endif

#ifdef __cplusplus
extern "C"
{
#endif

    /* ============================================================ */
    /* ff_jordan                                                    */
    /* ============================================================ */
    typedef struct
    {
        ff_mat_ext2_t *J;
        ff_mat_ext2_t *P;
        size_t block_count;
        size_t *block_sizes;
        ff_ext2_elem_t *block_eigenvalues;
    } ff_jordan_result_t;

    /* Jordan decomposition over an extension field. Eigenvalues must be supplied. */
    bool ff_jordan_ext2_with_eigenvalues(
        const ff_ext2_t *E,
        const ff_mat_ext2_t *A,
        const ff_ext2_elem_t *eigs,
        size_t eig_count,
        ff_jordan_result_t *out);

    /* Convenience wrapper: Converts a base-field matrix A into an extension-field matrix. */
    bool ff_jordan_base_with_eigenvalues_ext2(
        const ff_ext2_t *E,
        const ff_mat_t *A_base,
        const ff_ext2_elem_t *eigs,
        size_t eig_count,
        ff_jordan_result_t *out);

    void ff_jordan_result_clear(ff_jordan_result_t *res);

    /* ============================================================ */
    /* ff_staircase                                                 */
    /* ============================================================ */
    /* Deterministic finite-field staircase / Hessenberg similarity reduction. */
    bool ff_staircase_hessenberg_base(
        const ff_mat_t *A,
        ff_mat_t **P_out,
        ff_mat_t **H_out);

    /* ============================================================ */
    /* ff_ratrec                                                    */
    /* ============================================================ */
    /* Rational reconstruction using __int128. */
    bool ff_ratrec_from_u128(ff_u128_t x, ff_u128_t M, ff_i128_t *num, ff_i128_t *den);
    bool ff_ratrec_crt_i128(const ff_crt_t *crt, ff_i128_t *num, ff_i128_t *den);
    bool ff_ratrec_crt_i64(const ff_crt_t *crt, int64_t *num, int64_t *den);
    bool ff_ratrec_array_i128(
        const uint32_t *primes,
        const uint32_t *residues,
        size_t count,
        ff_i128_t *num,
        ff_i128_t *den);

#ifdef __cplusplus
}
#endif

#endif /* FF_DECOMP_H */