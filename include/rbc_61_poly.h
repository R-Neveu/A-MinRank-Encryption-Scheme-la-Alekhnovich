/**
 * \file rbc_61_poly.h
 * \brief Interface for polynomials over a finite field
 */

#ifndef RBC_61_POLY_H
#define RBC_61_POLY_H

#include "rbc_61.h"
#include "rbc_61_vspace.h"
#include "random_source.h"

void rbc_61_poly_init(rbc_61_poly* p, int32_t degree);

void rbc_61_poly_clear(rbc_61_poly p);

void rbc_61_poly_sparse_init(rbc_61_poly_sparse* p, uint32_t coeff_nb, uint32_t *coeffs);

void rbc_61_poly_sparse_clear(rbc_61_poly_sparse p);

void rbc_61_poly_set_zero(rbc_61_poly o, int32_t degree);

void rbc_61_poly_set(rbc_61_poly o, const rbc_61_poly p);

void rbc_61_poly_set_random(random_source* ctx, rbc_61_poly o, int32_t degree);

void rbc_61_poly_set_random_from_xof(rbc_61_poly o,
                                  int32_t degree,
                                  void (*xof)(uint8_t *, size_t, const uint8_t *, size_t),
                                  const uint8_t *xof_input,
                                  uint32_t xof_size);

void rbc_61_poly_set_random_full_rank(random_source* ctx, rbc_61_poly o, int32_t degree);

void rbc_61_poly_set_random_full_rank_with_one(random_source* ctx, rbc_61_poly o, int32_t degree);

void rbc_61_poly_set_random_from_support(random_source* ctx, rbc_61_poly o, int32_t degree, const rbc_61_vspace support, uint32_t support_size, uint8_t copy_flag);

void rbc_61_poly_set_random_pair_from_support(random_source* ctx, rbc_61_poly o1, rbc_61_poly o2, int32_t degree, const rbc_61_vspace support, uint32_t support_size, uint8_t copy_flag);

uint8_t rbc_61_poly_is_equal_to(const rbc_61_poly p1, const rbc_61_poly p2);

void rbc_61_poly_add(rbc_61_poly o, const rbc_61_poly p1, const rbc_61_poly p2);

void rbc_61_poly_mul(rbc_61_poly o, const rbc_61_poly p1, const rbc_61_poly p2);
void rbc_61_poly_mul2(rbc_61_poly o, const rbc_61_poly p1, const rbc_61_poly p2, int32_t p1_degree, int32_t p2_degree);

void rbc_61_poly_mulmod_sparse(rbc_61_poly o, const rbc_61_poly p1, const rbc_61_poly p2, const rbc_61_poly_sparse modulus);

void rbc_61_poly_inv(rbc_61_poly o, const rbc_61_poly p, const rbc_61_poly modulus);

void rbc_61_poly_to_string(uint8_t* str, const rbc_61_poly p);

void rbc_61_poly_from_string(rbc_61_poly p, const uint8_t* str);

void rbc_61_poly_print(const rbc_61_poly p);

void rbc_61_poly_sparse_print(const rbc_61_poly_sparse p);

void rbc_61_poly_set_monomial(rbc_61_poly o, uint32_t degree);

#endif

