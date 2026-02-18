/**
 * \file rbc_61_vec.h
 * \brief Interface for vectors of finite field elements
 */

#ifndef RBC_61_VEC_H
#define RBC_61_VEC_H

#include "rbc_61.h"
#include "rbc_61_elt.h"
#include "random_source.h"
#include "seedexpander_shake.h"

void rbc_61_vec_init(rbc_61_vec* v, uint32_t size);

void rbc_61_vec_clear(rbc_61_vec v);

void rbc_61_vec_set_zero(rbc_61_vec v, uint32_t size);

void rbc_61_vec_set(rbc_61_vec o, const rbc_61_vec v, uint32_t size);

void rbc_61_vec_set_random(random_source* ctx, rbc_61_vec o, uint32_t size);

void rbc_61_vec_set_random_tmp(seedexpander_shake_t* ctx, rbc_61_vec o, uint32_t size);

void rbc_61_vec_set_random_from_xof(rbc_61_vec o,
                                 uint32_t size,
                                 void (*xof)(uint8_t *, size_t, const uint8_t *, size_t),
                                 const uint8_t *xof_input,
                                 uint32_t xof_size);

void rbc_61_vec_set_random_fixed_rank(random_source* ctx, rbc_61_vec o, uint32_t size, int32_t rank);

void rbc_61_vec_set_random_full_rank(random_source*, rbc_61_vec o, uint32_t size);

void rbc_61_vec_set_random_full_rank_with_one(random_source* ctx, rbc_61_vec o, uint32_t size);

void rbc_61_vec_set_random_full_rank_with_one_tmp(seedexpander_shake_t* ctx, rbc_61_vec o, uint32_t size);

void rbc_61_vec_set_random_from_support(random_source* ctx, rbc_61_vec o, uint32_t size, const rbc_61_vec support, uint32_t support_size, uint8_t copy_flag);
void rbc_61_vec_set_random_from_support_tmp(seedexpander_shake_t* ctx, rbc_61_vec o, uint32_t size, const rbc_61_vec support, uint32_t support_size, uint8_t copy_flag);
void rbc_61_vec_set_random_pair_from_support(random_source* ctx, rbc_61_vec o1, rbc_61_vec o2, uint32_t size, const rbc_61_vec support, uint32_t support_size, uint8_t copy_flag);

uint32_t rbc_61_vec_gauss(rbc_61_vec v, uint32_t size, uint8_t reduced_flag, rbc_61_vec *other_matrices, uint32_t nMatrices);

uint32_t rbc_61_vec_get_rank(const rbc_61_vec v, uint32_t size);

uint32_t rbc_61_vec_echelonize(rbc_61_vec o, uint32_t size);

void rbc_61_vec_add(rbc_61_vec o, const rbc_61_vec v1, const rbc_61_vec v2, uint32_t size);

void rbc_61_vec_inner_product(rbc_61_elt o, const rbc_61_vec v1, const rbc_61_vec v2, uint32_t size);

void rbc_61_vec_scalar_mul(rbc_61_vec o, const rbc_61_vec v, const rbc_61_elt e, uint32_t size);

void rbc_61_vec_to_string(uint8_t* str, const rbc_61_vec v, uint32_t size);

void rbc_61_vec_from_string(rbc_61_vec v, uint32_t size, const uint8_t* str);

void rbc_61_vec_from_bytes(rbc_61_vec o, uint32_t size, uint8_t *random);

void rbc_61_vec_from_bytes_without_mask(rbc_61_vec o, uint32_t size, const uint8_t *random);

void rbc_61_vec_to_bytes( uint8_t *o, const rbc_61_vec v, uint32_t size);

void rbc_61_vec_print(const rbc_61_vec v, uint32_t size);

#endif

