/**
 * \file rbc_61_qre.h
 * \brief Interface for quotient ring elements 
 * */

#ifndef RBC_61_QRE_H
#define RBC_61_QRE_H

#include "rbc_61.h"
#include "random_source.h"

void rbc_61_qre_init_modulus(uint32_t degree);

void rbc_61_qre_clear_modulus(void);

rbc_61_poly_sparse rbc_61_qre_get_modulus(uint32_t degree);

void rbc_61_qre_init(rbc_61_qre* p);

void rbc_61_qre_clear(rbc_61_qre p);

void rbc_61_qre_set_zero(rbc_61_qre o);

void rbc_61_qre_set_random(random_source* ctx, rbc_61_qre o);

void rbc_61_qre_set_random_from_xof(rbc_61_qre o,
                                 void (*xof)(uint8_t *, size_t, const uint8_t *, size_t),
                                 const uint8_t *xof_input,
                                 uint32_t xof_size);

void rbc_61_qre_set_random_full_rank(random_source* ctx, rbc_61_qre o);

void rbc_61_qre_set_random_full_rank_with_one(random_source* ctx, rbc_61_qre o);

void rbc_61_qre_set_random_from_support(random_source* ctx, rbc_61_qre o, const rbc_61_vspace support, uint32_t support_size, uint8_t copy_flag);

void rbc_61_qre_set_random_pair_from_support(random_source* ctx, rbc_61_qre o1, rbc_61_qre o2, const rbc_61_vspace support, uint32_t support_size, uint8_t copy_flag);

uint8_t rbc_61_qre_is_equal_to(const rbc_61_qre p1, const rbc_61_qre p2);

void rbc_61_qre_add(rbc_61_qre o, const rbc_61_qre p1, const rbc_61_qre p2);

void rbc_61_qre_mul(rbc_61_qre o, const rbc_61_qre p1, const rbc_61_qre p2);

void rbc_61_qre_inv(rbc_61_qre o, const rbc_61_qre p);

void rbc_61_qre_to_string(uint8_t* str, const rbc_61_qre p);

void rbc_61_qre_from_string(rbc_61_qre o, const uint8_t* str);

void rbc_61_qre_print(const rbc_61_qre p);

void rbc_61_qre_set_monomial(rbc_61_qre o, uint32_t degree);

void rbc_61_qre_rot(rbc_61_qre o, rbc_61_qre e, uint32_t rot);

#endif

