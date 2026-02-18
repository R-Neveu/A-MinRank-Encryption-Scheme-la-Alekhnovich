/**
 * \file rbc_37_vspace.h
 * \brief Interface for subspaces of Fq^m
 */

#ifndef RBC_37_VSPACE_H
#define RBC_37_VSPACE_H

#include "rbc_37.h"
#include "random_source.h"
#include "seedexpander_shake.h"

void rbc_37_vspace_init(rbc_37_vspace* vs, uint32_t size);

void rbc_37_vspace_clear(rbc_37_vspace vs);

void rbc_37_vspace_set(rbc_37_vspace o, const rbc_37_vspace vs, uint32_t size);

void rbc_37_vspace_set_zero(rbc_37_vspace o, uint32_t size);

void rbc_37_vspace_set_random_full_rank(random_source* ctx, rbc_37_vspace o, uint32_t size);

void rbc_37_vspace_set_random_full_rank_with_one(random_source* ctx, rbc_37_vspace o, uint32_t size);

void rbc_37_vspace_set_random_full_rank_with_one_tmp(seedexpander_shake_t* ctx, rbc_37_vspace o, uint32_t size);

void rbc_37_vspace_direct_sum(rbc_37_vspace o, const rbc_37_vspace vs1, const rbc_37_vspace vs2, uint32_t vs1_size, uint32_t vs2_size);

uint32_t rbc_37_vspace_intersection(rbc_37_vspace o, const rbc_37_vspace vs1, const rbc_37_vspace vs2, uint32_t vs1_size, uint32_t vs2_size);

void rbc_37_vspace_product(rbc_37_vspace o, const rbc_37_vspace vs1, const rbc_37_vspace vs2, uint32_t vs1_size, uint32_t vs2_size);

void rbc_37_vspace_print(const rbc_37_vspace vs, uint32_t size);

#endif

