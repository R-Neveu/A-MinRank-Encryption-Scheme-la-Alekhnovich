/**
 * \file rbc_61_vspace.h
 * \brief Interface for subspaces of Fq^m
 */

#ifndef RBC_61_VSPACE_H
#define RBC_61_VSPACE_H

#include "rbc_61.h"
#include "random_source.h"
#include "seedexpander_shake.h"

void rbc_61_vspace_init(rbc_61_vspace* vs, uint32_t size);

void rbc_61_vspace_clear(rbc_61_vspace vs);

void rbc_61_vspace_set(rbc_61_vspace o, const rbc_61_vspace vs, uint32_t size);

void rbc_61_vspace_set_zero(rbc_61_vspace o, uint32_t size);

void rbc_61_vspace_set_random_full_rank(random_source* ctx, rbc_61_vspace o, uint32_t size);

void rbc_61_vspace_set_random_full_rank_with_one(random_source* ctx, rbc_61_vspace o, uint32_t size);

void rbc_61_vspace_set_random_full_rank_with_one_tmp(seedexpander_shake_t* ctx, rbc_61_vspace o, uint32_t size);

void rbc_61_vspace_direct_sum(rbc_61_vspace o, const rbc_61_vspace vs1, const rbc_61_vspace vs2, uint32_t vs1_size, uint32_t vs2_size);

uint32_t rbc_61_vspace_intersection(rbc_61_vspace o, const rbc_61_vspace vs1, const rbc_61_vspace vs2, uint32_t vs1_size, uint32_t vs2_size);

void rbc_61_vspace_product(rbc_61_vspace o, const rbc_61_vspace vs1, const rbc_61_vspace vs2, uint32_t vs1_size, uint32_t vs2_size);

void rbc_61_vspace_print(const rbc_61_vspace vs, uint32_t size);

#endif

