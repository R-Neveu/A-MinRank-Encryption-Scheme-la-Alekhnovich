/**
 * \file rbc_61_perm.h
 * \brief Interface for permutation matrices over Fq
 */

#ifndef RBC_61_PERM_H
#define RBC_61_PERM_H

#include "rbc_61.h"
#include "rbc_61_elt.h"
#include "rbc_61_vec.h"
#include "rbc_61_mat_fq.h"
#include "seedexpander.h"

void rbc_61_perm_init(rbc_61_perm* perm, uint32_t size);
void rbc_61_perm_clear(rbc_61_perm m);
void rbc_61_perm_set_zero(rbc_61_perm o, uint32_t size);
void rbc_61_perm_set(rbc_61_perm o, const rbc_61_perm perm, uint32_t size);
void rbc_61_perm_set_random(random_source* ctx, rbc_61_perm perm, uint32_t size);
void rbc_61_perm_set_random_from_xof(rbc_61_perm perm, uint32_t size, void (*xof)(uint8_t *, size_t, const uint8_t *, size_t), const uint8_t *xof_input, uint32_t xof_size);
void rbc_61_perm_apply(rbc_61_mat_fq o, const rbc_61_mat_fq m, const rbc_61_perm perm, uint32_t size);
#endif

