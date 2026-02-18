/**
 * \file rbc_79_perm.h
 * \brief Interface for permutation matrices over Fq
 */

#ifndef RBC_79_PERM_H
#define RBC_79_PERM_H

#include "rbc_79.h"
#include "rbc_79_elt.h"
#include "rbc_79_vec.h"
#include "rbc_79_mat_fq.h"
#include "seedexpander.h"

void rbc_79_perm_init(rbc_79_perm* perm, uint32_t size);
void rbc_79_perm_clear(rbc_79_perm m);
void rbc_79_perm_set_zero(rbc_79_perm o, uint32_t size);
void rbc_79_perm_set(rbc_79_perm o, const rbc_79_perm perm, uint32_t size);
void rbc_79_perm_set_random(random_source* ctx, rbc_79_perm perm, uint32_t size);
void rbc_79_perm_set_random_from_xof(rbc_79_perm perm, uint32_t size, void (*xof)(uint8_t *, size_t, const uint8_t *, size_t), const uint8_t *xof_input, uint32_t xof_size);
void rbc_79_perm_apply(rbc_79_mat_fq o, const rbc_79_mat_fq m, const rbc_79_perm perm, uint32_t size);
#endif

