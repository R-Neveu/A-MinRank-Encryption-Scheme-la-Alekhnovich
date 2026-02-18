/**
 * \file rbc_37_perm.h
 * \brief Interface for permutation matrices over Fq
 */

#ifndef RBC_37_PERM_H
#define RBC_37_PERM_H

#include "rbc_37.h"
#include "rbc_37_elt.h"
#include "rbc_37_vec.h"
#include "rbc_37_mat_fq.h"
#include "seedexpander.h"

void rbc_37_perm_init(rbc_37_perm* perm, uint32_t size);
void rbc_37_perm_clear(rbc_37_perm m);
void rbc_37_perm_set_zero(rbc_37_perm o, uint32_t size);
void rbc_37_perm_set(rbc_37_perm o, const rbc_37_perm perm, uint32_t size);
void rbc_37_perm_set_random(random_source* ctx, rbc_37_perm perm, uint32_t size);
void rbc_37_perm_set_random_from_xof(rbc_37_perm perm, uint32_t size, void (*xof)(uint8_t *, size_t, const uint8_t *, size_t), const uint8_t *xof_input, uint32_t xof_size);
void rbc_37_perm_apply(rbc_37_mat_fq o, const rbc_37_mat_fq m, const rbc_37_perm perm, uint32_t size);
#endif

