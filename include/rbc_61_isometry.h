/**
 * \file rbc_61_isometry.h
 * \brief Interface isometries over Fq^m
 */

#ifndef RBC_61_ISOMETRY_H
#define RBC_61_ISOMETRY_H

#include "rbc_61.h"
#include "rbc_61_elt.h"
#include "rbc_61_vec.h"
#include "rbc_61_mat_fq.h"
#include "seedexpander.h"

void rbc_61_isometry_init(rbc_61_isometry* o, uint32_t n);
void rbc_61_isometry_clear(rbc_61_isometry* o);
void rbc_61_isometry_set_zero(rbc_61_isometry *iso, uint32_t n);
void rbc_61_isometry_set(rbc_61_isometry *o, const rbc_61_mat_fq P, const rbc_61_mat_fq Q, uint32_t n);
void rbc_61_isometry_set_random(random_source* ctx, rbc_61_isometry *o, uint32_t n);
void rbc_61_isometry_set_random_from_xof(rbc_61_isometry *o,
                                      uint32_t n,
                                      void (*xof)(uint8_t *, size_t, const uint8_t *, size_t),
                                      const uint8_t *xof_input,
                                      uint32_t xof_size);
void rbc_61_isometry_apply(rbc_61_vec o, const rbc_61_vec v, uint32_t n, const rbc_61_isometry *iso);
void rbc_61_isometry_compose(rbc_61_isometry *o, const rbc_61_isometry *iso1, const rbc_61_isometry *iso2);
#endif

