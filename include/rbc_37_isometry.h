/**
 * \file rbc_37_isometry.h
 * \brief Interface isometries over Fq^m
 */

#ifndef RBC_37_ISOMETRY_H
#define RBC_37_ISOMETRY_H

#include "rbc_37.h"
#include "rbc_37_elt.h"
#include "rbc_37_vec.h"
#include "rbc_37_mat_fq.h"
#include "seedexpander.h"

void rbc_37_isometry_init(rbc_37_isometry* o, uint32_t n);
void rbc_37_isometry_clear(rbc_37_isometry* o);
void rbc_37_isometry_set_zero(rbc_37_isometry *iso, uint32_t n);
void rbc_37_isometry_set(rbc_37_isometry *o, const rbc_37_mat_fq P, const rbc_37_mat_fq Q, uint32_t n);
void rbc_37_isometry_set_random(random_source* ctx, rbc_37_isometry *o, uint32_t n);
void rbc_37_isometry_set_random_from_xof(rbc_37_isometry *o,
                                      uint32_t n,
                                      void (*xof)(uint8_t *, size_t, const uint8_t *, size_t),
                                      const uint8_t *xof_input,
                                      uint32_t xof_size);
void rbc_37_isometry_apply(rbc_37_vec o, const rbc_37_vec v, uint32_t n, const rbc_37_isometry *iso);
void rbc_37_isometry_compose(rbc_37_isometry *o, const rbc_37_isometry *iso1, const rbc_37_isometry *iso2);
#endif

