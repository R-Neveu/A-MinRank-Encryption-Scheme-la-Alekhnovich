/**
 * \file rbc_79_isometry.h
 * \brief Interface isometries over Fq^m
 */

#ifndef RBC_79_ISOMETRY_H
#define RBC_79_ISOMETRY_H

#include "rbc_79.h"
#include "rbc_79_elt.h"
#include "rbc_79_vec.h"
#include "rbc_79_mat_fq.h"
#include "seedexpander.h"

void rbc_79_isometry_init(rbc_79_isometry* o, uint32_t n);
void rbc_79_isometry_clear(rbc_79_isometry* o);
void rbc_79_isometry_set_zero(rbc_79_isometry *iso, uint32_t n);
void rbc_79_isometry_set(rbc_79_isometry *o, const rbc_79_mat_fq P, const rbc_79_mat_fq Q, uint32_t n);
void rbc_79_isometry_set_random(random_source* ctx, rbc_79_isometry *o, uint32_t n);
void rbc_79_isometry_set_random_from_xof(rbc_79_isometry *o,
                                      uint32_t n,
                                      void (*xof)(uint8_t *, size_t, const uint8_t *, size_t),
                                      const uint8_t *xof_input,
                                      uint32_t xof_size);
void rbc_79_isometry_apply(rbc_79_vec o, const rbc_79_vec v, uint32_t n, const rbc_79_isometry *iso);
void rbc_79_isometry_compose(rbc_79_isometry *o, const rbc_79_isometry *iso1, const rbc_79_isometry *iso2);
#endif

