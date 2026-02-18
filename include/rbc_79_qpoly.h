/** 
 * \file rbc_79_qpoly.h
 * \brief Functions to manipulate q-polynomials. 
 *
 * The rbc_79_qpolys are polynomials over a field of the form \f$ P(x) = \sum_{i=0}^{n} p_i \times x^{q^i} \f$ with \f$ q \f$ a rational prime.
 * See \cite ore:qpolynomials for a description of their main properties.
 *
 */

#ifndef RBC_79_QPOLY_H
#define RBC_79_QPOLY_H

#include "rbc_79_elt.h"
#include "rbc_79_poly.h"

/**
  * \typedef rbc_79_qpoly
  * \brief Structure of a q-polynomial
  */
typedef struct {
  rbc_79_vec values; /**< Coefficients of the q-polynomial */
  int32_t max_degree; /**< Maximum q-degree that the q-polynomial may have. This value never changes */
  int32_t degree; /**< Q-degree of the q-polynomial. This value is updated whenever necessary */
} rbc_79_qpoly_struct;
typedef rbc_79_qpoly_struct* rbc_79_qpoly;

void rbc_79_qpoly_init(rbc_79_qpoly* p, int32_t max_degree);

void rbc_79_qpoly_clear(rbc_79_qpoly o);

void rbc_79_qpoly_update_degree(rbc_79_qpoly p, int32_t position);

int8_t rbc_79_qpoly_set(rbc_79_qpoly o, const rbc_79_qpoly p);

int8_t rbc_79_qpoly_set_mask(rbc_79_qpoly o, const rbc_79_qpoly p1, const rbc_79_qpoly p2, uint32_t mask);

int8_t rbc_79_qpoly_set_zero(rbc_79_qpoly o);

int8_t rbc_79_qpoly_set_one(rbc_79_qpoly o);

int8_t rbc_79_qpoly_set_interpolate_vect_and_zero(rbc_79_qpoly o1, rbc_79_qpoly o2, const rbc_79_poly p1, const rbc_79_poly p2, int32_t size);

int8_t rbc_79_qpoly_annihilator(rbc_79_qpoly p, const rbc_79_vec v, int32_t size);

int8_t rbc_79_qpoly_is_zero(const rbc_79_qpoly p);

void rbc_79_qpoly_evaluate(rbc_79_elt o, const rbc_79_qpoly p, const rbc_79_elt e);

int8_t rbc_79_qpoly_scalar_mul(rbc_79_qpoly o, const rbc_79_qpoly p, const rbc_79_elt e);

int8_t rbc_79_qpoly_qexp(rbc_79_qpoly o, const rbc_79_qpoly p);

int8_t rbc_79_qpoly_add(rbc_79_qpoly o, const rbc_79_qpoly p1, const rbc_79_qpoly p2);

int8_t rbc_79_qpoly_mul(rbc_79_qpoly o, const rbc_79_qpoly p1, const rbc_79_qpoly p2);

int8_t rbc_79_qpoly_mul2(rbc_79_qpoly o, const rbc_79_qpoly p1, const rbc_79_qpoly p2, uint32_t p1_degree, uint32_t p2_degree);

int8_t rbc_79_qpoly_left_div(rbc_79_qpoly q, rbc_79_qpoly r, const rbc_79_qpoly a, const rbc_79_qpoly b);

int8_t rbc_79_qpoly_left_div2(rbc_79_qpoly q, rbc_79_qpoly r, const rbc_79_qpoly a, const rbc_79_qpoly b, uint32_t capacity, uint32_t k);

void rbc_79_qpoly_print(const rbc_79_qpoly p);

#endif

