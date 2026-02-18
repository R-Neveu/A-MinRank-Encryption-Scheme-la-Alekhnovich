/**
 * \file gabidulin.h
 * \brief Functions to encode and decode messages using Gabidulin codes
 *
 * The decoding algorithm provided is based on q_polynomials reconstruction, see \cite gabidulin:welch and \cite gabidulin:generalized for details.
 *
 */

#ifndef RBC_79_GABIDULIN_H
#define RBC_79_GABIDULIN_H

#include "rbc_79_vec.h"
#include "rbc_79_qre.h"


/**
  * \typedef rbc_79_gabidulin
  * \brief Structure of a gabidulin code
  */
typedef struct rbc_79_gabidulin {
  rbc_79_poly g; /**< Generator vector defining the code */
  uint32_t k; /**< Size of vectors representing messages */
  uint32_t n; /**< Size of vectors representing codewords */
} rbc_79_gabidulin;


void rbc_79_gabidulin_init(rbc_79_gabidulin* code, const rbc_79_poly g, uint32_t k, uint32_t n);

void rbc_79_gabidulin_encode(rbc_79_qre c, const rbc_79_gabidulin gc, const rbc_79_vec m);
void rbc_79_gabidulin_decode(rbc_79_vec m, const rbc_79_gabidulin gc, const rbc_79_qre y);

#endif

