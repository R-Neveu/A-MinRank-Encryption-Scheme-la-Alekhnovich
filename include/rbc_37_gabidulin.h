/**
 * \file gabidulin.h
 * \brief Functions to encode and decode messages using Gabidulin codes
 *
 * The decoding algorithm provided is based on q_polynomials reconstruction, see \cite gabidulin:welch and \cite gabidulin:generalized for details.
 *
 */

#ifndef RBC_37_GABIDULIN_H
#define RBC_37_GABIDULIN_H

#include "rbc_37_vec.h"
#include "rbc_37_qre.h"


/**
  * \typedef rbc_37_gabidulin
  * \brief Structure of a gabidulin code
  */
typedef struct rbc_37_gabidulin {
  rbc_37_poly g; /**< Generator vector defining the code */
  uint32_t k; /**< Size of vectors representing messages */
  uint32_t n; /**< Size of vectors representing codewords */
} rbc_37_gabidulin;


void rbc_37_gabidulin_init(rbc_37_gabidulin* code, const rbc_37_poly g, uint32_t k, uint32_t n);

void rbc_37_gabidulin_encode(rbc_37_qre c, const rbc_37_gabidulin gc, const rbc_37_vec m);
void rbc_37_gabidulin_decode(rbc_37_vec m, const rbc_37_gabidulin gc, const rbc_37_qre y);

#endif

