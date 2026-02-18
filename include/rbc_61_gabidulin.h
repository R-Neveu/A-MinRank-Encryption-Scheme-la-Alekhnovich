/**
 * \file gabidulin.h
 * \brief Functions to encode and decode messages using Gabidulin codes
 *
 * The decoding algorithm provided is based on q_polynomials reconstruction, see \cite gabidulin:welch and \cite gabidulin:generalized for details.
 *
 */

#ifndef RBC_61_GABIDULIN_H
#define RBC_61_GABIDULIN_H

#include "rbc_61_vec.h"
#include "rbc_61_qre.h"


/**
  * \typedef rbc_61_gabidulin
  * \brief Structure of a gabidulin code
  */
typedef struct rbc_61_gabidulin {
  rbc_61_poly g; /**< Generator vector defining the code */
  uint32_t k; /**< Size of vectors representing messages */
  uint32_t n; /**< Size of vectors representing codewords */
} rbc_61_gabidulin;


void rbc_61_gabidulin_init(rbc_61_gabidulin* code, const rbc_61_poly g, uint32_t k, uint32_t n);

void rbc_61_gabidulin_encode(rbc_61_qre c, const rbc_61_gabidulin gc, const rbc_61_vec m);
void rbc_61_gabidulin_decode(rbc_61_vec m, const rbc_61_gabidulin gc, const rbc_61_qre y);

#endif

