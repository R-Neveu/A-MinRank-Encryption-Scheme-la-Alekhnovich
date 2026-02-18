/**
 * \file rbc_37_qre.c
 * \brief Implementation of rbc_37_qre.h
 */

#include "rbc_37_elt.h"
#include "rbc_37_qre.h"
#include "rbc_37_poly.h"

static uint32_t RBC_37_QRE_DEG;
static int rbc_37_init_qre_modulus = 0;
static rbc_37_poly_sparse rbc_37_qre_modulus;
static rbc_37_poly rbc_37_qre_modulus_inv;




/**
 * \fn void rbc_37_qre_init_modulus(uint32_t degree) {
 * \brief This function init the quotient ring modulus corresponding to PARAM_N.
 *
 * \param[in] degree Degree of the modulus
 */
void rbc_37_qre_init_modulus(uint32_t degree) {
  if(rbc_37_init_qre_modulus == 0) {
    RBC_37_QRE_DEG = degree - 1;
    rbc_37_qre_modulus = rbc_37_qre_get_modulus(degree);
    rbc_37_poly_init(&rbc_37_qre_modulus_inv, RBC_37_QRE_DEG + 1);
    rbc_37_poly_set_zero(rbc_37_qre_modulus_inv, RBC_37_QRE_DEG + 1);
    for(uint32_t i = 0 ; i < rbc_37_qre_modulus->coeffs_nb ; i++) {
      rbc_37_elt_set_one(rbc_37_qre_modulus_inv->v[rbc_37_qre_modulus->coeffs[i]]);
    }
  }

  rbc_37_init_qre_modulus++;
}




/**
 * \fn void rbc_37_qre_clear_modulus(void) {
 * \brief This function clears the quotient ring modulus in use.
 */
void rbc_37_qre_clear_modulus(void) {
  rbc_37_init_qre_modulus--;

  if(rbc_37_init_qre_modulus == 0) {
    rbc_37_poly_sparse_clear(rbc_37_qre_modulus);
    rbc_37_poly_clear(rbc_37_qre_modulus_inv);
  }
}

/**
 * \fn rbc_37_poly_sparse rbc_37_qre_get_modulus(uint32_t degree) {
 * \brief This function return the sparse polynomial used as the quotient ring modulus for PARAM_N.
 *
 * \param[in] degree Degree of the modulus
 */
rbc_37_poly_sparse rbc_37_qre_get_modulus(uint32_t degree) {
  rbc_37_poly_sparse modulus = NULL;
  uint32_t* values;
  if(degree == 37) {
    values = (uint32_t*) malloc(5 * sizeof(uint32_t));
    values[0] = 0;
    values[1] = 1;
    values[2] = 4;
    values[3] = 6;
    values[4] = 37;
    rbc_37_poly_sparse_init(&modulus, 5, values);
    free(values);
  }
  else if(degree == 61) {
    values = (uint32_t*) malloc(5 * sizeof(uint32_t));
    values[0] = 0;
    values[1] = 1;
    values[2] = 2;
    values[3] = 5;
    values[4] = 61;
    rbc_37_poly_sparse_init(&modulus, 5, values);
    free(values);
  }
  else if(degree == 79) {
    values = (uint32_t*) malloc(3 * sizeof(uint32_t));
    values[0] = 0;
    values[1] = 9;
    values[2] = 79;
    rbc_37_poly_sparse_init(&modulus, 3, values);
    free(values);
  }

 return modulus;
}




/**
 * \fn void rbc_37_qre_init(rbc_37_qre* p)
 * \brief This function allocates the memory for a rbc_37_qre element.
 *
 * \param[out] p Pointer to the allocated rbc_37_qre
 */
void rbc_37_qre_init(rbc_37_qre* p) {
  rbc_37_poly_init(p, RBC_37_QRE_DEG);
}




/**
 * \fn void rbc_37_qre_clear(rbc_37_qre p)
 * \brief This function clears the memory allocated for a rbc_37_qre element.
 *
 * \param[in] p rbc_37_qre
 */
void rbc_37_qre_clear(rbc_37_qre p) {
  rbc_37_poly_clear(p);
}




/**
 * \fn void rbc_37_qre_set_zero(rbc_37_qre o)
 * \brief This functions sets a rbc_37_qre to zero.
 *
 * \param[in] o rbc_37_qre
 */
void rbc_37_qre_set_zero(rbc_37_qre o) {
  rbc_37_poly_set_zero(o, RBC_37_QRE_DEG);
}




/**
 * \fn void rbc_37_qre_set_random(random_source* ctx, rbc_37_qre o)
 * \brief This function sets a rbc_37_qre with random values using NIST seed expander.
 *
 * \param[out] ctx random source
 * \param[out] o rbc_37_qre
 */
void rbc_37_qre_set_random(random_source* ctx, rbc_37_qre o) {
  rbc_37_poly_set_random(ctx, o, RBC_37_QRE_DEG);
}




/**
 * \fn rbc_37_qre_set_random_from_xof(rbc_37_qre o, void (*xof)(uint8_t *, size_t, const uint8_t *, size_t), const uint8_t *xof_input, uint32_t xof_size)
 * \brief This function sets a rbc_37_qre with random values using XOF.
 *
 * \param[out] o rbc_37_qre
 * \param[in] xof XOF procedure
 * \param[in] xof_input XOF input byte array
 * \param[in] xof_size XOF input byte array length of xof_input
 */
void rbc_37_qre_set_random_from_xof(rbc_37_qre o,
                                 void (*xof)(uint8_t *, size_t, const uint8_t *, size_t),
                                 const uint8_t *xof_input,
                                 uint32_t xof_size) {
  rbc_37_poly_set_random_from_xof(o, RBC_37_QRE_DEG, xof, xof_input, xof_size);
}




/**
 * \fn void rbc_37_qre_set_random_full_rank(random_source* ctx, rbc_37_qre o)
 * \brief This function sets a rbc_37_qre with random values using the NIST seed expander. The returned rbc_37_qre has full rank.
 *
 * \param[out] ctx random source
 * \param[out] o rbc_37_qre
 */
void rbc_37_qre_set_random_full_rank(random_source* ctx, rbc_37_qre o) {
  rbc_37_poly_set_random_full_rank(ctx, o, RBC_37_QRE_DEG);
}




/**
 * \fn void rbc_37_qre_set_random_full_rank_with_one(random_source* ctx, rbc_37_qre o)
 * \brief This function sets a rbc_37_qre with random values using the NIST seed expander. The returned rbc_37_qre has full rank and contains one.
 *
 * \param[out] ctx random source
 * \param[out] o rbc_37_qre
 */
void rbc_37_qre_set_random_full_rank_with_one(random_source* ctx, rbc_37_qre o) {
  rbc_37_poly_set_random_full_rank_with_one(ctx, o, RBC_37_QRE_DEG);
}




/**
 * \fn void rbc_37_qre_set_random_from_support(random_source* ctx, rbc_37_qre o, const rbc_37_vspace support, uint32_t support_size)
 * \brief This function sets a rbc_37_qre with random values using the NIST seed expander. The support of the rbc_37_qre returned by this function is the one given in input.
 *
 * \param[out] ctx random source
 * \param[out] o rbc_37_qre
 * \param[in] support Support of <b>o</b>
 * \param[in] support_size Size of the support
 * \param[in] copy_flag If not 0, the support is copied into random coordinates of the resulting vector
 */
void rbc_37_qre_set_random_from_support(random_source* ctx, rbc_37_qre o, const rbc_37_vspace support, uint32_t support_size, uint8_t copy_flag) {
  rbc_37_poly_set_random_from_support(ctx, o, RBC_37_QRE_DEG, support, support_size, copy_flag);
}




/**
 * \fn void rbc_37_qre_set_random_pair_from_support(random_source* ctx, rbc_37_qre o1, rbc_37_qre o2, const rbc_37_vspace support, uint32_t support_size)
 * \brief This function sets a pair of rbc_37_qre with random values using the NIST seed expander. The support of (o1 || o2) is the one given in input.
 *
 * \param[out] ctx random_source
 * \param[out] o1 rbc_37_qre
 * \param[out] o2 rbc_37_qre
 * \param[in] support Support
 * \param[in] support_size Size of the support
 * \param[in] copy_flag If not 0, the support is copied into random coordinates of the resulting vector
 */
void rbc_37_qre_set_random_pair_from_support(random_source* ctx, rbc_37_qre o1, rbc_37_qre o2, const rbc_37_vspace support, uint32_t support_size, uint8_t copy_flag) {
  rbc_37_poly_set_random_pair_from_support(ctx, o1, o2, RBC_37_QRE_DEG, support, support_size, copy_flag);
}




/**
 * \fn uint8_t rbc_37_qre_is_equal_to(const rbc_37_qre p1, const rbc_37_qre p2)
 * \brief This function test if two rbc_37_qre are equal.
 *
 * \param[in] p1 rbc_37_qre
 * \param[in] p2 rbc_37_qre
 * \return 1 if the rbc_37_qre are equal, 0 otherwise
 */
uint8_t rbc_37_qre_is_equal_to(const rbc_37_qre p1, const rbc_37_qre p2) {
  return rbc_37_poly_is_equal_to(p1, p2);
}




/**
 * \fn void rbc_37_qre_add(rbc_37_qre o, const rbc_37_qre p1, const rbc_37_qre p2)
 * \brief This function adds two rbc_37_qre.
 *
 * \param[out] o Sum of <b>p1</b> and <b>p2</b>
 * \param[in] p1 rbc_37_qre
 * \param[in] p2 rbc_37_qre
 */
void rbc_37_qre_add(rbc_37_qre o, const rbc_37_qre p1, const rbc_37_qre p2) {
  rbc_37_poly_add(o, p1, p2);
}




/**
 * \fn void rbc_37_qre_mul(rbc_37_qre o, const rbc_37_qre p1, const rbc_37_qre p2)
 * \brief This function multiplies two rbc_37_qre.
 *
 * \param[out] o Product of <b>p1</b> and <b>p2</b>
 * \param[in] p1 rbc_37_qre
 * \param[in] p2 rbc_37_qre
 */
void rbc_37_qre_mul(rbc_37_qre o, const rbc_37_qre p1, const rbc_37_qre p2) {
  if(rbc_37_init_qre_modulus == 0) {
    printf("Call to rbc_37_qre_mul with uninitialized modulus\n");
    exit(1);
  }

  rbc_37_poly_mulmod_sparse(o, p1, p2, rbc_37_qre_modulus);
}




/**
 * \fn void rbc_37_qre_inv(rbc_37_qre o, const rbc_37_qre p)
 * \brief This function performs the extended euclidean algorithm to compute the inverse of p.
 *
 * \param[out] o Inverse of <b>e</b> modulo <b>modulus</b>
 * \param[in] p rbc_37_qre
 * \param[in] modulus Polynomial
 */
void rbc_37_qre_inv(rbc_37_qre o, const rbc_37_qre p) {
  if(rbc_37_init_qre_modulus == 0) {
    printf("Call to rbc_37_qre_inv with uninitialized modulus\n");
    exit(1);
  }

  rbc_37_poly_inv(o, p, rbc_37_qre_modulus_inv);
}




/**
 * \fn void rbc_37_qre_to_string(uint8_t* str, const rbc_37_qre p)
 * \brief This function parses a rbc_37_qre into a string.
 *
 * \param[out] str String
 * \param[in] e rbc_37_qre
 */
void rbc_37_qre_to_string(uint8_t* str, const rbc_37_qre p) {
  rbc_37_poly_to_string(str, p);
}




/**
 * \fn void rbc_37_qre_from_string(rbc_37_qre o, const uint8_t* str)
 * \brief This function parses a string into a rbc_37_qre.
 *
 * \param[out] o rbc_37_qre
 * \param[in] str String to parse
 */
void rbc_37_qre_from_string(rbc_37_qre o, const uint8_t* str) {
  rbc_37_poly_from_string(o, str);
}




/**
 * \fn void rbc_37_qre_print(const rbc_37_qre p)
 * \brief This function displays a rbc_37_qre.
 *
 * \param[in] p rbc_37_qre
 */
void rbc_37_qre_print(const rbc_37_qre p) {
  rbc_37_poly_print(p);
}




/**
 * \fn void rbc_37_poly_set_monial(rbc_37_qre o, uint32_t degree)
 * \brief This functions sets a polynomial to x^degree.
 *
 * \param[in] p Polynomial
 * \param[in] degree Degree of the monomial
 */
void rbc_37_qre_set_monomial(rbc_37_qre o, uint32_t degree) {
  rbc_37_poly_set_monomial(o, degree);
}




/**
 * \fn void rbc_37_qre_rot(rbc_37_qre v, uint32_t r)
 * \brief This functions rotates entries of a vector with coefficients in GF(q^m).
 *
 * \param[out] v rbc_37_qre with entries rotated by position
 * \param[in] position Integer between 1 and (size - 1) such that v is rotated-right r times
 */
void rbc_37_qre_rot(rbc_37_qre o, rbc_37_qre e, uint32_t rot) {

  if(rot > (uint32_t)(e->degree + 1)) {
    exit(EXIT_FAILURE);
  }

  if (rot == (uint32_t)(e->degree + 1)) {
    rbc_37_poly_set(o, e);
  }
  else {
    rbc_37_qre tmp;
    rbc_37_qre_init(&tmp);
    rbc_37_qre_set_monomial(tmp, rot);
    rbc_37_qre_mul(o, e, tmp);
    rbc_37_qre_clear(tmp);
  }
}


