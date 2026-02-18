/**
 * \file rbc_37_isometry.c
 * \brief Implementation of rbc_37_isometry.h
 */

#include "rbc_37_isometry.h"




/**
 * \fn void rbc_37_isometry_init(rbc_37_isometry* o, uint32_t n)
 * \brief This function allocates the memory for a rbc_37_isometry element.
 *
 * \param[out] o Pointer to the allocated rbc_37_isometry
 * \param[in] n Dimension size of the rbc_37_isometry
 */
void rbc_37_isometry_init(rbc_37_isometry* o, uint32_t n) {
  rbc_37_mat_fq_init(&(o->P), RBC_37_FIELD_M, RBC_37_FIELD_M);
  rbc_37_mat_fq_init(&(o->Q), n, n);
  o->n = n;
}




/**
 * \fn void rbc_37_isometry_clear(rbc_37_isometry* o)
 * \brief This function clear the memory for an rbc_37_isometry element.
 *
 * \param[in] o Pointer to a rbc_37_isometry
 */
void rbc_37_isometry_clear(rbc_37_isometry* o) {
  rbc_37_mat_fq_clear(o->P);
  rbc_37_mat_fq_clear(o->Q);
  o->n = 0;
}




/**
 * \fn void rbc_37_isometry_set_zero(rbc_37_isometry *o, uint32_t n)
 * \brief This function copies a pair matrices over GF(q) to an rbc_37_isometry.
 *
 * \param[out] o Pointer to the allocated rbc_37_isometry
 * \param[in] P rbc_37_mat_fq of dimension m
 * \param[in] Q rbc_37_mat_fq of dimension n
 * \param[in] n Dimension size of the rbc_37_isometry
 */
void rbc_37_isometry_set_zero(rbc_37_isometry *o, uint32_t n) {
  if (n != o->n) exit(EXIT_FAILURE);
  rbc_37_mat_fq_set_zero(o->P, RBC_37_FIELD_M, RBC_37_FIELD_M);
  rbc_37_mat_fq_set_zero(o->Q, n, n);
}




/**
 * \fn void rbc_37_isometry_set(rbc_37_isometry *o, const rbc_37_mat_fq P, const rbc_37_mat_fq Q, uint32_t n)
 * \brief This function copies a pair matrices over GF(q) to an rbc_37_isometry.
 *
 * \param[out] o Pointer to the allocated rbc_37_isometry
 * \param[in] P rbc_37_mat_fq of dimension m
 * \param[in] Q rbc_37_mat_fq of dimension n
 * \param[in] n Dimension size of the rbc_37_isometry
 */
void rbc_37_isometry_set(rbc_37_isometry *o, const rbc_37_mat_fq P, const rbc_37_mat_fq Q, uint32_t n) {
  if (n != o->n) exit(EXIT_FAILURE);
  rbc_37_mat_fq_set(o->P, P, RBC_37_FIELD_M, RBC_37_FIELD_M);
  rbc_37_mat_fq_set(o->Q, Q, n, n);
}




/**
 * \fn void rbc_37_isometry_set_random(random_source* ctx, rbc_37_isometry *o, uint32_t n)
 * \brief This function copies a pair matrices over GF(q) to an rbc_37_isometry.
 *
 * \param[out] ctx NIST seed expander
 * \param[out] o Pointer to the allocated rbc_37_isometry
 * \param[in] n Dimension size of the rbc_37_isometry
 */
void rbc_37_isometry_set_random(random_source* ctx, rbc_37_isometry *o, uint32_t n) {
  if (n != o->n) exit(EXIT_FAILURE);
  rbc_37_perm perm;
  rbc_37_perm_init(&perm, RBC_37_FIELD_M);
  rbc_37_mat_fq_gen(ctx, perm, o->P, RBC_37_FIELD_M);
  rbc_37_perm_clear(perm);
  rbc_37_perm_init(&perm, n);
  rbc_37_mat_fq_gen(ctx, perm, o->Q, n);
  rbc_37_perm_clear(perm);
}




/**
 * \fn void rbc_37_isometry_set_random_from_xof(rbc_37_isometry *o, uint32_t n, void (*xof)(uint8_t *, size_t, const uint8_t *, size_t), const uint8_t *xof_input, uint32_t xof_size)
 * \brief This function sets a pair of rbc_37_isometry with random values using XOF.
 *
 * \param[out] o Pointer to the allocated rbc_37_isometry
 * \param[in] n Dimension size of the rbc_37_isometry
 * \param[in] xof XOF procedure
 * \param[in] xof_input XOF input byte array
 * \param[in] xof_size XOF input byte array length of xof_input
 */
void rbc_37_isometry_set_random_from_xof(rbc_37_isometry *o, uint32_t n, void (*xof)(uint8_t *, size_t, const uint8_t *, size_t), const uint8_t *xof_input, uint32_t xof_size) {
  if (n != o->n) exit(EXIT_FAILURE);
  rbc_37_perm perm;
  rbc_37_perm_init(&perm, RBC_37_FIELD_M);
  rbc_37_mat_fq_gen_from_xof(perm, o->P, RBC_37_FIELD_M, xof, xof_input, xof_size);
  rbc_37_perm_clear(perm);
  rbc_37_perm_init(&perm, n);
  rbc_37_mat_fq_gen_from_xof(perm, o->Q, n, xof, &xof_input[xof_size], xof_size);
  rbc_37_perm_clear(perm);
}




/**
 * \fn void rbc_37_isometry_apply(rbc_37_vec o, const rbc_37_vec v, uint32_t n, const rbc_37_isometry *iso)
 * \brief This function applies a rbc_37_isometry element on a rbc_37_vec element v. Here, v^T is required as in theory being an mxn matrix over GF(q).
 *
 * \param[in] o rbc_37_vector equals \f$ Q^T \times v \times P^T \f$
 * \param[in] v rbc_37_vector
 * \param[in] iso Pointer to a rbc_37_isometry
 */
void rbc_37_isometry_apply(rbc_37_vec o, const rbc_37_vec v, uint32_t n, const rbc_37_isometry *iso) {
  if (n != iso->n) exit(EXIT_FAILURE);
  rbc_37_vec t;
  rbc_37_vec_init(&t, n);
  rbc_37_tr_mat_fq_vec_mul(t, iso->Q, v, n);
  rbc_37_vec_mat_fq_tr_mul(o, t, iso->P, n);
  rbc_37_vec_clear(t);
}




/**
 * \fn void rbc_37_isometry_compose(rbc_37_isometry *o, const rbc_37_isometry *iso1, const rbc_37_isometry *iso2)
 * \brief This function compose two rbc_37_isometry elements.
 *
 * \param[out] o Pointer to the rbc_37_isometry equals \f$ iso1 \circ iso2 \f$
 * \param[in] iso1 Pointer to a rbc_37_isometry
 * \param[in] iso2 Pointer to a rbc_37_isometry
 */
void rbc_37_isometry_compose(rbc_37_isometry *o, const rbc_37_isometry *iso1, const rbc_37_isometry *iso2) {
  if (iso1->n != iso2->n) exit(EXIT_FAILURE);
  uint32_t n = iso1->n;
  rbc_37_mat_fq t1, t2;
  rbc_37_mat_fq_init(&t1, RBC_37_FIELD_M, RBC_37_FIELD_M);
  rbc_37_mat_fq_init(&t2, RBC_37_FIELD_M, RBC_37_FIELD_M);
  rbc_37_mat_fq_set(t1, iso1->P, RBC_37_FIELD_M, RBC_37_FIELD_M);
  rbc_37_mat_fq_set(t2, iso2->P, RBC_37_FIELD_M, RBC_37_FIELD_M);
  rbc_37_mat_fq_mul(o->P, t2, t1, RBC_37_FIELD_M, RBC_37_FIELD_M, RBC_37_FIELD_M);
  rbc_37_mat_fq_clear(t1);
  rbc_37_mat_fq_clear(t2);
  rbc_37_mat_fq_init(&t1, n, n);
  rbc_37_mat_fq_init(&t2, n, n);
  rbc_37_mat_fq_set(t1, iso1->Q, n, n);
  rbc_37_mat_fq_set(t2, iso2->Q, n, n);
  rbc_37_mat_fq_mul(o->Q, t1, t2, n, n, n);
  rbc_37_mat_fq_clear(t1);
  rbc_37_mat_fq_clear(t2);
}

