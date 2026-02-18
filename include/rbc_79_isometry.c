/**
 * \file rbc_79_isometry.c
 * \brief Implementation of rbc_79_isometry.h
 */

#include "rbc_79_isometry.h"




/**
 * \fn void rbc_79_isometry_init(rbc_79_isometry* o, uint32_t n)
 * \brief This function allocates the memory for a rbc_79_isometry element.
 *
 * \param[out] o Pointer to the allocated rbc_79_isometry
 * \param[in] n Dimension size of the rbc_79_isometry
 */
void rbc_79_isometry_init(rbc_79_isometry* o, uint32_t n) {
  rbc_79_mat_fq_init(&(o->P), RBC_79_FIELD_M, RBC_79_FIELD_M);
  rbc_79_mat_fq_init(&(o->Q), n, n);
  o->n = n;
}




/**
 * \fn void rbc_79_isometry_clear(rbc_79_isometry* o)
 * \brief This function clear the memory for an rbc_79_isometry element.
 *
 * \param[in] o Pointer to a rbc_79_isometry
 */
void rbc_79_isometry_clear(rbc_79_isometry* o) {
  rbc_79_mat_fq_clear(o->P);
  rbc_79_mat_fq_clear(o->Q);
  o->n = 0;
}




/**
 * \fn void rbc_79_isometry_set_zero(rbc_79_isometry *o, uint32_t n)
 * \brief This function copies a pair matrices over GF(q) to an rbc_79_isometry.
 *
 * \param[out] o Pointer to the allocated rbc_79_isometry
 * \param[in] P rbc_79_mat_fq of dimension m
 * \param[in] Q rbc_79_mat_fq of dimension n
 * \param[in] n Dimension size of the rbc_79_isometry
 */
void rbc_79_isometry_set_zero(rbc_79_isometry *o, uint32_t n) {
  if (n != o->n) exit(EXIT_FAILURE);
  rbc_79_mat_fq_set_zero(o->P, RBC_79_FIELD_M, RBC_79_FIELD_M);
  rbc_79_mat_fq_set_zero(o->Q, n, n);
}




/**
 * \fn void rbc_79_isometry_set(rbc_79_isometry *o, const rbc_79_mat_fq P, const rbc_79_mat_fq Q, uint32_t n)
 * \brief This function copies a pair matrices over GF(q) to an rbc_79_isometry.
 *
 * \param[out] o Pointer to the allocated rbc_79_isometry
 * \param[in] P rbc_79_mat_fq of dimension m
 * \param[in] Q rbc_79_mat_fq of dimension n
 * \param[in] n Dimension size of the rbc_79_isometry
 */
void rbc_79_isometry_set(rbc_79_isometry *o, const rbc_79_mat_fq P, const rbc_79_mat_fq Q, uint32_t n) {
  if (n != o->n) exit(EXIT_FAILURE);
  rbc_79_mat_fq_set(o->P, P, RBC_79_FIELD_M, RBC_79_FIELD_M);
  rbc_79_mat_fq_set(o->Q, Q, n, n);
}




/**
 * \fn void rbc_79_isometry_set_random(random_source* ctx, rbc_79_isometry *o, uint32_t n)
 * \brief This function copies a pair matrices over GF(q) to an rbc_79_isometry.
 *
 * \param[out] ctx NIST seed expander
 * \param[out] o Pointer to the allocated rbc_79_isometry
 * \param[in] n Dimension size of the rbc_79_isometry
 */
void rbc_79_isometry_set_random(random_source* ctx, rbc_79_isometry *o, uint32_t n) {
  if (n != o->n) exit(EXIT_FAILURE);
  rbc_79_perm perm;
  rbc_79_perm_init(&perm, RBC_79_FIELD_M);
  rbc_79_mat_fq_gen(ctx, perm, o->P, RBC_79_FIELD_M);
  rbc_79_perm_clear(perm);
  rbc_79_perm_init(&perm, n);
  rbc_79_mat_fq_gen(ctx, perm, o->Q, n);
  rbc_79_perm_clear(perm);
}




/**
 * \fn void rbc_79_isometry_set_random_from_xof(rbc_79_isometry *o, uint32_t n, void (*xof)(uint8_t *, size_t, const uint8_t *, size_t), const uint8_t *xof_input, uint32_t xof_size)
 * \brief This function sets a pair of rbc_79_isometry with random values using XOF.
 *
 * \param[out] o Pointer to the allocated rbc_79_isometry
 * \param[in] n Dimension size of the rbc_79_isometry
 * \param[in] xof XOF procedure
 * \param[in] xof_input XOF input byte array
 * \param[in] xof_size XOF input byte array length of xof_input
 */
void rbc_79_isometry_set_random_from_xof(rbc_79_isometry *o, uint32_t n, void (*xof)(uint8_t *, size_t, const uint8_t *, size_t), const uint8_t *xof_input, uint32_t xof_size) {
  if (n != o->n) exit(EXIT_FAILURE);
  rbc_79_perm perm;
  rbc_79_perm_init(&perm, RBC_79_FIELD_M);
  rbc_79_mat_fq_gen_from_xof(perm, o->P, RBC_79_FIELD_M, xof, xof_input, xof_size);
  rbc_79_perm_clear(perm);
  rbc_79_perm_init(&perm, n);
  rbc_79_mat_fq_gen_from_xof(perm, o->Q, n, xof, &xof_input[xof_size], xof_size);
  rbc_79_perm_clear(perm);
}




/**
 * \fn void rbc_79_isometry_apply(rbc_79_vec o, const rbc_79_vec v, uint32_t n, const rbc_79_isometry *iso)
 * \brief This function applies a rbc_79_isometry element on a rbc_79_vec element v. Here, v^T is required as in theory being an mxn matrix over GF(q).
 *
 * \param[in] o rbc_79_vector equals \f$ Q^T \times v \times P^T \f$
 * \param[in] v rbc_79_vector
 * \param[in] iso Pointer to a rbc_79_isometry
 */
void rbc_79_isometry_apply(rbc_79_vec o, const rbc_79_vec v, uint32_t n, const rbc_79_isometry *iso) {
  if (n != iso->n) exit(EXIT_FAILURE);
  rbc_79_vec t;
  rbc_79_vec_init(&t, n);
  rbc_79_tr_mat_fq_vec_mul(t, iso->Q, v, n);
  rbc_79_vec_mat_fq_tr_mul(o, t, iso->P, n);
  rbc_79_vec_clear(t);
}




/**
 * \fn void rbc_79_isometry_compose(rbc_79_isometry *o, const rbc_79_isometry *iso1, const rbc_79_isometry *iso2)
 * \brief This function compose two rbc_79_isometry elements.
 *
 * \param[out] o Pointer to the rbc_79_isometry equals \f$ iso1 \circ iso2 \f$
 * \param[in] iso1 Pointer to a rbc_79_isometry
 * \param[in] iso2 Pointer to a rbc_79_isometry
 */
void rbc_79_isometry_compose(rbc_79_isometry *o, const rbc_79_isometry *iso1, const rbc_79_isometry *iso2) {
  if (iso1->n != iso2->n) exit(EXIT_FAILURE);
  uint32_t n = iso1->n;
  rbc_79_mat_fq t1, t2;
  rbc_79_mat_fq_init(&t1, RBC_79_FIELD_M, RBC_79_FIELD_M);
  rbc_79_mat_fq_init(&t2, RBC_79_FIELD_M, RBC_79_FIELD_M);
  rbc_79_mat_fq_set(t1, iso1->P, RBC_79_FIELD_M, RBC_79_FIELD_M);
  rbc_79_mat_fq_set(t2, iso2->P, RBC_79_FIELD_M, RBC_79_FIELD_M);
  rbc_79_mat_fq_mul(o->P, t2, t1, RBC_79_FIELD_M, RBC_79_FIELD_M, RBC_79_FIELD_M);
  rbc_79_mat_fq_clear(t1);
  rbc_79_mat_fq_clear(t2);
  rbc_79_mat_fq_init(&t1, n, n);
  rbc_79_mat_fq_init(&t2, n, n);
  rbc_79_mat_fq_set(t1, iso1->Q, n, n);
  rbc_79_mat_fq_set(t2, iso2->Q, n, n);
  rbc_79_mat_fq_mul(o->Q, t1, t2, n, n, n);
  rbc_79_mat_fq_clear(t1);
  rbc_79_mat_fq_clear(t2);
}

