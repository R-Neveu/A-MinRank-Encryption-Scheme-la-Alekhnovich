/**
 * \file rbc_79_mat.h
 * \brief Interface for matrices over Fq^m
 */

#ifndef RBC_79_MAT_H
#define RBC_79_MAT_H

#include "rbc_79.h"
#include "rbc_79_elt.h"
#include "rbc_79_vec.h"
#include "seedexpander.h"
#include "seedexpander_shake.h"

void rbc_79_mat_init(rbc_79_mat* m, uint32_t rows, uint32_t columns);
void rbc_79_mat_clear(rbc_79_mat m);
void rbc_79_mat_set_zero(rbc_79_mat m, uint32_t rows, uint32_t columns);
void rbc_79_mat_set(rbc_79_mat o, const rbc_79_mat m, uint32_t rows, uint32_t columns);
void rbc_79_mat_set_random(random_source* ctx, rbc_79_mat o, uint32_t rows, uint32_t columns);
void rbc_79_mat_set_random_tmp(seedexpander_shake_t* ctx, rbc_79_mat o, uint32_t rows, uint32_t columns);
void rbc_79_mat_add(rbc_79_mat o, const rbc_79_mat m1, const rbc_79_mat m2, uint32_t rows, uint32_t columns);
void rbc_79_mat_mul(rbc_79_mat o, const rbc_79_mat m1, const rbc_79_mat m2, uint32_t rows1, uint32_t columns1_rows2, uint32_t columns2);
void rbc_79_mat_vec_mul(rbc_79_vec o, const rbc_79_mat m, const rbc_79_vec v, uint32_t rows, uint32_t columns);
void rbc_79_vec_mat_mul(rbc_79_vec o, const rbc_79_mat m, const rbc_79_vec v, uint32_t rows, uint32_t columns);
void rbc_79_mat_to_string(uint8_t* str, const rbc_79_mat m, uint32_t rows, uint32_t columns);
void rbc_79_mat_from_string(rbc_79_mat m, uint32_t rows, uint32_t columns, const uint8_t* str);
void rbc_79_mat_print(const rbc_79_mat m, uint32_t rows, uint32_t columns);
#endif

