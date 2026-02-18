#ifndef RBC_37_MAT_FQ1_H
#define RBC_37_MAT_FQ1_H

#include <rbc_37_mat_fq.h>
#include <rbc_37_qre.h>
#include <time.h>


void rbc_37_mat_fq_from_qre(rbc_37_mat_fq o, rbc_37_qre m);

void rbc_37_mat_fq_to_qre(rbc_37_qre o, rbc_37_mat_fq m);

void rbc_37_mat_fq_add(rbc_37_mat_fq o, rbc_37_mat_fq m1, rbc_37_mat_fq m2, size_t rows, size_t columns);

#endif


#ifndef RBC_61_MAT_FQ1_H
#define RBC_61_MAT_FQ1_H

#include <rbc_61_mat_fq.h>
#include <rbc_61_qre.h>

void rbc_61_mat_fq_from_qre(rbc_61_mat_fq o, rbc_61_qre m);

void rbc_61_mat_fq_to_qre(rbc_61_qre o, rbc_61_mat_fq m);

void rbc_61_mat_fq_add(rbc_61_mat_fq o, rbc_61_mat_fq m1, rbc_61_mat_fq m2, size_t rows, size_t columns);

#endif

#ifndef RBC_79_MAT_FQ1_H
#define RBC_79_MAT_FQ1_H

#include <rbc_79_mat_fq.h>
#include <rbc_79_qre.h>

void rbc_79_mat_fq_from_qre(rbc_79_mat_fq o, rbc_79_qre m);

void rbc_79_mat_fq_to_qre(rbc_79_qre o, rbc_79_mat_fq m);

void rbc_79_mat_fq_add(rbc_79_mat_fq o, rbc_79_mat_fq m1, rbc_79_mat_fq m2, size_t rows, size_t columns);

#endif

void vectorization(rbc_79_mat_fq o, rbc_79_mat_fq *m,uint32_t nbmat, uint32_t rows, uint32_t columns);
void addmat(rbc_79_mat_fq o, rbc_79_mat_fq m1,rbc_79_mat_fq m2, uint32_t rows, uint32_t columns);
void concatene2(rbc_79_mat_fq o, rbc_79_mat_fq m1, uint32_t rows1, uint32_t columns1, uint32_t columns2);
void concatene(rbc_79_mat_fq o, rbc_79_mat_fq m1, uint32_t rows1, uint32_t columns1, uint32_t columns2);
void truncation(rbc_79_mat_fq o1,rbc_79_mat_fq o2, rbc_79_mat_fq m1, uint32_t rows1, uint32_t columns1, uint32_t columns2);
void multiplication(rbc_79_mat_fq o, const rbc_79_mat_fq m1, const rbc_79_mat_fq m2, uint32_t rows1, uint32_t columns1_rows2, uint32_t columns2);

static long long cpucycles() {

  unsigned long long result;

  __asm__ volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax" : "=a" (result) ::  "%rdx");

  return result;

}
 double now_ms() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1.0e6;

}
