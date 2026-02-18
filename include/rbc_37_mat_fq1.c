#include <rbc_37_mat_fq.h>
#include <rbc_37_qre.h>

#include <rbc_61_mat_fq.h>
#include <rbc_61_qre.h>

#include <rbc_79_mat_fq.h>
#include <rbc_79_qre.h>
#include <time.h>

void rbc_37_mat_fq_from_qre(rbc_37_mat_fq o, rbc_37_qre m){
    uint64_t m_ij;
    rbc_37_mat_fq_set_zero(o, RBC_37_FIELD_M, m->max_degree);
    for(int i = 0 ; i < RBC_37_FIELD_M; i++){
        for(int j = 0 ; j <= m->max_degree; j++){
            m_ij = (m->v[i][j / 64] >> (j % 64)) & 1;
            o[i][j / 64] ^= m_ij << (j % 64);
        }
    }
}

void rbc_37_mat_fq_to_qre(rbc_37_qre o, rbc_37_mat_fq m){
    uint64_t m_ij;
    rbc_37_qre_set_zero(o);
    for(int i = 0 ; i < RBC_37_FIELD_M; i++){
        for(int j = 0; j <= o->max_degree; j++){
            m_ij = (m[i][j / 64] >> (j % 64)) & 1;
            o->v[i][j / 64] ^= m_ij << (j % 64);
        }
    }
}

//don't use with o = m1 or m2
void rbc_37_mat_fq_add(rbc_37_mat_fq o, rbc_37_mat_fq m1, rbc_37_mat_fq m2, size_t rows, size_t columns){
    uint32_t words = (columns + 63) / 64;
    rbc_37_mat_fq_set_zero(o, rows, columns);
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < words; j++){
            o[i][j] = m1[i][j] ^ m2[i][j];
            //printf("o  = %lx\nm1 = %lx\nm2 = %lx\n", o[i][j], m1[i][j], m2[i][j]);
        }
    }
}


void rbc_61_mat_fq_from_qre(rbc_61_mat_fq o, rbc_61_qre m){
    uint64_t m_ij;
    rbc_61_mat_fq_set_zero(o, RBC_61_FIELD_M, m->max_degree);
    for(int i = 0 ; i < RBC_61_FIELD_M; i++){
        for(int j = 0 ; j <= m->max_degree; j++){
            m_ij = (m->v[i][j / 64] >> (j % 64)) & 1;
            o[i][j / 64] ^= m_ij << (j % 64);
        }
    }
}

void rbc_61_mat_fq_to_qre(rbc_61_qre o, rbc_61_mat_fq m){
    uint64_t m_ij;
    rbc_61_qre_set_zero(o);
    for(int i = 0 ; i < RBC_61_FIELD_M; i++){
        for(int j = 0; j <= o->max_degree; j++){
            m_ij = (m[i][j / 64] >> (j % 64)) & 1;
            o->v[i][j / 64] ^= m_ij << (j % 64);
        }
    }
}

//don't use with o = m1 or m2
void rbc_61_mat_fq_add(rbc_61_mat_fq o, rbc_61_mat_fq m1, rbc_61_mat_fq m2, size_t rows, size_t columns){
    uint32_t words = (columns + 63) / 64;
    rbc_61_mat_fq_set_zero(o, rows, columns);
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < words; j++){
            o[i][j] = m1[i][j] ^ m2[i][j];
            //printf("o  = %lx\nm1 = %lx\nm2 = %lx\n", o[i][j], m1[i][j], m2[i][j]);
        }
    }
}



void rbc_79_mat_fq_from_qre(rbc_79_mat_fq o, rbc_79_qre m){
    uint64_t m_ij;
    rbc_79_mat_fq_set_zero(o, RBC_79_FIELD_M, m->max_degree);
    for(int i = 0 ; i < RBC_79_FIELD_M; i++){
        for(int j = 0 ; j <= m->max_degree; j++){
            m_ij = (m->v[i][j / 64] >> (j % 64)) & 1;
            o[i][j / 64] ^= m_ij << (j % 64);
        }
    }
}

void rbc_79_mat_fq_to_qre(rbc_79_qre o, rbc_79_mat_fq m){
    uint64_t m_ij;
    rbc_79_qre_set_zero(o);
    for(int i = 0 ; i < RBC_79_FIELD_M; i++){
        for(int j = 0; j <= o->max_degree; j++){
            m_ij = (m[i][j / 64] >> (j % 64)) & 1;
            o->v[i][j / 64] ^= m_ij << (j % 64);
        }
    }
}

//don't use with o = m1 or m2
void rbc_79_mat_fq_add(rbc_79_mat_fq o, rbc_79_mat_fq m1, rbc_79_mat_fq m2, size_t rows, size_t columns){
    uint32_t words = (columns + 63) / 64;
    rbc_79_mat_fq_set_zero(o, rows, columns);
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < words; j++){
            o[i][j] = m1[i][j] ^ m2[i][j];
            //printf("o  = %lx\nm1 = %lx\nm2 = %lx\n", o[i][j], m1[i][j], m2[i][j]);
        }
    }
}



//Given a matrix m1 and a matrix m2, multiply m1 with transpose(m2). This is done to obtain a more efficient matrix multiplication.
// Otherwise, to multiply m1 and m2, this would require accessing the bits in a much less practical way.

 void multiplication(rbc_79_mat_fq o, const rbc_79_mat_fq m1, const rbc_79_mat_fq m2, uint32_t rows1, uint32_t columns1_rows2, uint32_t columns2) {
  uint64_t acc;
 
  rbc_79_mat_fq_set_zero(o, rows1, columns2);
  for(size_t i = 0 ; i < rows1 ; ++i) {
    for(size_t j = 0 ; j < columns2 ; ++j) {
      acc = 0;
      for(size_t k = 0 ; k < ((columns1_rows2 + 63) / 64) ; ++k) {
        acc ^= (m1[i][k] & m2[j][k]);
      }
      o[i][j / 64] ^= ((uint64_t)__builtin_popcountll(acc) & 0x1) << (j % 64);
    }
  }
}


// Given nbmat matrices, flatten them and then put them row-wise in a matrix o.
void vectorization(rbc_79_mat_fq o, rbc_79_mat_fq *m,uint32_t nbmat, uint32_t rows, uint32_t columns) {
  uint64_t m_ij;
  rbc_79_mat_fq_set_zero(o, nbmat, rows*columns);
  for(size_t ell = 0 ; ell < nbmat ; ++ell) {
  for(size_t i = 0 ; i < rows ; ++i) {
    for(size_t j= 0 ; j<columns ; ++j) {
        m_ij = (uint8_t)((m[ell][i][j / 64] >> (j % 64)) & 1);
        o[ell][(i*columns+j)/64] ^= m_ij<<((i*columns+j)%64);
    }
}
}
}

// Addition of two matrices
void addmat(rbc_79_mat_fq o, rbc_79_mat_fq m1,rbc_79_mat_fq m2, uint32_t rows, uint32_t columns) {
uint64_t acc, m1_ij, m2_ij;
  rbc_79_mat_fq_set_zero(o, rows,columns);
uint32_t words = (columns + 63) / 64;
for(size_t i = 0 ; i < rows ; ++i) {
    for(size_t j = 0 ; j < words ; ++j) {
        m1_ij = m1[i][j];
        m2_ij = m2[i][j];
        o[i][j] = (uint64_t)(m1_ij^m2_ij);
  }
}
}


// Given two matrices o and m1, recover, in o, the matrix [m1|0].
void concatene(rbc_79_mat_fq o, rbc_79_mat_fq m1, uint32_t rows1, uint32_t columns1, uint32_t columns2) {

uint64_t acc, m1_ij, m2_ij;
  rbc_79_mat_fq_set_zero(o, rows1,columns1+columns2);

  for(size_t i = 0 ; i < rows1 ; ++i) {
    for(size_t j = 0 ; j < columns1 ; ++j) {
        m1_ij = (uint8_t)((m1[i][j / 64] >> (j % 64)) & 1);
        o[i][j/64] ^=  (uint64_t)(m1_ij^m2_ij) << (j % 64);
  }
}



}
// Given two matrices o and m1, recover, in o, the matrix [m1|0]. Same as previous, but better
void concatene2(rbc_79_mat_fq o, rbc_79_mat_fq m1, uint32_t rows1, uint32_t columns1, uint32_t columns2) {

uint64_t acc, m1_ij, m2_ij;
  rbc_79_mat_fq_set_zero(o, rows1,columns1+columns2);
uint32_t words = (columns1 + 63) / 64;
  for(size_t i = 0 ; i < rows1 ; ++i) {
    for(size_t j = 0 ; j < words ; ++j) {
        m1_ij = m1[i][j];
        o[i][j] ^=  m1_ij;
  }
}
}



// Given one matrices m1, recover o1 and o2 such that m1 = [o1|o2]

void truncation(rbc_79_mat_fq o1,rbc_79_mat_fq o2, rbc_79_mat_fq m1, uint32_t rows1, uint32_t columns1, uint32_t columns2) {//
uint64_t acc, m1_ij, m2_ij;

rbc_79_mat_fq_set_zero(o1,rows1,columns1);

rbc_79_mat_fq_set_zero(o2,rows1,columns2-columns1);
uint64_t mask;
 uint64_t coefficient;
size_t to_where = columns1/64;
size_t adjust = columns1%64;
size_t words = (columns2)/64;

  for(size_t i = 0 ; i < rows1 ; ++i) {
    for(size_t j = 0 ; j < to_where ; ++j) {
        m1_ij = m1[i][j];
        o1[i][j] ^= m1_ij;
  }
 coefficient = m1[i][to_where];
 mask = (1ULL << adjust) - 1;
o1[i][to_where] ^= coefficient & mask;


}
if(to_where == words){
  for(size_t i = 0 ; i < rows1 ; ++i) {
 
for(size_t j = columns1 ; j<(columns2) ; ++j) {
        uint64_t coefficient = m1[i][j/64];
        m1_ij = (( coefficient >> (j % 64)) & 1);
        o2[i][(j-columns1)/64] ^=  m1_ij<<((j-columns1)%64);
    }
  }
}
else{
for(size_t i = 0 ; i < rows1 ; ++i) {
for(size_t j = to_where ; j<(words-1) ; ++j) {
    coefficient = m1[i][j] >> adjust;
    o2[i][j-to_where] ^= coefficient;
    coefficient = m1[i][j+1] << (64-adjust);
    o2[i][j-to_where] ^= coefficient;
}

 for(size_t j = (words-1)*64+adjust ; j<(columns2) ; ++j) {
         uint64_t coefficient = m1[i][j/64];
         m1_ij = (( coefficient >> (j % 64)) & 1);
        o2[i][(j-columns1)/64] ^=  m1_ij<<((j-columns1)%64);
    }

}
}
}