#include <rbc.h>
#include <rbc_79.h>

#include <rbc_79_mat_fq.h>

#include <random_source.h>
#include <rbc_79_elt.h>
#include <rbc_79_vec.h>
#include <rbc_79_qre.h>

#include <rbc_37_mat_fq1.h>
 
#include <stdio.h>
#include <stdlib.h>
 
#define R 6

#define D 6

#define N 115

#define M 115

#define MN M*N

#define K 6613

#define ELL1 79

#define NUM_OF_TEST 300


int main(){
   rbc_79_mat_fq S,SE,H,V,pk,Vp,E,pkD,A,ct1,ct2,E_k,Ht,A1,A2,W,E1,E2,Spoub,pkD1,U,Wf,EFt;

    rbc_79_mat_fq *P = malloc(ELL1 * sizeof(rbc_79_mat_fq));
    rbc_79_mat_fq *Pa = malloc(ELL1 * sizeof(rbc_79_mat_fq));
    rbc_79_mat_fq *Ep = malloc(ELL1 * sizeof(rbc_79_mat_fq));
    rbc_79_mat_fq *Ap = malloc(ELL1 * sizeof(rbc_79_mat_fq));
 rbc_79_gabidulin code,code_place;
  rbc_79_mat_fq T, C, V1;
  rbc_79_qre c, g,c1,g_place;
  rbc_79_vec m,m1;

  for (int i = 0; i < ELL1; i++) {
        rbc_79_mat_fq_init(&P[i], N, R);  
        rbc_79_mat_fq_init(&Pa[i], M, D);   
        rbc_79_mat_fq_init(&Ep[i], M, N);  
        rbc_79_mat_fq_init(&Ap[i], M, N);  
    }
    random_source prng;

    random_source sk_seedexpander;
    random_source pk_seedexpander;
    random_source enc_seedexpander;

    random_source_init(&sk_seedexpander, RANDOM_SOURCE_SEEDEXP);
    random_source_init(&pk_seedexpander, RANDOM_SOURCE_SEEDEXP);
    random_source_init(&enc_seedexpander, RANDOM_SOURCE_SEEDEXP);

    uint8_t sk_seed[40] = {0};
    uint8_t pk_seed[40] = {0};
    uint8_t enc_seed[40] = {0};

    long long t0 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5=0, t01=0, t02=0,t_H=0;
 rbc_79_field_init();

  rbc_79_qre_init_modulus(ELL1);
  rbc_79_mat_fq_init(&T, ELL1, ELL1);
  rbc_79_mat_fq_init(&C, ELL1, ELL1);
  rbc_79_mat_fq_init(&V1, ELL1, ELL1);
  rbc_79_qre_init(&g);
    rbc_79_qre_init(&g_place);

  rbc_79_qre_init(&c);
  rbc_79_qre_init(&c1);
  rbc_79_vec_init(&m, 3);
  rbc_79_vec_init(&m1, 3);


    rbc_79_mat_fq_init(&H, MN-K, K);

    rbc_79_mat_fq_init(&V, M, R);

    rbc_79_mat_fq_init(&Ht, K, MN-K);

  
    rbc_79_mat_fq_init(&S, ELL1, K);

    rbc_79_mat_fq_init(&SE, ELL1, K);
    
    rbc_79_mat_fq_init(&pk, ELL1, MN);

    rbc_79_mat_fq_init(&Spoub, ELL1, MN-K);

   rbc_79_mat_fq_init(&E_k, ELL1, K);

    rbc_79_mat_fq_init(&pkD,ELL1,MN-K);

    rbc_79_mat_fq_init(&pkD1,ELL1,MN-K);

    rbc_79_mat_fq_init(&Vp, N, D);

    rbc_79_mat_fq_init(&E, ELL1, MN);

    rbc_79_mat_fq_init(&E1, ELL1, MN-K);

    rbc_79_mat_fq_init(&E2, ELL1, K);

    rbc_79_mat_fq_init(&A, ELL1, MN);

    rbc_79_mat_fq_init(&A1, ELL1, MN-K);

    rbc_79_mat_fq_init(&A2, ELL1, K);

    rbc_79_mat_fq_init(&ct1, ELL1,K);

    rbc_79_mat_fq_init(&U, ELL1,K);

    rbc_79_mat_fq_init(&ct2, ELL1, ELL1);

    rbc_79_mat_fq_init(&W, ELL1, ELL1);

    rbc_79_mat_fq_init(&Wf, ELL1, ELL1);

    rbc_79_mat_fq_init(&EFt, ELL1, ELL1);


    random_source_init(&prng, RANDOM_SOURCE_PRNG);

double start,end;
  double dif,dif1;

    for(int i=0; i < NUM_OF_TEST; i++){
      t0 = now_ms();

        random_source_get_bytes(&prng, sk_seed, 40);
        random_source_get_bytes(&prng, pk_seed, 40); 
        random_source_seed(&sk_seedexpander, sk_seed);
        random_source_seed(&pk_seedexpander, pk_seed);

        rbc_79_mat_fq_set_random(&pk_seedexpander, H, MN-K, K); // Sample H a random MN-K times K matrix (we assume standard form)
        
        t_H += now_ms() -t0;

        t0 = now_ms();

        rbc_79_mat_fq_set_random(&sk_seedexpander, V, M, R); // Sample V a random M times R matrix (this is the column support)

        for (int j = 0; j < ELL1; j++) {
        rbc_79_mat_fq_set_random(&sk_seedexpander, P[j], N, R);         // Sample ELL1 R times N matrices P_j

        }
        rbc_79_qre_set_random(&pk_seedexpander, g);
        rbc_79_gabidulin_init(&code, g, 3, ELL1);

        t1 += now_ms() - t0;

        t0 = now_ms();

        for (int j = 0; j < ELL1; j++) {

         multiplication(Ep[j], V, P[j], M, R, N); // Compute ELL1 M times N matrices Ep_j, of rank R, of column support V

        }

        vectorization(E,Ep,ELL1,M,N); // Set E as the matrix whose ith line is Ep_i (written as a line)
    
        truncation(E1,E2,E,ELL1,MN-K,MN); 

       multiplication(pkD, E2,H,ELL1,K,MN-K); 

        addmat(pkD1,E1,pkD,ELL1,MN-K); 
 
        t2 += now_ms()-t0;

        t0 = now_ms();

        rbc_79_mat_fq_set_transpose(Ht,H,MN-K,K); // Transpose H. In practice, this would be generated directly in the step just above instead of transposed (as Ht would have already been generated thanks to the seed in its transposed form).
        
        t01 += now_ms()-t0;

        t0 = now_ms();

        random_source_get_bytes(&prng, enc_seed, 40);
  
        random_source_seed(&enc_seedexpander, enc_seed);

        for (int j = 0; j < ELL1; j++) {
        rbc_79_mat_fq_set_random(&enc_seedexpander, Pa[j], M, D); 
        }

        t02 += now_ms()-t0;

        t0= now_ms();

        
        rbc_79_vec_set_random(&enc_seedexpander, m, 3);
        rbc_79_gabidulin_encode(c, code, m);

        rbc_79_mat_fq_from_qre(C, c);
        rbc_79_mat_fq_set_random(&enc_seedexpander, Vp, N, D); // Generate a matrix Vp of size M times D

        for (int j = 0; j < ELL1; j++) {
        multiplication(Ap[j], Pa[j], Vp, M, D, N); // Compute Ap_j of size M times N of rank d, row support Vp
        }

        vectorization(A,Ap,ELL1,M,N); // Set A as the matrix whose ith line is Transpose(Ap_i) (written as a line)
        truncation(A1,A2,A,ELL1,MN-K,MN); // Write A = [A_1 | A_2 ] where A_1 is ELL2 times K and A_2 is ELL2 times MN-K

        multiplication(ct1, A1, Ht, ELL1, MN-K, K); // Compute A_1 * Transpose(H)

        addmat(U,ct1,A2,ELL1,K); // Compute ct1 = A_1 * Transpose(H)+ A_2 = A*Transpose(H)

        concatene2(pk,pkD1,ELL1,MN-K,K); // Compute pk = (pkD1, 0)

        multiplication(ct2,A,pk,ELL1,MN,ELL1); // Compute ct2 = A*Transpose(pk)
        addmat(V1,ct2,C,ELL1,ELL1);


        t3+=now_ms()-t0;
    

        t0=now_ms();



        concatene2(pk,pkD1,ELL1,MN-K,K); // Compute pk = (pkD, 0)

        truncation(Spoub,S,pk,ELL1,MN-K,MN);

        truncation(Spoub,E_k,E,ELL1,MN-K,MN);

        addmat(SE,E_k,S,ELL1,K);

        multiplication(W,U,SE,ELL1,K,ELL1);

        addmat(Wf,W,V1,ELL1,ELL1);


        rbc_79_mat_fq_to_qre(c1, Wf);
        rbc_79_gabidulin_decode(m1, code, c1);


        t4 += now_ms()-t0;
       // printf("\nDecryption done, decrypted\n");
       // rbc_79_vec_print(m1,3);
       // printf("\n Message encrypted was \n");
       // rbc_79_vec_print(m,3);
        
    }

    printf("Average time to generate public key H in KeyGen: %f ms\n", t_H/(double)(NUM_OF_TEST));

    printf("Average time to generate the other randomness in KeyGen: %f ms\n", t1/(double)(NUM_OF_TEST));

    printf("Average computation time in KeyGen: %f ms\n", t2/(double)(NUM_OF_TEST));

    printf("Average KeyGen time: %f ms\n", (t_H+t1+t2)/(double)(NUM_OF_TEST));

    printf("Average time to generate randomness in encryption: %f ms\n", t02/(double)(NUM_OF_TEST));

    printf("Average computation time in encryption: %f ms\n", t3/(double)(NUM_OF_TEST));

    printf("Average encryption time: %f ms\n", (t3+t02)/(double)(NUM_OF_TEST));

    printf("Average Decryption time: %f ms\n", (t4)/(double)(NUM_OF_TEST));

    printf("Average computation of transpose(H) (added to the KeyGen): %f ms\n", (t01)/(double)(NUM_OF_TEST));

    random_source_clear(&prng);

    rbc_79_mat_fq_clear(S);
    rbc_79_mat_fq_clear(E);
    rbc_79_mat_fq_clear(V);
    rbc_79_mat_fq_clear(E1);
    rbc_79_mat_fq_clear(E2);

    rbc_79_mat_fq_clear(A);
    rbc_79_mat_fq_clear(Vp);
    rbc_79_mat_fq_clear(H);
    rbc_79_mat_fq_clear(ct1);
    rbc_79_mat_fq_clear(ct2);
    rbc_79_mat_fq_clear(Ht);

    rbc_79_mat_fq_clear(pk);

    rbc_79_mat_fq_clear(SE);
    rbc_79_mat_fq_clear(pkD);
    rbc_79_mat_fq_clear(E_k);
    rbc_79_mat_fq_clear(A1);
    rbc_79_mat_fq_clear(A2);
    rbc_79_mat_fq_clear(W);
    rbc_79_mat_fq_clear(Spoub);
    rbc_79_mat_fq_clear(pkD1);
    rbc_79_mat_fq_clear(U);
    rbc_79_mat_fq_clear(Wf);
    rbc_79_mat_fq_clear(EFt);
   rbc_79_mat_fq_clear(C);
   rbc_79_mat_fq_clear(V1);
   rbc_79_vec_clear(m);
   rbc_79_vec_clear(m1);
   rbc_79_qre_clear(c);
   rbc_79_qre_clear(g);
   rbc_79_qre_clear(c1);
   rbc_79_qre_clear(g_place);
    for (int i = 0; i < ELL1; i++) {
        rbc_79_mat_fq_clear(P[i]); 
        rbc_79_mat_fq_clear(Pa[i]);
        rbc_79_mat_fq_clear(Ep[i]);
        rbc_79_mat_fq_clear(Ap[i]);
    }
    free(P);
    return 0;

}

