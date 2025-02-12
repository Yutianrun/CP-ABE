#pragma once

#include "circuit.h"
#include "matrix.h"
#include "sampling.h"
#include "cprf.h"
#include "common.h"
#include "abe.h"
/*
Represents a pair of public / secret keys used in SH_ABE
  - public key : A = [A0, A1, ..., Ak] where Ai in Zq^{n * l}
  - secret key : Tf in Z^{l * l} the trap used to compute A0
*/
typedef struct {
    matrix* A;
    signed_matrix Tf;
} SH_ABE_keys;

typedef struct {
  ABE_msk msk;
  int user_index;
} SH_ABE_sk;

typedef struct {
    matrix B;      // The trapdoor matrix
    matrix* A;     // Array of matrices A0,...,Ak
    matrix A_big;     // Array of matrices A0,...,Ak
    // matrix v;      // Random vector v
} SH_ABE_pp;

typedef struct {
  matrix B_tau;      // The trapdoor matrix Bτ−01
  int64_t sigma;     // σ
} SH_ABE_msk;

typedef struct {
  SH_ABE_pp pp;
} SH_ABE_setup_result;

typedef struct {
  int64_t sk_f_int;     // The main ciphertext component
  bool* sk_f_bool;     // The auxiliary components
  matrix u0;     // First auxiliary component    // Second auxiliary component
  int64_t u1;
  matrix u2;     // Third auxiliary component
} SH_ABE_ct;

typedef struct {
  matrix u0;     // First auxiliary component    // Second auxiliary component
  int64_t u1;     // Third auxiliary component
}SH_ABE_ct_one;

typedef struct {
  ABE_pp pp;
  matrix B;
  matrix v;
  int user_index;
}SH_ABE_pk;

typedef struct {
  SH_ABE_sk sk;
  SH_ABE_pk pk;
  int user_index;
} SH_ABE_key_pair;


/*
Given a circuit f (and a sampler s) computes a new pair of SH_ABE_keys
    - A = [A0, A1, ..., Ak] where Ai in Zq^{n * l}
    - Tf in Z^{l * l} the trap used to compute A0
*/
SH_ABE_key_pair SH_ABE_KeyGen(SH_ABE_pp pp, int32_t user_index);
SH_ABE_pp  SH_ABE_Setup(void);

SH_ABE_ct SH_ABE_Enc(SH_ABE_pk pk, bool u, ClauseF* clauses, int num_clauses);
bool SH_Decj_1(SH_ABE_sk sk, SH_ABE_ct_one ct, SH_ABE_pp pp, SH_ABE_pk pk);
bool SH_Decj_2(SH_ABE_sk sk, SH_ABE_ct ct, SH_ABE_pk pk) ;
