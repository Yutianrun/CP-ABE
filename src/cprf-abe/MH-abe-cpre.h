#pragma once

#include "circuit.h"
#include "matrix.h"
#include "sampling.h"
#include "cprf.h"
#include "common.h"
#include "abe.h"
/*
Represents a pair of public / secret keys used in MH_ABE
  - public key : A = [A0, A1, ..., Ak] where Ai in Zq^{n * l}
  - secret key : Tf in Z^{l * l} the trap used to compute A0
*/
typedef struct {
    matrix* A;
    signed_matrix Tf;
} MH_ABE_keys;

typedef struct {
  ABE_msk msk;
  int user_index;
} MH_ABE_sk;

typedef struct {    // The trapdoor matrix
    matrix* A;     // Array of matrices A0,...,Ak
    matrix A_big;     // Array of matrices A0,...,Ak
    // matrix v;      // Random vector v
} MH_ABE_pp;

typedef struct {
  matrix B_tau;      // The trapdoor matrix Bτ−01
  int64_t sigma;     // σ
} MH_ABE_msk;

typedef struct {
  MH_ABE_pp pp;
} MH_ABE_setup_result;

typedef struct {
  int64_t sk_f_int;     // The main ciphertext component
  bool* sk_f_bool;     // The auxiliary components
  matrix u0;     // First auxiliary component    // Second auxiliary component
  int64_t u1;
  matrix u2;     // Third auxiliary component
  circuit*** constrain;
} MH_ABE_ct;

typedef struct {
  matrix u0;     // First auxiliary component    // Second auxiliary component
  int64_t u1;     // Third auxiliary component
}MH_ABE_ct_one;

typedef struct {
  ABE_pp pp;
  matrix B;
  matrix v;
  int user_index;
}MH_ABE_pk;

typedef struct {
  MH_ABE_sk sk;
  MH_ABE_pk pk;
  int user_index;
} MH_ABE_key_pair;

typedef struct{
    bool * r;
    matrix Z;
    circuit ** eval;
}MH_ABE_ReEnc_key ;


/*
Given a circuit f (and a sampler s) computes a new pair of MH_ABE_keys
    - A = [A0, A1, ..., Ak] where Ai in Zq^{n * l}
    - Tf in Z^{l * l} the trap used to compute A0
*/
MH_ABE_key_pair MH_ABE_KeyGen(MH_ABE_pp pp, int32_t user_index);
MH_ABE_pp  MH_ABE_Setup(void);
MH_ABE_ct_one MH_ABE_Enc_step1(MH_ABE_pk pk, bool u);
MH_ABE_ct MH_ABE_Enc(MH_ABE_pk pk, bool u, ClauseF* clauses, int num_clauses);
bool MH_Decj_1(MH_ABE_sk sk, MH_ABE_ct_one ct, MH_ABE_pp pp, MH_ABE_pk pk);
bool MH_Decj_2(MH_ABE_sk sk, MH_ABE_ct ct, MH_ABE_pk pk) ;
MH_ABE_ReEnc_key MH_ABE_Re_KeyGen(ClauseT* clauses, int num_clauses, MH_ABE_sk sk, MH_ABE_pk pk, int64_t x);
MH_ABE_ct MH_ReEnc(MH_ABE_ReEnc_key rkx, MH_ABE_ct ct, MH_ABE_pp pp, MH_ABE_pk pk);