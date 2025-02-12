#pragma once

#include "circuit.h"
#include "matrix.h"
#include "sampling.h"
#include "cprf.h"
#include "common.h"
/*
Represents a pair of public / secret keys used in ABE
  - public key : A = [A0, A1, ..., Ak] where Ai in Zq^{n * l}
  - secret key : Tf in Z^{l * l} the trap used to compute A0
*/
typedef struct {
    matrix* A;
    signed_matrix Tf;
} ABE_keys;

typedef struct {
    bool * r;
    matrix k;
} ABE_skx;

typedef struct {
    matrix B;      // The trapdoor matrix
    matrix* A;     // Array of matrices A0,...,Ak
    matrix A_big;     // Array of matrices A0,...,Ak
    matrix v;      // Random vector v
} ABE_pp;

typedef struct {
  matrix B_tau;      // The trapdoor matrix Bτ−01
  int64_t sigma;     // σ
} ABE_msk;

typedef struct {
  ABE_pp pp;
  ABE_msk msk_out;
} ABE_setup_result;

typedef struct {
  int64_t sk_f_int;     // The main ciphertext component
  bool* sk_f_bool;     // The auxiliary components
  matrix u0;     // First auxiliary component
  matrix u1;     // Second auxiliary component
  int64_t u2;     // Third auxiliary component
} ABE_ct;

/*
Given a circuit f (and a sampler s) computes a new pair of ABE_keys
    - A = [A0, A1, ..., Ak] where Ai in Zq^{n * l}
    - Tf in Z^{l * l} the trap used to compute A0
*/
ABE_skx ABE_KeyGen(ClauseT* clauses, int num_clauses,  ABE_msk abe_msk, ABE_pp pp, int32_t x);
ABE_setup_result  ABE_Setup(uint64_t msk);

ABE_ct ABE_Enc(ClauseF* clauses, int num_clauses, ABE_msk msk , ABE_pp pp, bool flag_u);
bool ABE_dec(ABE_ct ct,  ABE_skx abe_sk, ABE_setup_result setup_res);
/*
Given a bit to encrypt u, a sampler s and
A = [A0, A1, ..., Ak] where Ai in Zq^{n * l}
returns CTf = [C0, C1-0, C1-1, ..., Ck-0, Ck-1]
*/
matrix* ABE_OfflineEnc(matrix* A, bool u);
