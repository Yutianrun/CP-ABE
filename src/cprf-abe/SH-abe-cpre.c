#include "SH-abe-cpre.h"
#include "abe.h"

SH_ABE_pp  SH_ABE_Setup(void) {
    // Allocating new matrixes

    init_params_default();
    
    matrix* A = new_matrixes(PARAMS.Att_num + 1, PARAMS.N, PARAMS.L);

    matrix BIG = new_matrix(PARAMS.N, PARAMS.L * PARAMS.Att_num);
    for (int i = 1; i < PARAMS.Att_num+1; i++) {
        matrix ti = copy_matrix(A[i]); 
        for (int j = 0; j < PARAMS.N; j++)
            for (int k = 0; k < PARAMS.L; k++)
                matrix_element(BIG, j, (i-1) * PARAMS.L + k) = matrix_element(ti, j, k);
        // free_matrix(ti);
    }

    // Generate A1, ..., Ak uniformely over Zq^{n * l}
    for (int i = 0; i < PARAMS.Att_num; i++) sample_Zq_uniform_matrix(A[i + 1]); 
    // Create SH_ABE_keys struct to return both pp and msk
    SH_ABE_pp pp = {
        .A = A,
        .A_big = BIG,
    };
    
    return pp;
}

SH_ABE_key_pair SH_ABE_KeyGen(SH_ABE_pp pp, int32_t user_index) {

    init_G();
    matrix B = new_matrix(PARAMS.N, PARAMS.MBAR + PARAMS.L);
    matrix R = new_matrix(PARAMS.MBAR, PARAMS.L);
    matrix v = new_matrix(PARAMS.N, 1);
    sample_Zq_uniform_matrix(v);

    single_TrapGen(B, R);
   ABE_setup_result abe_setup_res = ABE_Setup(user_index);
    SH_ABE_key_pair key_pair = {
        .pk = {
            .pp = abe_setup_res.pp,
            .B= B,
            .v = v,
            .user_index = user_index
        },
        .sk = {
            .msk = abe_setup_res.msk_out,
            .user_index = user_index
        },
        .user_index = user_index
    };
    return key_pair;

}


SH_ABE_ct_one SH_ABE_Enc_step1(SH_ABE_pk pk, bool u) {    
    matrix s_vector_trans = new_matrix(1,pk.pp.B.rows);
    signed_matrix e_vector_signed_trans = new_signed_matrix(1,PARAMS.MBAR);

    sample_Zq_uniform_matrix(s_vector_trans);
    sample_Z_centered_matrix(e_vector_signed_trans);
    int64_t e1 = uniform_mod_n_64(PARAMS.Q);

    // Transpose s_vector so that we can perform a row-vector multiplication.

    // Compute u0 = s^T * Bα + (e_vector_signed)^T
    matrix tmp_u0 = new_matrix(s_vector_trans.rows, pk.pp.B.columns);
    mul_matrix(s_vector_trans, pk.pp.B, tmp_u0);
    matrix u0 = new_matrix(tmp_u0.rows, tmp_u0.columns);
    add_matrix_error(tmp_u0, e_vector_signed_trans, u0);

    // Compute u1_temp = s^T * vα
    matrix tmp_u1 = new_matrix(1,1);
    mul_matrix(s_vector_trans, pk.pp.v, tmp_u1);

    // Compute u2 = (s^T * vα) + e1 + (u ? ⌊q/2⌉ : 0)
    int64_t u1 = matrix_element(tmp_u1, 0, 0) + e1;
    if (u)
        u1 += PARAMS.Q / 2;

    // Assemble ciphertext (additional components may be set as needed)
    SH_ABE_ct_one ct;
    ct.u0 = u0;
    ct.u1 = u1;
    return ct;
}



SH_ABE_ct SH_ABE_Enc(SH_ABE_pk pk, bool u, ClauseF* clauses, int num_clauses)
{
   //STEP 1: ENC 1
   matrix s_vector_trans = new_matrix(1,pk.pp.B.rows);
   signed_matrix e_vector_signed_trans = new_signed_matrix(1,PARAMS.MBAR);

   sample_Zq_uniform_matrix(s_vector_trans);
   sample_Z_centered_matrix(e_vector_signed_trans);
   int64_t e1 = uniform_mod_n_64(PARAMS.Q);

   // Transpose s_vector so that we can perform a row-vector multiplication.

   // Compute u0 = s^T * Bα + (e_vector_signed)^T
   matrix tmp_u0 = new_matrix(s_vector_trans.rows, pk.pp.B.columns);
   mul_matrix(s_vector_trans, pk.pp.B, tmp_u0);
   matrix u0 = new_matrix(tmp_u0.rows, tmp_u0.columns);
   add_matrix_error(tmp_u0, e_vector_signed_trans, u0);

   // Compute u1_temp = s^T * vα
   matrix tmp_u1 = new_matrix(1,1);
   mul_matrix(s_vector_trans, pk.pp.v, tmp_u1);

   // Compute u2 = (s^T * vα) + e1 + (u ? ⌊q/2⌉ : 0)
   int64_t u1 = matrix_element(tmp_u1, 0, 0) + e1;
   if (u)
       u1 += PARAMS.Q / 2;

   // Assemble ciphertext (additional components may be set as needed)


// Enc Step2

    int prf_k = PRF_K;
    circuit** msk_circuits = (circuit**)malloc(prf_k * sizeof(circuit*));
    for (int i = 0; i < prf_k; i++) {
        msk_circuits[i] = gen_leaf(i + 1, true);
    }

    int S_len;  
    Pair* S = (Pair*)malloc(S_len * sizeof(Pair));
    circuit*** sk_f = build_constrain_circuit(clauses, num_clauses, S, &S_len, prf_k, msk_circuits);

    printf("param S_len: %d\n", S_len);
    printf("param PARAMS.M: %d\n", PARAMS.M);
    bool* skf_bool = (bool*)malloc(2 * prf_k * S_len * sizeof(bool));
        printf(">KEY-ENC STEP 1: Computing sk_f values.\n");
    for (int i = 0; i < S_len; i++) {
        uint64_t input = rand();
        for (attribute j = 0; j < 2 * PRF_K; j++) {
            // skf_bool[i * 2 * PRF_K + j] = compute_f(*sk_f[i][j], input);
            skf_bool[i * 2 * PRF_K + j] =(bool)(rand() % 2);
        }
        printf("key-enc Completed sk_f computation for attribute index %d.\n", i);
    }

    int count = 0;
    for (int i = 2 * PRF_K * S_len - 1; i >= 0; i--) {
        printf("%d", skf_bool[i] ? 1 : 0);
        count++;
        if (count % 8 == 0)
            printf(" ");
    }
    printf("\n");
    int64_t skf_int = 0;
    for (int i = 0; i < 2 * prf_k * S_len; i++) {
        if (skf_bool[i]) {
            skf_int |= ((int64_t)1 << i);
        }
    }
    printf("Converted skf int64: %lld\n", (long long)skf_int);



    matrix* u1_matrixes = new_matrixes(S_len * prf_k, PARAMS.N, 1);
    matrix bigAf_concat = new_matrix(PARAMS.N, S_len * 2 * prf_k * PARAMS.M);
    int col_offset = 0;
    for (int i = 0; i < S_len; i++) {
        for (int j = 0; j < 2 * prf_k; j++) {
            matrix currAf = compute_Af(pk.pp.A, *sk_f[i][j]);
            if (skf_bool[i * 2 * PRF_K + j])
                sub_matrix(currAf, G, currAf);
            for (unsigned int r = 0; r < PARAMS.N; r++) {
                for (unsigned int c = 0; c < PARAMS.M; c++) {
                    matrix_element(bigAf_concat, r, col_offset + c) = matrix_element(currAf, r, c);
                }
            }
            col_offset += PARAMS.M;
        }
    }
    matrix u2 = new_matrix(PARAMS.M * 2 * prf_k * S_len, 1);
    matrix temp = new_matrix(PARAMS.N, PARAMS.L);
    printf("dimension of bigAf_concat: %d x %d\n", bigAf_concat.rows, bigAf_concat.columns);
    printf("dimension of G: %d x %d\n", G.rows, G.columns);
    matrix u2T = new_matrix(1, PARAMS.M * 2 * prf_k * S_len);
    mul_matrix(s_vector_trans, bigAf_concat, u2T);
    matrix e1T = new_matrix(1, PARAMS.M * 2 * prf_k * S_len);

    signed_matrix e2_vector = new_signed_matrix(S_len * 2 * prf_k * PARAMS.M, 1);
    sample_Z_centered_matrix(e2_vector);
    add_matrix_error(u2, e2_vector, u2);
    transpose_matrix(u2T, u2);


    SH_ABE_ct ct = {
        .sk_f_bool = skf_bool,
        .sk_f_int = skf_int,
        .u0 = u0,
        .u1 = u1,  // Using first matrix from array
        .u2 = u2T
    };


    return ct;

}


bool simple_function(bool* input) {
    return input[0] ^ input[1];
}

bool simple_function_clasuse2(bool* input) {
    return input[0] & input[1];
}


bool SH_Decj_1(SH_ABE_sk sk, SH_ABE_ct_one ct, SH_ABE_pp pp, SH_ABE_pk pk) {

    // Parse ciphertext components from ct
    matrix parsed_u0 = ct.u0;
    int64_t parsed_u1 = ct.u1;

    // Sample kα using preimage sampling: kα ← SamplePre(Bα, TBα, vα, δ)
    
    // Output 1 (true) if |u| > q/4, else 0 (false).
    int64_t u ;

    signed_matrix k = new_signed_matrix(pk.B.columns, 1);
    
    //samplePre
    SampleD(pk.pp.B, pk.v, sk.msk.B_tau, PARAMS.SIGMA, k);

    signed_matrix temp = new_signed_matrix(1, 1);
    mul_matrix_signed( parsed_u0, k, temp);
    u =parsed_u1- matrix_element(temp,0,0);
    return (llabs(u) > (PARAMS.Q / 4));
}

bool SH_Decj_2(SH_ABE_sk sk, SH_ABE_ct ct, SH_ABE_pk pk) {

    // Parse ciphertext components from ct
    matrix parsed_u0 = ct.u0;
    int64_t parsed_u1 = ct.u1;

    // Sample kα using preimage sampling: kα ← SamplePre(Bα, TBα, vα, δ)
    
    // Output 1 (true) if |u| > q/4, else 0 (false).
    int64_t u ;

    signed_matrix k = new_signed_matrix(pk.B.columns, 1);
    
    //samplePre
    SampleD(pk.pp.B, pk.v, sk.msk.B_tau, PARAMS.SIGMA, k);

    signed_matrix temp = new_signed_matrix(1, 1);
    mul_matrix_signed( parsed_u0, k, temp);
    u =parsed_u1- matrix_element(temp,0,0);
    return (llabs(u) > (PARAMS.Q / 4));
}