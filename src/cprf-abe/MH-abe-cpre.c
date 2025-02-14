#include "MH-abe-cpre.h"
#include "SH-abe-cpre.h"
#include "abe.h"
#include "cprf.h"
#include "sampling.h"
#include "common.h"


MH_ABE_pp  MH_ABE_Setup(void) {
    // Allocating new matrixes

    if(PARAMS.N == 0)
    {
        init_params_default();
    }

    matrix* A = new_matrixes(PARAMS.Att_num + 1, PARAMS.N, PARAMS.L);
    for (int i = 0; i < PARAMS.Att_num; i++) sample_Zq_uniform_matrix(A[i + 1]); 


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
    // Create MH_ABE_keys struct to return both pp and msk
    MH_ABE_pp pp = {
        .A = A,
        .A_big = BIG,
    };
    
    return pp;
}

MH_ABE_key_pair MH_ABE_KeyGen(MH_ABE_pp pp, int32_t user_index) {

    init_G();
    matrix B = new_matrix(PARAMS.N, PARAMS.MBAR + PARAMS.L);
    matrix R = new_matrix(PARAMS.MBAR, PARAMS.L);
    matrix v = new_matrix(PARAMS.N, 1);
    sample_Zq_uniform_matrix(v);

    single_TrapGen(B, R);
   ABE_setup_result abe_setup_res = ABE_Setup(user_index);
    MH_ABE_key_pair key_pair = {
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


MH_ABE_ct_one MH_ABE_Enc_step1(MH_ABE_pk pk, bool u) {    
    matrix s_vector_trans = new_matrix(1,pk.pp.B.rows);
    signed_matrix e_vector_signed_trans = new_signed_matrix(1,PARAMS.MBAR+PARAMS.M);

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
    MH_ABE_ct_one ct;
    ct.u0 = u0;
    ct.u1 = u1;
    return ct;
}



MH_ABE_ct MH_ABE_Enc(MH_ABE_pk pk, bool u, ClauseF* clauses, int num_clauses)
{
   //STEP 1: ENC 1
   matrix s_vector_trans = new_matrix(1,pk.pp.B.rows);
   signed_matrix e_vector_signed_trans = new_signed_matrix(1,PARAMS.MBAR+PARAMS.M);

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


    MH_ABE_ct ct = {
        .sk_f_bool = skf_bool,
        .sk_f_int = skf_int,
        .u0 = u0,
        .u1 = u1,  // Using first matrix from array
        .u2 = u2T,
        .constrain = sk_f
    };


    return ct;

}


static bool simple_function3(bool* input) {
    return input[0] ^ input[1];
}

static bool simple_function_clasuse4(bool* input) {
    return input[0] & input[1];
}

bool MH_Decj_1(MH_ABE_sk sk, MH_ABE_ct_one ct, MH_ABE_pp pp, MH_ABE_pk pk) {

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

bool MH_Decj_2(MH_ABE_sk sk, MH_ABE_ct ct, MH_ABE_pk pk) {

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




MH_ABE_ReEnc_key MH_ABE_Re_KeyGen(ClauseT* clauses, int num_clauses, MH_ABE_sk sk, MH_ABE_pk pk, int64_t x){

    real start, end;
    int prf_k = 8; // 比特宽度，可以根据需要调整

    // KEY-GEN STEP 1: 构造 clause 以及 PRF 电路
    printf(">KEY-GEN STEP 1: 构造 PRF 电路.\n");

    // 构建 PRF 电路
    circuit** x_cir = (circuit**)malloc(prf_k * sizeof(circuit*));
    circuit** msk_cir = (circuit**)malloc(prf_k * sizeof(circuit*));
    for (int i = 0; i < prf_k; i++) {
        x_cir[i] = gen_leaf(i + 1, true);
        msk_cir[i] = gen_leaf(i + 1 + prf_k, true);
    }
    circuit** prf_output = build_eval_circuit(prf_k, clauses, num_clauses, msk_cir, x_cir);

    matrix bigAx_concat = new_matrix(PARAMS.N, prf_k * PARAMS.M);
    matrix bigAH_concat = new_matrix(PARAMS.L * PARAMS.Att_num, prf_k * PARAMS.M);

    printf("initial BIG finished\n");
    printf("dimension of A: %d x %d\n", pk.pp.A[0].rows, pk.pp.A[0].columns);
    printf("dimension of BIG: %d x %d\n", pk.pp.A_big.rows, pk.pp.A_big.columns);
    printf("dimension of G: %d x %d\n", G.rows, G.columns);
    printf("dimension of bigAf_concat: %d x %d\n", bigAx_concat.rows, bigAx_concat.columns);
    printf("dimension of bigAH_concat: %d x %d\n", bigAH_concat.rows, bigAH_concat.columns);
    
    int col_offset = 0;
    for (int j = 0; j < prf_k; j++) {
        printf("Progress: computing Af for index %d out of %d...\n", j + 1, prf_k);
        matrix currAx = compute_Af(pk.pp.A, *prf_output[j]);
        printf("dimension of currAx: %d x %d \n", currAx.rows, currAx.columns);
        for (unsigned int r = 0; r < PARAMS.N; r++) {
            for (unsigned int c = 0; c < PARAMS.M; c++) {
                matrix_element(bigAx_concat, r, col_offset + c) = matrix_element(currAx, r, c);
            }
        }
        col_offset += PARAMS.M;
    }

    // KEY-GEN STEP 2: 计算 r 值.
    printf(">KEY-GEN STEP 2: 计算 r 值.\n");
    bool *r = (bool*)malloc(prf_k * sizeof(bool));
    uint32_t input = (sk.msk.sigma << prf_k) | x; // 组合前k位和 msk.sigma
    printf("Input (binary): ");
    for (int bit = sizeof(uint32_t)*8 - 1; bit >= 0; bit--) {
        printf("%d", (input >> bit) & 1);
        if (bit % 8 == 0 && bit != 0)
            printf(" ");
    }
    printf("\n");

    for (int j = 0; j < PRF_K; j++) {
        r[j] = compute_f(*prf_output[j], input);
    }
    printf("Reversed binary r : ");
    for (int i = PRF_K - 1; i >= 0; i--) {
        printf("%d ", r[i] ? 1 : 0);
    }
    printf("\n");

    
    // Build an equality circuit for key verification.
    circuit* bit_check[PRF_K];
    for (int i = 0; i < PRF_K; i++) {
        circuit* bit = gen_leaf(i + 1, true);
        if (!r[i]) {
            bit_check[i] = circuit_not(bit);
        } else {
            bit_check[i] = bit;
        }
    }

    circuit* equality_circuit = circuit_consecutive_and(bit_check, PRF_K);
    matrix Arf = compute_Af(pk.pp.A, *equality_circuit);
    printf("dimension of Arf: %d x %d\n", Arf.rows, Arf.columns);

    // 使用 r 值构造 key 仿真矩阵 BIG_eight.
    matrix BIG_eight = new_matrix(PARAMS.N, PARAMS.L * prf_k);
    for (int i = 0; i < prf_k ; i++) {
        matrix ti = copy_matrix(pk.pp.A[i]);
        for (int j = 0; j < PARAMS.N; j++)
            for (int k = 0; k < PARAMS.L; k++)
                matrix_element(BIG_eight, j, i * PARAMS.L + k) = matrix_element(ti, j, k);
        free_matrix(ti);
    }

    matrix Hr = compute_H_from_A_Af(&BIG_eight, &Arf);
    printf("dimension of Hr: %d x %d\n", Hr.rows, Hr.columns);

    matrix Axr = new_matrix(bigAx_concat.rows, Hr.columns);
    mul_matrix(bigAx_concat, Hr, Axr);
    printf("dimension of Axr: %d x %d\n", Axr.rows, Axr.columns);

    //sample left
    matrix result_d = new_matrix(Axr.columns+pk.pp.B.columns, 1);
    printf("dimension of pk.v: %d x %d\n", pk.v.rows, pk.v.columns);
    SampleLeft(pk.pp.B, sk.msk.B_tau, Axr, pk.v, PARAMS.SIGMA, result_d );

    // Choosing randomly R1 ∈ χ(m′+m)⌈log q⌉×n, R2 ∈ χ(m′+m)⌈log q⌉×m′ and r3 ∈ χ(m′+m)⌈log q⌉
    matrix R1 = new_matrix((PARAMS.MBAR + PARAMS.M+PARAMS.M)*PARAMS.K, PARAMS.N);
    matrix R2 = new_matrix((PARAMS.MBAR + PARAMS.M+PARAMS.M)*PARAMS.K, PARAMS.MBAR+PARAMS.M);
    matrix r3 = new_matrix((PARAMS.MBAR + PARAMS.M+PARAMS.M)*PARAMS.K, 1);

    sample_Zq_uniform_matrix(R1);
    sample_Zq_uniform_matrix(R2);
    sample_Zq_uniform_matrix(r3);

    // Construct matrix Z as:
    // Z = [ R1 * B_beta + R2,   R1 * v_beta + r3 - P2_d ]
    //     [      0_{1×m'}   ,              1         ]

    // Compute the left block: R1*B_beta + R2
    matrix temp1 = new_matrix(R1.rows, pk.B.columns);
    mul_matrix(R1, pk.B, temp1);
    matrix left_block = new_matrix(R1.rows, R2.columns);
    add_matrix(temp1, R2, left_block);
    free_matrix(temp1);

    // Compute the right block: R1*v_beta + r3 - P2_d
    matrix temp2 = new_matrix(R1.rows, pk.v.columns);
    mul_matrix(R1, pk.v, temp2);
    matrix temp3 = new_matrix(R1.rows, 1);
    add_matrix(temp2, r3, temp3);
    matrix right_block = new_matrix(R1.rows, 1);
    sub_matrix(temp3, P2(result_d), right_block);
    free_matrix(temp2);
    free_matrix(temp3);

    // Concatenate left_block and right_block horizontally to form the top part of Z
    int top_rows = left_block.rows;
    int top_cols = left_block.columns + right_block.columns;
    matrix Z_top = new_matrix(top_rows, top_cols);
    for (unsigned int i = 0; i < left_block.rows; i++) {
        for (unsigned int j = 0; j < left_block.columns; j++) {
            matrix_element(Z_top, i, j) = matrix_element(left_block, i, j);
        }
        for (unsigned int j = 0; j < right_block.columns; j++) {
            matrix_element(Z_top, i, left_block.columns + j) = matrix_element(right_block, i, j);
        }
    }
    free_matrix(left_block);
    free_matrix(right_block);

    // Create the bottom row: a zero row of length (top_cols - 1) with a trailing 1
    matrix Z_bottom = new_matrix(1, top_cols);
    for (unsigned int j = 0; j < top_cols - 1; j++) {
        matrix_element(Z_bottom, 0, j) = 0;
    }
    matrix_element(Z_bottom, 0, top_cols - 1) = 1;

    // Vertically concatenate Z_top and Z_bottom to build the final Z
    matrix Z = new_matrix(Z_top.rows + 1, top_cols);
    for (unsigned int i = 0; i < Z_top.rows; i++) {
        for (unsigned int j = 0; j < top_cols; j++) {
            matrix_element(Z, i, j) = matrix_element(Z_top, i, j);
        }
    }
    for (unsigned int j = 0; j < top_cols; j++) {
        matrix_element(Z, Z_top.rows, j) = matrix_element(Z_bottom, 0, j);
    }

//    int prf_k = PRF_K;
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
    matrix bigAf_concat = new_matrix(PARAMS.N, S_len * 2 * prf_k * PARAMS.M);
    // int col_offset = 0;
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
    mul_matrix(R1, bigAf_concat, u2T);
    matrix e1T = new_matrix(1, PARAMS.M * 2 * prf_k * S_len);

    signed_matrix e3_vector = new_signed_matrix(S_len * 2 * prf_k * PARAMS.M, 1);
    sample_Z_centered_matrix(e3_vector);
    add_matrix_error(u2, e3_vector, u2);


    // Z is now constructed as required.
    
    // Set the re-encryption ciphertext
    MH_ABE_ReEnc_key re_enc_key = {
        .r = r,
        .Z = Z,
        .eval = prf_output
    };

    return re_enc_key;
}


static bool simple_function(bool* input) {
    if (!input) {
        fprintf(stderr, "Error: simple_function received NULL pointer\n");
        exit(EXIT_FAILURE);
    }
    return input[0] ^ input[1];
}

static bool simple_function_clasuse2(bool* input) {
    if (!input) {
        fprintf(stderr, "Error: simple_function_clasuse2 received NULL pointer\n");
        exit(EXIT_FAILURE);
    }
    return input[0] & input[1];
}

MH_ABE_ct MH_ReEnc(MH_ABE_ReEnc_key rkx, MH_ABE_ct ct, MH_ABE_pp pp, MH_ABE_pk pk) {

   int prf_k = PRF_K; // 比特宽度，可以根据需要调整
    // 解析密钥 rkx = (r, k)
    bool* r = rkx.r;

    
    // Parse ciphertext components from ct
    matrix parsed_u0 = ct.u0;
    matrix parsed_u1 = ct.u2;
    int64_t parsed_u2 = ct.u1;
    int64_t parsed_sf = ct.sk_f_int;
    bool* parsed_sf_bool = ct.sk_f_bool;

    int count = 0;
    // for (int i = 2 * PRF_K * 3 - 1; i >= 0; i--) {
    //     printf("%d", ct.sk_f_bool[i] ? 1 : 0);
    //     count++;
    //     if (count % 8 == 0)
    //         printf(" ");
    // }
    printf("\n");

    printf("Decapsulation: Parsed ciphertext components:\n");
    printf("  - int sf: %lld\n", parsed_sf);
    printf("  - u0 dimensions: %d x %d\n", parsed_u0.rows, parsed_u0.columns);
    printf("  - u1 dimensions: %d x %d\n", parsed_u1.rows, parsed_u1.columns);
    printf("  - u2 value: %lld\n", (long long)parsed_u2);

    int num_clausesT = 2;
    int num_clausesF = 3;

    ClauseT* clauses = (ClauseT*)malloc(num_clausesT * sizeof(ClauseT));
    clauses[0].T = (int*)malloc(2 * sizeof(int));
    clauses[0].t_len = 2;
    clauses[0].T[0] = 0;
    clauses[0].T[1] = 2;

    clauses[1].T = (int*)malloc(2 * sizeof(int));
    clauses[1].t_len = 2;
    clauses[1].T[0] = 2;
    clauses[1].T[1] = 3;
    // 构建 PRF 电路

    ClauseF* clausesf = (ClauseF*)malloc(num_clausesT * sizeof(ClauseF));

    clausesf[0].f = simple_function;
    clausesf[1].f = simple_function_clasuse2;
    clausesf[0].t_len = clauses[0].t_len;
    clausesf[0].T = (int*)malloc(clauses[0].t_len * sizeof(int));
    for (int i = 0; i < clauses[0].t_len; i++) {
        clausesf[0].T[i] = clauses[0].T[i];
    }
    clausesf[1].t_len = clauses[1].t_len;
    clausesf[1].T = (int*)malloc(clauses[1].t_len * sizeof(int));
    for (int i = 0; i < clauses[1].t_len; i++) {
        clausesf[1].T[i] = clauses[1].T[i];
    }

    circuit** x = (circuit**)malloc((prf_k) * sizeof(circuit*));
    circuit** msk = (circuit**)malloc((prf_k) * sizeof(circuit*));
    circuit** x_cir = (circuit**)malloc((prf_k) * sizeof(circuit*));
    circuit** msk_cir = (circuit**)malloc((prf_k) * sizeof(circuit*));
    for(int i = 0;i< prf_k; i++){
        x[i] = gen_leaf(i+1, true);
        msk[i] = gen_leaf(i+1, true);
        x_cir[i] = gen_leaf(i+1, true);
        msk_cir[i] = gen_leaf(i+1+prf_k, true);
    }
    printf("initial circuit start\n");

    int S_len =100; // 示例长度
    Pair* S = (Pair*)malloc(S_len * sizeof(Pair));


    circuit** eval = build_eval_circuit(prf_k, clauses, num_clausesT, msk_cir, x_cir);
    // circuit** eval = build_eval_circuit_fixed_x(prf_k, clauses, num_clausesT, msk_cir, x_cir, x_bool);

    circuit*** sk_f = malloc(num_clausesF * sizeof(circuit**));
    for (int i = 0; i <  num_clausesF; i++) {
        sk_f[i] = malloc(2* prf_k * sizeof(circuit*));
        for(int j = 0; j < 2* prf_k; j++) {
            sk_f[i][j] = gen_leaf(2*prf_k*i + j +1, true);
        }
    }

    bool xvalue[] = {0, 0, 1, 1, 0, 0, 0, 0};
    bool* x_bool = xvalue;

    // circuit** constrain_eval  = build_constrain_eval_circuit(clauses, num_clausesF, num_clausesT, prf_k, sk_f, x);
    circuit** constrain_eval  = build_constrain_eval_circuit_fixed_x(clauses, num_clausesF, num_clausesT, prf_k, sk_f, x_bool);
    
    circuit*** constrain = build_constrain_circuit(clausesf, num_clausesT, S, &S_len, prf_k, msk);


    printf("initial circuit finished\n");
   bool r_prime[PRF_K];

    uint64_t input = parsed_sf ;
    printf("Decryption: Computing r_prime values:%lld \n", input);
    
    printf("binary r_prime: ");
    for(int i = 0; i < PRF_K; i++) {
        r_prime[i] = compute_f(*constrain_eval[0], input);
        printf("%d ", r_prime[i] ? 1 : 0);
    }
    printf("\n");

    printf("Decryption: abe sk r values:\n");
    for(int i =0;i<PRF_K;i++){
        printf("%d ", rkx.r[i] ? 1 : 0);
    }
    printf("\n");
    // // 如果 r 等于 r′ 则中止
    bool all_equal = true;
    for (int i = 0; i < PRF_K; i++) {
        if (rkx.r[i] != r_prime[i]) {
            all_equal = false;
            break;
        }
    }
    if (all_equal) {
        printf("Decryption aborted: r is completely equal to r'.\n");
        abort();
    }

    int r_prime_int = 0;
    int col_offset = 0;
    for (int i = 0; i < PRF_K; i++) {
        if (r_prime[i])
            r_prime_int |= (1 << i);
    }
    printf("Combined r_prime into int: %d\n", r_prime_int);


    printf("dimension of A: %d x %d\n", pk.pp.A[0].rows,pk.pp.A[0].columns);
    printf("dimension of BIG: %d x %d\n", pk.pp.A_big.rows,pk.pp.A_big.columns);
    printf("dimension of G: %d x %d\n", G.rows,G.columns);

    // print_matrix(pk.pp.A[1]);
    // 计算  Ax，依据  KeyGen 的实现
    matrix bigAf_concat = new_matrix(PARAMS.N, S_len * 2 * prf_k * PARAMS.M);
    col_offset = 0;
    for (int i = 0; i < S_len; i++) {
        for (int j = 0; j < 2 * prf_k; j++) {
            matrix currAf = compute_Af(pk.pp.A, *constrain[0][0]);
            for (unsigned int r = 0; r < PARAMS.N; r++) {
                for (unsigned int c = 0; c < PARAMS.M; c++) {
                    matrix_element(bigAf_concat, r, col_offset + c) = matrix_element(currAf, r, c);
                }
            }
            col_offset += PARAMS.M;
        }
    }

    printf("dimension of bigAf_concat: %d x %d\n", bigAf_concat.rows,bigAf_concat.columns);
    col_offset = 0;
    printf("S_len: %d\n", S_len);

    matrix bigAx_concat = new_matrix(PARAMS.N,  prf_k * PARAMS.M);
    printf("dimension of bigAx_concat: %d x %d\n", bigAx_concat.rows,bigAx_concat.columns);
    for (int j = 0; j <  prf_k; j++) {
        // printf("Progress: computing Ax for index %d out of %d...\n", j + 1, prf_k);
        matrix currAfx = compute_Af(pk.pp.A, *eval[0]);
        // printf("computed AF finished\n");
        // printf(" dimension of currAx: %d x %d \n", currAfx.rows,currAfx.columns);

        for (unsigned int r = 0; r < PARAMS.N; r++) {
            for (unsigned int c = 0; c < PARAMS.M; c++) {
                matrix_element(bigAx_concat, r, col_offset + c) = matrix_element(currAfx, r, c);
            }
        }
        col_offset += PARAMS.M;
    }

    printf("dimension of bigAx_concat: %d x %d\n", bigAx_concat.rows,bigAx_concat.columns);

     //Ir circuit
    circuit* bit_check[PRF_K];
    for (int i = 0; i < PRF_K; i++) {
        circuit* bit = gen_leaf(i + 1, true);
        if (!r[i]) {
            bit_check[i] = circuit_not(bit);
        } else {
            bit_check[i] = bit;
        }
    }

    circuit* equality_circuit = circuit_consecutive_and(bit_check, PRF_K);


    printf("start to compute Hr\n");
    matrix bigHconstrain_eval = new_matrix(PARAMS.L * PARAMS.Att_num,  prf_k  * PARAMS.M);

    printf("alloc matrix constrain_eval_Af\n");
    matrix* constrain_eval_Af = new_matrixes(PARAMS.Att_num+1, PARAMS.N, PARAMS.M);
    for (int i = 1; i <= PARAMS.Att_num; i++) {
        for (unsigned int r = 0; r < PARAMS.N; r++) {
            for (unsigned int c = 0; c < PARAMS.M; c++) {
                matrix_element(constrain_eval_Af[i], r, c) = matrix_element(bigAf_concat, r, (i - 1) * PARAMS.M + c);
            }
        }
    }

    printf("alloc matrix eval_Ax\n");
    matrix* eval_Ax = new_matrixes(prf_k+1, PARAMS.N, PARAMS.M);
    for(int i = 1; i<prf_k+1; i++){
        for (unsigned int r = 0; r < PARAMS.N; r++) {
            for (unsigned int c = 0; c < PARAMS.M; c++) {
                matrix_element(eval_Ax[i], r, c) = matrix_element(bigAx_concat, r, (i - 1) * PARAMS.M + c);
            }
        }
    }

    // for(int i =0; i<prf_k; i++){
    //     matrix Hconstrain_eval = compute_H(constrain_eval_Af, *constrain_eval[i], parsed_sf);
    //     printf("Progress: computing Hr for index %d out of %d...\n", i + 1, prf_k);
    //     printf("dimension of Hconstrain_eval: %d x %d\n", Hconstrain_eval.rows,Hconstrain_eval.columns);
    //     int col_offset = i * Hconstrain_eval.columns;
    //     for (unsigned int r = 0; r < Hconstrain_eval.rows; r++) {
    //         for (unsigned int c = 0; c < Hconstrain_eval.columns; c++) {
    //             matrix_element(bigHconstrain_eval, r, col_offset + c) = matrix_element(Hconstrain_eval, r, c);
    //         }
    //     }
    // } 
    omp_set_num_threads(8);
    #pragma omp parallel for
    for (int i = 0; i < prf_k; i++) {
        matrix Hconstrain_eval = compute_H(constrain_eval_Af, *constrain_eval[0], parsed_sf);
        // 使用一次临界区输出调试信息
        #pragma omp critical
        {
            printf("Progress: computing Hr for index %d out of %d...\n", i + 1, prf_k);
            printf("dimension of Hconstrain_eval: %d x %d\n", Hconstrain_eval.rows, Hconstrain_eval.columns);
        }
        int col_offset = i * Hconstrain_eval.columns;
        // 创建线程本地数组(假设bigHconstrain_eval的写入操作需要同步)
        // 这里仅作为示例展示如何减少多次临界区调用，实际需要根据matrix结构调整逻辑
        for (unsigned int r = 0; r < Hconstrain_eval.rows; r++) {
            for (unsigned int c = 0; c < Hconstrain_eval.columns; c++) {
                int value = matrix_element(Hconstrain_eval, r, c);
                // 合并写入：用一次临界区
                #pragma omp critical
                {
                    matrix_element(bigHconstrain_eval, r, col_offset + c) = value;
                }
            }
        }
    }
    printf("\n");

    printf("dimension of big Hconstrain_eval: %d x %d\n", bigHconstrain_eval.rows,bigHconstrain_eval.columns);

    // for(int i =0;i<prf_k;i++){
    matrix Hidentity= compute_H_prfk(eval_Ax, *equality_circuit,r_prime_int);
    printf("dimension of Hidentity: %d x %d\n", Hidentity.rows,Hidentity.columns);


    matrix HH = new_matrix(bigHconstrain_eval.rows, Hidentity.columns);
    
    mul_matrix(bigHconstrain_eval, Hidentity, HH);
    
    matrix uHH = new_matrix(parsed_u1.rows, HH.columns);

    mul_matrix(parsed_u1, HH, uHH);

}

