#include "abe.h"
#include "common.h"
#include "cprf.h"

#define PRF_K 8

ABE_setup_result  ABE_Setup(uint64_t msk) {
    // Allocating new matrixes

    // int32_t K = 63;
    // uint64_t Q = 9223372036854775807ULL;
    // int32_t N = 1;
    // real SIGMA = 7.00;
    // int32_t P = 30; // cp trap size M=32
    // init_params(N, Q, K, P, SIGMA);

    init_params_default();

    const unsigned int n = PARAMS.N;
    const unsigned int m_bar = PARAMS.MBAR;
    const unsigned int w = PARAMS.L;
    const unsigned int m = m_bar + w;
    
    printf("Dimensions: n=%u, m_bar=%u, w=%u, m=%u\n", n, m_bar, w, m);
    init_G();
    // 2. 分配矩阵内存
    printf("m,n: %d,%d\n", m, n);
    matrix B = new_matrix(n, m);
    matrix R = new_matrix(m_bar, w);

    single_TrapGen(B, R);
    // Generate vector v uniformly over Zq^n
    matrix v = new_matrix(PARAMS.N, 1);
    sample_Zq_uniform_matrix(v);

    matrix* A = new_matrixes(PARAMS.Att_num + 1, PARAMS.N, PARAMS.L);
    // Generate A1, ..., Ak uniformely over Zq^{n * l}
    for (int i = 0; i < PARAMS.Att_num; i++) sample_Zq_uniform_matrix(A[i + 1]); 
    // Create ABE_keys struct to return both pp and msk
    ABE_pp pp = {
        .B = B,
        .A = A,
        .v = v
    };
    
    ABE_msk msk_out = {
        .B_tau = R,
        .sigma = msk
    };

    // Assign to output
    ABE_setup_result result = {pp, msk_out};// 新建一个结构体
    return result;
}

ABE_ct ABE_KeyEnc(circuit*** f, ABE_msk msk, int S_len, ABE_pp pp, bool flag_u) {
    // printf("key-enc Step 1: Allocating new matrices A.\n");
    // matrix* A = new_matrixes(PARAMS.K + 1, PARAMS.N, PARAMS.L);
    
    printf("param S_len: %d\n", S_len);
    printf("param PARAMS.M: %d\n", PARAMS.M);
    int prf_k = PRF_K;
    bool sk_f[2 * PRF_K * S_len];
    printf("key-enc Step 2: Computing sk_f values.\n");
    for (int i = 0; i < S_len; i++) {
        uint64_t input = msk.sigma;
        for (attribute j = 0; j < 2 * PRF_K; j++) {
            sk_f[i * 2 * PRF_K + j] = compute_f(*f[i][j], input) ;
        }
        printf("key-enc   Completed sk_f computation for attribute index %d.\n", i);
    }

    printf("key-enc Step 3: Generating s_vector.\n");
    matrix s_vector = new_matrix(PARAMS.N, 1);
    sample_Zq_uniform_matrix(s_vector);

    printf("key-enc Step 4: Generating e_vector.\n");
    signed_matrix e_vector = new_signed_matrix(PARAMS.M+PARAMS.MBAR, 1);
    sample_Z_centered_matrix(e_vector);

    printf("key-enc Step 5: Allocating and sampling e1_vector.\n");
    signed_matrix e1_vector = new_signed_matrix(S_len * 2 * prf_k * PARAMS.M, 1);
    sample_Z_centered_matrix(e1_vector);
    
    printf("key-enc Step 6: Sampling e3.\n");
    signed e3 = sample_Z_centered();

    printf("key-enc Step 7: Computing u0 matrix.\n");
    

    matrix u0 = new_matrix(PARAMS.M+PARAMS.MBAR, 1);
    matrix u0T = new_matrix(1, PARAMS.M+PARAMS.MBAR);

    printf("dimension of pp.B: %d x %d\n", pp.B.rows,pp.B.columns);
    printf("dimension of s_vector: %d x %d\n", s_vector.rows,s_vector.columns);
    mul_transpose_matrix(s_vector, pp.B, u0T);
    transpose_matrix(u0T, u0);
    add_matrix_error(u0, e_vector, u0);

    free_matrix(u0T);
    printf("key-enc Step 8: Computing u1 matrices.\n");
    matrix* u1_matrixes = new_matrixes(S_len * prf_k, PARAMS.N, 1);
    

    matrix bigAf_concat = new_matrix(PARAMS.N, S_len * 2 * prf_k * PARAMS.M);

    printf("dimension of bigAf_concat: %d x %d\n", bigAf_concat.rows,bigAf_concat.columns);
    int col_offset = 0;
    for (int i = 0; i < S_len; i++) {
        for (int j = 0; j < 2 * prf_k; j++) {
            matrix currAf = compute_Af(pp.A, *f[i][j]);

            if (sk_f[i * 2 * PRF_K + j])
                sub_matrix(currAf, G, currAf);

            for (unsigned int r = 0; r < PARAMS.N; r++) {
                for (unsigned int c = 0; c < PARAMS.M; c++) {
                    matrix_element(bigAf_concat, r, col_offset + c) = matrix_element(currAf, r, c);
                }
            }
            col_offset += PARAMS.M;
            free_matrix(currAf);
        }
    }
    matrix u1 = new_matrix(PARAMS.M*2*prf_k*S_len, 1);
    matrix temp = new_matrix(PARAMS.N, PARAMS.L);
    printf("dimension of bigAf_concat: %d x %d\n", bigAf_concat.rows,bigAf_concat.columns);
    printf("dimension of G: %d x %d\n", G.rows, G.columns);

    matrix u1T = new_matrix(1, PARAMS.M*2*prf_k*S_len);
    mul_transpose_matrix(s_vector, bigAf_concat, u1T);
    transpose_matrix(u1T, u1);
    // free_matrix(u1T);
    add_matrix_error(u1, e1_vector, u1);

    printf("key-enc Step 9: Computing u2 value.\n");
    matrix temp_result_matrix = new_matrix(1, 1);
    mul_matrix_transpose(s_vector, pp.v, temp_result_matrix);
    int64_t temp_result = matrix_element(temp_result_matrix, 0, 0);
    free_matrix(temp_result_matrix);
    int64_t u2 = temp_result + e3;
    if (flag_u)
        u2 += (PARAMS.Q / 2);

    printf("key-enc Step 10: Finished computing all components. Returning ciphertext.\n");
    ABE_ct ct = {
        .sf = sk_f,
        .u0 = u0,
        .u1 = u1,  // Using first matrix from array
        .u2 = u2
    };
    return ct;
}



ABE_keys ABE_KeyGen(ABE_msk abe_msk, ABE_pp pp, int32_t x) {
    // Allocating new matrixes
    
    real start, end;
    int prf_k = 8; // 比特宽度，可以根据需要调整

    int num_clauses = 2;
    ClauseT* clauses = (ClauseT*)malloc(num_clauses * sizeof(ClauseT));
    clauses[0].T = (int*)malloc(2 * sizeof(int));
    clauses[0].t_len = 2;
    clauses[0].T[0] = 0;
    clauses[0].T[1] = 2;

    clauses[1].t_len = 2;
    clauses[1].T = (int*)malloc(2 * sizeof(int));
    clauses[1].T[0] = 2;
    clauses[1].T[1] = 3;
       // 构建 PRF 电路

    circuit** x_cir = (circuit**)malloc((prf_k) * sizeof(circuit*));
    circuit** msk_cir = (circuit**)malloc((prf_k) * sizeof(circuit*));
    for(int i = 0;i< prf_k; i++){
        x_cir[i] = gen_leaf(i+1, false);
        msk_cir[i] = gen_leaf(i+1+prf_k, false);
    }
    circuit** prf_output = build_eval_circuit(prf_k, clauses, num_clauses, msk_cir, x_cir);
    
    matrix bigAf_concat = new_matrix(PARAMS.N, prf_k * PARAMS.M);
    matrix bigAH_concat = new_matrix(PARAMS.L * PARAMS.Att_num,  prf_k * PARAMS.M);


    matrix BIG = new_matrix(PARAMS.N, PARAMS.L * PARAMS.Att_num);
    for (int i = 1; i < PARAMS.Att_num+1; i++) {
        matrix ti = copy_matrix(pp.A[i]); 
        for (int j = 0; j < PARAMS.N; j++)
            for (int k = 0; k < PARAMS.L; k++)
                matrix_element(BIG, j, (i-1) * PARAMS.L + k) = matrix_element(ti, j, k);
        // free_matrix(ti);
    }

    printf("initial BIG finished\n");

    printf("dimension of A: %d x %d\n", pp.A[0].rows,pp.A[0].columns);
    printf("dimension of BIG: %d x %d\n", BIG.rows,BIG.columns);
    printf("dimension of G: %d x %d\n", G.rows,G.columns);
    printf("dimension of bigAf_concat: %d x %d\n", bigAf_concat.rows,bigAf_concat.columns);
    printf("dimension of bigAH_concat: %d x %d\n", bigAH_concat.rows,bigAH_concat.columns);
    int col_offset = 0;

    for (int j = 0; j <  prf_k; j++) {
        printf("Progress: computing Af for index %d out of %d...\n", j + 1, prf_k);
        matrix currAf = compute_Af(pp.A, *prf_output[j]);
        // printf("computed AF finished\n");
        printf(" dimension of currAf: %d x %d ", currAf.rows,currAf.columns);
        matrix currAH = compute_H_from_A_Af(&BIG, &currAf);
        // matrix currAH = new_matrix(PARAMS.L * PARAMS.Att_num, PARAMS.M);
        // printf("computed Ah finished\n");
        printf(" dimension of currAH: %d x %d \n", currAH.rows,currAH.columns);

        for (unsigned int r = 0; r < PARAMS.N; r++) {
            for (unsigned int c = 0; c < PARAMS.M; c++) {
                matrix_element(bigAf_concat, r, col_offset + c) = matrix_element(currAf, r, c);
            }
        }
        for (unsigned int r = 0; r < bigAH_concat.rows; r++) {
            for (unsigned int c = 0; c < PARAMS.M; c++) {
                matrix_element(bigAH_concat, r, col_offset + c) = matrix_element(currAH, r, c);
            }
        }
        col_offset += PARAMS.M;
        // free_matrix(currAf);
        // free_matrix(currAH);
    }

    // /* Compute the product of bigAf_concat and bigAH_concat */
    matrix product = new_matrix(BIG.rows, bigAH_concat.columns);
    mul_matrix(BIG, bigAH_concat, product);



    /* Compute the transpose of BIG */
    matrix BIG_trans = new_matrix(BIG.columns, BIG.rows);
    transpose_matrix(BIG, BIG_trans);

    /* Verify element–wise equality between product and BIG_trans */
    bool equal = true;
    for (unsigned int i = 0; i < product.rows; i++) {
        for (unsigned int j = 0; j < product.columns; j++) {
            if (matrix_element(product, i, j) != matrix_element(bigAf_concat, i, j)) {
                equal = false;
                break;
            }
        }
        if (!equal)
            break;
    }

    if (equal)
        printf("Verification passed: big * bigAh equals big At.\n");
    else
    {
        printf("Verification failed: big * bigAh does not equal big At.\n");

    }
        
    free_matrix(product);

    // bool r[PRF_K];
    // printf("key-gen Step 2: Computing r values.\n");
    // uint32_t input = (x << prf_k) | abe_msk.sigma; // 组合前k位和后k位
    // for (int j = 0; j <  PRF_K; j++) {
    //     r[j] = compute_f(*prf_output[j], input) ;
    // }

    // Build an equality circuit that checks if the input’s first 8 bits equal r[0..7].
    // For each bit the circuit output is defined as follows:
    //    if r[i] == 1 then output the i-th input bit,
    //    if r[i] == 0 then output the negation of the i-th input bit.
    // Finally, take the “AND” (implemented via NAND gates) of all these 8 bits so that 
    // the overall circuit outputs 1 exactly when all bits match.
    // circuit* bit_check[PRF_K];
    // for (int i = 0; i < PRF_K; i++) {
    //     // Obtain the i-th input bit from the public attribute: a leaf circuit reading bit (i+1)
    //     circuit* bit = gen_leaf(i + 1, true);
    //     if (!r[i]) {
    //         // Negate the bit:  NOT(x) = NAND(x, x)
    //         bit_check[i] = circuit_not(bit);
    //     } else {
    //         bit_check[i] = bit;
    //     }
    // }

    // circuit* equality_circuit = circuit_consecutive_and(bit_check, PRF_K);

    // matrix Arf = compute_Af(pp.A, *equality_circuit);
    // printf("dimension of Arf: %d x %d\n", Arf.rows,Arf.columns);

    // matrix BIG_eight = new_matrix(PARAMS.N, PARAMS.L * prf_k);
    // for (int i = 1; i < prf_k + 1; i++) {
    //     matrix ti = copy_matrix(pp.A[i]); 
    //     if (r[i-1]) add_matrix(ti, G, ti);
    //     for (int j = 0; j < PARAMS.N; j++)
    //         for (int k = 0; k < PARAMS.L; k++)
    //             matrix_element(BIG_eight, j, (i-1) * PARAMS.L + k) = matrix_element(ti, j, k);
    //     free_matrix(ti);
    // }

    // matrix Hr = compute_H_from_A_Af(&BIG_eight, &Arf);
    // printf("dimension of Hr: %d x %d\n", Hr.rows,Hr.columns);

    // matrix Axr = new_matrix(bigAf_concat.rows, Hr.columns);
    // mul_matrix(bigAf_concat, Hr, Axr);

    // printf("dimension of Axr: %d x %d\n", Axr.rows,Axr.columns);

    // printf("dimension of pp.B %d x %d\n", pp.B.rows,pp.B.columns);
    // signed_matrix k1 = new_signed_matrix(PARAMS.M+PARAMS.MBAR, 1);
    // matrix target_u = new_matrix(PARAMS.N, 1);
    // matrix k2 = new_matrix(Axr.columns, 1);
    // sample_Zq_uniform_matrix(k2);
    
    // matrix temp = new_matrix(1, 1); 

    // printf("dimension of k2 %d x %d\n", k2.rows,k2.columns);
    // mul_matrix(Axr, k2, temp);

    // sub_matrix(pp.v, temp, target_u);
    // SampleD(pp.B, target_u, abe_msk.B_tau, PARAMS.SIGMA, k1);

    // matrix computed_u = new_matrix(pp.B.rows, 1);
    // mul_matrix_trap(pp.B, k1, computed_u);

    // bool equal_bk1 = true;
    // for (unsigned int i = 0; i < computed_u.rows; i++) {
    //     if (matrix_element(computed_u, i, 0) != matrix_element(target_u, i, 0)) {
    //         equal_bk1 = false;
    //         break;
    //     }
    // }
    // if (equal_bk1)
    //     printf("Verification passed: B*k1 equals target_u.\n");
    // else
    //     printf("Verification failed: B*k1 does not equal target_u.\n");

    // free_matrix(computed_u);

    // matrix k_vector = new_matrix(k1.rows + k2.rows, 1);
    // for (unsigned int i = 0; i < k1.rows; i++) {
    //     matrix_element(k_vector, i, 0) = matrix_element(k1, i, 0);
    // }
    // for (unsigned int i = 0; i < k2.rows; i++) {
    //     matrix_element(k_vector, k1.rows + i, 0) = matrix_element(k2, i, 0);
    // }

    // printf("dimension of k_vector %d x %d\n", k_vector.rows,k_vector.columns);


    // // Construct the concatenated matrix [B | Axr]
    // matrix concatenated = new_matrix(pp.B.rows, pp.B.columns + Axr.columns);
    // for (unsigned int i = 0; i < pp.B.rows; i++) {
    //     for (unsigned int j = 0; j < pp.B.columns; j++)
    //         matrix_element(concatenated, i, j) = matrix_element(pp.B, i, j);
    //     for (unsigned int j = 0; j < Axr.columns; j++)
    //         matrix_element(concatenated, i, pp.B.columns + j) = matrix_element(Axr, i, j);
    // }

    // // Multiply concatenated matrix with k_vector
    // matrix product = new_matrix(concatenated.rows, 1);
    // mul_matrix(concatenated, k_vector, product);

    // // Check if product equals v (pp.v)
    // bool equal = true;
    // for (unsigned int i = 0; i < product.rows; i++) {
    //     if (matrix_element(product, i, 0) != matrix_element(pp.v, i, 0)) {
    //         equal = false;
    //         break;
    //     }
    // }
    // if (equal)
    //     printf("Verification passed: concatenated matrix * k equals v.\n");
    // else
    //     printf("Verification failed: concatenated matrix * k does not equal v.\n");

    // free_matrix(concatenated);
    // free_matrix(product);
    // // Compute Af
    // matrix Af = compute_Af(A, f);

    // // Generate Tf
    // sample_Z_centered_matrix(Tf);

    // // Compute A0
    // mul_matrix_trap(Af, Tf, A[0]);
    // free_matrix(Af);

    // ABE_keys key = {A, Tf};
    // return key;
}
