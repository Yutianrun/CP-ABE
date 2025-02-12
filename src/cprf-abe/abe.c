#include "abe.h"


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
    
    init_G();
    matrix B = new_matrix(n, m);
    matrix R = new_matrix(m_bar, w);

    single_TrapGen(B, R);
    // Generate vector v uniformly over Zq^n
    matrix v = new_matrix(PARAMS.N, 1);
    sample_Zq_uniform_matrix(v);

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
    // Create ABE_keys struct to return both pp and msk
    ABE_pp pp = {
        .B = B,
        .A = A,
        .A_big = BIG,
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

ABE_ct ABE_Enc(ClauseF* clauses, int num_clauses, ABE_msk msk , ABE_pp pp, bool flag_u) {    
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
        uint64_t input = msk.sigma;
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
    printf("skf (binary): ");
    for (int bit = sizeof(int64_t)*8 - 1; bit >= 0; bit--) {
        printf("%d", (int)((skf_int >> bit) & 1));
        if (bit % 8 == 0 && bit != 0) {
            printf(" ");
        }
    }

    printf("\n");
    printf(">KEY-ENC STEP 2: Generating s_vector.\n");
    matrix s_vector = new_matrix(PARAMS.N, 1);
    sample_Zq_uniform_matrix(s_vector);

    printf(">KEY-ENC STEP 3: Generating e_vector.\n");
    signed_matrix e_vector = new_signed_matrix(PARAMS.M+PARAMS.MBAR, 1);
    sample_Z_centered_matrix(e_vector);

    printf(">KEY-ENC STEP 4: Allocating and sampling e1_vector.\n");
    signed_matrix e1_vector = new_signed_matrix(S_len * 2 * prf_k * PARAMS.M, 1);
    sample_Z_centered_matrix(e1_vector);
    
    printf(">KEY-ENC STEP 5: Sampling e3.\n");
    signed e3 = sample_Z_centered();

    printf(">KEY-ENC STEP 6: Computing u0 matrix.\n");
    matrix u0 = new_matrix(PARAMS.M+PARAMS.MBAR, 1);
    matrix u0T = new_matrix(1, PARAMS.M+PARAMS.MBAR);
    printf("dimension of pp.B: %d x %d\n", pp.B.rows, pp.B.columns);
    printf("dimension of s_vector: %d x %d\n", s_vector.rows, s_vector.columns);
    mul_transpose_matrix(s_vector, pp.B, u0T);
    transpose_matrix(u0T, u0);
    add_matrix_error(u0, e_vector, u0);
    transpose_matrix(u0, u0T);
    // free_matrix(u0T);

    printf(">KEY-ENC STEP 7: Computing u1 matrices.\n");
    matrix* u1_matrixes = new_matrixes(S_len * prf_k, PARAMS.N, 1);
    matrix bigAf_concat = new_matrix(PARAMS.N, S_len * 2 * prf_k * PARAMS.M);
    int col_offset = 0;
    for (int i = 0; i < S_len; i++) {
        for (int j = 0; j < 2 * prf_k; j++) {
            matrix currAf = compute_Af(pp.A, *sk_f[i][j]);
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
    matrix u1 = new_matrix(PARAMS.M * 2 * prf_k * S_len, 1);
    matrix temp = new_matrix(PARAMS.N, PARAMS.L);
    printf("dimension of bigAf_concat: %d x %d\n", bigAf_concat.rows, bigAf_concat.columns);
    printf("dimension of G: %d x %d\n", G.rows, G.columns);
    matrix u1T = new_matrix(1, PARAMS.M * 2 * prf_k * S_len);
    mul_transpose_matrix(s_vector, bigAf_concat, u1T);
    matrix e1T = new_matrix(1, PARAMS.M * 2 * prf_k * S_len);
    add_matrix_error(u1, e1_vector, u1);
    transpose_matrix(u1T, u1);

    printf(">KEY-ENC STEP 8: Computing u2 value.\n");
    matrix temp_result_matrix = new_matrix(1, 1);
    mul_matrix_transpose(s_vector, pp.v, temp_result_matrix);
    int64_t temp_result = matrix_element(temp_result_matrix, 0, 0);
    free_matrix(temp_result_matrix);
    int64_t u2 = temp_result + e3;
    if (flag_u)
        u2 += (PARAMS.Q / 2);

    printf(">KEY-ENC STEP 9: Finished computing all components. Returning ciphertext.\n");


    ABE_ct ct = {
        .sk_f_bool = skf_bool,
        .sk_f_int = skf_int,
        .u0 = u0T,
        .u1 = u1T,  // Using first matrix from array
        .u2 = u2
    };


    return ct;
}



ABE_skx ABE_KeyGen(ClauseT* clauses, int num_clauses,  ABE_msk abe_msk, ABE_pp pp, int32_t x) {
    // Allocating new matrixes
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

    matrix bigAf_concat = new_matrix(PARAMS.N, prf_k * PARAMS.M);
    matrix bigAH_concat = new_matrix(PARAMS.L * PARAMS.Att_num, prf_k * PARAMS.M);

    printf("initial BIG finished\n");
    printf("dimension of A: %d x %d\n", pp.A[0].rows, pp.A[0].columns);
    printf("dimension of BIG: %d x %d\n", pp.A_big.rows, pp.A_big.columns);
    printf("dimension of G: %d x %d\n", G.rows, G.columns);
    printf("dimension of bigAf_concat: %d x %d\n", bigAf_concat.rows, bigAf_concat.columns);
    printf("dimension of bigAH_concat: %d x %d\n", bigAH_concat.rows, bigAH_concat.columns);
    
    // int col_offset = 0;
    // for (int j = 0; j < prf_k; j++) {
    //     printf("Progress: computing Af for index %d out of %d...\n", j + 1, prf_k);
    //     matrix currAf = compute_Af(pp.A, *prf_output[j]);
    //     printf("dimension of currAf: %d x %d \n", currAf.rows, currAf.columns);
    //     matrix currAH = compute_H_from_A_Af(&pp.A_big, &currAf);
    //     printf(" dimension of currAH: %d x %d \n", currAH.rows, currAH.columns);
    //     for (unsigned int r = 0; r < PARAMS.N; r++) {
    //         for (unsigned int c = 0; c < PARAMS.M; c++) {
    //             matrix_element(bigAf_concat, r, col_offset + c) = matrix_element(currAf, r, c);
    //         }
    //     }
    //     for (unsigned int r = 0; r < bigAH_concat.rows; r++) {
    //         for (unsigned int c = 0; c < PARAMS.M; c++) {
    //             matrix_element(bigAH_concat, r, col_offset + c) = matrix_element(currAH, r, c);
    //         }
    //     }
    //     col_offset += PARAMS.M;
    // }

    // 验证 bigAf_concat 与 bigAH_concat 的运算结果
    matrix product = new_matrix(pp.A_big.rows, bigAH_concat.columns);
    mul_matrix(pp.A_big, bigAH_concat, product);
    matrix BIG_trans = new_matrix(pp.A_big.columns, pp.A_big.rows);
    transpose_matrix(pp.A_big, BIG_trans);
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
        printf("Verification failed: big * bigAh does not equal big At.\n");

    // KEY-GEN STEP 2: 计算 r 值.
    printf(">KEY-GEN STEP 2: 计算 r 值.\n");
    bool *r = (bool*)malloc(prf_k * sizeof(bool));
    uint32_t input = (abe_msk.sigma << prf_k) | x; // 组合前k位和 msk.sigma
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
    matrix Arf = compute_Af(pp.A, *equality_circuit);
    printf("dimension of Arf: %d x %d\n", Arf.rows, Arf.columns);

    // 使用 r 值构造 key 仿真矩阵 BIG_eight.
    matrix BIG_eight = new_matrix(PARAMS.N, PARAMS.L * prf_k);
    for (int i = 0; i < prf_k ; i++) {
        matrix ti = copy_matrix(pp.A[i]);
        for (int j = 0; j < PARAMS.N; j++)
            for (int k = 0; k < PARAMS.L; k++)
                matrix_element(BIG_eight, j, i * PARAMS.L + k) = matrix_element(ti, j, k);
        free_matrix(ti);
    }

    matrix Hr = compute_H_from_A_Af(&BIG_eight, &Arf);
    printf("dimension of Hr: %d x %d\n", Hr.rows, Hr.columns);

    matrix Axr = new_matrix(bigAf_concat.rows, Hr.columns);
    mul_matrix(bigAf_concat, Hr, Axr);
    printf("dimension of Axr: %d x %d\n", Axr.rows, Axr.columns);

    // KEY-GEN STEP 3: 生成密钥向量并进行验证.
    printf(">KEY-GEN STEP 3: 生成密钥向量并进行验证.\n");
    printf("dimension of pp.B %d x %d\n", pp.B.rows, pp.B.columns);
    signed_matrix k1 = new_signed_matrix(PARAMS.M + PARAMS.MBAR, 1);
    matrix target_u = new_matrix(PARAMS.N, 1);
    matrix k2 = new_matrix(Axr.columns, 1);
    sample_Zq_uniform_matrix(k2);

    matrix temp = new_matrix(1, 1);
    printf("dimension of k2 %d x %d\n", k2.rows, k2.columns);
    mul_matrix(Axr, k2, temp);
    sub_matrix(pp.v, temp, target_u);
    SampleD(pp.B, target_u, abe_msk.B_tau, PARAMS.SIGMA, k1);

    matrix computed_u = new_matrix(pp.B.rows, 1);
    mul_matrix_trap(pp.B, k1, computed_u);
    bool equal_bk1 = true;
    for (unsigned int i = 0; i < computed_u.rows; i++) {
        if (matrix_element(computed_u, i, 0) != matrix_element(target_u, i, 0)) {
            equal_bk1 = false;
            break;
        }
    }
    if (equal_bk1)
        printf("Verification passed: B*k1 equals target_u.\n");
    else
        printf("Verification failed: B*k1 does not equal target_u.\n");
    free_matrix(computed_u);

    matrix k_vector = new_matrix(k1.rows + k2.rows, 1);
    for (unsigned int i = 0; i < k1.rows; i++) {
        matrix_element(k_vector, i, 0) = matrix_element(k1, i, 0);
    }
    for (unsigned int i = 0; i < k2.rows; i++) {
        matrix_element(k_vector, k1.rows + i, 0) = matrix_element(k2, i, 0);
    }
    printf("dimension of k_vector %d x %d\n", k_vector.rows, k_vector.columns);


    // Concatenate pp.B and Axr horizontally into a new matrix called concatenated
    matrix concatenated = new_matrix(pp.B.rows, pp.B.columns + Axr.columns);
    for (unsigned int i = 0; i < pp.B.rows; i++) {
        for (unsigned int j = 0; j < pp.B.columns; j++) {
            matrix_element(concatenated, i, j) = matrix_element(pp.B, i, j);
        }
        for (unsigned int j = 0; j < Axr.columns; j++) {
            matrix_element(concatenated, i, pp.B.columns + j) = matrix_element(Axr, i, j);
        }
    }
    printf("Concatenated matrix dimensions: %d x %d\n", concatenated.rows, concatenated.columns);
    matrix prod = new_matrix(concatenated.rows, 1);
    mul_matrix(concatenated, k_vector, prod);
    bool match = true;
    for (unsigned int i = 0; i < prod.rows; i++) {
        if (matrix_element(prod, i, 0) != matrix_element(pp.v, i, 0)) {
            match = false;
            break;
        }
    }
    if (match)
        printf("Verification passed: B||Axr * k equals v.\n");
    else
        printf("Verification failed: B||Axr * k does not equal v.\n");
    free_matrix(prod);
    // 最后构造私钥结构体 ABE_skx 并返回.
    ABE_skx key = {r, k_vector };
    return key;
}

bool simple_function(bool* input) {
    return input[0] ^ input[1];
}

bool simple_function_clasuse2(bool* input) {
    return input[0] & input[1];
}

// abe_dec 函数实现
bool ABE_dec(ABE_ct ct,  ABE_skx abe_sk, ABE_setup_result setup_res) {

    int prf_k = PRF_K; // 比特宽度，可以根据需要调整
    // 解析密钥 abe_sk = (r, k)
    bool* r = abe_sk.r;
    matrix k = abe_sk.k;
    
    // Parse ciphertext components from ct
    matrix parsed_u0 = ct.u0;
    matrix parsed_u1 = ct.u1;
    int64_t parsed_u2 = ct.u2;
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
    printf("  - k dimension %d x %d\n", k.rows, k.columns);

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
    clausesf[0].clauseT.t_len = clauses[0].t_len;
    clausesf[0].clauseT.T = (int*)malloc(clauses[0].t_len * sizeof(int));
    for (int i = 0; i < clauses[0].t_len; i++) {
        clausesf[0].clauseT.T[i] = clauses[0].T[i];
    }
    clausesf[1].clauseT.t_len = clauses[1].t_len;
    clausesf[1].clauseT.T = (int*)malloc(clauses[1].t_len * sizeof(int));
    for (int i = 0; i < clauses[1].t_len; i++) {
        clausesf[1].clauseT.T[i] = clauses[1].T[i];
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

    uint64_t input = (parsed_sf << prf_k) | 12;
    printf("Decryption: Computing r_prime values:%lld \n", input);
    
    printf("binary r_prime: ");
    for(int i = 0; i < PRF_K; i++) {
        r_prime[i] = compute_f(*constrain_eval[i], input);
        printf("%d ", r_prime[i] ? 1 : 0);
    }
    printf("\n");

    // input =  12;    
    // printf("Decryption: Computing constrain values:%lld \n", input);
    // for(int i=0;i<S_len;i++)
    // {
    //     printf("Slen: %d ", i);
    //     for(int j=0 ;j< 2*prf_k;j++){
    //         uint64_t input =  12;
    //         r_prime[i] = compute_f(*constrain[i][j], input);
    //         printf("%d", r_prime[i] ? 1 : 0);
    //         if((j+1)%8 == 0 )
    //             printf(" ");
    //     }
    //     printf("\n");
    // }


    printf("Decryption: abe sk r values:\n");
    for(int i =0;i<PRF_K;i++){
        printf("%d ", abe_sk.r[i] ? 1 : 0);
    }
    printf("\n");
    // // 如果 r 等于 r′ 则中止
    bool all_equal = true;
    for (int i = 0; i < PRF_K; i++) {
        if (abe_sk.r[i] != r_prime[i]) {
            all_equal = false;
            break;
        }
    }
    if (all_equal) {
        printf("Decryption aborted: r is completely equal to r'.\n");
        return false;
    }

    int r_prime_int = 0;
    int col_offset = 0;
    for (int i = 0; i < PRF_K; i++) {
        if (r_prime[i])
            r_prime_int |= (1 << i);
    }
    printf("Combined r_prime into int: %d\n", r_prime_int);
    matrix BIG = new_matrix(PARAMS.N, PARAMS.L * PARAMS.Att_num);
    for (int i = 1; i < PARAMS.Att_num+1; i++) {
        matrix ti = copy_matrix(setup_res.pp.A[i]); 
        for (int j = 0; j < PARAMS.N; j++)
            for (int k = 0; k < PARAMS.L; k++)
                matrix_element(BIG, j, (i-1) * PARAMS.L + k) = matrix_element(ti, j, k);
    }
    printf("initial BIG finished\n");



    // 计算  Ax，依据  KeyGen 的实现
    // matrix bigAf_concat = new_matrix(PARAMS.N, prf_k * PARAMS.M);

    printf("dimension of A: %d x %d\n", setup_res.pp.A[0].rows,setup_res.pp.A[0].columns);
    printf("dimension of BIG: %d x %d\n", BIG.rows,BIG.columns);
    printf("dimension of G: %d x %d\n", G.rows,G.columns);
    // printf("dimension of bigAf_concat: %d x %d\n", bigAf_concat.rows,bigAf_concat.columns);
    col_offset = 0;
    printf("S_len: %d\n", S_len);

    matrix bigAx_concat = new_matrix(PARAMS.N, S_len * 2 * prf_k * PARAMS.M);
    printf("dimension of bigAx_concat: %d x %d\n", bigAx_concat.rows,bigAx_concat.columns);
    for (int j = 0; j <  prf_k; j++) {
        printf("Progress: computing Ax for index %d out of %d...\n", j + 1, prf_k);
        matrix currAfx = compute_Af(setup_res.pp.A, *eval[j]);
        // printf("computed AF finished\n");
        printf(" dimension of currAx: %d x %d \n", currAfx.rows,currAfx.columns);

        for (unsigned int r = 0; r < PARAMS.N; r++) {
            for (unsigned int c = 0; c < PARAMS.M; c++) {
                matrix_element(bigAx_concat, r, col_offset + c) = matrix_element(currAfx, r, c);
            }
        }
        col_offset += PARAMS.M;
    }

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

    for(int i =0; i<prf_k; i++){
        r_prime[i] = compute_f(*constrain_eval[i], input);

        printf("constrain eval: r prime value:%d ", r_prime[i] ? 1 : 0);
        matrix Hconstrain_eval = compute_H(setup_res.pp.A, *constrain_eval[i], parsed_sf);
        printf("Progress: computing Hr for index %d out of %d...\n", i + 1, prf_k);
        printf("dimension of Hconstrain_eval: %d x %d\n", Hconstrain_eval.rows,Hconstrain_eval.columns);
        int col_offset = i * Hconstrain_eval.columns;
        for (unsigned int r = 0; r < Hconstrain_eval.rows; r++) {
            for (unsigned int c = 0; c < Hconstrain_eval.columns; c++) {
                matrix_element(bigHconstrain_eval, r, col_offset + c) = matrix_element(Hconstrain_eval, r, c);
            }
        }
    } 
    printf("\n");

    printf("dimension of big Hconstrain_eval: %d x %d\n", bigHconstrain_eval.rows,bigHconstrain_eval.columns);

    // for(int i =0;i<prf_k;i++){
    matrix Hidentity= compute_H_prfk(setup_res.pp.A, *equality_circuit,r_prime_int);
    printf("dimension of Hidentity: %d x %d\n", Hidentity.rows,Hidentity.columns);


    matrix HH = new_matrix(bigHconstrain_eval.rows, Hidentity.columns);
    
    mul_matrix(bigHconstrain_eval, Hidentity, HH);
    
    matrix uHH = new_matrix(parsed_u1.rows, HH.columns);

    mul_matrix(parsed_u1, HH, uHH);
    unsigned int cols_u0 = parsed_u0.columns;
    unsigned int cols_uhh = uHH.columns;
    unsigned int total_cols = cols_u0 + cols_uhh;
    unsigned int rows = parsed_u0.rows; // assume parsed_u0 and uHH have the same number of rows
    matrix concat = new_matrix(rows, total_cols);
    
    for (unsigned int i = 0; i < rows; i++) {
        for (unsigned int j = 0; j < cols_u0; j++) {
            matrix_element(concat, i, j) = matrix_element(parsed_u0, i, j);
        }
        for (unsigned int j = 0; j < cols_uhh; j++) {
            matrix_element(concat, i, cols_u0 + j) = matrix_element(uHH, i, j);
        }
    }
    
    printf("dimension of concat: %d x %d\n", concat.rows, concat.columns);
    // Compute dot product by treating concat as a flattened vector in row-major order.
    // It is assumed that the secret key vector *k is organized in the corresponding order.


    matrix dot_product_matrix = new_matrix(1,1);

    mul_matrix(concat, k, dot_product_matrix);

    int64_t dot_product = matrix_element(dot_product_matrix, 0, 0);
    
    int64_t u = parsed_u2 - dot_product;
    free_matrix(concat);
    
    if (llabs(u) >= (PARAMS.Q / 4)) {
        printf("Decryption output: 1\n");
        return true;
    } else {
        printf("Decryption output: 0\n");
        return false;
    }
    
}



