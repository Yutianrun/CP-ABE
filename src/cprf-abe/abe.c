#include "abe.h"

#include "common.h"
#include "cprf.h"
#define PRF_K 4


matrix* ABE_OfflineEnc(matrix* A, bool u, circuit* constrain, circuit* constrain_eval, ClauseF* f_clauses, int num_clauses, uint8_t msk_int) {
    int prf_k = PRF_K;
    matrix Af_constrain = compute_Af(A, *constrain);
    matrix* CTf = new_matrixes(2 * PARAMS.K + 1, PARAMS.M, PARAMS.L);
    circuit** msk = (circuit**)malloc(prf_k * sizeof(circuit*));
    int S_len = 0; // 初始化长度
    Pair* S = (Pair*)malloc(S_len * sizeof(Pair));
    circuit*** sk_f = build_constrain_circuit(f_clauses, num_clauses, S, &S_len, prf_k, msk);

    // 生成基于 msk_int 的叶子电路
    for (int i = 0; i < msk_int; i++) {
        msk[i] = gen_leaf(i, false);
    }

    uint32_t input = msk_int;

    // 分配 sk_f
    sk_f = (int **)malloc(num_clauses * sizeof(int *));
    for(int i = 0; i < num_clauses; i++) {
        sk_f[i] = (int *)malloc(PARAMS.K * sizeof(int));
        // 初始化 sk_f[i] 根据具体需求
    }

    // 生成 LWE 均匀秘密
    matrix S_matrix = new_matrix(PARAMS.M, PARAMS.N);
    sample_Zq_uniform_matrix(S_matrix);

    // 短高斯误差向量
    signed_matrix E = new_signed_matrix(PARAMS.M, PARAMS.L);
    sample_Z_centered_matrix(E);

    signed_matrix E_1 = new_signed_matrix(PARAMS.M, PARAMS.L);
    sample_Z_centered_matrix(E_1);

    signed_scalar e_2 = sample_Z_centered();

    matrix u_0 = new_matrix(PARAMS.M, PARAMS.L);
    mul_matrix(S_matrix, G, u_0);

    // 实例化 CTf 的第一个部分 (`u` in p.12 或 C0 in p.13)
    if (u)
        sample_Zq_uniform_matrix(CTf[0]);
    else
        mul_matrix_trap(S_matrix, Af_constrain, CTf[0]);

    // 计算 S * G 仅一次
    matrix SG = new_matrix(PARAMS.M, PARAMS.L);
    mul_matrix(S_matrix, G, SG);

    for (int i = 0; i < PARAMS.K; i++) {
        mul_matrix(SG, A[i + 1], CTf[1 + i]);
        add_matrix(CTf[1 + i], E, CTf[1 + i]);
    }

    // 处理 ABE_OfflineEnc 的剩余部分
    for (int i = 0; i < num_clauses; i++) {
        for (int b = 0; b < 2; b++) {
            sample_Z_centered_matrix(E);
            mul_matrix(S_matrix, A[1 + i], CTf[1 + 2 * i + b]);
            if (b) add_matrix(CTf[1 + 2 * i + b], SG, CTf[1 + 2 * i + b]);
            add_matrix_error(CTf[1 + 2 * i + b], E, CTf[1 + 2 * i + b]);
        }
    }

    free_matrix(S_matrix);
    free_signed_matrix(E);
    free_signed_matrix(E_1);
    free_matrix(SG);
    free_sk_tv(sk_f, num_clauses);
    free(S);
    free(msk);

    return CTf;
}


ABE_keys ABE_KeyGen(uint8_t x, uint8_t msk) {
    // Allocating new matrixes
    matrix* A = new_matrixes(PARAMS.K + 1, PARAMS.N, PARAMS.L);
    signed_matrix Tf = new_signed_matrix(PARAMS.L, PARAMS.L);

    // Generate A1, ..., Ak uniformely over Zq^{n * l}
    for (int i = 0; i < PARAMS.K; i++) sample_Zq_uniform_matrix(A[i + 1]);

    // Compute Af


    // Generate Tf
    sample_Z_centered_matrix(Tf);

    // Compute A0
    mul_matrix_trap(Af, Tf, A[0]);
    free_matrix(Af);

    ABE_keys key = {A, Tf};
    return key;
}

matrix* ABE_OfflineEnc(matrix* A, bool u, circuit* constrain, circuit* constrain_eval,  ClauseF* f_clasuses, int num_clauses, uint8_t msk_int) {

    int prf_k = PRF_K;
    matrix Af_constrain = compute_Af(A, *constrain);
    matrix* CTf = new_matrixes(2 * PARAMS.K + 1, PARAMS.M, PARAMS.L);
    circuit** msk = (circuit**)malloc((prf_k) * sizeof(circuit*));
        int S_len ; // 示例长度
    Pair* S = (Pair*)malloc(S_len * sizeof(Pair));
    circuit*** sk_f = build_constrain_circuit(f_clasuses, num_clauses, S, &S_len, prf_k, msk);

    // Generate leaf circuits based on msk_int
    for (int i = 0; i < msk_int; i++) {
        sk_f[i] = create_leaf_circuit(i+1, true);
    }

    uint32_t input = msk_int;

    int **sk_f;

    sk_f = (int **)malloc(num_clauses * sizeof(int *));

    for(int i =0;i<num_clauses;i++){
        sk_f[i] = (int *)malloc(prf_k * sizeof(int));
    }

    for(int i =0;i<num_clauses;i++){
        for(int j = 0;j<prf_k;j++){
            sk_f[i][j] = compute_f(*sk_f[i][j], input);
        }
    }

    // Generating LWE uniform secret
    matrix S_matrix = new_matrix(PARAMS.M, PARAMS.N);
    sample_Zq_uniform_matrix(S_matrix);

    // Short gaussian error vector
    signed_matrix E = new_signed_matrix(PARAMS.M, PARAMS.L);
    sample_Z_centered_matrix(E);

    signed_matrix E_1 = new_signed_matrix(PARAMS.M, PARAMS.L);
    sample_Z_centered_matrix(E_1);

    signed_scalar e_2 = sample_Z_centered();

    matrix u_0 =  new_matrix(PARAMS.M, PARAMS.L);

    mul_matrix(S_matrix, B, u_0);
    // Instantiating first term of CTf (`u` in p.12 or C0 in p.13)
    if (u)
        sample_Zq_uniform_matrix(CTf[0]);
    else {
        mul_matrix(S, A[0], CTf[0]);
        add_matrix_error(CTf[0], E, CTf[0]);
    }

    // Computing S * G only once
    matrix SG = new_matrix(PARAMS.M, PARAMS.L);
    mul_matrix(S, G, SG);

    for (int i = 0; i < PARAMS.K; i++) {
        for (int b = 0; b < 2; b++) {
            sample_Z_centered_matrix(E);
            mul_matrix(S, A[1 + i], CTf[1 + 2 * i + b]);
            if (b) add_matrix(CTf[1 + 2 * i + b], SG, CTf[1 + 2 * i + b]);
            add_matrix_error(CTf[1 + 2 * i + b], E, CTf[1 + 2 * i + b]);
        }
    }

    free_matrix(S);
    free_signed_matrix(E);
    free_matrix(SG);

    return CTf;
}



ABE_keys ABE_KeyGen(uint8_t x, uint8_t msk) {
    // Allocating new matrices
    matrix* A = new_matrixes(PARAMS.K + 1, PARAMS.N, PARAMS.L);
    signed_matrix Tf = new_signed_matrix(PARAMS.L, PARAMS.L);

    // Generate A1, ..., Ak uniformly over Zq^{n * l}
    for (int i = 0; i < PARAMS.K; i++) {
        sample_Zq_uniform_matrix(A[i + 1]);
    }

    // Compute Af
    signed_matrix Af = new_signed_matrix(PARAMS.N, PARAMS.L);
    // 假设有一个函数 compute_Af 用于计算 Af
    compute_Af(A, Tf, Af);

    // Generate Tf
    sample_Z_centered_matrix(Tf);

    // Compute A0 = Af * Tf
    mul_matrix_trap(Af, Tf, A[0]);
    free_matrix(Af);

    ABE_keys key = {A, Tf};
    return key;
}


ABE_keys ABE_KeyGen(uint8_t x, uint8_t msk) {
    // Allocating new matrices
    matrix* A = new_matrixes(PARAMS.K + 1, PARAMS.N, PARAMS.L);
    signed_matrix Tf = new_signed_matrix(PARAMS.L, PARAMS.L);

    // Generate A1, ..., Ak uniformly over Zq^{n * l}
    for (int i = 0; i < PARAMS.K; i++) {
        sample_Zq_uniform_matrix(A[i + 1]);
    }

    // Compute Af
    signed_matrix Af = new_signed_matrix(PARAMS.N, PARAMS.L);
    // 假设有一个函数 compute_Af 用于计算 Af
    compute_Af(A, Tf, Af);

    // Generate Tf
    sample_Z_centered_matrix(Tf);

    // Compute A0 = Af * Tf
    mul_matrix_trap(Af, Tf, A[0]);
    free_matrix(Af);

    ABE_keys key = {A, Tf};
    return key;
}



matrix* ABE_OfflineEnc(matrix* A, bool u, circuit* constrain, circuit* constrain_eval, ClauseF* f_clauses, int num_clauses, uint8_t msk_int) {
    int prf_k = PRF_K;
    matrix Af_constrain = compute_Af(A, *constrain);
    matrix* CTf = new_matrixes(2 * PARAMS.K + 1, PARAMS.M, PARAMS.L);
    circuit** msk = (circuit**)malloc(prf_k * sizeof(circuit*));
    int S_len = 0; // 初始化长度
    Pair* S = (Pair*)malloc(S_len * sizeof(Pair));
    circuit*** sk_f = build_constrain_circuit(f_clauses, num_clauses, S, &S_len, prf_k, msk);

    // 生成基于 msk_int 的叶子电路
    for (int i = 0; i < msk_int; i++) {
        msk[i] = gen_leaf(i, false);
    }

    uint32_t input = msk_int;

    // 分配 sk_f
    sk_f = (int **)malloc(num_clauses * sizeof(int *));
    for(int i = 0; i < num_clauses; i++) {
        sk_f[i] = (int *)malloc(PARAMS.K * sizeof(int));
        // 初始化 sk_f[i] 根据具体需求
    }

    // 生成 LWE 均匀秘密
    matrix S_matrix = new_matrix(PARAMS.M, PARAMS.N);
    sample_Zq_uniform_matrix(S_matrix);

    // 短高斯误差向量
    signed_matrix E = new_signed_matrix(PARAMS.M, PARAMS.L);
    sample_Z_centered_matrix(E);

    signed_matrix E_1 = new_signed_matrix(PARAMS.M, PARAMS.L);
    sample_Z_centered_matrix(E_1);

    signed_scalar e_2 = sample_Z_centered();

    matrix u_0 = new_matrix(PARAMS.M, PARAMS.L);
    mul_matrix(S_matrix, G, u_0);

    // 实例化 CTf 的第一个部分 (`u` in p.12 或 C0 in p.13)
    if (u)
        sample_Zq_uniform_matrix(CTf[0]);
    else
        mul_matrix_trap(S_matrix, Af_constrain, CTf[0]);

    // 计算 S * G 仅一次
    matrix SG = new_matrix(PARAMS.M, PARAMS.L);
    mul_matrix(S_matrix, G, SG);

    for (int i = 0; i < PARAMS.K; i++) {
        mul_matrix(SG, A[i + 1], CTf[1 + i]);
        add_matrix(CTf[1 + i], E, CTf[1 + i]);
    }

    // 处理 ABE_OfflineEnc 的剩余部分
    for (int i = 0; i < num_clauses; i++) {
        for (int b = 0; b < 2; b++) {
            sample_Z_centered_matrix(E);
            mul_matrix(S_matrix, A[1 + i], CTf[1 + 2 * i + b]);
            if (b) add_matrix(CTf[1 + 2 * i + b], SG, CTf[1 + 2 * i + b]);
            add_matrix_error(CTf[1 + 2 * i + b], E, CTf[1 + 2 * i + b]);
        }
    }

    free_matrix(S_matrix);
    free_signed_matrix(E);
    free_signed_matrix(E_1);
    free_matrix(SG);
    free_sk_tv(sk_f, num_clauses);
    free(S);
    free(msk);

    return CTf;
}


uint8_t ABE_Decrypt(ABE_params params, ABE_keys key, matrix* CTf, int user_id) {
    // 1. 提取密文组件
    matrix u0 = CTf[0];
    matrix* u = CTf + 1; // 假设密文从索引1开始

    // 2. 生成用户秘密密钥 sk_x
    matrix sk_x = new_matrix(params.m, 1);
    mul_matrix(params.MSK, key.A[user_id], sk_x); // 根据 user_id 选择对应的 A

    // 3. 计算 y = u1 - sk_x^T * u0 (模 q)
    matrix skx_T = transpose_matrix(sk_x);
    matrix temp = new_matrix(1, params.n);
    mul_matrix(skx_T, u0, temp);
    matrix y = new_matrix(1, 1);
    sub_matrix_mod(u[0], temp, y, params.q);

    // 4. 恢复消息 mu
    uint8_t mu;
    if (y->data[0][0] > params.q / 2) {
        mu = 1;
    } else {
        mu = 0;
    }

    // 5. 释放资源
    free_matrix(sk_x);
    free_matrix(skx_T);
    free_matrix(temp);
    free_matrix(y);

    return mu;
}
// ...existing code...