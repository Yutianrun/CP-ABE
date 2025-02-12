// Language: c
#include <assert.h>
#include <assert.h>
#include <stdio.h>
#include "abe.h"
#include "cprf.h"

// Test ABE_Setup
ABE_setup_result  test_ABE_Setup() {
    uint64_t test_msk = 103;
    real start, end;
    // Separator line indicating the start of the setup process
    printf("------Start ABE_Setup process------\n");
    ABE_setup_result res = ABE_Setup(test_msk);
    CHRONO("Generated A in %fs\n", {
        for (int i = 0; i < PARAMS.K + 1; i++) sample_Zq_uniform_matrix(res.pp.A[i]);
    });
        
    // Check public parameters
    assert(res.pp.B.data != NULL);
    assert(res.pp.A[0].data != NULL);
    assert(res.pp.v.data != NULL);

    // Check master secret key
    assert(res.msk_out.B_tau.data != NULL);
    assert(res.msk_out.sigma == test_msk);
    printf("ABE_Setup: sigma = %llu\n", (unsigned long long)res.msk_out.sigma);
    printf("binary of sigma = ");
    {
        int bit_count = sizeof(res.msk_out.sigma) * 8;
        int printed = 0;
        for (int i = bit_count - 1; i >= 0; i--) {
            printf("%d", (int)((res.msk_out.sigma >> i) & 1));
            printed++;
            if (printed % 8 == 0 && i > 0)
                printf(" ");
        }
        printf("\n");
    }
    printf("dimension of B: %d x %d\n", res.pp.B.rows, res.pp.B.columns);
    printf("dimension of B_tau: %d x %d\n", res.msk_out.B_tau.rows, res.msk_out.B_tau.columns);
    printf("num of A: %d\n", PARAMS.Att_num);
    printf("dimension of A[0]: %d x %d\n", res.pp.A[0].rows, res.pp.A[0].columns);
    printf("dimension of A_big: %d x %d\n", res.pp.A_big.rows, res.pp.A_big.columns);
    printf("dimension of v: %d x %d\n", res.pp.v.rows, res.pp.v.columns);
    printf("ABE_Setup test passed\n");


    // 构造[R|I]矩阵
    matrix RI = new_matrix(PARAMS.MBAR+PARAMS.M, PARAMS.M);
    // 复制R部分
    for(unsigned int i = 0; i < PARAMS.MBAR; i++) {
        for(unsigned int j = 0; j < PARAMS.M; j++) {
            matrix_element(RI, i, j) = matrix_element(res.msk_out.B_tau, i, j);
        }
    }
    // 添加单位矩阵I部分
    for(unsigned int i = 0; i < PARAMS.M; i++) {
        for(unsigned int j = 0; j < PARAMS.M; j++) {
            matrix_element(RI, i + PARAMS.MBAR, j) = (i == j) ? 1 : 0;
        }
    }
    printf("dimension of RI: %d %d\n", RI.rows, RI.columns);
    // 计算G_p = B * [R|I]
    // printf("Computing G_p = B * [R|I]...\n");
    matrix G_p = new_matrix(PARAMS.N, PARAMS.M);
    mul_matrix(res.pp.B, RI, G_p);

    // 验证G_p == G
    // printf("Verifying G_p == G...\n");
    bool valid = true;
    for(unsigned int i = 0; i < PARAMS.N && valid; i++) {
        for(unsigned int j = 0; j < PARAMS.M && valid; j++) {
            if(matrix_element(G_p, i, j) != matrix_element(G, i, j)) {
                valid = false;
                printf("Mismatch at (%u,%u): G_p=%lu, G=%lu\n", 
                    i, j,
                    matrix_element(G_p, i, j),
                    matrix_element(G, i, j));
            }
        }
    }

    if(valid) {
        printf("✓ Trapgen test passed!\n");
    } else {
        printf("✗ Trapgen test does not pass!\n");
    }
    
    // Separator line indicating the end of the setup process
    printf("------End ABE_Setup process------\n");

    // free_matrix(RI);
    // free_matrix(G_p);
    return res;
}

bool simple_function(bool* input) {
    return input[0] ^ input[1];
}

bool simple_function_clasuse2(bool* input) {
    return input[0] & input[1];
}

// 重写后的测试 ABE_KeyEnc 函数
ABE_ct test_ABE_Enc(ABE_setup_result setup_res) {
    real start, end;
    int prf_k = PRF_K; // PRF 门电路比特宽度
    int num_clauses = 2;
    // 分配并初始化 dummy clauses
    ClauseF* clauses = (ClauseF*)malloc(num_clauses * sizeof(ClauseF));
    // Clause 0
    clauses[0].f = simple_function;
    clauses[0].clauseT.t_len = 2;
    clauses[0].clauseT.T = (int*)malloc(2 * sizeof(int));
    clauses[0].clauseT.T[0] = 0;
    clauses[0].clauseT.T[1] = 2;
    // Clause 1
    clauses[1].f = simple_function_clasuse2;
    clauses[1].clauseT.t_len = 2;
    clauses[1].clauseT.T = (int*)malloc(2 * sizeof(int));
    clauses[1].clauseT.T[0] = 2;
    clauses[1].clauseT.T[1] = 3;

    // Separator line indicating the start of the encryption process
    printf("------start ABE_Enc process------\n");

    // 测试 flag u = true
    ABE_ct ct;
    CHRONO("Finished ABE_Enc computation in %fs\n", {
        ct = ABE_Enc(clauses, num_clauses, setup_res.msk_out, setup_res.pp, true);
    });

    // 打印密文信息
    printf("Ciphertext field: sf\n");
    printf("sf = %lld\n", (long long)ct.sk_f_int);

    // Separator line indicating success of encryption test
    printf("------End ABE_Enc process------\n");

    // 释放 clauses 内部分配的数组
    for (int i = 0; i < num_clauses; i++) {
        free(clauses[i].clauseT.T);
    }
    free(clauses);


    int count = 0;
    for (int i = 2 * PRF_K * 3 - 1; i >= 0; i--) {
        printf("%d", ct.sk_f_bool[i] ? 1 : 0);
        count++;
        if (count % 8 == 0)
            printf(" ");
    }
    printf("\n");

    return ct;
}

ABE_skx test_ABE_KeyGen(ABE_setup_result setup_res) {
    real start, end;

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
    printf("-----start ABE_KeyGen process-----\n");

    ABE_skx keys;
    CHRONO("Finished ABE_KeyGen computation in %fs\n", {
         keys = ABE_KeyGen(clauses, num_clauses, setup_res.msk_out, setup_res.pp, 12);
    });
    // Validate that the produced key components are allocated.
    assert(keys.r != NULL);
    assert(keys.k.data != NULL);

    printf("-------End ABE_KeyGen process------\n");
    return keys;
}

void test_ABE_Dec(ABE_ct ct, ABE_skx skx, ABE_setup_result setup_res) {
    real start, end;
    printf("------start ABE_Dec process------\n");
    // 解密
    bool result;
    CHRONO("Finished ABE_Dec computation in %fs\n", {
         result = ABE_dec(ct, skx, setup_res);
    });
    printf("Decryption result: %s\n", result ? "true" : "false");
    printf("------End ABE_Dec process------\n");
}

int main(void) {

    ABE_setup_result setup_res = test_ABE_Setup();
    ABE_ct ct = test_ABE_Enc(setup_res);

    int count = 0;
    for (int i = 2 * PRF_K * 3 - 1; i >= 0; i--) {
        printf("%d", ct.sk_f_bool[i] ? 1 : 0);
        count++;
        if (count % 8 == 0)
            printf(" ");
    }
    printf("\n");
    ABE_skx skx = test_ABE_KeyGen(setup_res);

    printf("dimension of skx %d x %d \n", skx.k.rows, skx.k.columns);
    test_ABE_Dec(ct, skx, setup_res);

    return 0;
}
 