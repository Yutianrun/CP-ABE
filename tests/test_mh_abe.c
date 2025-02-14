#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "MH-abe-cpre.h"
// #include "SH-abe-cpre.h"
#include "cprf.h"

static bool simple_function(bool* input) {
    return input[0] ^ input[1];
}

static bool simple_function_clasuse2(bool* input) {
    return input[0] & input[1];
}


// 测试 MH_ABE_Setup
MH_ABE_pp test_MH_ABE_Setup(void) {
    printf("------Start MH_ABE_Setup process------\n");
    MH_ABE_pp pp = MH_ABE_Setup();
    // 简单检查部分矩阵是否分配成功
    assert(pp.A != NULL);
    printf("Dimension of A[0]: %d x %d\n", pp.A[0].rows, pp.A[0].columns);
    printf("Dimension of A_big: %d x %d\n", pp.A_big.rows, pp.A_big.columns);
    printf("------End MH_ABE_Setup process------\n");
    return pp;
}

// 测试 MH_ABE_KeyGen
MH_ABE_key_pair test_MH_ABE_KeyGen(MH_ABE_pp pp, int user_index) {
    printf("------Start MH_ABE_KeyGen process------\n");
    MH_ABE_key_pair key_pair = MH_ABE_KeyGen(pp, user_index);
    printf("Generated key pair for user_index: %d\n", key_pair.user_index);
    printf("------End MH_ABE_KeyGen process------\n");
    return key_pair;
}

MH_ABE_ct_one test_MH_ABE_Enc_step1(MH_ABE_pk pk) {
    printf("------Start MH_ABE_Enc_step1 process------\n");
    bool u = true;  // 使用true作为加密标志，可根据需要修改
    MH_ABE_ct_one ct_one = MH_ABE_Enc_step1(pk, u);
    // 根据需要打印ct_one内部重要的组件=
    printf("Step1 encryption: main integer component: %lld\n", (long long)ct_one.u1);
    printf("------End MH_ABE_Enc_step1 process------\n");
    return ct_one;
}

// 测试 MH_ABE_Enc
MH_ABE_ct test_MH_ABE_Enc(MH_ABE_pk pk) {
    printf("------Start MH_ABE_Enc process------\n");
    // 构造两个 dummy clause 用于加密
    int num_clauses = 2;
    ClauseF* clauses = (ClauseF*)malloc(num_clauses * sizeof(ClauseF));
    // Clause 0
    clauses[0].f = simple_function;
    clauses[0].t_len = 2;
    clauses[0].T = (int*)malloc(2 * sizeof(int));
    clauses[0].T[0] = 0;
    clauses[0].T[1] = 2;
    // Clause 1
    clauses[1].f = simple_function_clasuse2;
    clauses[1].t_len = 2;
    clauses[1].T = (int*)malloc(2 * sizeof(int));
    clauses[1].T[0] = 2;
    clauses[1].T[1] = 3;
    
    // 使用 true 作为加密标志
    MH_ABE_ct ct = MH_ABE_Enc(pk, true, clauses, num_clauses);
    
    printf("Ciphertext main integer component: %lld\n", (long long)ct.sk_f_int);
    printf("Dimension of u0: %d x %d\n", ct.u0.rows, ct.u0.columns);
    
    // 释放 clauses 内部动态分配的内存

    printf("------End MH_ABE_Enc process------\n");
    return ct;
}

// 测试 MH_ABE_Re_KeyGen：生成重加密密钥
MH_ABE_ReEnc_key test_MH_ABE_ReKeyGen(MH_ABE_sk sk, MH_ABE_pk pk) {
    printf("------Start MH_ABE_Re_KeyGen process------\n");
    int num_clauses = 2;
    ClauseT* clauses = (ClauseT*)malloc(num_clauses * sizeof(ClauseT));
    // Clause 0
    clauses[0].t_len = 2;
    clauses[0].T = (int*)malloc(2 * sizeof(int));
    clauses[0].T[0] = 0;
    clauses[0].T[1] = 2;
    // Clause 1
    clauses[1].t_len = 2;
    clauses[1].T = (int*)malloc(2 * sizeof(int));
    clauses[1].T[0] = 2;
    clauses[1].T[1] = 3;
    
    // 选择一个测试参数 x
    int64_t x = 15;
    MH_ABE_ReEnc_key rekey = MH_ABE_Re_KeyGen(clauses, num_clauses, sk, pk, x);

    // 释放分配的 clause 内存


    printf("Re-KeyGen: 重加密密钥生成成功.\n");
    printf("------End MH_ABE_Re_KeyGen process------\n");
    return rekey;
}
MH_ABE_ct test_enc2(MH_ABE_pk pk) {
    printf("------Start MH_ABE_Enc process (test_enc2)------\n");

    int num_clauses = 2;
    ClauseF* clauses = (ClauseF*)malloc(num_clauses * sizeof(ClauseF));

    // 初始化 Clause 0
    clauses[0].f = simple_function;
    clauses[0].t_len = 2;
    clauses[0].T = (int*)malloc(2 * sizeof(int));
    clauses[0].T[0] = 0;
    clauses[0].T[1] = 2;

    // 初始化 Clause 1
    clauses[1].f = simple_function_clasuse2;
    clauses[1].t_len = 2;
    clauses[1].T = (int*)malloc(2 * sizeof(int));
    clauses[1].T[0] = 2;
    clauses[1].T[1] = 3;
    
    MH_ABE_ct ct = MH_ABE_Enc(pk, true, clauses, num_clauses);
    printf("Ciphertext main integer component: %lld\n", (long long)ct.sk_f_int);
    printf("Dimension of u0: %d x %d\n", ct.u0.rows, ct.u0.columns);
    
    // 释放动态分配的 clauses 内存


    printf("------End MH_ABE_Enc process (test_enc2)------\n");
    return ct;
}


// Testing decryption (Dec)
bool test_dec(MH_ABE_key_pair key_pair, MH_ABE_ct_one ct_one, MH_ABE_pp pp) {
    printf("------Start MH_Decj_1 process (test_dec)------\n");
    bool dec_result = MH_Decj_1(key_pair.sk, ct_one, pp, key_pair.pk);
    printf("Decryption result: %s\n", dec_result ? "true" : "false");
    printf("------End MH_Decj_1 process (test_dec)------\n");
    return dec_result;
}
// 测试 MH_ReEnc：重加密操作
MH_ABE_ct test_MH_ReEnc(MH_ABE_ReEnc_key rekey, MH_ABE_ct ct, MH_ABE_pp pp, MH_ABE_pk pk) {
    printf("------Start MH_ReEnc process------\n");
    MH_ABE_ct new_ct = MH_ReEnc(rekey, ct, pp, pk);
    printf("ReEnc: 重加密后密文主组件: %lld\n", (long long)new_ct.sk_f_int);
    printf("------End MH_ReEnc process------\n");
    return new_ct;
}

int main(void) {
    // 测试 Setup
    // 测试 Setup
    // Define a list of parameter values to test.
    // 这里我们用参数数组来改变某些函数中的参数（例如user_index和x）
    // int param_values[] = {5, 10, 15, 20, 25};
    int param_values[] = {5};
    int num_params = sizeof(param_values) / sizeof(param_values[0]);

    for (int j = 0; j < num_params; j++) {
        int current_param = param_values[j];
        printf("\n===== Testing Parameter Value: %d =====\n", current_param);

        clock_t start, end;
        double setup_time, keygen_time, enc1_time, enc_time, dec_time, rekeygen_time, reenc_time;

        // 测试 Setup
        start = clock();
        MH_ABE_pp pp = test_MH_ABE_Setup();
        end = clock();
        setup_time = (double)(end - start) / CLOCKS_PER_SEC;

        // 使用current_param作为user_index
        start = clock();
        MH_ABE_key_pair key_pair = test_MH_ABE_KeyGen(pp, current_param);
        end = clock();
        keygen_time = (double)(end - start) / CLOCKS_PER_SEC;

        // 测试 Encryption Step1
        start = clock();
        MH_ABE_ct_one ct_one = test_MH_ABE_Enc_step1(key_pair.pk);
        end = clock();
        enc1_time = (double)(end - start) / CLOCKS_PER_SEC;

        // 测试 Encryption
        start = clock();
        MH_ABE_ct ct = test_MH_ABE_Enc(key_pair.pk);
        end = clock();
        enc_time = (double)(end - start) / CLOCKS_PER_SEC;

        // 测试 Decryption
        start = clock();
        bool dec_status = test_dec(key_pair, ct_one, pp);
        end = clock();
        dec_time = (double)(end - start) / CLOCKS_PER_SEC;

        // 测试 Re-KeyGen，利用current_param作为x参数
        start = clock();
        {
            int num_clauses = 2;
            ClauseT* clauses = (ClauseT*)malloc(num_clauses * sizeof(ClauseT));
            // Clause 0
            clauses[0].t_len = 2;
            clauses[0].T = (int*)malloc(2 * sizeof(int));
            clauses[0].T[0] = 0;
            clauses[0].T[1] = 2;
            // Clause 1
            clauses[1].t_len = 2;
            clauses[1].T = (int*)malloc(2 * sizeof(int));
            clauses[1].T[0] = 2;
            clauses[1].T[1] = 3;
            int64_t x = current_param; // 使用current_param作为重加密的参数x
            // 生成重加密密钥
            MH_ABE_ReEnc_key rekey = MH_ABE_Re_KeyGen(clauses, num_clauses, key_pair.sk, key_pair.pk, x);

            end = clock();
            rekeygen_time = (double)(end - start) / CLOCKS_PER_SEC;

            // 测试 ReEnc
            start = clock();
            MH_ABE_ct new_ct = test_MH_ReEnc(rekey, ct, pp, key_pair.pk);
            end = clock();
            reenc_time = (double)(end - start) / CLOCKS_PER_SEC;
        }

        // 打印每个参数对应的运行时间
        printf("\n------ Timing Results for Parameter Value: %d ------\n", current_param);
        printf("| %-15s | %-15s |\n", "Step", "Time (seconds)");
        printf("-----------------------------------------\n");
        printf("| %-15s | %-15f |\n", "Setup", setup_time);
        printf("| %-15s | %-15f |\n", "KeyGen", keygen_time);
        printf("| %-15s | %-15f |\n", "Enc_step1", enc1_time);
        printf("| %-15s | %-15f |\n", "Enc", enc_time);
        printf("| %-15s | %-15f |\n", "Dec", dec_time);
        printf("| %-15s | %-15f |\n", "Re-KeyGen", rekeygen_time);
        printf("| %-15s | %-15f |\n", "ReEnc", reenc_time);
        printf("------------------------------------------------\n");
    }
    return 0;
}