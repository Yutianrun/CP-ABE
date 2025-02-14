#define _POSIX_C_SOURCE 199309L
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "SH-abe-cpre.h"
#include "cprf.h"
#include <stdbool.h>



// 测试 SH_ABE_Setup
SH_ABE_pp test_SH_ABE_Setup(void) {
    printf("------Start SH_ABE_Setup process------\n");
    SH_ABE_pp pp = SH_ABE_Setup();
    // 简单检查部分矩阵是否分配成功
    assert(pp.A != NULL);
    printf("Dimension of A[0]: %d x %d\n", pp.A[0].rows, pp.A[0].columns);
    printf("Dimension of A_big: %d x %d\n", pp.A_big.rows, pp.A_big.columns);
    printf("------End SH_ABE_Setup process------\n");
    return pp;
}

// 测试 SH_ABE_KeyGen
SH_ABE_key_pair test_SH_ABE_KeyGen(SH_ABE_pp pp, int user_index) {
    printf("------Start SH_ABE_KeyGen process------\n");
    SH_ABE_key_pair key_pair = SH_ABE_KeyGen(pp, user_index);
    printf("Generated key pair for user_index: %d\n", key_pair.user_index);
    printf("------End SH_ABE_KeyGen process------\n");
    return key_pair;
}

SH_ABE_ct_one test_SH_ABE_Enc_step1(SH_ABE_pk pk) {
    printf("------Start SH_ABE_Enc_step1 process------\n");
    bool u = true;  // 使用true作为加密标志，可根据需要修改
    SH_ABE_ct_one ct_one = SH_ABE_Enc_step1(pk, u);
    // 根据需要打印ct_one内部重要的组件=
    printf("Step1 encryption: main integer component: %lld\n", (long long)ct_one.u1);
    printf("------End SH_ABE_Enc_step1 process------\n");
    return ct_one;
}

// 测试 SH_ABE_Enc
SH_ABE_ct test_SH_ABE_Enc(SH_ABE_pk pk) {
    printf("------Start SH_ABE_Enc process------\n");
    // 构造两个 dummy clause 用于加密

    // 使用 true 作为加密标志
    SH_ABE_ct ct = SH_ABE_Enc(pk, true);
    
    printf("Ciphertext main integer component: %lld\n", (long long)ct.sk_f_int);
    printf("Dimension of u0: %d x %d\n", ct.u0.rows, ct.u0.columns);
    

    printf("------End SH_ABE_Enc process------\n");
    return ct;
}

// 测试 SH_ABE_Re_KeyGen：生成重加密密钥
SH_ABE_ReEnc_key test_SH_ABE_ReKeyGen(SH_ABE_sk sk, SH_ABE_pk pk) {
    printf("------Start SH_ABE_Re_KeyGen process------\n");
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
    SH_ABE_ReEnc_key rekey = SH_ABE_Re_KeyGen(clauses, num_clauses, sk, pk, x);

    // 释放分配的 clause 内存


    printf("Re-KeyGen: 重加密密钥生成成功.\n");
    printf("------End SH_ABE_Re_KeyGen process------\n");
    return rekey;
}


SH_ABE_ct test_enc2(SH_ABE_pk pk) {
    printf("------Start SH_ABE_Enc process (test_enc2)------\n");


    
    SH_ABE_ct ct = SH_ABE_Enc(pk, true);
    printf("Ciphertext main integer component: %lld\n", (long long)ct.sk_f_int);
    printf("Dimension of u0: %d x %d\n", ct.u0.rows, ct.u0.columns);
    


    printf("------End SH_ABE_Enc process (test_enc2)------\n");
    return ct;
}


// Testing decryption (Dec)
bool test_dec(SH_ABE_key_pair key_pair, SH_ABE_ct_one ct_one, SH_ABE_pp pp) {
    printf("------Start SH_Decj_1 process (test_dec)------\n");
    bool dec_result = SH_Decj_1(key_pair.sk, ct_one, pp, key_pair.pk);
    printf("Decryption result: %s\n", dec_result ? "true" : "false");
    printf("------End SH_Decj_1 process (test_dec)------\n");
    return dec_result;
}
// 测试 SH_ReEnc：重加密操作
SH_ABE_ct test_SH_ReEnc(SH_ABE_ReEnc_key rekey, SH_ABE_ct ct, SH_ABE_pp pp, SH_ABE_pk pk) {
    printf("------Start SH_ReEnc process------\n");
    SH_ABE_ct new_ct = SH_ReEnc(rekey, ct, pp, pk);
    printf("ReEnc: 重加密后密文主组件: %lld\n", (long long)new_ct.sk_f_int);
    printf("------End SH_ReEnc process------\n");
    return new_ct;
}

static inline unsigned long long read_cycles(void) {
    unsigned int lo, hi;
    __asm__ __volatile__("rdtsc" : "=a"(lo), "=d"(hi));
    return ((unsigned long long)hi << 32) | lo;
}

void set_params(int N, int q, int k, int sigma, int attnum) {
// 根据实际情况更新全局参数 PARAMS
    PARAMS.N = N;
    PARAMS.Q = q;
    PARAMS.K = k;
    PARAMS.SIGMA = sigma;
    PARAMS.Att_num = attnum;
}
int main(void) {
    // 测试 Setup
    // 测试 Setup
    // Define a list of parameter values to test.
    // 这里我们用参数数组来改变某些函数中的参数（例如user_index和x）
    // int param_values[] = {5, 10, 15, 20, 25};
    // Helper: Reads the CPU cycle counter (requires an x86-compatible processor)

typedef struct {
    int N;
    int q;      // q 为素数，且须满足 ceil(log2(q)) == k
    int k;      // 分别取 10、15、20
    int sigma;  // sigma 的值在 20~30
    int attnum; // 对应 P，取值不小于 48
} TestParams;

TestParams test_params[] = {
    // k = 10, q = 523
    {64,   523, 10, 20, 48},
    {64,   523, 10, 21, 49},
    {128,  523, 10, 22, 50},
    {128,  523, 10, 23, 51},
    {256,  523, 10, 24, 52},
    {256,  523, 10, 25, 53},
    {64,   523, 10, 26, 54},
    {128,  523, 10, 27, 55},
    {256,  523, 10, 28, 56},
    {64,   523, 10, 29, 57},
    // k = 15, q = 16387
    {64,   16387, 15, 20, 48},
    {64,   16387, 15, 21, 49},
    {128,  16387, 15, 22, 50},
    {128,  16387, 15, 23, 51},
    {256,  16387, 15, 24, 52},
    {256,  16387, 15, 25, 53},
    {64,   16387, 15, 26, 54},
    {128,  16387, 15, 27, 55},
    {256,  16387, 15, 28, 56},
    {64,   16387, 15, 29, 57},
    // k = 20, q = 1048573
    {64,   1048573, 20, 20, 48},
    {64,   1048573, 20, 21, 49},
    {128,  1048573, 20, 22, 50},
    {128,  1048573, 20, 23, 51},
    {256,  1048573, 20, 24, 52},
    {256,  1048573, 20, 25, 53},
    {64,   1048573, 20, 26, 54},
    {128,  1048573, 20, 27, 55},
    {256,  1048573, 20, 28, 56},
    {64,   1048573, 20, 29, 57}
};

int num_params = sizeof(test_params) / sizeof(test_params[0]);

    // Optionally, if your implementation supports updating parameters,
    // loop through each test group and set the system parameters accordingly.
    // For example:
    // for (int i = 0; i < num_params; i++) {
    //     set_params(test_params[i].N, test_params[i].q, test_params[i].k,
    //                test_params[i].sigma, test_params[i].attnum);
    // }

    for (int i = 0; i < num_params; i++) {
        printf("\n测试组 %d: N=%d, q=%d, k=%d, sigma=%d, attnum=%d\n", i, test_params[i].N, test_params[i].q, test_params[i].k, test_params[i].sigma, test_params[i].attnum);
        // 更新参数
        // set_params(test_params[i].N, test_params[i].q, test_params[i].k,
        //            test_params[i].sigma, test_params[i].attnum);
        init_params(test_params[i].N, test_params[i].q, test_params[i].k, 1, test_params[i].sigma, test_params[i].attnum);

        printf("\n===== Testing Parameter Value: N=%d, q=%d, k=%d, sigma=%d, attnum=%d =====\n", test_params[i].N, test_params[i].q, test_params[i].k, test_params[i].sigma, test_params[i].attnum);
        print_params();
        struct timespec start, end;
        unsigned long long start_cycles, end_cycles;
        double setup_time, keygen_time, enc1_time, enc_time, dec_time, rekeygen_time, reenc_time;
        unsigned long long setup_cycles, keygen_cycles, enc1_cycles, enc_cycles, dec_cycles, rekeygen_cycles, reenc_cycles;

        // 测试 Setup
        clock_gettime(CLOCK_MONOTONIC, &start);
        start_cycles = read_cycles();
        SH_ABE_pp pp = test_SH_ABE_Setup();
        clock_gettime(CLOCK_MONOTONIC, &end);
        end_cycles = read_cycles();
        setup_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        setup_cycles = end_cycles - start_cycles;

        // 测试 KeyGen（使用 current_param 作为 user_index）
        clock_gettime(CLOCK_MONOTONIC, &start);
        start_cycles = read_cycles();
        SH_ABE_key_pair key_pair = test_SH_ABE_KeyGen(pp, i);
        clock_gettime(CLOCK_MONOTONIC, &end);
        end_cycles = read_cycles();
        keygen_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        keygen_cycles = end_cycles - start_cycles;

        // 测试 Encryption Step1
        clock_gettime(CLOCK_MONOTONIC, &start);
        start_cycles = read_cycles();
        SH_ABE_ct_one ct_one = test_SH_ABE_Enc_step1(key_pair.pk);
        clock_gettime(CLOCK_MONOTONIC, &end);
        end_cycles = read_cycles();
        enc1_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        enc1_cycles = end_cycles - start_cycles;

        // 测试 Encryption
        clock_gettime(CLOCK_MONOTONIC, &start);
        start_cycles = read_cycles();
        SH_ABE_ct ct = test_SH_ABE_Enc(key_pair.pk);
        clock_gettime(CLOCK_MONOTONIC, &end);
        end_cycles = read_cycles();
        enc_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        enc_cycles = end_cycles - start_cycles;

        // 测试 Decryption
        clock_gettime(CLOCK_MONOTONIC, &start);
        start_cycles = read_cycles();
        bool dec_status = test_dec(key_pair, ct_one, pp);
        clock_gettime(CLOCK_MONOTONIC, &end);
        end_cycles = read_cycles();
        dec_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
        dec_cycles = end_cycles - start_cycles;

        // 测试 Re-KeyGen 和 ReEnc（current_param 用作 x 参数）
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
            int64_t x = rand()%PARAMS.Q; // 使用 current_param 作为重加密的参数 x

            clock_gettime(CLOCK_MONOTONIC, &start);
            start_cycles = read_cycles();
            SH_ABE_ReEnc_key rekey = SH_ABE_Re_KeyGen(clauses, num_clauses, key_pair.sk, key_pair.pk, x);
            clock_gettime(CLOCK_MONOTONIC, &end);
            end_cycles = read_cycles();
            rekeygen_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
            rekeygen_cycles = end_cycles - start_cycles;

            clock_gettime(CLOCK_MONOTONIC, &start);
            start_cycles = read_cycles();
            SH_ABE_ct new_ct = test_SH_ReEnc(rekey, ct, pp, key_pair.pk);
            clock_gettime(CLOCK_MONOTONIC, &end);
            end_cycles = read_cycles();
            reenc_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
            reenc_cycles = end_cycles - start_cycles;
        }

        // 打印各个功能的物理时间和逻辑 cycle 数
        printf("\n------ Timing Results for Parameter Value: N=%d, q=%d, k=%d, sigma=%d, attnum=%d =====\n",i, test_params[i].N, test_params[i].q, test_params[i].k, test_params[i].sigma, test_params[i].attnum);
        printf("| %-15s | %-30s |\n", "Step", "Time (sec) / Cycles");
        printf("---------------------------------------------------------------\n");
        printf("| %-15s | %-15f sec / %-15llu           |\n", "Setup", setup_time, setup_cycles);
        printf("| %-15s | %-15f sec / %-15llu           |\n", "KeyGen", keygen_time, keygen_cycles);
        printf("| %-15s | %-15f sec / %-15llu           |\n", "Enc_step1", enc1_time, enc1_cycles);
        printf("| %-15s | %-15f sec / %-15llu           |\n", "Enc", enc_time, enc_cycles);
        printf("| %-15s | %-15f sec / %-15llu           |\n", "Dec", dec_time, dec_cycles);
        printf("| %-15s | %-15f sec / %-15llu           |\n", "Re-KeyGen", rekeygen_time, rekeygen_cycles);
        printf("| %-15s | %-15f sec / %-15llu(8 cores)     |\n", "ReEnc", reenc_time, reenc_cycles);
        printf("---------------------------------------------------------------\n");
    }
    return 0;
}