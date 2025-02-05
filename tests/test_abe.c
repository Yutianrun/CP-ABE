// Language: c
#include <assert.h>
#include <assert.h>
#include <stdio.h>
#include "abe.h"
#include "cprf.h"
#define PRF_K 8

// Test ABE_Setup
ABE_setup_result  test_ABE_Setup() {
    uint64_t test_msk = 103;
    ABE_setup_result res = ABE_Setup(test_msk);
        
    // Check public parameters
    assert(res.pp.B.data != NULL);
    assert(res.pp.A[0].data != NULL);
    assert(res.pp.v.data != NULL);

    // Check master secret key
    assert(res.msk_out.B_tau.data != NULL);
    assert(res.msk_out.sigma == test_msk);
    
    printf("ABE_Setup: sigma = %llu\n", (unsigned long long)res.msk_out.sigma);

    // Free allocated resources
    // free_matrix(res.pp.B);
    // free_matrixes(res.pp.A, PARAMS.K + 1);
    // free_matrix(res.pp.v);
    // free_matrix(res.msk_out.B_tau);


    printf("ABE_Setup test passed\n");
    return res;
}

bool simple_function(bool* input) {
    return input[0] ^ input[1];
}

bool simple_function_clasuse2(bool* input) {
    return input[0] & input[1];
}

// 重写后的测试 ABE_KeyEnc 函数
void test_ABE_KeyEnc(ABE_setup_result setup_res) {
    uint64_t test_msk = 103;
    int prf_k = PRF_K; // PRF 门电路比特宽度
    int num_clauses = 2;

    // 分配并初始化 dummy clauses
    Clause* clauses = (Clause*)malloc(num_clauses * sizeof(Clause));
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

    // 构造 dummy PRF leaf circuits 用于主密钥
    circuit** msk_circuits = (circuit**)malloc(prf_k * sizeof(circuit*));
    for (int i = 0; i < prf_k; i++) {
        msk_circuits[i] = gen_leaf(i + 1, true);
    }

    // 分配 dummy S 并构造约束电路树
    int S_len = prf_k;  // 为测试起见，假设 S 的长度等于 prf_k
    Pair* S = (Pair*)malloc(S_len * sizeof(Pair));
    circuit*** sk_f = build_constrain_circuit(clauses, num_clauses, S, &S_len, prf_k, msk_circuits);

    // 获取 ABA_Setup 的输出
    // ABE_setup_result setup_res = ABE_Setup(test_msk);

    // 测试 flag u = true
    ABE_ct ct = ABE_KeyEnc(sk_f, setup_res.msk_out, S_len, setup_res.pp, true);
    // 测试 flag u = true
    printf("Ciphertext field: sf\n");
    for (int i = 0; i < S_len; i++) {
        printf("sf[%d] =", i);
        for(int j = 0; j < 2 * prf_k; j++) {
            printf(" %d",ct.sf[i * 2 * prf_k + j]);
        }
        printf("\n");
    }

    printf("Ciphertext field: u0\n");
    print_matrix(ct.u0);

    printf("Ciphertext field: u1\n");
    print_matrix(ct.u1);

    printf("Ciphertext field: u2 = %lld\n", (long long)ct.u2);
    

    // 释放分配的资源
    // free_matrix(setup_res.pp.B);
    // free_matrixes(setup_res.pp.A, PARAMS.K + 1);
    // free_matrix(setup_res.pp.v);
    // free_matrix(setup_res.msk_out.B_tau);
    
    // 释放 clauses 内部分配的数组
    for (int i = 0; i < num_clauses; i++) {
        free(clauses[i].T);
    }
    free(clauses);
    free(msk_circuits);
    free(S);
    // 根据 build_constrain_circuit 的实现，还可能需要释放 sk_f 内部结构，此处假设有对应的释放函数：
    // free_sk_f(sk_f, num_clauses, prf_k);
}

void test_ABE_KeyGen(ABE_setup_result setup_res) {
    printf("Testing ABE_KeyGen...\n");

    // Generate keys using master secret and public parameters from setup
    ABE_keys keys = ABE_KeyGen(setup_res.msk_out, setup_res.pp, 48);

    // Validate that the produced key components are allocated.
    assert(keys.A != NULL);
    assert(keys.Tf.data != NULL);

    // Print the secret key matrix Tf and first public key component as an example.
    // printf("Secret key Tf matrix:\n");
    // print_matrix(keys.Tf);
    
    // printf("Public key component A[0]:\n");
    // print_matrix(keys.A[0]);

    printf("ABE_KeyGen test passed\n");

    // Free keys resources if necessary (assuming corresponding free functions exist)
    // free_matrix(keys.Tf);
    // free_matrixes(keys.A, PARAMS.K + 1);
}

int main(void) {

    ABE_setup_result setup_res = test_ABE_Setup();
    test_ABE_KeyEnc(setup_res);
    test_ABE_KeyGen(setup_res);
    return 0;
}
 