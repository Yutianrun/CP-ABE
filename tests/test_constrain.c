#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

#include "attribute.h"
#include "circuit.h"
#include "common.h"
#include "gen_circuit.h"
#include "matrix.h"
#include "sampling.h"
#include "cprf.h"

#define PRF_K 8


bool simple_function(bool* input) {
    return input[0] ^ input[1];
}

bool simple_function_clasuse2(bool* input) {
    return input[0] & input[1];
}

int main() {
    init_params_default();
    // Printing parameters

    bool test_input = {0,0};
    printf("test simple fucntion : %d\n", simple_function(&test_input));
    printf("Testing circuit with parameters\n");
    print_params();

    // 使用辅助函数简化主电路
    int prf_k = PRF_K; // 比特宽度，可以根据需要调整

    int num_clauses = 2;
    Clause* clauses = (Clause*)malloc(num_clauses * sizeof(Clause));

    clauses[0].f = simple_function;
    clauses[1].f = simple_function_clasuse2;
    clauses[0].T = (int*)malloc(2 * sizeof(int));
    clauses[0].t_len = 2;
    clauses[0].T[0] = 0;
    clauses[0].T[1] = 2;
    clauses[1].T = (int*)malloc(2 * sizeof(int));
    clauses[1].t_len = 2;
    clauses[1].t_len = 2;
    clauses[1].T = (int*)malloc(2 * sizeof(int));
    clauses[1].T[0] = 2;
    clauses[1].T[1] = 3;

    // 构建 PRF 电路
    

    // circuit** x = (circuit**)malloc((prf_k) * sizeof(circuit*));
    circuit** msk = (circuit**)malloc((prf_k) * sizeof(circuit*));
    for(int i = 0;i< prf_k; i++){
        // x[i] = gen_leaf(i+1, false);
        msk[i] = gen_leaf(i+1, true);
    }
    int S_len ; // 示例长度
    Pair* S = (Pair*)malloc(S_len * sizeof(Pair));

    circuit*** sk_f = build_constrain_circuit(clauses, num_clauses, S, &S_len, prf_k, msk);

    printf("S_len = %d\n", S_len);
    
    printf("Circuit : ");
    print_circuit(*sk_f[0][3]);
    printf("\n");

    int x_max = 104;
    // for (int i = 0; i < prf_k/2; i++) x_max *= 2;


    // printf("\n%d\n", compute_f(*eval[0], 300));
    // uint32_t mask = rand() % (1 << prf_k); // 随机生成一个kbit的掩码

    for (attribute x = 103; x < x_max; x++) {
        for(int u = 0; u < S_len; u++) { // 前k位遍历0-x_max
        // uint32_t input = (x << prf_k) | mask; // 组合前k位和后k位
        uint32_t input = x;
        char concatenated_output[256] = "";
        for (attribute i = 0; i < 2 * prf_k; i++) {
            // input = 0;
            // printf("u = %d, i = %d\n", u, i);
            sprintf(concatenated_output + strlen(concatenated_output), "%d", compute_f(*sk_f[u][i], input));
            if (i == prf_k - 1) {
            sprintf(concatenated_output + strlen(concatenated_output), " ");
            }
            // sprintf(concatenated_output + strlen(concatenated_output), "%d", compute_f(*msk[i], input));
        }
        char binary_input[prf_k+1];
        for (int j = 0; j < prf_k; j++) {
            binary_input[prf_k - 1 - j] = (input & (1 << j)) ? '1' : '0';
        }
        binary_input[prf_k] = '\0';
        printf("constrain f(%d)=%s (binary:%.*s)\n", input, concatenated_output, prf_k, binary_input);
    }
}

}