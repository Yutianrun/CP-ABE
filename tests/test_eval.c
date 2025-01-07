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

#define PRF_K 128


int main() {
    init_params_default();
    // Printing parameters
    printf("Testing circuit with parameters\n");
    print_params();

    // 使用辅助函数简化主电路
    int prf_k = 8; // 比特宽度，可以根据需要调整

    int num_clauses = 2;
    Clause* clauses = (Clause*)malloc(num_clauses * sizeof(Clause));
    clauses[0].T = (int*)malloc(2 * sizeof(int));
    clauses[0].t_len = 2;
    clauses[0].T[0] = 0;
    clauses[0].T[1] = 2;
    clauses[1].t_len = 2;
    clauses[1].T = (int*)malloc(2 * sizeof(int));
    clauses[1].T[0] = 2;
    clauses[1].T[1] = 3;

    // 构建 PRF 电路

    circuit** x = (circuit**)malloc((prf_k) * sizeof(circuit*));
    circuit** msk = (circuit**)malloc((prf_k) * sizeof(circuit*));
    for(int i = 0;i< prf_k; i++){
        x[i] = gen_leaf(i+1, false);
        msk[i] = gen_leaf(i+1+prf_k, false);
    }

    circuit** eval = build_eval_circuit(prf_k, clauses, num_clauses, msk, x);

    // printf("Circuit : ");
    // print_circuit(*eval[3]);
    // printf("\n");
    // PARAMS.K = p r f;
    int x_max = 1;
    for (int i = 0; i < prf_k/2; i++) x_max *= 2;


    // printf("\n%d\n", compute_f(*eval[0], 300));
    uint32_t mask = rand() % (1 << prf_k); // 随机生成一个kbit的掩码
    for (attribute x = 0; x < x_max; x++) { // 前k位遍历0-x_max
        uint32_t input = (x << prf_k) | mask; // 组合前k位和后k位
        // uint8_t input = x;
        char concatenated_output[256] = "";
        for (attribute i = 0; i < prf_k; i++) {
            // input = 0;
            sprintf(concatenated_output + strlen(concatenated_output), "%d", compute_f(*eval[i], input));
        }
        char binary_input[2*prf_k+1];
        for (int j = 0; j < 2*prf_k; j++) {
            binary_input[2*prf_k - 1 - j] = (input & (1 << j)) ? '1' : '0';
        }
        binary_input[2*prf_k] = '\0';
        printf("eval: f(%d)=%s (binary: %.*s %.*s)\n", input, concatenated_output, prf_k, binary_input, prf_k, binary_input + prf_k);
    }
    

}