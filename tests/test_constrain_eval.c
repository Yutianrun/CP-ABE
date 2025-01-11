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

#define PRF_K 16



int main() {
    init_params_default();
    // Printing parameters

    printf("Testing circuit with parameters\n");
    print_params();

    // 使用辅助函数简化主电路
    int prf_k = 8; // 比特宽度，可以根据需要调整

    int num_clausesT = 2;
    int num_clausesF = 3;

    ClauseT* clauses = (ClauseT*)malloc(num_clausesT * sizeof(ClauseT));
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
    for(int i = 0;i< prf_k; i++){
        x[i] = gen_leaf(i+1, false);
    }

    circuit*** sk_f = malloc(num_clausesF * sizeof(circuit**));
    for (int i = 0; i <  num_clausesF; i++) {
        sk_f[i] = malloc(2* prf_k * sizeof(circuit*));
        for(int j = 0; j < 2* prf_k; j++) {
            sk_f[i][j] = gen_leaf(prf_k*i+j+1, true);
        }
    }
    int S_len ; // 示例长度

    circuit** eval  =build_constrain_eval_circuit(clauses, num_clausesF, num_clausesT, prf_k, sk_f, x);

    // printf("S_len = %d\n", S_len);
    
    // printf("Circuit : ");
    // print_circuit(*sk_f[0][3]);
    // printf("\n");

    int x_max = 4;
    // for (int i = 0; i < prf_k/2; i++) x_max *= 2;


    // printf("\n%d\n", compute_f(*eval[0], 300));
    // uint32_t mask = rand() % (1 << prf_k); // 随机生成一个kbit的掩码

    for (attribute x = 0; x < x_max; x++) {
        // printf("x = %d\n", x);
        // 前k位遍历0-x_max
        // uint32_t input = (x << prf_k) | mask; // 组合前k位和后k位
        uint32_t input = x;
        char concatenated_output[256] = "";
        for (attribute i = 0; i <  prf_k; i++) {
            // input = 0;
            // printf("u = %d, i = %d\n", u, i);
            sprintf(concatenated_output + strlen(concatenated_output), "%d", compute_f(*eval[i], input));
            if (i == prf_k - 1) {
            sprintf(concatenated_output + strlen(concatenated_output), " ");
            }
            // sprintf(concatenated_output + strlen(concatenated_output), "%d", compute_f(*msk[i], input));
        }
        int total_bits = prf_k * 2 * num_clausesF;
        int spaces = total_bits / 8;
        char binary_input[total_bits + spaces + 1];
        int space_interval = 8;
        int index = 0;
        for (int j = total_bits - 1; j >= 0; j--) {
            binary_input[index++] = (input & ((uint64_t)1 << j)) ? '1' : '0';
            if ((j) % space_interval == 0 && j != 0) {
                binary_input[index++] = ' ';
            }
        }
        binary_input[index] = '\0';
        printf("constrain eval f(%d)=%s (binary:%s)\n", input, concatenated_output, binary_input);
    }
}