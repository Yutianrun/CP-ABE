#include <assert.h>
#include <stdio.h>
#include <time.h>
#include <string.h>

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
    int prf_k = 16; // 比特宽度，可以根据需要调整

    // 构建 PRF 电路
    circuit** prf_output = build_prp_circuit(prf_k);

    // printf("Circuit : ");
    // print_circuit(**prf_output);
    // printf("\n");
    // PARAMS.K = k;
    int x_max = 1;
    for (int i = 0; i < prf_k; i++) x_max *= 2;

    uint32_t mask = rand() % (1 << prf_k); // 随机生成一个kbit的掩码
    for (attribute x = 0; x < x_max; x++) { // 前k位遍历0-x_max
        uint32_t input = (x << prf_k) | mask; // 组合前k位和后k位
        // uint8_t input = x;
        char concatenated_output[256] = "";
        for (attribute i = 0; i < prf_k; i++) {
            sprintf(concatenated_output + strlen(concatenated_output), "%d", compute_f(*prf_output[i], input));
        }
        char binary_input[2*prf_k+1];
        for (int j = 0; j < 2*prf_k; j++) {
            binary_input[2*prf_k - 1 - j] = (input & (1 << j)) ? '1' : '0';
        }
        binary_input[2*prf_k] = '\0';
        printf("prf_output: f(%d)=%s (binary: %.*s %.*s)\n", input, concatenated_output, prf_k, binary_input, prf_k, binary_input + prf_k);
    }
    

}