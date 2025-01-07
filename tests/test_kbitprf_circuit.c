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
    init_sampler();
    init_params_default();
    init_G();
    real start, end;

    // Printing parameters
    printf("Testing circuit with parameters\n");
    print_params();

    // Generating A
    matrix* A = new_matrixes(PARAMS.K + 1, PARAMS.N, PARAMS.L);
    CHRONO("Generated A in %fs\n", {
        for (int i = 0; i < PARAMS.K + 1; i++) sample_Zq_uniform_matrix(A[i]);
    });

    // Testing G * G^-1(A) = A
    matrix inv = new_matrix(PARAMS.L, PARAMS.L);
    matrix res = new_matrix(PARAMS.N, PARAMS.L);
    CHRONO("Checked G * G^-1(A) = A in %fs\n", {
        inv_G(A[0], inv);
        mul_matrix(G, inv, res);
        assert(equals(A[0], res));
    });
    free_matrix(inv);
    free_matrix(res);

    // 使用辅助函数简化主电路
    int prf_k = 8; // 比特宽度，可以根据需要调整
    // 构建 PRF 电路
    circuit** prf_output = build_prp_circuit(prf_k);

    // printf("Circuit : ");
    // print_circuit(**prf_output);
    // printf("\n");

    char output[80];

    matrix T = new_matrix(PARAMS.N, PARAMS.L);
    matrix BIG = new_matrix(PARAMS.N, PARAMS.L * PARAMS.K);

    int x_max = 1;
    for (int i = 0; i < prf_k; i++) x_max *= 2;

    uint32_t mask = rand() % (1 << prf_k); // 随机生成一个kbit的掩码
    for (attribute x = 0; x < x_max; x++) { // 前k位遍历0-x_max
        uint32_t input = (x << prf_k) | mask; // 组合前k位和后k位
        char concatenated_output[256] = "";
        for (attribute i = 0; i < prf_k; i++) {
            sprintf(concatenated_output + strlen(concatenated_output), "%d", compute_f(*prf_output[i], input));
        }
        char binary_input[2 * prf_k + 1];
        for (int j = 0; j < 2 * prf_k; j++) {
            binary_input[2 * prf_k - 1 - j] = (input & (1 << j)) ? '1' : '0';
        }
        binary_input[2 * prf_k] = '\0';
        printf("prf_output: f(%d)=%s (binary: %.*s %.*s)\n", input, concatenated_output, prf_k, binary_input, prf_k, binary_input + prf_k);

        for (int i = 0; i < prf_k; i++) {
            sprintf(output, "BIG * H = Af + f(x)G for x = %d, bit %d = %d : done in %%fs\n", input, i, compute_f(*prf_output[i], input));
            CHRONO(output, {
                matrix Af = compute_Af(A, *prf_output[i]);
                matrix H = compute_H(A, *prf_output[i], input);
                printf("Norm H : %f\n", norm(H));
                matrix R = copy_matrix(Af);
                if (compute_f(*prf_output[i], input)) add_matrix(R, G, R);

                for (int j = 1; j < PARAMS.K + 1; j++) {
                    matrix ti = copy_matrix(A[j]);
                    if (get_xn(input, j)) add_matrix(ti, G, ti);
                    for (int m = 0; m < PARAMS.N; m++)      // ti.rows = PARAMS.N
                    for (int n = 0; n < PARAMS.L; n++)  // ti.columns = PARAMS.L
                        matrix_element(BIG, m, (j - 1) * PARAMS.L + n) =
                        matrix_element(ti, m, n);
                    free_matrix(ti);
                }
                mul_matrix(BIG, H, T);

                assert(equals(R, T));
                free_matrix(H);
                free_matrix(Af);
                free_matrix(R);
            });
        }
    }


    free_matrix(BIG);
    free_matrix(T);

    free_matrixes(A, PARAMS.K + 1);
    free_matrix(G);
}