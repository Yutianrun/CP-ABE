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
#include <stdlib.h>
// #define PRF_K 128



int main() {
    init_sampler();
    init_params_default();
    init_G();
    real start, end;
    int prf_k = 8; // 比特宽度，可以根据需要调整

    // Printing parameters
    printf("Testing eval circuit with parameters\n");
    print_params();

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


       // 构建 PRF 电路

    circuit** x = (circuit**)malloc((prf_k) * sizeof(circuit*));
    circuit** msk = (circuit**)malloc((prf_k) * sizeof(circuit*));
    for(int i = 0;i< prf_k; i++){
        x[i] = gen_leaf(i+1, false);
        msk[i] = gen_leaf(i+1+prf_k, false);
    }
    circuit** prf_output = build_eval_circuit(prf_k, clauses, num_clauses, msk, x);

    // printf("compute Circuit : %d\n", compute_f(*prf_output[0], 0));


       // Generating A
    matrix* A = new_matrixes(PARAMS.K + 1, PARAMS.N, PARAMS.L);
    CHRONO("Generated A in %fs\n", {
        for (int i = 0; i < PARAMS.K + 1; i++) sample_Zq_uniform_matrix(A[i]);
    });

    matrix currAf = compute_Af(A, *prf_output[0]);

    // Testing G * G^-1(A) = A
    // matrix inv = new_matrix(PARAMS.L, PARAMS.L);
    // matrix res = new_matrix(PARAMS.N, PARAMS.L);
    // CHRONO("Checked G * G^-1(A) = A in %fs\n", {
    //     inv_G(A[0], inv);
    //     mul_matrix(G, inv, res);
    //     assert(equals(A[0], res));
    // });
    // free_matrix(inv);
    // free_matrix(res);


    // printf("Circuit : ");
    // print_circuit(**prf_output);
    // printf("\n");
    int x_max = 1;
    for (int i = 0; i < prf_k/2; i++) x_max *= 2;

    char output[80];
    matrix T = new_matrix(PARAMS.N, PARAMS.L);
    matrix BIG = new_matrix(PARAMS.N, PARAMS.L * PARAMS.K);


    uint32_t mask = rand() % (1 << prf_k); // 随机生成一个kbit的掩码
    for (attribute x = 0; x < x_max; x++) { // 前k位遍历0-x_max
        uint32_t input = (x << prf_k) | mask; // 组合前k位和后k位
        char concatenated_output[256] = "";
        for (attribute i = 0; i < prf_k; i++) {
            sprintf(concatenated_output + strlen(concatenated_output), "%d", compute_f(*prf_output[i], input));
        }

        // printf("computed over\n");
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