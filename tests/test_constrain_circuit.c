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


bool simple_function(bool* input) {
    return input[0] ^ input[1];
}

bool simple_function_clasuse2(bool* input) {
    return input[0] & input[1];
}

int main() {
    init_sampler();
    init_params_default();
    init_G();
    real start, end;
    int prf_k = 8; // 比特宽度，可以根据需要调整

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
    // Printing parameters
    printf("Testing circuit with parameters\n");
    print_params();


    circuit** msk = (circuit**)malloc((prf_k) * sizeof(circuit*));
    for(int i = 0;i< prf_k; i++){
        // x[i] = gen_leaf(i+1, false);
        msk[i] = gen_leaf(i+1, true);
    }
    int S_len ; // 示例长度
    Pair* S = (Pair*)malloc(S_len * sizeof(Pair));

    circuit*** prf_output = build_constrain_circuit(clauses, num_clauses, S, &S_len, prf_k, msk);

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
       // 构建 PRF 电路

    // printf("Circuit : ");
    // print_circuit(**prf_output);
    printf("\n");
    int x_max = 1;
    for (int i = 0; i < prf_k/2; i++) x_max *= 2;

    char output[80];
    matrix T = new_matrix(PARAMS.N, PARAMS.L);
    matrix BIG = new_matrix(PARAMS.N, PARAMS.L * PARAMS.K);


    uint32_t mask = rand() % (1 << prf_k); // 随机生成一个kbit的掩码
    for (attribute x = 0; x < x_max; x++) { // 前k位遍历0-x_max
        for(int u = 0; u < S_len; u++) {
            uint32_t input = x;
            char concatenated_output[256] = "";
            for (attribute i = 0; i < prf_k; i++) {
                sprintf(concatenated_output + strlen(concatenated_output), "%d", compute_f(*prf_output[u][i], input));
            }
            char binary_input[prf_k+1];
            for (int j = 0; j < prf_k; j++) {
                binary_input[prf_k - 1 - j] = (input & (1 << j)) ? '1' : '0';
            }
            binary_input[prf_k] = '\0';
            printf("constrain f(%d)=%s (binary:%.*s)\n", input, concatenated_output, prf_k, binary_input);

            for (int i = 0; i < prf_k; i++) {
                sprintf(output, "BIG * H = Af + f(x)G for x = %d, bit %d = %d : done in %%fs\n", input, i, compute_f(*prf_output[u][i], input));
                CHRONO(output, {
                    matrix Af = compute_Af(A, *prf_output[u][i]);
                    matrix H = compute_H(A, *prf_output[u][i], input);
                    printf("Norm H : %f\n", norm(H));
                    matrix R = copy_matrix(Af);
                    if (compute_f(*prf_output[u][i], input)) add_matrix(R, G, R);

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
    }


    free_matrix(BIG);
    free_matrix(T);

    free_matrixes(A, PARAMS.K + 1);
    free_matrix(G);
}