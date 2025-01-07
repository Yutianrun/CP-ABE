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
// #define PRF_K 128



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
    int k = 4; // 比特宽度，可以根据需要调整

    // 构建 PRF 电路
    // 定义字句集合
    int num_clauses = 2;
    Clause* clauses = (Clause*)malloc(num_clauses * sizeof(Clause));
    if (clauses == NULL) {
        fprintf(stderr, "内存分配失败\n");
        exit(1);
    }

    // 示例字句：(T1, v1), (T2, v2)
    // 例如，T1 = {1, 2}, T2 = {2, 3}
    for(int i = 0; i < num_clauses; i++) {
        clauses[i].t_size = 2;
        clauses[i].T = (int*)malloc(clauses[i].t_size * sizeof(int));
        if (clauses[i].T == NULL) {
            fprintf(stderr, "内存分配失败\n");
            exit(1);
        }
    }

    // 初始化字句
    clauses[0].T[0] = 1;
    clauses[0].T[1] = 2;
    clauses[1].T[0] = 2;
    clauses[1].T[1] = 3;

    // 构建复杂 PRF 电路
    circuit* eval = build_eval_prf(3, clauses, num_clauses);

    printf("Circuit : ");
    print_circuit(*eval);
    printf("\n");


    circuit* prf_output[k];
    for (int i = 0; i < k; i++) {
        prf_output[i] = build_eval_prf(3, clauses, num_clauses);
    }

    int x_max = 1;
    for (int i = 0; i < PARAMS.K; i++) x_max *= 2;
    char output[80];

    matrix T = new_matrix(PARAMS.N, PARAMS.L);
    matrix BIG = new_matrix(PARAMS.N, PARAMS.L * PARAMS.K);
    for (attribute x = 0; x < x_max; x++) {
        char concatenated_output[256] = "";
        for (attribute i = 0; i < k; i++) {
            sprintf(concatenated_output + strlen(concatenated_output), "%d", compute_f(*prf_output[i], x));
        }
        printf("prf_output: f(%d)=%s\n", x, concatenated_output);
        
        for (int i = 0; i < k; i++) {
            sprintf(output, "BIG * H = Af + f(x)G for x = %d, bit %d = %d : done in %%fs\n", x, i, compute_f(*prf_output[i], x));
            CHRONO(output, {

            matrix Af = compute_Af(A, *prf_output[i]);
            matrix H = compute_H(A, *prf_output[i], x);
            printf("Norm H : %f\n", norm(H));
            matrix R = copy_matrix(Af);
            if (compute_f(*prf_output[i], x)) add_matrix(R, G, R);

            for (int j = 1; j < PARAMS.K + 1; j++) {
                matrix ti = copy_matrix(A[j]);
                if (get_xn(x, j)) add_matrix(ti, G, ti);
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

    for (int i = 0; i < k; i++) {
        free_circuit(prf_output[i]);
    }

    free_matrixes(A, PARAMS.K + 1);
    free_matrix(G);
}