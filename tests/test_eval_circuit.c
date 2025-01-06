#include <assert.h>
#include <stdio.h>
#include <time.h>

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
    int k = 8; // 比特宽度，可以根据需要调整

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

    // printf("Circuit : ");
    // print_circuit(*eval);
    // printf("\n");
    // matrix Af = compute_Af(A, *eval);

    // int x_max = 1;
    // for (int i = 0; i < PARAMS.K; i++) x_max *= 2;
    // char output[80];

    // matrix T = new_matrix(PARAMS.N, PARAMS.L);
    // matrix BIG = new_matrix(PARAMS.N, PARAMS.L * PARAMS.K);
    // for (attribute x = 0; x < x_max; x++) {
    //     printf("f(%d)=%d\n", x, compute_f(*eval, x));
    //     sprintf(output, "BIG * H = Af + f(x)G for x = %d : done in %%fs\n", x);
    //     CHRONO(output, {
    //         matrix H = compute_H(A, *eval, x);
    //         printf("Norm H : %f\n", norm(H));
    //         matrix R = copy_matrix(Af);
    //         if (compute_f(*eval, x)) add_matrix(R, G, R);

    //         for (int i = 1; i < PARAMS.K + 1; i++) {
    //             matrix ti = copy_matrix(A[i]);
    //             if (get_xn(x, i)) add_matrix(ti, G, ti);
    //             for (int j = 0; j < PARAMS.N; j++)      // ti.rows = PARAMS.N
    //                 for (int k = 0; k < PARAMS.L; k++)  // ti.columns = PARAMS.L
    //                     matrix_element(BIG, j, (i - 1) * PARAMS.L + k) =
    //                         matrix_element(ti, j, k);
    //             free_matrix(ti);
    //         }
    //         mul_matrix(BIG, H, T);

    //         assert(equals(R, T));
    //         free_matrix(H);
    //         free_matrix(R);
    //     });
    // }



    // circuit** constrain = build__prf(k);

    // printf("Circuit : ");
    // print_circuit(*eval);
    // printf("\n");
    // matrix Af = compute_Af(A, *eval);

    // int x_max = 1;
    // for (int i = 0; i < PARAMS.K; i++) x_max *= 2;
    // char output[80];



    // free_matrix(BIG);
    // free_matrix(T);

    // free_matrix(Af);

    // free_matrixes(A, PARAMS.K + 1);
    // free_matrix(G);
}