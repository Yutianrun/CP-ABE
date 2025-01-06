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
    circuit** prf_output = build_prf_circuit(k);

    printf("Circuit : ");
    print_circuit(**prf_output);
    printf("\n");
    matrix Af = compute_Af(A, **prf_output);

    int x_max = 1;
    for (int i = 0; i < PARAMS.K; i++) x_max *= 2;
    char output[80];

    matrix T = new_matrix(PARAMS.N, PARAMS.L);
    matrix BIG = new_matrix(PARAMS.N, PARAMS.L * PARAMS.K);
    for (attribute x = 0; x < x_max; x++) {
        printf("f(%d)=%d\n", x, compute_f(**prf_output, x));
        sprintf(output, "BIG * H = Af + f(x)G for x = %d : done in %%fs\n", x);
        CHRONO(output, {
            matrix H = compute_H(A, **prf_output, x);
            printf("Norm H : %f\n", norm(H));
            matrix R = copy_matrix(Af);
            if (compute_f(**prf_output, x)) add_matrix(R, G, R);

            for (int i = 1; i < PARAMS.K + 1; i++) {
                matrix ti = copy_matrix(A[i]);
                if (get_xn(x, i)) add_matrix(ti, G, ti);
                for (int j = 0; j < PARAMS.N; j++)      // ti.rows = PARAMS.N
                    for (int k = 0; k < PARAMS.L; k++)  // ti.columns = PARAMS.L
                        matrix_element(BIG, j, (i - 1) * PARAMS.L + k) =
                            matrix_element(ti, j, k);
                free_matrix(ti);
            }
            mul_matrix(BIG, H, T);

            assert(equals(R, T));
            free_matrix(H);
            free_matrix(R);
        });
    }


    free_matrix(BIG);
    free_matrix(T);

    free_matrix(Af);

    free_matrixes(A, PARAMS.K + 1);
    free_matrix(G);
}