#include <assert.h>
#include <stdio.h>
#include <time.h>

#include "attribute.h"
#include "circuit.h"
#include "common.h"
#include "gen_circuit.h"
#include "matrix.h"
#include "sampling.h"

int main() {
    init_sampler();
    int32_t N = 1;
    // int32_t K = 56;
    int32_t K = 20;
    // int32_t K = 4;
    // int64_t Q = 72057594037927936;
    int64_t Q = 870367;
    // int64_t Q = 16;


    int32_t P = 1;
    real SIGMA = 20;
    int att_num = 56;
    init_params(N, Q, K, P, SIGMA, att_num);
    // init_params_default();
    init_G();
    real start, end;

    // Printing parameters
    printf("Testing circuit with parameters\n");
    print_params();

    // Generating A
    matrix* A = new_matrixes(PARAMS.Att_num + 1, PARAMS.N, PARAMS.L);
    CHRONO("Generated A in %fs\n", {
        for (int i = 0; i < PARAMS.K + 1; i++) sample_Zq_uniform_matrix_64(A[i]);
    });

    // Testing G * G^-1(A) = A
    matrix inv = new_matrix(PARAMS.L, PARAMS.L);
    matrix res = new_matrix(PARAMS.N, PARAMS.L);


    CHRONO("Checked G * G^-1(A) = A in %fs\n", {
        inv_G(A[0], inv);
        mul_matrix(G, inv, res);
    });
    free_matrix(inv);
    free_matrix(res);

    circuit* f = circuit_xor(gen_leaf(1,true), gen_leaf(2,true));

    circuit** list = (circuit**)malloc(4 * sizeof(circuit*));
    for (size_t i = 0; i < 4; i++)
    {
        list[i] = gen_leaf(i+1, true);
    }

    printf("Circuit : ");
    print_circuit(*f);
    printf("\n");

    matrix Af = compute_Af(A, *f);

    matrix T = new_matrix(PARAMS.N, PARAMS.L);
    matrix BIG = new_matrix(PARAMS.N, PARAMS.L * PARAMS.K);
    for (int i = 1; i < PARAMS.K + 1; i++) {
    matrix ti = copy_matrix(A[i]); 
    for (int j = 0; j < PARAMS.N; j++)
        for (int k = 0; k < PARAMS.L; k++)
            matrix_element(BIG, j, (i-1) * PARAMS.L + k) = matrix_element(ti, j, k);
    free_matrix(ti);
    }


    int x_max = 16;
    // for (int i = 0; i < PARAMS.K; i++) x_max *= 2;
    char output[80];

    attribute x = 3;

    matrix H = compute_H(A, *f, x);
    // matrix H_not_hat = compute_not_hat_H(A, *f);
    matrix H_from_A_Af = compute_H_from_A_Af(&BIG, &Af);


     // 添加调试信息
    printf("BIG dimensions: %d x %d\n", BIG.rows, BIG.columns);
    printf("H_from_A_Af dimensions: %d x %d\n", H_from_A_Af.rows, H_from_A_Af.columns);
    printf("Af dimensions: %d x %d\n", Af.rows, Af.columns);
    // Test that A * H_from_A_Af == Af
    matrix test_result = new_matrix(PARAMS.N, PARAMS.L);
    mul_matrix(BIG, H_from_A_Af, test_result);
    assert(equals(Af, test_result));
    printf("Test A * H_from_A_Af == Af passed\n");
    free_matrix(test_result);
    // for (attribute x = 0; x < x_max; x++) {
    //     printf("f(%d)=%d\n", x, compute_f(*f, x));
    //     sprintf(output, "BIG * H = Af + f(x)G for x = %d : done in %%fs\n", x);
    //     // CHRONO(output, {
    //         matrix H = compute_H(A, *f, x);
    //         matrix H_not_hat = compute_not_hat_H(A, *f);
    //         printf("Norm H : %f\n", norm(H));
    //         matrix R = copy_matrix(Af);
    //         if (compute_f(*f, x)) add_matrix(R, G, R);

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
    //     // });
    // }
    free_matrix(H);
     
     
    free_matrix(H_from_A_Af);

    free_matrix(BIG);
    free_matrix(T);

    free_matrix(Af);

    free_matrixes(A, PARAMS.K + 1);
    free_matrix(G);
}