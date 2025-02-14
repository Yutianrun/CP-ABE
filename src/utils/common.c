#include "common.h"

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

cp_params PARAMS;

static bool initialized = false;

void init_params(int32_t N, int64_t Q, int32_t K, int32_t P, real SIGMA, int32_t att_num) {
    // We can initialize parameters only once
    // assert(!initialized);

    // Checking parameters
    // assert(0 < N && N <= UINT32_MAX);
    // assert(0 < Q && Q <= UINT32_MAX);
    // assert(0 < K && K <= 32);
    assert(0 < N && N <= UINT64_MAX);
    assert(0 < Q && Q <= INT64_MAX);
    assert(0 < K && K <= 64);
    assert((uint64_t)ceil(log2(Q)) == (uint64_t)K);
    assert(0 < P && P <= 30);
    assert(0 < SIGMA);

    // Assigning parameters
    PARAMS.N = N;
    PARAMS.Q = Q;
    PARAMS.K = K;
    PARAMS.L = N * K;
    PARAMS.P = P;
    PARAMS.M = N * K;
    PARAMS.SIGMA = SIGMA;
    PARAMS.MBAR = PARAMS.L / 8;
    // Empirical formula from is_short tests
    real SHORT_THRESHOLD = (real)N * (Q / 2) * pow(P, 1.45) * pow(SIGMA, 2);
    PARAMS.Att_num = att_num;
    // Arbitrary factor to compensate reality so it works better in practice
    // ie for a wider range of parameters
    SHORT_THRESHOLD /= 10;
    PARAMS.SHORT_THRESHOLD = SHORT_THRESHOLD;

    // Not possible to initialize it again
    initialized = true;
}

void init_params_default() {
    // int32_t N = 1;
    // real SIGMA = 20;
    // uint32_t Q = 1073707009;
    // int32_t K = 30;
    // int32_t K = 10;
    // uint32_t Q = 1024;

    // int32_t K = 16;
    // uint32_t Q = 65536;
    // int32_t K = 64;
    // uint64_t Q = 18446744073709551616;
    int32_t N = 1;
    int64_t Q = 870367;
    int32_t K = 20;
    // int32_t K = 56;
    // int32_t K = 4;
    // int64_t Q = 13;
    // int64_t Q = 72057594037927936;

    int32_t P = 1;
    real SIGMA = 20;
    // int32_t att_num = 56;
     int32_t att_num = 48;
    init_params(N, Q, K, P, SIGMA, att_num);
}

void print_params() {
    printf("Printing parameters...\n");
    printf("\tN = %d\n", PARAMS.N);
    printf("\tQ = %ld\n", PARAMS.Q);
    printf("\tK = %d\n", PARAMS.K);
    printf("\tL = %d\n", PARAMS.L);
    printf("\tP = %d\n", PARAMS.P);
    printf("\tM = P + 2 = %d\n", PARAMS.M);
    printf("\tSIGMA = %.2f\n", PARAMS.SIGMA);
    printf("\tSHORT_THRESHOLD = %.2g\n", PARAMS.SHORT_THRESHOLD);
}
