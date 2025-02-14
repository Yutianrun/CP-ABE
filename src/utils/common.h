#pragma once

#include <stdint.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
// Defining alias types
typedef uint64_t scalar;
typedef int64_t signed_scalar;
typedef double real;

// Struct containing all parameters of the CP-ABE
typedef struct _cp_params {
    int32_t N;   // should be 1 for now - degree of polynomials
    int64_t Q;  // Q <= 2^K - modulus
    int32_t K;   // K = ceil(logQ) - attribute length and boolean circuits arity
    int32_t L;   // L = (N * K) - KP-ABE matrices dimension
    int32_t P;   // CP trap size
    int32_t M;   // M = (P + 2) - artificial CP trap size for computation
    int32_t Att_num; // number of attributes
    real SIGMA;  // used for discrete gaussian distribution
    real SHORT_THRESHOLD;  // threshold used in is_short
    int32_t MBAR;  // m_bar = N/8
} cp_params;

// Global params globally available
extern cp_params PARAMS;

void init_params(int32_t N, int64_t Q, int32_t K, int32_t P, real SIGMA, int32_t att_num);

// Init params with default values, check common.c for values
void init_params_default(void);

// Nicely print current parameters
void print_params(void);

// void* safe_malloc(size_t size) {
//     void* ptr = malloc(size);
//     if(ptr == NULL) {
//         fprintf(stderr, "Error: Memory allocation failed for size %zu\n", size);
//         exit(EXIT_FAILURE);
//     }
//     return ptr;
// }
/*
double start and double end need to be defined before !
comment need to include %f specifier to include duration
*/
// #define CHRONO(comment, code)                        \
//     do {                                             \
//         start = (real)clock() / CLOCKS_PER_SEC;      \
//         {code} end = (real)clock() / CLOCKS_PER_SEC; \
//         printf(comment, end - start);                \
//     } while (0)
#define CHRONO(comment, code)                                 \
    do {                                                    \
        struct timespec start, end;                         \
        clock_gettime(CLOCK_REALTIME, &start);              \
        { code }                                           \
        clock_gettime(CLOCK_REALTIME, &end);                \
        double elapsed = (end.tv_sec - start.tv_sec) +      \
            (end.tv_nsec - start.tv_nsec) / 1e9;            \
        printf(comment, elapsed);                           \
    } while (0)
// Debug parameters
// #define DEBUG_NORM
