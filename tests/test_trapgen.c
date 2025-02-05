#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "matrix.h"
#include "sampling.h"
#include "common.h"

void test_single_TrapGen() {
    // init_params_default();
    printf("Starting test_single_TrapGen...\n");
    
    // 1. 初始化维度参数
    const unsigned int n = PARAMS.N;
    const unsigned int m_bar = PARAMS.MBAR;
    const unsigned int w = PARAMS.L;
    const unsigned int m = m_bar + w;
    
    printf("Dimensions: n=%u, m_bar=%u, w=%u, m=%u\n", n, m_bar, w, m);
    
    // 2. 分配矩阵内存
    printf("m,n: %d,%d\n", m, n);
    matrix B = new_matrix(n, m);
    matrix R = new_matrix(m_bar, w);
    
    printf("Matrices allocated\n");
    
    // 3. 执行TrapGen
    printf("Executing single_TrapGen...\n");
    single_TrapGen(B, R);

    printf("TrapGen completed\n");
    printf("dimension of B: %d %d\n", B.rows, B.columns);
    printf("dimension of R: %d %d\n", R.rows, R.columns);

    matrix G = new_matrix(PARAMS.N, PARAMS.L);
    for (int i = 0; i < PARAMS.N; i++)
        for (int k = 0; k < PARAMS.K; k++)
            matrix_element(G, i, i * PARAMS.K + k) = (1 << k) % PARAMS.Q;

    // 构造[R|I]矩阵
    matrix RI = new_matrix(m, w);
    // 复制R部分
    for(unsigned int i = 0; i < m_bar; i++) {
        for(unsigned int j = 0; j < w; j++) {
            matrix_element(RI, i, j) = matrix_element(R, i, j);
        }
    }
    // 添加单位矩阵I部分
    for(unsigned int i = 0; i < w; i++) {
        for(unsigned int j = 0; j < w; j++) {
            matrix_element(RI, i + m_bar, j) = (i == j) ? 1 : 0;
        }
    }
    printf("dimension of RI: %d %d\n", RI.rows, RI.columns);
    // 计算G_p = B * [R|I]
    // printf("Computing G_p = B * [R|I]...\n");
    matrix G_p = new_matrix(n, w);
    mul_matrix(B, RI, G_p);

    // 验证G_p == G
    // printf("Verifying G_p == G...\n");
    bool valid = true;
    for(unsigned int i = 0; i < n && valid; i++) {
        for(unsigned int j = 0; j < w && valid; j++) {
            if(matrix_element(G_p, i, j) != matrix_element(G, i, j)) {
                valid = false;
                printf("Mismatch at (%u,%u): G_p=%lu, G=%lu\n", 
                    i, j,
                    matrix_element(G_p, i, j),
                    matrix_element(G, i, j));
            }
        }
    }

    if(valid) {
        printf("✓ G_p equals G!\n");
    } else {
        printf("✗ G_p does not equal G!\n");
    }

    // Print G matrix
    // printf("G matrix:\n");
    // print_matrix(G);

    // printf("Gp matrix:\n");
    // print_matrix(G_p);

    // 5. 清理内存
    printf("Cleaning up...\n");
    // free_matrix(B);
    // free_matrix(R);
    // free_matrix(G);
    // free_matrix(G_p);
    // free_matrix(RI);

    printf("Test completed\n");
}

int main() {
    printf("Initializing...\n");
    init_params_default();
    init_sampler();
    
    test_single_TrapGen();
    
    return 0;
}