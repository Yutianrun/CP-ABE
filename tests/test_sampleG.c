#include <stdio.h>
#include <stdlib.h>
#include "sampling.h"
#include "common.h"
#include "matrix.h"

void test_SampleG() {

    printf("Testing SampleG function...\n");
    
    // 测试参数
//     n = 64
    // int64_t q = 870367;  // 从PARAMS获取q值
    // int64_t u = 1213;      // 测试输入值
    // double s = 20;        // sigma参数
    // int k = 20;     // 二进制长度
    
    int64_t q = PARAMS.Q;  // 从PARAMS获取q值
    double s = PARAMS.SIGMA;        // sigma参数
    int k = PARAMS.K;     // 二进制长度

    // 创建g向量 (1,2,4,8,...)
    int64_t* g_v = (int64_t*)malloc(k * sizeof(int64_t));
    int64_t power = 1;
    for(int i = 0; i < k; i++) {
        g_v[i] = power;
        power = (power * 2) % q;
    }
    
    // 运行SampleG
        // 运行SampleG
    for (int test = 0; test < 1000; test++) {
        // 确保q不为0
        if (q == 0) {
            printf("Error: q cannot be zero\n");
            break;
        }

        // Generate random u for each test (确保在有效范围内)
        int64_t random_u = (rand() % (q-1)) + 1;
        
        int* result = (int*)malloc(k * sizeof(int));
        if (!result) {
            printf("Memory allocation failed\n");
            break;
        }

        // 添加调试信息
        printf("Debug - Test %d: q=%ld, u=%ld, s=%f\n", test+1, q, random_u, s);
        
        SampleG(q, random_u, s, k, result);
        
        // 计算g_v * x mod q (使用临时变量避免溢出)
        int64_t dot_product = 0;
        for(int i = 0; i < k; i++) {
            int64_t product = ((int64_t)g_v[i] * result[i]);
            int64_t temp = product % q;
            if(temp < 0) temp += q;
            dot_product = (dot_product + temp) % q;
            if(dot_product < 0) dot_product += q;
        }

        // free(result);
        // ...existing code...
        
        // 验证结果
        printf("Test %d - Input u: %d, ", test + 1, random_u);
        printf("g_v * x mod q = %d, ", dot_product);
        printf("Result: %s\n", (dot_product == random_u) ? "PASSED" : "FAILED");
        
        if (dot_product != random_u) {
            printf("g_v vector: [");
            for(int i = 0; i < k; i++) {
                printf("%ld%s", g_v[i], (i < k-1) ? ", " : "");
            }
            printf("]\n");
            printf("x vector: [");
            for(int i = 0; i < k; i++) {
                printf("%d%s", result[i], (i < k-1) ? ", " : "");
            }
            printf("]\n");
            
            return;
        }

        // Optionally print sampled vector (commented out to reduce output)
        /*printf("Sampled vector x: [");
        for(int i = 0; i < k; i++) {
            printf("%d%s", result[i], (i < k-1) ? ", " : "");
        }
        printf("]\n");*/
        
        free(result);
    }
    
    free(g_v);
}

void test_SampleD() {
    printf("Starting test_SampleD...\n");
    
    // 1. 初始化参数
    const int n = PARAMS.N;
    const int m_bar = PARAMS.MBAR;
    const int w = PARAMS.L;
    const int m = m_bar + w;
    
    printf("Dimensions: n=%u, m_bar=%u, w=%u, m=%u\n", n, m_bar, w, m);
    // 2. 创建并生成B和R

    matrix B = new_matrix(n, m);
    matrix R = new_matrix(m_bar, w);
    single_TrapGen(B, R);
    
    printf("Matrix B:\n");
    print_matrix(B);
    printf("Matrix R:\n");
    print_matrix(R);
    // printf("TrapGen completed\n");
    // 3. 创建随机目标向量u
    matrix u = new_matrix(n, 1);
    sample_Zq_uniform_matrix(u);
    

    // 4. 创建结果向量x和临时向量p
    // signed_matrix x = new_signed_matrix(w, 1);
    // matrix x = new_matrix(m, 1);
    signed_matrix x = new_signed_matrix(m, 1);
    
    SampleD(B, u, R, PARAMS.SIGMA, x);
 
 
    matrix test = new_matrix(n, 1);
    // mul_matrix_trap(B, x_final, test);
    mul_matrix_trap(B, x, test);
    
    

    bool valid = equals(test, u);

    // printf(" vector x:\n");
    // print_signed_matrix(x);
    printf("Target vector u:\n");
    print_matrix(u);
    printf("Result vector B*x:\n");
    print_matrix(test);
    printf(valid ? "✓ SampleD test passed: B * x = u\n" 
                 : "✗ SampleD test failed: B * x != u\n");
    
    // // 10. 释放内存
    // free_matrix(B);
    // free_matrix(R);
    // free_matrix(u);
    // free_signed_matrix(x);
    // free_signed_matrix(x_final);
    // free_signed_matrix(p);
    // free_matrix(Bp);
    // free_matrix(v);
    // free_matrix(test);

    // free(g_samples);
    
    // printf("Test completed\n");
}

int main() {
    int64_t q = 17;  // 从PARAMS获取q值
    double s = 2;        // sigma参数
    int k = 5;     // 二进制长度
    int n = 1;
    int p= 1;

    // init_params_default();
    init_params(n, q, k, p , s, 1);
    PARAMS.MBAR=1;
    init_sampler();
    // test_SampleG();
    test_SampleD();
    return 0;
}