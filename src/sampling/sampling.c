#include "sampling.h"
#include "common.h"
#include "random.h"
#include <assert.h>

void init_sampler() { random_bytes_init(); }

void TrapGen(matrix* B, signed_matrix T) {
    // Pre trap computation
    signed_matrix rho = new_signed_matrix(PARAMS.P, 1);
    signed_matrix mu = new_signed_matrix(PARAMS.P, 1);
    sample_Z_centered_matrix(rho);
    sample_Z_centered_matrix(mu);

    // Trap computation : T = [rho | -g + mu | Ip] : size P * M
    for (int i = 0; i < PARAMS.P; i++)
        matrix_element(T, i, 0) = matrix_element(rho, i, 0);

    signed_scalar g = 1;
    for (int i = 0; i < PARAMS.P; i++) {
        matrix_element(T, i, 1) = -g + matrix_element(mu, i, 0);
        g *= 2;
    }

    for (int i = 0; i < PARAMS.P; i++) matrix_element(T, i, i + 2) = 1;

    // B matrixes computation
    // For each matrix Bk
    for (int k = 0; k < 2 * PARAMS.K + 1; k++) {
        // We treat each column independantly
        // For each column c of matrix Bk
        for (int c = 0; c < PARAMS.N; c++) {
            scalar a = uniform_mod_n(PARAMS.Q);

            // Two first terms by hand
            matrix_element(B[k], 0, c) = a;
            matrix_element(B[k], 1, c) = 1;

            // Then use the formula for the PARAMS.P last terms
            scalar g = 1;
            for (int i = 0; i < PARAMS.P; i++) {
                signed_scalar t = g;
                t -= a * matrix_element(rho, i, 0);
                t -= matrix_element(mu, i, 0);
                t %= PARAMS.Q;
                if (t < 0) t += PARAMS.Q;
                matrix_element(B[k], i + 2, c) = (scalar)t;
                g *= 2;
            }
        }
    }

    free_signed_matrix(rho);
    free_signed_matrix(mu);
}

signed_matrix TrapSamp(matrix* B, signed_matrix T, attribute x) {
    signed_matrix Tx = new_signed_matrix(PARAMS.P, PARAMS.M);
    // We simply copy T for now...
    // We give full knowledge instead of specific knowledge related to x
    for (int i = 0; i < PARAMS.P; i++) {
        for (int j = 0; j < PARAMS.M; j++) {
            matrix_element(Tx, i, j) = matrix_element(T, i, j);
        }
    }
    return Tx;
}

/* -------------------- */
/* Functions for matrix */
/* -------------------- */

void sample_Zq_uniform_matrix(matrix A) {
    for (int i = 0; i < A.rows; i++)
        for (int j = 0; j < A.columns; j++)
            matrix_element(A, i, j) = uniform_mod_n(PARAMS.Q);
}


void sample_Zq_uniform_matrix_64(matrix A) {
    for (int i = 0; i < A.rows; i++)
        for (int j = 0; j < A.columns; j++)
            matrix_element(A, i, j) = uniform_mod_n_64(PARAMS.Q);
}

signed_scalar sample_Z_centered() { return algorithmF(0, PARAMS.SIGMA); }

void sample_Z_centered_matrix(signed_matrix A) {
    for (int i = 0; i < A.rows; i++)
        for (int j = 0; j < A.columns; j++)
            matrix_element(A, i, j) = sample_Z_centered();
}

void single_TrapGen(matrix B, matrix R) 
{
    // 获取维度参数

    unsigned int n = PARAMS.N;
    unsigned int m_bar = PARAMS.MBAR;  
    unsigned int w = PARAMS.L;

    // print_matrix(B);
    // 生成随机矩阵R
    // matrix R = new_matrix(m_bar, w);
    sample_Zq_uniform_matrix(R);

    // 生成随机矩阵A_bar
    matrix A_bar = new_matrix(n, m_bar);
    sample_Zq_uniform_matrix(A_bar);

    // 生成工具矩阵G - 修正维度参数
    // matrix G = new_matrix(n, w);
    // for (int i = 0; i < n; i++)
    //     for (int k = 0; k < w; k++)
    //         matrix_element(G, i, k) = (1 << k) % PARAMS.Q;
    

    matrix G = new_matrix(n, w);
    for (int i = 0; i < n; i++)
        for (int k = 0; k < PARAMS.K; k++)
            matrix_element(G, i, i * PARAMS.K + k) = (1 << k) % PARAMS.Q;
    // printf("inside G\n");
    // print_matrix(G);
    // 计算a = G - A_bar * R
    matrix temp = new_matrix(n, w);
    matrix a = new_matrix(n, w);
    mul_matrix(A_bar, R, temp);
    sub_matrix(G, temp, a);
    
    // 构造公钥矩阵B = [A_bar | a]
    for(unsigned int i = 0; i < n; i++) {
        for(unsigned int j = 0; j < m_bar; j++) {
            matrix_element(B, i, j) = matrix_element(A_bar, i, j);
        }
        for(unsigned int j = 0; j < w; j++) {
            matrix_element(B, i, m_bar + j) = matrix_element(a, i, j);
        }
    }
    

    
    // 释放临时矩阵
    free_matrix(A_bar);
    free_matrix(G);
    free_matrix(temp);
    free_matrix(a);
}


// 二进制表示转换函数
void bin_rp(scalar q, int* result, int k) {
    for(int i = 0; i < k; i++) {
        result[i] = q & 1;
        q >>= 1;
    }
}

// 离散高斯分布采样 O_1
int sample_O1(double s, double c) {
    // 使用一个更好的随机源生成两个随机数
    scalar r1 = uniform_mod_n(UINT32_MAX);
    scalar r2 = uniform_mod_n(UINT32_MAX);
    double u1 = (double)r1 / UINT32_MAX;
    double u2 = (double)r2 / UINT32_MAX;
    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return (int)round(z * (s/2) + c);
}

// 离散高斯分布采样 O_2
int sample_O2(double s, int v) {
    double u1 = (double)rand() / RAND_MAX;
    double u2 = (double)rand() / RAND_MAX;
    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
    return 2 * (int)round(z * s) + v;
}


// SampleG 主函数
void SampleG(scalar q, scalar u, double s, int k, int* result) {
    int* u_v = (int*)malloc(k * sizeof(int));
    int* q_v = (int*)malloc(k * sizeof(int));
    int* v_v = (int*)malloc(k * sizeof(int));
    
    // 转换为二进制表示
    bin_rp(u, u_v, k);
    // Debug: print u_v
    // printf("u_v: ");
    // for(int i = 0; i < k; i++) {
    //     printf("%d ", u_v[i]);
    // }
    // printf("\n");


    bin_rp(q, q_v, k);
    // Debug: print q_v
    // printf("q_v: ");
    // for(int i = 0; i < k; i++) {
    //     printf("%d ", q_v[i]);
    // }
    // printf("\n");
    // 计算 c 和 y
    double c = -((double)u/q);
    // printf("c: %f\n", c);
    int y = sample_O1(s, c);
    
    // printf("y: %d\n", y);
    // y=4;
    // 计算 v_v
    for(int i = 0; i < k; i++) {
        v_v[i] = u_v[i] + y * q_v[i];
    }
    
    // printf("v_v: ");
    // for(int i = 0; i < k; i++) {
    //     printf("%d ", v_v[i]);
    // }
    // printf("\n");
    // // 采样过程

    // int x[] = {30, -20, -35, 93, -43, -30, -2, 58, 29, -38, -73, 65, -43, 48, -64, -41, 0, 73, -11,};
    for(int i = 0; i < k-1; i++) {
        result[i] = sample_O2(s, v_v[i]);
        // result[i]=x[i];
        v_v[i+1] = v_v[i+1] + (v_v[i] - result[i])/2;
    }

    // printf("v_v: ");
    // for(int i = 0; i < k; i++) {
    //     printf("%d ", v_v[i]);
    // }
    // printf("\n");
    result[k-1] = v_v[k-1];
    

    free(u_v);
    free(q_v);
    free(v_v);
}

void SampleD(matrix A, matrix u, matrix R, double sigma, signed_matrix result) {
    // 创建临时向量和矩阵

    // printf("dimension of A: %d %d\n", A.rows, A.columns);
    // printf("dimension of u: %d %d\n", u.rows, u.columns);
    // printf("dimension of R: %d %d\n", R.rows, R.columns);
    // printf("dimension of result: %d %d\n", result.rows, result.columns);

    // 4. 创建结果向量x和临时向量p
    signed_matrix p = new_signed_matrix(PARAMS.M+PARAMS.MBAR, 1);
    sample_Z_centered_matrix(p);

    // // 6. 计算v = u - A * p
    matrix Bp = new_matrix(PARAMS.N, 1);
    mul_matrix_trap(A, p, Bp);
    
    matrix v = new_matrix(PARAMS.N, 1);
    sub_matrix(u, Bp, v);
    // // 7. 对每个v[i]进行SampleG采样

    signed_matrix x_temp = new_signed_matrix(PARAMS.K, PARAMS.N);
    int* g_samples = (int*)malloc((PARAMS.K) * sizeof(int));
    for(int i = 0; i < PARAMS.N; i++) {
        SampleG(PARAMS.Q, matrix_element(v,i,0), sqrt(2), PARAMS.K, g_samples);
        // printf(" \ngsamples:");
        // for(int j = 0; j < PARAMS.K; j++) {
        //     printf("%d ", g_samples[j]);
        // }
        for(int j = 0; j < PARAMS.K; j++) {
            matrix_element(x_temp,  j, i) = g_samples[j];
        }
    }

    // printf("v: \n");
    // print_matrix(v);
    // printf("dimension of x_temp: %d %d\n", x_temp.rows, x_temp.columns);
    // print_signed_matrix(x_temp);

     // 8. 构造[R|I]矩阵
    matrix RI = new_matrix(PARAMS.M+PARAMS.MBAR, PARAMS.M);
    for(int i = 0; i < PARAMS.MBAR; i++) {
        for(int j = 0; j < PARAMS.M; j++) {
            matrix_element(RI, i, j) = matrix_element(R, i, j);
        }
    }
    for(int i = 0; i < PARAMS.M; i++) {
        matrix_element(RI, i + PARAMS.MBAR, i) = 1;
    }

    // signed_matrix x_final = new_signed_matrix(PARAMS.M+PARAMS.MBAR, 1);
    mul_matrix_signed(RI, x_temp, result);
    // Add p and x_final
    for(int i = 0; i < PARAMS.MBAR+PARAMS.M; i++) {
        matrix_element(result, i, 0) += matrix_element(p, i, 0);
    }
}


void SamplePre(matrix A, matrix R, matrix u, double sigma, signed_matrix result){  
    SampleD(A, u, R, sigma, result);
}


void SampleLeft(matrix A, matrix R, matrix B, matrix u, double sigma, matrix result){  

    signed_matrix k1 = new_signed_matrix(A.columns, 1);
    matrix k2 = new_matrix(B.columns, 1);
    matrix target_u = new_matrix(u.rows, u.columns);
    sample_Zq_uniform_matrix(k2);
    matrix temp = new_matrix(u.rows, u.columns);
    // printf("dimension of k2 %d x %d\n", k2.rows, k2.columns);
    mul_matrix(B, k2, temp);
    sub_matrix(u, temp, target_u);
    SampleD(A, target_u, R, PARAMS.SIGMA, k1);

    // matrix k_vector = new_matrix(k1.rows + k2.rows, 1);

    assert(result.rows == k1.rows + k2.rows);
    for (unsigned int i = 0; i < k1.rows; i++) {
        matrix_element(result, i, 0) = matrix_element(k1, i, 0);
    }
    for (unsigned int i = 0; i < k2.rows; i++) {
        matrix_element(result, k1.rows + i, 0) = matrix_element(k2, i, 0);
    }
}
