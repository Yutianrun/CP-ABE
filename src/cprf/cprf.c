// prf_circuit.c
#define PRF_K 128
#include "cprf.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>


// 构建单个 PRF 电路的函数（已修正 n 的索引）
circuit** build_prf_circuit(int k) {
    // 创建一个数组存储输出电路
    circuit** y = (circuit**)malloc(k * sizeof(circuit*));
    if (y == NULL) {
        fprintf(stderr, "内存分配失败\n");
        exit(1);
    }

    for(int i = 0; i < k; i++) {
        // 使用1基索引，确保不会传递 n <= 0
        circuit* x_i = gen_leaf(i + 1, false); // 输入 x[i] -> n = i + 1
        circuit* x_ip1 = gen_leaf((i + 1) % k + 1, false); // 输入 x[i+1] -> n = (i + 1) % k + 1
        circuit* msk_i = gen_leaf(i + 1, true); // 密钥 msk[i] -> n = i + 1
        circuit* msk_ip1 = gen_leaf((i + 1) % k + 1, true); // 密钥 msk[i+1] -> n = (i + 1) % k + 1
        circuit* msk_ip2 = gen_leaf((i + 2) % k + 1, true); // 密钥 msk[i+2] -> n = (i + 2) % k + 1

        // 第一层逻辑门组合
        circuit* a = circuit_and(x_i, msk_i); // a = x[i] AND msk[i]
        circuit* b = circuit_or(x_ip1, msk_ip1); // b = x[i+1] OR msk[i+1]
        circuit* c = circuit_not(msk_ip2); // c = NOT msk[i+2]

        // 第二层逻辑门组合
        circuit* d = circuit_xor(a, b); // d = a XOR b
        circuit* e = circuit_and(d, c); // e = d AND c
        circuit* f = circuit_or(e, a); // f = e OR a

        // 第三层逻辑门组合
        circuit* g = circuit_not(f); // g = NOT f
        circuit* h = circuit_xor(g, d); // h = g XOR d

        // 将结果作为输出位 y[i]
        y[i] = h;
    }

    return y;
}

// 辅助函数：提取 v 比特
circuit* extract_v(int k, int* T, int t_size, circuit** x) {
    // 假设 v 的长度为 t_size，依次提取对应的 x 位
    // 比如 T = (1, 2)，则 v = x1, x2
    // 将 v 通过 XOR 组合起来作为唯一输入
    circuit* v = NULL;
    for(int i = 0; i < t_size; i++) {
        circuit* bit = x[T[i] - 1]; // T 是1基
        if (v == NULL) {
            v = bit;
        } else {
            v = circuit_or(v, bit); // 你可以根据需要选择合适的逻辑门组合
        }
    }
    return v;
}

circuit* build_eval_prf(int k, Clause* clauses, int num_clauses) {
    // 生成 x 和 msk 的叶节点
    circuit** x = (circuit**)malloc(k * sizeof(circuit*));
    circuit** msk = (circuit**)malloc(k * sizeof(circuit*));
    if (x == NULL || msk == NULL) {
        fprintf(stderr, "内存分配失败\n");
        exit(1);
    }
    for(int i = 0; i < k; i++) {
        x[i] = gen_leaf(i + 1, false); // x1, x2, ..., xk
        msk[i] = gen_leaf(k + i + 1, true); // msk1, msk2, ..., mskk
    }

    // 存储所有 sk_(T,v)
    circuit** sk = (circuit**)malloc(num_clauses * sizeof(circuit*));
    if (sk == NULL) {
        fprintf(stderr, "内存分配失败\n");
        exit(1);
    }

    // 第一步：计算 sk_(T,v) = PRF(msk, (T,v)) for each clause
    for(int i = 0; i < num_clauses; i++) {
        Clause current_clause = clauses[i];
        // 提取 v 比特
        circuit* v = extract_v(k, current_clause.T, current_clause.t_size, x);
        if (v == NULL) {
            fprintf(stderr, "提取 v 失败\n");
            exit(1);
        }

        // 构建 (T, v) 的输入电路
        circuit* T_v = NULL;
        for(int j = 0; j < current_clause.t_size; j++) {
            if (T_v == NULL) {
                T_v = circuit_not(v);
            } else {
                T_v = circuit_xor(T_v, v);
            }
            if (T_v == NULL) {
                fprintf(stderr, "构建 T_v 失败\n");
                exit(1);
            }
        }

        // 构建 PRF 输入电路：将 msk 和 T_v 组合
        circuit* prf_input = circuit_xor(v, x[0]); // 仅示例，实际需要组合所有相关位
        if (prf_input == NULL) {
            fprintf(stderr, "构建 PRF 输入失败\n");
            exit(1);
        }

        // 生成 PRF 电路
        circuit** prf = build_prf_circuit(k);
        if (prf == NULL || prf[0] == NULL) {
            fprintf(stderr, "生成 PRF 电路失败\n");
            exit(1);
        }

        sk[i] = prf[0];
    }

    // 第二步：计算 PRF(sk_(T,v), x) for each clause
    circuit** prf_outputs = (circuit**)malloc(num_clauses * sizeof(circuit*));
    if (prf_outputs == NULL) {
        fprintf(stderr, "内存分配失败\n");
        exit(1);
    }

    for(int i = 0; i < num_clauses; i++) {
        circuit** prf = build_prf_circuit(k);
        if (prf == NULL || prf[0] == NULL) {
            fprintf(stderr, "生成 PRF 输出失败\n");
            exit(1);
        }
        prf_outputs[i] = prf[0];
    }

    // 将所有 PRF 输出异或
    circuit* final_y = NULL;
    for(int i = 0; i < num_clauses; i++) {
        if (final_y == NULL) {
            final_y = prf_outputs[i];
        } else {
            final_y = circuit_xor(final_y, prf_outputs[i]);
            if (final_y == NULL) {
                fprintf(stderr, "异或 PRF 输出失败\n");
                exit(1);
            }
        }
    }

    // 清理中间电路
    // for(int i = 0; i < num_clauses; i++) {
    //     free_circuit(sk[i]);
    //     // 移除对 prf_outputs[i] 的释放
    //     // free_circuit(prf_outputs[i]);
    // }
    free(x);
    free(msk);
    free(sk);
    free(prf_outputs); // 仅释放数组指针，不释放其中的电路节点

    return final_y;
}


// 新增函数：compute_policy_prf
#include "cprf.h"
#include <string.h>

// 定义政策函数f，基于t-CNF
int policy(int* T, int t_size) {
    // 示例：简单的t-CNF判断，实际实现根据具体政策
    // 返回1表示满足，0表示不满足
    // 此处假设所有bit为1时满足
    for(int i = 0; i < t_size; i++) {
        if(T[i] != 1) {
            return 0;
        }
    }
    return 1;
}

circuit* compute_policy_prf(int k, int* T, int t_size, Clause* clauses, int num_clauses) {
    // 第一步：计算所有满足policy的(T, v)对
    // 这里假设v为T的复制
    if(!policy(T, t_size)) {
        return NULL;
    }

    // 第二步：计算PRF值
    circuit** prf_values = (circuit**)malloc(num_clauses * sizeof(circuit*));
    if(prf_values == NULL) {
        fprintf(stderr, "内存分配失败\n");
        exit(1);
    }

    for(int i = 0; i < num_clauses; i++) {
        // 使用现有的build_prf_circuit函数
        circuit** prf = build_prf_circuit(k);
        prf_values[i] = prf[0];
    }

    // 第三步：异或所有PRF值
    circuit* final_prf = NULL;
    for(int i = 0; i < num_clauses; i++) {
        if(final_prf == NULL) {
            final_prf = prf_values[i];
        } else {
            final_prf = circuit_xor(final_prf, prf_values[i]);
        }
    }

    // 清理
    free(prf_values);

    return final_prf;
}
