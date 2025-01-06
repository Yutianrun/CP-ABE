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

// 构建复杂 PRF 电路，包含多个字句
circuit* build_eval_prf(int k, Clause* clauses, int num_clauses) {
    // 假设前 k 位是 x，接下来的 k 位是 msk
    // 每个字句包含一个 T (长度 t_size) 和相应的 v

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

        // 构建 (T, v) 的输入电路
        // 这里假设将 T 和 v 合并为一个电路，例如通过 OR 或其他逻辑门
        // 具体实现取决于如何表示 (T, v)
        // 这里简单地通过 XOR 组合 T 和 v
        circuit* T_v = NULL;
        for(int j = 0; j < current_clause.t_size; j++) {
            // 生成 T[j] 的电路
            // 这里 T[j] 是一个索引，作为电路的输入
            // 假设 T 本身也是一些特定位的输入
            // 具体取决于 (T, v) 的表示方式
            // 这里假设 T 只是 v 的提取位置，不需要额外电路
            // 所以直接使用 v 作 PRF 的输入
            if (T_v == NULL) {
                T_v = circuit_not(v); // 举例操作
            } else {
                T_v = circuit_xor(T_v, v); // 举例操作
            }
        }

        // 生成 PRF(msk, (T,v))，即 PRF 的输入是 T_v，密钥是 msk
        // 假设 PRF(msk, input) 已通过 build_prf_circuit 实现
        // 这里需要将 msk 和 T_v 作为 PRF 的输入
        // 假设 PRF 以某种方式将 msk 与输入组合
        // 具体实现需要根据 PRF 的定义调整

        // 这里简化为直接调用 PRF
        // 需要一个函数将输入电路和密钥电路传递给 PRF
        // 参考 build_prf_circuit 函数

        // 假设 build_prf_circuit 接受 x 和 msk
        // 这里需要构建 PRF(meshk, T_v)
        // 需要将 T_v 作为 x 的一部分

        // 为简化示例，假设每个 PRF 只以 T_v 为输入，并共享 msk
        // 这里需要更多的逻辑门组合来将 msk 与 (T,v) 结合
        // 具体实现可能需要扩展 PRF 电路的输入

        // 这里假设存在一个函数 PRF(msk, T_v)，返回 sk_(T,v)
        // 具体细节需要根据你的 PRF 实现调整
        // 因此，此处仅提供示例代码

        // 示例：sk_(T,v) = PRF(msk, T_v)
        // 假设 PRF 函数已能处理这样的输入
        // 这里使用 build_prf_circuit 生成 sk_(T,v)
        // 将 (T,v) 作为 PRF 的输入，需要将其与 msk 结合

        // 构建 PRF 输入电路：将 msk 和 T_v 组合
        // 这里假设组合方式为 OR, XOR 等，具体取决于 PRF 的要求
        // 示例中将 msk 和 T_v 通过 XOR 组合
        circuit* prf_input = circuit_xor(v, x[0]); // 仅示例，实际需要组合所有相关位

        // 生成 PRF 电路
        circuit** prf = build_prf_circuit(k);
        // 在实际实现中，应该将 prf_input 作为 PRF 的输入
        // 这里简化为直接使用 prf[0] 作为 sk_(T,v)
        sk[i] = prf[0]; // 示例，实际需要正确映射
    }

    // 第二步：计算 PRF(sk_(T,v), x) for each clause
    circuit** prf_outputs = (circuit**)malloc(num_clauses * sizeof(circuit*));
    if (prf_outputs == NULL) {
        fprintf(stderr, "内存分配失败\n");
        exit(1);
    }

    for(int i = 0; i < num_clauses; i++) {
        // 计算 PRF(sk_(T,v), x)
        // 类似于第一步，需要将 sk[i] 和 x 作为输入
        // 这里简化为直接调用 PRF 电路
        // 实际实现需要根据 PRF 的定义调整

        // 这里假设 PRF 接受 sk[i] 作为密钥，并以 x 作为输入
        // 因此需要将 sk[i] 作为 PRF 的密钥，x 作为输入
        // 示例中使用 build_prf_circuit 生成 PRF 输出
        circuit** prf = build_prf_circuit(k);
        prf_outputs[i] = prf[0]; // 示例，实际需要正确映射
    }

    // 将所有 PRF 输出异或
    circuit* final_y = NULL;
    for(int i = 0; i < num_clauses; i++) {
        if (final_y == NULL) {
            final_y = prf_outputs[i];
        } else {
            final_y = circuit_xor(final_y, prf_outputs[i]);
        }
    }

    // 清理中间电路
    for(int i = 0; i < k; i++) {
        // 假设不需要释放 x 和 msk，因为它们可能被多次引用
    }
    for(int i = 0; i < num_clauses; i++) {
        free_circuit(sk[i]);
        free_circuit(prf_outputs[i]);
    }
    free(x);
    free(msk);
    free(sk);
    free(prf_outputs);

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
