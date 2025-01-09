// prf_circuit.c
#include "cprf.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#define ROUND 2

// 构建单个 PRF 电路的函数（已修正 n 的索引）
// 修改后的构建单个 PRP 电路的函数，使用Feistel结构
circuit** build_prp_circuit(int k) {
    int rounds = ROUND;
    // 创建一个数组存储输出电路
    circuit** y = (circuit**)malloc(k * sizeof(circuit*));
    if (y == NULL) {
        fprintf(stderr, "内存分配失败\n");
        exit(1);
    }

    for(int i = 0; i < k/2; i++) {
        // 分割输入为左右两部分
        circuit* left = gen_leaf(i + 1, false); // 左半部分 x[i]
        circuit* right = gen_leaf((i + 1) +k/2, false); // 右半部分 x[i+1]
        circuit* msk = gen_leaf(i + 1 + k, false); // 密钥 msk[i]
        circuit* msk2 = gen_leaf(i + 1 + k + k/2, false); // 密钥 msk[i]

        // circuit* msk3 = gen_leaf(i + 1 + k, false); // 密钥 msk[i]

        // circuit* msk4 = gen_leaf(i + 1 + k, false); // 密钥 msk[i]
        // circuit* round_keys[rounds];
        // for (int r = 0; r < rounds; r++) {
        //     round_keys[r] = gen_leaf(i + 1 + k + r * k / rounds, false);
        // }

        circuit* L = left;
        circuit* R = right;

        for(int r = 0; r < rounds; r++) {
            circuit* F = new_circuit();
            circuit* new_L = R;
            // 轮函数
            if (r % 2 == 0) {
            F = circuit_xor(R, msk); // 偶数轮使用 msk
            F->left_double_pointed = true;
            } else {
            F = circuit_xor(R, msk2); // 奇数轮使用 msk2
            F->left_double_pointed = true;
            }
              
            circuit* new_R = circuit_xor(L, F); // 新右 = L XOR 


            L = new_L;
            R = new_R;

        }

        // 合并左右部分作为输出
        y[i] = L; // 输出左部分
        y[(i + k/2) % k] = R; // 输出右部分
    }

    return y;
}


circuit** initial_prp_circuit(int k, circuit** x, circuit** msk) {
    int rounds = ROUND;
    // 创建一个数组存储输出电路
    circuit** y = (circuit**)malloc(k * sizeof(circuit*));
    if (y == NULL) {
        fprintf(stderr, "内存分配失败\n");
        exit(1);
    }

    for(int i = 0; i < k/2; i++) {
        // 使用传入的 x 和 msk 数组
        circuit* left = x[i]; // 左半部分 x[i]
        circuit* right = x[i + k/2]; // 右半部分 x[i + k/2]
        circuit* msk1 = msk[i]; // 密钥 msk[i]
        circuit* msk2 = msk[i + k/2]; // 密钥 msk[i + k/2]

        circuit* L = left;
        circuit* R = right;

        for(int r = 0; r < rounds; r++) {
            // 轮函数
            circuit* F;
            if (r % 2 == 0) {
                F = circuit_xor(R, msk1); // 偶数轮使用 msk1
                F->left_double_pointed = true;
            } else {
                F = circuit_xor(R, msk2); // 奇数轮使用 msk2
                F->left_double_pointed = true;  
            }
            // Feistel交换
            circuit* new_L = R;
            circuit* new_R = circuit_xor(L, F); // 新右 = L XOR F
            L = new_L;
            R = new_R;
        }

        // 合并左右部分作为输出
        y[i] = L; // 输出左部分
        y[(i + k/2) % k] = R; // 输出右部分
    }

    return y;
}


// typedef struct {
//     int* T;       // T 的位置数组
//     int t_len;    // T 的长度
//     int* v;       // 对应位置上的值数组
//     int v_len;    // v 的长度（应与 t_len 相同）
// } ClauseT;

circuit* circuit_set_bit(bool bit_value, circuit* leaf_circuit) {
    if (bit_value) {
        // 如果 bit_value 是 1，创建一个恒定返回 1 的电路
        circuit* not_one = circuit_not(leaf_circuit); // 取反，得到 0
        return circuit_or(leaf_circuit, not_one); // 1 OR 0 恒定返回 1
    } else {
        // 如果 bit_value 是 0，创建一个恒定返回 0 的电路
        circuit* not_zero = circuit_not(leaf_circuit); // 取反，得到 1
        return circuit_and(leaf_circuit, not_zero); // 0 AND 1 恒定返回 0
    }
}


circuit** encode_Tx(int k, ClauseT clause, circuit** x) {
    // 初始化一个 k 位全0的电路
    circuit** encoded = (circuit**)malloc(sizeof(circuit*) * k*clause.t_len);
    if (encoded == NULL) {
        fprintf(stderr, "内存分配失败\n");
        exit(1);
    }
    for(int i =0;i<clause.t_len;i++){
        // printf("T[%d] = %d\n",i,clause.T[i]);
        encoded[clause.T[i]] = circuit_set_bit(1, x[0]);
        encoded[i+k/2] = x[clause.T[i]];
    }

    for(int i = 0; i < k; i++) {
        if (encoded[i] == NULL) {
            // printf("encoded[%d] is NULL\n", i);
            encoded[i] = circuit_set_bit(0, x[0]);
        }
    }

    return encoded;
}

circuit** encode_pair(int k, Pair pair, circuit* arbitray) {
    // 初始化一个 k 位全0的电路
    circuit** encoded = (circuit**)malloc(sizeof(circuit*) * k);
    if (encoded == NULL) {
        fprintf(stderr, "内存分配失败\n");
        exit(1);
    }
    for(int i =0;i<pair.t_len;i++){
        // printf("T[%d] = %d\n",i,clause.T[i]);
        encoded[pair.T[i]] = circuit_set_bit(1, arbitray);
        encoded[i+k/2] = circuit_set_bit(pair.v[i], arbitray);
    }

    for(int i = 0; i < k; i++) {
        if (encoded[i] == NULL) {
            // printf("encoded[%d] is NULL\n", i);
            encoded[i] = circuit_set_bit(0, arbitray);
        }
    }
    return encoded;
}


// 分配 sk_tv 的函数
circuit*** allocate_sk_tv(int num_clauses, int k) {
    circuit*** sk_tv = (circuit***)malloc(num_clauses * sizeof(circuit**));
    if (sk_tv == NULL) {
        fprintf(stderr, "内存分配失败用于 sk_tv\n");
        exit(1);
    }

    for(int i = 0; i < num_clauses; i++) {
        sk_tv[i] = (circuit**)malloc(k * sizeof(circuit*));
        if (sk_tv[i] == NULL) {
            fprintf(stderr, "内存分配失败用于 sk_tv[%d]\n", i);
            for(int j = 0; j < i; j++) {
                free(sk_tv[j]);
            }
            free(sk_tv);
            exit(1);
        }
    }

    return sk_tv;
}

// 释放 sk_tv 的函数
void free_sk_tv(circuit*** sk_tv, int num_clauses) {
    for(int i = 0; i < num_clauses; i++) {
        free(sk_tv[i]);
    }
    free(sk_tv);
}

// 构建 sk_tv 的函数
circuit*** build_sk_tv(int k, ClauseT* clauses, int num_clauses, circuit** msk, circuit** x) {
    circuit*** sk_tv = allocate_sk_tv(num_clauses, k);

    for(int i = 0; i < num_clauses; i++) {
        // 编码 (T, v) 为8位
        circuit** encoded_Tv = encode_Tx(k, clauses[i], x);

        // 构建 PRF 的输入：前8位为 encoded_Tv，后8位为 msk
        // 这里 k=8
        sk_tv[i] = initial_prp_circuit(k, encoded_Tv, msk);
        if (sk_tv[i] == NULL) {
            fprintf(stderr, "合并输入失败用于 sk_tv[%d]\n", i);
            free_sk_tv(sk_tv, num_clauses);
            exit(1);
        }
    }
    return sk_tv;
}


circuit** final_prf(int k, circuit*** sk_tv, int num_clauses, circuit** x) {
    circuit** prf_outputs = (circuit**)malloc(k * sizeof(circuit*));
    if (prf_outputs == NULL) {
        fprintf(stderr, "内存分配失败用于 prf_outputs\n");
        exit(1);
    }

    for (int i = 0; i < num_clauses; i++) {
        // 使用 sk_tv[i] 作为密钥，x 作为输入，计算 PRF 电路
        circuit** prf_output = initial_prp_circuit(k, x, sk_tv[i]);
        if (prf_output == NULL) {
            fprintf(stderr, "PRF 计算失败用于 prf_output[%d]\n", i);
            free(prf_outputs);
            exit(1);
        }

        if (i==0){
            for(int j = 0; j < k; j++) {
                prf_outputs[j] = prf_output[j];
            }
        }
        else{
            for(int j = 0; j < k; j++) {
                prf_outputs[j] = circuit_xor(prf_outputs[j], prf_output[j]);
            }
        }
    }

    return prf_outputs;

}

circuit** build_eval_circuit(int k, ClauseT* clauses, int num_clauses, circuit** msk, circuit** x) {
    // 构建 sk_tv
    circuit*** sk_tv = build_sk_tv(k, clauses, num_clauses, msk, x);

    // 计算最终的 PRF 输出
    circuit** prf_outputs = final_prf(k, sk_tv, num_clauses, x);

    return prf_outputs;

    
}


// 计算 S^f_i
Pair* compute_Sf_i(Clause clause, Pair* S, int S_len, int* result_len) {
    Pair* Sf_i = (Pair*)malloc(S_len * sizeof(Pair)); // 最多 S_len 个结果
    *result_len = 0;

    for (int i = 0; i < (1 << clause.t_len); i++) {
        // printf("i = %d\n",i);
        // Generate v array based on the current value of i
        bool v[clause.t_len];
        for (int bit = 0; bit < clause.t_len; bit++) {
            v[bit] = (i & (1 << bit)) ? true : false;
        }

        // for(int i = 0;i<clause.t_len;i++){
        //     printf("v[%d] = %d ",i,v[i]);
        // }
        // printf("\n");
        // Check if f(v) == 1


        if (clause.f(v)) {
            // Create a new Pair with t and v
            // printf("here test_right!\n");
            Pair new_pair;
            // Initialize new_pair.t based on clause.t
            // Assuming Clause has a member 't' which is an array
            new_pair.T = malloc(clause.t_len * sizeof(int));
            new_pair.v = malloc(clause.t_len * sizeof(int));
            new_pair.t_len = clause.t_len;
            // printf("new_pair.t_len = %d\n", new_pair.t_len);
            for (int t_idx = 0; t_idx < clause.t_len; t_idx++) {
                new_pair.T[t_idx] = clause.T[t_idx];
                new_pair.v[t_idx] = v[t_idx] ? 1 : 0;
            }

            // Add the new_pair to Sf_i
            Sf_i[*result_len] = new_pair;
            (*result_len)++;
        }
    }

    for(int i = 0; i < *result_len; i++) {
        printf("Sf_i[%d] : ", i);
        for (int j = 0; j < clause.t_len; j++) {
            printf("T[%d] = %d, v[%d] = %d ", j, Sf_i[i].T[j], j, Sf_i[i].v[j]);

        }
        printf("\n");
    }
    return Sf_i;
}

// 计算 S^f_rest
Pair* compute_Sf_rest(Clause* clauses, int num_clauses, Pair* S, int S_len, int* result_len) {
    Pair* Sf_rest = (Pair*)malloc(S_len * sizeof(Pair)); // 最多 S_len 个结果
    if (Sf_rest == NULL) {
        fprintf(stderr, "Memory allocation failed for Sf_rest\n");
        exit(1);
    }
    *result_len = 0;

    for (int i = 0; i < S_len; i++) {
        Pair current = S[i];
        bool has_clause = false;
        // printf("here_debug\n");
        // 检查是否有任何子句与当前 (T, v) 不一致
        for (int j = 0; j < num_clauses; j++) {
            if (clauses[j].f(current.v) == 1) {

                has_clause = true;
                break;
            }
        }

        if (!has_clause) {
            Sf_rest[*result_len] = current; // 添加到结果中
            (*result_len)++;
        }
    }

    return Sf_rest;
}


// 计算 S^f
Pair* compute_Sf(Clause* clauses, int num_clauses, Pair* S, int S_len, int* result_len) {

    int total_len = 0;
    Pair* Sf = (Pair*)malloc(S_len * sizeof(Pair)); // 最多 S_len 个结果

    // 计算 S^f_rest
    // int rest_len;

    // Pair* Sf_rest = compute_Sf_rest(clauses, num_clauses, S, S_len, &rest_len);

    // for (int i = 0; i < rest_len; i++) {
    //     Sf[total_len++] = Sf_rest[i];
    // }
    // free(Sf_rest);

    // 计算每个 S^f_i
    for (int i = 0; i < num_clauses; i++) {
        int sf_i_len;

        Pair* Sf_i = compute_Sf_i(clauses[i], S, S_len, &sf_i_len);
        for (int j = 0; j < sf_i_len; j++) {
            Sf[total_len++] = Sf_i[j];
        }
        // free(Sf_i);
    }

    *result_len = total_len;
    printf("result_len :%d\n", *result_len);
    return Sf;
}


// // 计算 sk_f
// circuit** compute_sk_f(Pair* Sf, int Sf_len, circuit** msk) {
//     circuit** sk_f = (circuit**)malloc(Sf_len * sizeof(circuit*));

//     for (int i = 0; i < Sf_len; i++) {
//         Pair current = Sf[i];
//         // 计算 sk_{T,v} = P.Eval_{msk}(T || v)
//         // 假设有一个函数 P_Eval 来计算
//         sk_f[i] = initial_prp_circuit(k, encode_Tx(k, )); // 伪代码，具体实现取决于电路库
//     }

//     return sk_f;
// }


circuit*** compute_sk_f(Pair* Sf, int Sf_len, circuit** msk, int k) {
    circuit*** sk_f = (circuit***)malloc(Sf_len * k * sizeof(circuit*));
    if (sk_f == NULL) {
        fprintf(stderr, "内存分配失败用于 sk_f\n");
        exit(1);
    }

    for (int i = 0; i < Sf_len; i++) {
        Pair current = Sf[i];
        // 编码 (T, v) 为 k 位
        circuit** encoded_Tv = encode_pair(k, current, msk[0]);

        // 构建 PRF 的输入：前k位为 encoded_Tv，后k位为 msk
        // 生成 PRF 电路： sk_{T,v} = PRF(msk, encoded_Tv)
        circuit** prf_output = initial_prp_circuit(k, encoded_Tv, msk);
        // circuit** prf_output = encoded_Tv;
        if (prf_output == NULL) {
            fprintf(stderr, "PRF 构建失败用于 sk_f[%d]\n", i);
            exit(1);
        }
        sk_f[i] = prf_output; // 示例：仅存储第一个比特，实际应按需求调整
        // sk_f[i] = msk; 
    }
    return sk_f;
}

// 综合构建电路的主函数
circuit*** build_constrain_circuit(Clause* clauses, int num_clauses, Pair* S, int* result_len, int k, circuit** msk) {

    int S_len = (1 << k) * num_clauses; // 示例长度

    // printf("debug here\n");
    Pair* Sf = compute_Sf(clauses, num_clauses, S, S_len, result_len);
    for(int i =0;i<*result_len;i++){
        printf("Sf[%d] : ", i);
        // printf("Sf t_len :%d \n", Sf[i].t_len);
        for (int j = 0; j < Sf[i].t_len; j++) {
            printf("T[%d] = %d, v[%d] = %d ", j, Sf[i].T[j], j, Sf[i].v[j]);

        }
        printf("\n");
    }
    // printf("Sf_len = %d\n", result_len);

    circuit*** sk_f = compute_sk_f(Sf, S_len, msk, k);

    // 释放 Sf
    free(Sf);

    return sk_f; // 返回计算出的 sk_f
}

