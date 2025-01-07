#include "../bgg/circuit.h"
#include "../bgg/gen_circuit.h"
#ifndef CPRF_H
#define CPRF_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

// ClauseT结构体定义
// typedef struct {
//     int* T;       // T 的位置数组
//     int t_len;    // T 的长度
//     int* v;       // 对应位置上的值数组
//     int v_len;    // v 的长度（应与 t_len 相同）
// } ClauseT;

typedef struct {
    int* T;       // T 的位置数组
    int t_len;    // T 的长度
} ClauseT;

typedef struct {
    int* T;       // T 的位置数组
    int t_len;    // T 的长度
    circuit *f; // 子句函数
} ClauseF;


typedef struct {
    int* T;       // T 的位置数组
    int t_len;    // T 的长度
    int (*f)(int*); // 函数 f 的指针
} Clause;

typedef struct {
    int* T;       // T 的位置数组
    int t_len; 
    int* v;       // 对应位置上的值数组
} Pair;

// 函数声明
circuit** build_prp_circuit(int k);
circuit** initial_prp_circuit(int k, circuit** x, circuit** msk);
circuit** encode_Tx(int k, ClauseT clause, circuit** x);
circuit*** allocate_sk_tv(int num_clauses, int k);
void free_sk_tv(circuit*** sk_tv, int num_clauses);
circuit*** build_sk_tv(int k, ClauseT* clauses, int num_clauses, circuit** msk, circuit** x);
circuit** final_prf(int k, circuit*** sk_tv, int num_clauses, circuit** x);
circuit** build_eval_circuit(int k, ClauseT* clauses, int num_clauses, circuit** msk, circuit** x);

#ifdef __cplusplus
}
#endif

#endif // CPRF_H