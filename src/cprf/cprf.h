#ifndef CPRF_H
#define CPRF_H

#include "../bgg/circuit.h"
#include "../bgg/gen_circuit.h"
// #include "params.h" // 假设 Params 和 Clause 定义在 params.h 中

#ifdef __cplusplus
extern "C" {
#endif

#define PRF_K 128

// Clause 结构体定义
typedef struct {
    int* T;
    int t_size;
} Clause;

// 外部参数
// extern Params PARAMS;

// 函数声明
circuit** build_prf_circuit(int k);
circuit* extract_v(int k, int* T, int t_size, circuit** x);
circuit* build_eval_prf(int k, Clause* clauses, int num_clauses);

#ifdef __cplusplus
}
#endif

#endif // CPRF_H