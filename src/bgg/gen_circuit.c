#include "gen_circuit.h"

#include "common.h"

circuit* gen_leaf(int n, bool xn) {
    // n is still 1-indexed
    circuit* f = new_circuit();
    f->left = f->right = NULL;
    f->n = n;
    if (xn) return f;
    return circuit_not(f);
}

circuit* circuit_not(circuit* f) {
    // NOT(A) = NAND(A, A)
    circuit* g = new_circuit();
    g->left = g->right = f;
    g->right_double_pointed = true;
    return g;
}

circuit* circuit_and(circuit* f, circuit* g) {
    // AND(A, B) = NOT(C) where C = NAND(A, B)
    circuit* c = new_circuit();
    c->left = f;
    c->right = g;
    return circuit_not(c);
}


circuit* circuit_recurssive_and(circuit** f, int k) {
    // AND(A, B) = NOT(C) where C = NAND(A, B)

    if (k==2)
    {
        return circuit_and(*f,*(f+1));/* code */  
    }
    
    return circuit_and(circuit_recurssive_and(f,k/2), circuit_recurssive_and(f+k/2, k/2));
}


circuit* circuit_consecutive_and(circuit** f, int k) {
    // AND(A, B) = NOT(C) where C = NAND(A, B)

    if (k == 2) {
        return circuit_and(*f, *(f + 1)); // 处理两个电路的 AND
    }

    // 递归处理 k - 1 个电路的 AND
    circuit* partial_result = circuit_consecutive_and(f + 1, k - 1);
    return circuit_and(*f, partial_result); // 处理当前电路与递归结果的 AND
}



circuit* circuit_or(circuit* f, circuit* g) {
    // OR(A, B) = NAND(NOT(A), NOT(B))
    circuit* o = new_circuit();
    o->left = circuit_not(f);
    o->right = circuit_not(g);
    return o;
}

// 辅助函数实现 XOR
circuit* circuit_xor(circuit* a, circuit* b) {
    // XOR(A, B) = (A AND NOT(B)) OR (NOT(A) AND B)
    // a AND NOT(b)
    circuit* not_b = circuit_not(b);

    circuit* a_and_not_b = new_circuit();
    a_and_not_b->left = a;
    a_and_not_b->right = not_b;

    // NOT(a) AND b
    circuit* not_a = circuit_not(a);
    not_a->left_double_pointed = true;

    circuit* not_a_and_b = new_circuit();
    not_a_and_b->left = not_a;
    not_a_and_b->right = b;
    not_a_and_b->right_double_pointed = true;
    

    // OR(a AND NOT(b), NOT(a) AND b)
    circuit* xor_result = new_circuit();
    xor_result->left = a_and_not_b;
    xor_result->right = not_a_and_b;


    return xor_result;
}
// circuit* circuit_eq(circuit* a, circuit* b) {
//     // XOR(a, b) using existing circuit_xor
//     circuit* xor_ab = circuit_xor(a, b);
//     // EQ(a, b) = NOT(XOR(a, b))
//     circuit* eq = circuit_not(xor_ab);
//     return eq;
// }

//  ((((1 ^ (2 ^ 2)) ^ (1 ^ (2 ^ 2))) ^ ((1 ^ (2 ^ 2)) ^ (1 ^ (2 ^ 2)))) ^ ((((1 ^ 1) ^ 2) ^ ((1 ^ 1) ^ 2)) ^ (((1 ^ 1) ^ 2) ^ ((1 ^ 1) ^ 2))))
circuit* gen_circuit(attribute x) {
    /*
    We denote Gi the logic gate associated to xi.
    If xi = 1 : Gi is a leaf with Gi.n = i
    Else (xi = 0) : Gi = not (i) ie a circuit with
    two children such Gil.n = Gir.n = i (not gate using nand gates).
    Then the final circuit is just NOT( AND(Gi) forall i).
    Beware of NOT gate !
    */
    circuit* f;

    // Forall AND(Gi)
    circuit* G1 = gen_leaf(1, get_xn(x, 1));
    circuit* G2 = gen_leaf(2, get_xn(x, 2));
    f = circuit_and(G1, G2);
    for (int n = 3; n < PARAMS.K + 1; n++) {
        circuit* Gn = gen_leaf(n, get_xn(x, n));
        f = circuit_and(f, Gn);
    }

    // NOT gate
    circuit* g = new_circuit();
    g->left = g->right = f;
    return g;
}