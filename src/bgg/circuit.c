#include "circuit.h"

#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include <cprf.h>
#include "common.h"
#include "matrix.h"

/***************/
/* Computing G */
/***************/
#define EPSILON 1e-10
matrix G;

void init_G() {
    // We need that L = N * K
    G = new_matrix(PARAMS.N, PARAMS.L);
    for (int i = 0; i < PARAMS.N; i++)
        for (int k = 0; k < PARAMS.K; k++)
            matrix_element(G, i, i * PARAMS.K + k) = ((int64_t)1 << k) % PARAMS.Q;
}

// R <- G^-1(A)
void inv_G(matrix A, matrix R) {
    for (int n = 0; n < PARAMS.N; n++) {
        for (int l = 0; l < PARAMS.L; l++) {
            for (int k = 0; k < PARAMS.K; k++) {
                scalar Anl = matrix_element(A, n, l);
                matrix_element(R, n * PARAMS.K + k, l) = (Anl >> k) & 1;
            }
        }
    }
}

bool is_invertible(matrix A) {
    if (A.rows != A.columns) {
        return false;
    }
    
    unsigned int n = A.rows;
    matrix L = new_matrix(n, n);
    
    // 复制原矩阵用于LU分解
    matrix U = copy_matrix(A);
    
    // LU分解
    for(unsigned int i = 0; i < n; i++) {
        // 对角线元素设为1
        matrix_element(L, i, i) = 1;
        
        for(unsigned int j = i + 1; j < n; j++) {
            if(fabs(matrix_element(U, i, i)) < EPSILON) {
                free_matrix(L);
                free_matrix(U);
                return false;
            }
            
            scalar factor = matrix_element(U, j, i) / matrix_element(U, i, i);
            matrix_element(L, j, i) = factor;
            
            for(unsigned int k = i; k < n; k++) {
                matrix_element(U, j, k) -= factor * matrix_element(U, i, k);
            }
        }
    }
    
    // 计算行列式(对角线元素乘积)
    scalar det = 1;
    for(unsigned int i = 0; i < n; i++) {
        det *= matrix_element(U, i, i);
    }
    
    free_matrix(L);
    free_matrix(U);
    
    return fabs(det) >= EPSILON;
}   


matrix matrix_inverse(matrix A) {
    if (A.rows != A.columns) {
        fprintf(stderr, "Error: Matrix must be square for inversion\n");
        return new_matrix(0, 0);
    }
    
    unsigned int n = A.rows;
    matrix augmented = new_matrix(n, 2*n);
    
    // 填充增广矩阵
    for(unsigned int i = 0; i < n; i++) {
        for(unsigned int j = 0; j < n; j++) {
            matrix_element(augmented, i, j) = matrix_element(A, i, j);
            matrix_element(augmented, i, j+n) = (i == j) ? 1 : 0;
        }
    }
    
    // 高斯-约旦消元
    for(unsigned int i = 0; i < n; i++) {
        // 查找主元
        scalar pivot = matrix_element(augmented, i, i);
        if(fabs(pivot) < EPSILON) {
            fprintf(stderr, "Error: Matrix is singular\n");
            free_matrix(augmented);
            return new_matrix(0, 0);
        }
        
        // 归一化当前行
        for(unsigned int j = 0; j < 2*n; j++) {
            matrix_element(augmented, i, j) /= pivot;
        }
        
        // 消元
        for(unsigned int k = 0; k < n; k++) {
            if(k != i) {
                scalar factor = matrix_element(augmented, k, i);
                for(unsigned int j = 0; j < 2*n; j++) {
                    matrix_element(augmented, k, j) -= factor * matrix_element(augmented, i, j);
                }
            }
        }
    }
    
    // 提取逆矩阵
    matrix inverse = new_matrix(n, n);
    for(unsigned int i = 0; i < n; i++) {
        for(unsigned int j = 0; j < n; j++) {
            matrix_element(inverse, i, j) = matrix_element(augmented, i, j+n);
        }
    }
    
    free_matrix(augmented);
    return inverse;
}
/*****************/
/* Circuit utils */
/*****************/

circuit* new_circuit() {
    circuit* f = calloc(1, sizeof(circuit));
    f->left = NULL;
    f->right = NULL;
    f->n = 0;
    return f;
}


void circuit_to_string(circuit f, char* buffer, size_t size) {
    if (!f.left && !f.right) {
        snprintf(buffer + strlen(buffer), size - strlen(buffer), "%d", f.n);
        return;
    }
    snprintf(buffer + strlen(buffer), size - strlen(buffer), "(");
    circuit_to_string(*f.left, buffer, size);
    snprintf(buffer + strlen(buffer), size - strlen(buffer), " ^ ");
    circuit_to_string(*f.right, buffer, size);
    snprintf(buffer + strlen(buffer), size - strlen(buffer), ")");
}

void free_circuit(circuit* f) {

    if (!f->left && !f->right) {
        free(f);
        return;
    }
    if(!f->left_double_pointed)
    {
            free_circuit(f->left);
    }
    if(!f->right_double_pointed)
    {
            free_circuit(f->right);
    }

    free(f);
}

void print_circuit(circuit f) {
    assert((f.left && f.right) || (!f.left && !f.right));

    if (!f.left && !f.right) {
        printf("%d", f.n);
        return;
    }

    printf("(");
    print_circuit(*f.left);
    printf(" ^ ");
    print_circuit(*f.right);
    printf(")");
}

unsigned depth(circuit f) {
    assert((f.left && f.right) || (!f.left && !f.right));

    if (!f.left && !f.right) return 1;

    if (f.left == f.right) return 1 + depth(*f.left);

    unsigned dl = depth(*f.left);
    unsigned dr = depth(*f.right);
    unsigned m = (dl < dr) ? dr : dl;
    return m + 1;
}

/*
Returns A * G^-1(B) - G
where G is the gadget matrix and G^-1
the inverse function defined as follows in 2020-191
We define the inverse function G−1 : Zq^{n×m} → {0, 1}^{N×m}
which expands each entry a ∈ Zq of the input matrix into
a column of size K = ceil(log q) consisting of the bits of the
binary representation of a. We have the property that
for any matrix A ∈ Zq^{n×l}, it holds that G*G−1(A) = A.
*/
matrix nand(matrix A, matrix B) {
    // We need to heap-allocate for the matrix to survive
    matrix R = new_matrix(PARAMS.N, PARAMS.L);
    matrix temp = new_matrix(PARAMS.L, PARAMS.L);
    inv_G(B, temp);          // temp <- G^-1(B)
    mul_matrix(A, temp, R);  // R <- A * temp = A * G^-1(B)
    sub_matrix(R, G, R);     // R <- R - G = A * G^-1(B) - G
    free_matrix(temp);
    return R;
}

// Returns Af = f(A) = f(A1, ..., Ak)
matrix compute_Af(matrix* A, circuit f) {
    assert((f.left && f.right) || (!f.left && !f.right));

    if (!f.left && !f.right) return copy_matrix(A[f.n]);  // A fresh copy of An

    // f.left may be equal to f.right (exact same circuit cause exact same
    // pointer) then we only need to compute recursively on *f.left but if
    // f.left != f.right we need to recursively compute on *f.right too
    matrix R_left, R_right;
    R_left = R_right = compute_Af(A, *f.left);
    if (f.left != f.right) R_right = compute_Af(A, *f.right);
    matrix R = nand(R_left, R_right);

    free_matrix(R_left);
    if (f.left != f.right) free_matrix(R_right);

    return R;
}

/***************/
/* Computing H */
/***************/

typedef struct H_triplet {
    matrix A;
    bool x;
    matrix H;
} H_triplet;

H_triplet new_H_triplet() {
    matrix A = new_matrix(PARAMS.N, PARAMS.L);
    matrix H = new_matrix(PARAMS.Att_num * PARAMS.L, PARAMS.L);
    H_triplet t;
    t.A = A;
    t.x = 0;
    t.H = H;
    return t;
}

H_triplet new_H_triplet_prfk() {
    matrix A = new_matrix(PARAMS.N, PARAMS.L);
    matrix H = new_matrix(PRF_K * PARAMS.L, PARAMS.L);
    H_triplet t;
    t.A = A;
    t.x = 0;
    t.H = H;
    return t;
}

void free_H_triplet(H_triplet* t) {
    free_matrix(t->A);
    free_matrix(t->H);
}

H_triplet leaf(matrix* A, attribute x, int n) {
    H_triplet t = new_H_triplet();
    t.A = copy_matrix(A[n]);
    t.x = get_xn(x, n);
    // H seen as a column is empty except
    // in n-th position which is the identity
    for (int i = 0; i < PARAMS.L; i++)
        matrix_element(t.H, (n - 1) * PARAMS.L + i, i) = 1;
    return t;
}

H_triplet leaf_prfk(matrix* A, attribute x, int n) {
    H_triplet t = new_H_triplet_prfk();
    t.A = copy_matrix(A[n]);
    t.x = get_xn(x, n);
    // H seen as a column is empty except
    // in n-th position which is the identity
    for (int i = 0; i < PARAMS.L; i++)
        matrix_element(t.H, (n - 1) * PARAMS.L + i, i) = 1;
    return t;
}

H_triplet compute_H_triplet(matrix* A, circuit f, attribute x) {
    assert((f.left && f.right) || (!f.left && !f.right));

    if (!f.left && !f.right) return leaf(A, x, f.n);

    // f.left may be equal to f.right (exact same circuit cause exact same
    // pointer) then we only need to compute recursively on *f.left but if
    // f.left != f.right we need to recursively compute on *f.right too
    H_triplet tl, tr;
    tl = tr = compute_H_triplet(A, *f.left, x);
    if (f.left != f.right) tr = compute_H_triplet(A, *f.right, x);

    H_triplet t = new_H_triplet();

    matrix inv = new_matrix(PARAMS.L, PARAMS.L);
    matrix tempA = new_matrix(PARAMS.N, PARAMS.L);
    matrix tempH = new_matrix(PARAMS.Att_num * PARAMS.L, PARAMS.L);

    // Computing new A = Al * G^-1(Ar) - G
    inv_G(tr.A, inv);
    mul_matrix(tl.A, inv, tempA);
    sub_matrix(tempA, G, tempA);
    free_matrix(t.A);
    t.A = copy_matrix(tempA);

    // Computing new x = 1 - xl * xr
    t.x = 1 - tl.x * tr.x;

    // Computing new H = Hl * G^-1(Ar) - xl * Hr
    mul_matrix(tl.H, inv, tempH);
    if (tl.x) sub_matrix(tempH, tr.H, tempH);
    free_matrix(t.H);
    t.H = copy_matrix(tempH);

    // Free time !
    free_matrix(inv);
    free_matrix(tempA);
    free_matrix(tempH);
    free_H_triplet(&tl);
    if (f.left != f.right) free_H_triplet(&tr);

    return t;
}

H_triplet compute_H_triplet_prfk(matrix* A, circuit f, attribute x) {
    assert((f.left && f.right) || (!f.left && !f.right));

    if (!f.left && !f.right) return leaf_prfk(A, x, f.n);

    // f.left may be equal to f.right (exact same circuit cause exact same
    // pointer) then we only need to compute recursively on *f.left but if
    // f.left != f.right we need to recursively compute on *f.right too
    H_triplet tl, tr;
    tl = tr = compute_H_triplet_prfk(A, *f.left, x);
    if (f.left != f.right) tr = compute_H_triplet_prfk(A, *f.right, x);

    H_triplet t = new_H_triplet_prfk();

    matrix inv = new_matrix(PARAMS.L, PARAMS.L);
    matrix tempA = new_matrix(PARAMS.N, PARAMS.L);
    matrix tempH = new_matrix(PRF_K * PARAMS.L, PARAMS.L);

    // Computing new A = Al * G^-1(Ar) - G
    inv_G(tr.A, inv);
    mul_matrix(tl.A, inv, tempA);
    sub_matrix(tempA, G, tempA);
    free_matrix(t.A);
    t.A = copy_matrix(tempA);

    // Computing new x = 1 - xl * xr
    t.x = 1 - tl.x * tr.x;

    // Computing new H = Hl * G^-1(Ar) - xl * Hr
    mul_matrix(tl.H, inv, tempH);
    if (tl.x) sub_matrix(tempH, tr.H, tempH);
    free_matrix(t.H);
    t.H = copy_matrix(tempH);

    // Free time !
    free_matrix(inv);
    free_matrix(tempA);
    free_matrix(tempH);
    free_H_triplet(&tl);
    if (f.left != f.right) free_H_triplet(&tr);

    return t;
}



matrix compute_H(matrix* A, circuit f, attribute x) {
    H_triplet t = compute_H_triplet(A, f, x);
    matrix H = copy_matrix(t.H);
    // free_H_triplet(&t);
    return H;
}

matrix compute_H_prfk(matrix* A, circuit f, attribute x) {
    H_triplet t = compute_H_triplet_prfk(A, f, x);
    matrix H = copy_matrix(t.H);
    // free_H_triplet(&t);
    return H;
}
// 计算模Q下的逆元
scalar mod_inverse(scalar a, scalar m) {
    scalar m0 = m;
    scalar y = 0, x = 1;
    if (m == 1) return 0;
    while (a > 1) {
        scalar q = a / m;
        scalar t = m;
        m = a % m;
        a = t;
        t = y;
        y = x - q * y;
        x = t;
    }
    if (x < 0) x += m0;
    return x;
}

matrix matrix_inverse_mod_q(matrix* A) {
    assert(A->rows == A->columns);
    int n = A->rows;
    
    // 创建增广矩阵 [A|I]
    matrix aug = new_matrix(n, 2 * n);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            matrix_element(aug, i, j) = matrix_element(*A, i, j);
        }
        matrix_element(aug, i, i + n) = 1;
    }
    
    // 高斯-约旦消元
    for(int i = 0; i < n; i++) {
        // 找主元
        scalar pivot = matrix_element(aug, i, i);
        scalar pivot_inv = mod_inverse(pivot, PARAMS.Q);
        
        // 将主对角线元素化为1
        for(int j = 0; j < 2 * n; j++) {
            matrix_element(aug, i, j) = 
                (matrix_element(aug, i, j) * pivot_inv) % PARAMS.Q;
        }
        
        // 将其他行的这一列化为0
        for(int k = 0; k < n; k++) {
            if(k != i) {
                scalar factor = matrix_element(aug, k, i);
                for(int j = 0; j < 2 * n; j++) {
                    scalar temp = (matrix_element(aug, i, j) * factor) % PARAMS.Q;
                    matrix_element(aug, k, j) = 
                        (PARAMS.Q + matrix_element(aug, k, j) - temp) % PARAMS.Q;
                }
            }
        }
    }
    
    // 提取右半部分作为逆矩阵
    matrix inv = new_matrix(n, n);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            matrix_element(inv, i, j) = matrix_element(aug, i, j + n);
        }
    }
    
    free_matrix(aug);
    return inv;
}

matrix compute_H_from_A_Af(matrix* A, matrix* Af) {
    // A: 1×400, Af: 1×20, 目标 H: 400×20
    // 1. 计算 A^T (400×1)
    matrix AT = new_matrix(A->columns, A->rows);
    for (int i = 0; i < A->columns; i++) {
        for (int j = 0; j < A->rows; j++) {
            matrix_element(AT, i, j) = matrix_element(*A, j, i);
        }
    }
    
    // 2. 计算 ATA = A^T × A (400×400)
    matrix ATA = new_matrix(AT.rows, A->columns);
    mul_matrix(AT, *A, ATA);
    
    // 3. 计算 (ATA)^(-1) mod Q (400×400)
    matrix ATA_inv = matrix_inverse_mod_q(&ATA);
    
    // 4. 计算 temp = ATA_inv × A^T (400×1)
    matrix temp = new_matrix(ATA_inv.rows, AT.columns);
    mul_matrix(ATA_inv, AT, temp);
    
    // 5. 计算 H = temp × Af
    // temp: (400×1) 与 Af: (1×20) 得到 H: (400×20)
    matrix H = new_matrix(temp.rows, Af->columns);
    mul_matrix(temp, *Af, H);
    
    // 清理临时矩阵
    free_matrix(AT);
    free_matrix(ATA);
    free_matrix(ATA_inv);
    free_matrix(temp);
    
    return H;
}

bool compute_f(circuit f, attribute x) {
    assert((f.left && f.right) || (!f.left && !f.right));

    if (!f.left && !f.right) return get_xn(x, f.n);

    // f.left may be equal to f.right (exact same circuit cause exact same
    // pointer) then we only need to compute recursively on *f.left but if
    // f.left != f.right we need to recursively compute on *f.right too
    bool xl, xr;
    xl = xr = compute_f(*f.left, x);
    if (f.left != f.right) xr = compute_f(*f.right, x);
    return 1 - xl * xr;
}
