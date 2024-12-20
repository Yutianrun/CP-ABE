#pragma once

#include <stdbool.h>

#include "common.h"

/*
A matrix is a struct _matrix
We can modify it in place because the data
is stored outside the struct (only the pointer
to it get copied every time)
*/
typedef struct _matrix {
    unsigned int rows;
    unsigned int columns;
    scalar* data;
} matrix;

/*
A signed_matrix is a struct _signed_matrix
We can modify it in place because the data
is stored outside the struct (only the pointer
to it get copied every time)
*/
typedef struct _signed_matrix {
    unsigned int rows;
    unsigned int columns;
    signed_scalar* data;
} signed_matrix;

// Accessing (signed_)scalar Mi,j where M has type (signed_)matrix
#define matrix_element(M, i, j) ((M).data[(i) * (M).columns + (j)])

// Heap-allocate data filled with 0s and returns the corresponding matrix
matrix new_matrix(unsigned int rows, unsigned int columns);

// Heap-allocate data filled with 0s and returns the corresponding signed_matrix
signed_matrix new_signed_matrix(unsigned int rows, unsigned int columns);

// Returns a heap-allocated array of matrixes filled with 0s
matrix* new_matrixes(int n, unsigned int rows, unsigned int columns);

// Clean free of a matrix : freeing data
void free_matrix(matrix M);

// Clean free of a signed_matrix : freeing data
void free_signed_matrix(signed_matrix M);

// Clean free of a matrixes array
void free_matrixes(matrix* A, int n);

void print_matrix(matrix M);

void print_signed_matrix(signed_matrix M);

// Returns a fresh heap-allocated copy of A's data
matrix copy_matrix(matrix M);

// R <- A + B
void add_matrix(matrix A, matrix B, matrix R);

// R <- A + E
void add_matrix_error(matrix A, signed_matrix E, matrix R);

// R <- A - B
void sub_matrix(matrix A, matrix B, matrix R);

// R <- A * B
void mul_matrix(matrix A, matrix B, matrix R);

// R <- A * Tf where A in Zq^{d1 * d2} and Tf in Z^{d2 * d2}
void mul_matrix_trap(matrix A, signed_matrix Tf, matrix R);

// R <- T * B where T in Z^{d1 * d2} and B in Zq^{d2 * d3}
void mul_matrix_trap_left(signed_matrix T, matrix B, matrix R);

// R <- x * A where x is a scalar
void mul_matrix_scalar(scalar x, matrix A, matrix R);

bool equals(matrix A, matrix B);

/*
Returns "norm" of matrix A = sqrt(sum |Tij|^2) where
Tij = min(Aij, Q - Aij) is the minimal distance to 0 in Zq !!!
*/
real norm(matrix A);

// Returns euclidian norm of matrix A = sqrt(sum |Aij|^2)
real norm_signed(signed_matrix A);

/*
Returns true only if A is short
Check SHORT_THRESHOLD in common.h to see
the (computational) definition used
*/
bool is_short(matrix A);
