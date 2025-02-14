#pragma once

#include "attribute.h"
#include "matrix.h"
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
/*
Computes B = [B0 | B1,0 | B1,1 | ... | Bk,1]
and its trap T such TBi,b = 0 where Bi,b in Zq^{m * n} and T in Z^{m * m}.
Warning, some dimensions are hardcoded in TrapGen and we use the gadget vector
(with g=2).
*/
void TrapGen(matrix* B, signed_matrix T);
void single_TrapGen(matrix B, matrix T);

/*
Given B = [B0 | B1,0 | B1,1 | ... | Bk,1], its trap T such as TB0=TBi,b=0,
and an attribute x returns trap Tx (in Z^{P * M}) for the attribute :
- TxB0 = 0
- TxBi,xk = 0
- TxBi,(1-xk) is uniformely distributed
Warning dimensions are hardcoded in TrapGen and we use the gadget vector (g=2).
*/
signed_matrix TrapSamp(matrix* B, signed_matrix T, attribute x);

/* -------------------- */
/* Functions for matrix */
/* -------------------- */

// Init sampling functions
void init_sampler(void);

// Assigns uniformely scalars in Zq over a random matrix
void sample_Zq_uniform_matrix(matrix A);
void sample_Zq_uniform_matrix_64(matrix A);
// Samples an entire matrix from D_{Zq,sigma}
void sample_Z_centered_matrix(signed_matrix A);
void SampleD(matrix A, matrix u, matrix R, double sigma, signed_matrix result);
void sample_Zq_invertible_matrix(matrix A);
void SampleG(scalar q, scalar u, double s, int k, int* result);
void SamplePre(matrix A, matrix R, matrix u, double sigma, signed_matrix result);
void SampleLeft(matrix A, matrix R, matrix B, matrix u, double sigma, matrix result);