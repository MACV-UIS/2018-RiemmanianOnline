#ifndef MANAGEMATOPS_H
#define MANAGEMATOPS_H

#include "../NR3/nr3.h"

// Converts float array to MatDoub
void float2Mat(float **fmat, MatDoub &M, int rows, int cols);
// Converts MatDoub to float array
void mat2float(MatDoub &M, float **fmat, int rows, int cols);

void vec2float(VecDoub &M, float *fmat, int rows);

void displayMat(MatDoub M);
void displayMat(MatDoub M, char* x);


// Calculates inverse of A using gaussian elimination and stores in A_inverse
void inverse(MatDoub &A, MatDoub &A_inverse);

// Stores transpose of A in B
void transpose_of_any_matrix(MatDoub &A, MatDoub &B);

void matrix_multiplication(MatDoub &A, MatDoub &B, MatDoub &C);

// Performs Cholesky Decomposition on matrix A
// Output: L, I
void choleskyd(MatDoub &A, MatDoub &L, MatDoub &I);

void choleskyd(MatDoub &A, MatDoub &L);
void  restart_mat(MatDoub &mat_iter);


#endif // MANAGEMATOPS_H


