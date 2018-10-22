#include <iostream>
#include <math.h>

#include "../NR3/nr3.h"
#include "../NR3/cholesky.h"
#include "../NR3/gaussj.h"
#include "ManageMat_OPS.h"
using namespace std;


// Converts float array to MatDoub
void float2Mat(float **fmat, MatDoub &M, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            M[i][j] = fmat[i][j];
        }
    }
}

// Converts MatDoub to float array
void mat2float(MatDoub &M, float **fmat, int rows, int cols)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            fmat[i][j] = M[i][j];
        }
    }
}

void vec2float(VecDoub &M, float *fmat, int rows)
{
    for (int i = 0; i < rows; i++)
    {
        fmat[i] = M[i];
    }
}




void displayMat(MatDoub M)
{
    cout << "[\n";
    for (int i = 0; i < M.nrows(); i++)
    {
        for (int j = 0; j < M.ncols(); j++)
        {
            cout << setw(20) << M[i][j];
        }
        cout << endl;
    }
    cout <<"]" << endl;
}

void displayMat(MatDoub M, char* x)
{
    cout << x << " =\n[" << endl;
    for (int i = 0; i < M.nrows(); i++)
    {
        for (int j = 0; j < M.ncols(); j++)
        {
            cout << setw(20) << M[i][j];
        }
        cout << endl;
    }
    cout << "]\n";
}


// Calculates inverse of A using gaussian elimination and stores in A_inverse
void inverse(MatDoub &A, MatDoub &A_inverse)
{
    MatDoub C;
    C = A;
    gaussj(C);
    A_inverse = C;
}



// Stores transpose of A in B
void transpose_of_any_matrix(MatDoub &A, MatDoub &B)
{
    int rows = A.nrows();
    int cols = A.ncols();
    MatDoub Transpose(cols, rows);
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            Transpose[j][i] = A[i][j];
        }
    }
    B = Transpose;
}

// matrix_multiplication
void matrix_multiplication(MatDoub &A, MatDoub &B, MatDoub &C)
{
    MatDoub result(C.nrows(), C.ncols());
    for(int i=0; i < C.nrows(); i++)
    {
        for(int j=0; j < C.ncols(); j++)
        {
            double temp = 0;
            for(int a=0; a < C.ncols(); a++)
            {
                temp += A[i][a] * B[a][j];
            }
            result[i][j] = temp;
        }
    }
    C = result;

}




// Performs Cholesky Decomposition on matrix A
// Output: L, I
void choleskyd(MatDoub &A, MatDoub &L, MatDoub &I)
{
    Cholesky achol(A);
    L = achol.el;
    achol.inverse(I);
}

void choleskyd(MatDoub &A, MatDoub &L)
{
    Cholesky achol(A);
    L = achol.el;
}

void  restart_mat(MatDoub &mat_iter)
{

    int NFEATURES = mat_iter.nrows();

    for(int i=0; i<NFEATURES; i++)
        for(int j=0; j<NFEATURES; j++)
            mat_iter[i][j] = 0;


}
