#include "LU.h"
#include "Tools.h"

matrix SubstractMatrix(matrix A, matrix B, DIM,prime);
void ChangeSign(matrix A,DIM,prime);
matrix matmul(matrix a, matrix b, DIM,prime);
matrix SubMatrix(matrix A,DIM,int d);
matrix assembly( matrix A11, matrix A12, matrix A21, matrix A22,DIM);
matrix InvStrassen(matrix A,int s,prime);
void displayStrassen(char title[], matrix m,DIM);
