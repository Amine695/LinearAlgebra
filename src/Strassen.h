#include "LU.h"
#include "Tools.h"

matrix SubstractMatrix(matrix A, matrix B, DIM,prime);
matrix ChangeSign(matrix A,DIM,prime);
matrix AddMatrix(matrix A, matrix B,DIM,prime);
matrix Copy(matrix m,DIM);
matrix SubMatrix(matrix A,DIM,int d);
matrix assembly( matrix A11, matrix A12, matrix A21, matrix A22,DIM);
matrix InvStrassen(matrix A,DIM,prime);
matrix StrassenMult(matrix A, matrix B, DIM,prime);
matrix InvStrassenMultStrassen(matrix A,matrix I,matrix L,matrix U,matrix P,matrix Z,DIM,prime);
void displayStrassen(char title[], matrix m,DIM);
void RunStrassenInv(matrix A,DIM,prime);
void RunStrassenMult(matrix A,matrix B,DIM,prime);