#include "LU.h"

int computeZ(matrix L, matrix I,matrix Z,int col, int row,DIM,prime);
int computeInverse(matrix inverse, matrix U,matrix Z, int col, int row, DIM,prime);
matrix InverseMatrix(matrix I,matrix L, matrix U,matrix P, matrix Z, DIM,prime);
void correctionInv(matrix A, matrix s, DIM, prime);
void RunNaiveInverse(matrix A,matrix L,matrix U,matrix P,DIM,prime);
