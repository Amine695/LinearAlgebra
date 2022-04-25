#include "LU.h"
#define SIZE 5


void InvLU(matrix A,DIM,prime);
void InverseMatrix(matrix A,int s[][SIZE],DIM, prime);
void DisplayInv(int s[][SIZE], DIM);
void correctionInv(matrix A, int s[][SIZE], DIM, prime);