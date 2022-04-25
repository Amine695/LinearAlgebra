#include"LU.h"


void LinearSystem(matrix L, matrix U, int Z[], int B[], int X[],DIM,prime);
void generateSolution(int B[],DIM,prime);
void InitVector(int Z[], int X[], DIM);
void DisplaySolutions(int B[],int X[],int Z[] ,DIM);
void correctionLS(matrix L, matrix U, int Z[], int B[], int X[],DIM, prime);
