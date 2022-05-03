#include "LU.h"
#include "Tools.h"

void BenchmarkNaiveInverse(matrix A,matrix L,matrix U,matrix P,DIM,prime, FILE *fp);
void BenchmarkStrassenInv(matrix A,DIM,prime, FILE *fp);
void BenchmarksNaiveMul(matrix A,matrix B,DIM,prime, FILE *fp);
void BenchmarkStrassenMult(matrix A, matrix B, DIM, prime, FILE *fp);
void BenchmarkStrassenInvAndMult(matrix A, DIM, prime, FILE *fp);