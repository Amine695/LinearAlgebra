#include <time.h>

// some macro to go faster
#define foreach(a,b,c) for(int a = b; a<c ; a++) 
#define forneg(a,b,c) for(int a = b; a>=c ; a--) 
#define for_(a,n) foreach(a,0,n) 

#define DIM int n   // the size of the matrix
#define prime int p // the prime number
typedef int **matrix; // the matrix

// prototypes functions 
void null_matrix(matrix m, DIM);
matrix init_matrix(int s);
matrix loadMatrix(void *s, DIM);
matrix Identity(DIM);
matrix Copy(matrix m,DIM);
matrix generate(DIM,prime);
void freeMatrix(matrix m);
void displayMatrix(char title[], matrix m,DIM);
matrix mat_mult(matrix a, matrix b, DIM, prime);
matrix matmul(matrix a, matrix b, DIM , prime);
void pivot(matrix a, matrix p, DIM);
void LU(matrix A, matrix L, matrix U, matrix P, DIM,prime);
void correction(matrix A, matrix L, matrix U, matrix P, DIM, prime);
void RunLU(matrix A, matrix P, matrix L, matrix U, DIM, prime);