#include <time.h>
#define foreach(a,b,c) for(int a = b; a<c ; a++) 
#define forneg(a,b,c) for(int a = b; a>=c ; a--) 
#define for_(a,n) foreach(a,0,n) 
#define DIM int n 
#define prime int p
typedef int **matrix;


void null_matrix(matrix m, DIM);
matrix init_matrix(int s);
matrix loadMatrix(void *s, DIM);
matrix generate(DIM,prime);
void freeMatrix(matrix m);
void displayMatrix(char title[], matrix m,DIM);
matrix mat_mult(matrix a, matrix b, DIM, prime);
void pivot(matrix a, matrix p, DIM);
void LU(matrix A, matrix L, matrix U, matrix P, DIM,prime);
void correction(matrix A, matrix L, matrix U, matrix P, DIM, prime);
void RunLU(matrix A, matrix P, matrix L, matrix U, DIM, prime);