/**
 * @file Benchmarks.c
 * @author Amine Berbagui
 * @brief Launch the benchmarks of our different algorithms.
 * @date 2022-05-05
 * @copyright Copyright (c) 2022
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include "Benchmarks.h"
#include "Strassen.h"
#include "MatrixInv.h"


/**
 * @brief Launch the benchmarks for the Naive Inversion matrix
*/
void BenchmarkNaiveInverse(matrix A,matrix L,matrix U,matrix P,DIM,prime, FILE *fp)
{
    
    matrix identity = Identity(n);
    matrix Z = init_matrix(n);
    double start,end,time_elapsed;
   
    LU(A,L,U,P,n,p);
    start = clock();
    matrix Inv = InverseMatrix(identity,L,U,P,Z,n,p);
    end = clock();
    time_elapsed = ((double)end - start) / CLOCKS_PER_SEC;
    printf("Naive Inverse              n = %d : %f s\n",n,time_elapsed);
    fprintf(fp,"\n%-4.1d\t\t%4.5f",n,time_elapsed);

    freeMatrix(Inv);
    freeMatrix(identity);
    freeMatrix(Z);
}

/**
 * @brief Launch the benchmarks for the Strassen's Inversion matrix
*/
void BenchmarkStrassenInv(matrix A,DIM,prime, FILE *fp)
{
    double start,end,time_elapsed;
    start = clock();
    matrix inv = InvStrassen(A,n,p);
    end = clock();
    time_elapsed = ((double)end - start) / CLOCKS_PER_SEC;

    printf("Strassen's inverse         n = %d : %f s\n",n,time_elapsed);
    fprintf(fp,"\n%-4.1d\t\t%4.5f",n,time_elapsed);

    freeMatrix(inv);
  
}

/**
 * @brief Launch the benchmarks for the Naive multiplication matrix
*/
void BenchmarksNaiveMul(matrix A,matrix B,DIM,prime, FILE *fp)
{
    double start,end,time_elapsed;
    start = clock();
    matrix res = matmul(A,B,n,p);
    end = clock();
    time_elapsed = ((double)end - start) / CLOCKS_PER_SEC;

    printf("Naive multiplication       n = %d : %f s\n",n,time_elapsed);
    fprintf(fp,"\n%-4.1d\t\t%4.5f",n,time_elapsed);

    freeMatrix(res);
}

/**
 * @brief Launch the benchmarks for the Strassen's multiplication matrix
*/
void BenchmarkStrassenMult(matrix A, matrix B, DIM, prime, FILE *fp)
{
    double start,end,time_elapsed;
    start = clock();
    matrix res = StrassenMult(A,B,n,p);
    end = clock();
    time_elapsed = ((double)end - start) / CLOCKS_PER_SEC;

    printf("Strassen's multiplication  n = %d : %f s\n",n,time_elapsed);
    fprintf(fp,"\n%-4.1d\t\t%4.5f",n,time_elapsed);


    freeMatrix(res);
}

/**
 * @brief Launch the benchmarks for the Strassen's inversion matrix using strassen's multiplication algorithm
 * use the naive inversion for time improvements
*/
void BenchmarkStrassenInvAndMult(matrix A, DIM, prime, FILE *fp)
{
    matrix I = Identity(n);
    matrix Z = init_matrix(n);
    matrix L = init_matrix(n);
    matrix U = init_matrix(n);
    matrix P = init_matrix(n);
    LU(A,L,U,P,n,p);
    double start,end,time_elapsed;
    start = clock();
    matrix res = InvStrassenMultStrassen(A,I,L,U,P,Z,n,p);
    end = clock();
    time_elapsed = ((double)end - start) / CLOCKS_PER_SEC;

    printf("Strassen's inv and mult    n = %d : %f s\n",n,time_elapsed);
    fprintf(fp,"\n%-4.1d\t\t%4.5f",n,time_elapsed);

    freeMatrix(I);
    freeMatrix(Z);
    freeMatrix(L);
    freeMatrix(U);
    freeMatrix(P);
    freeMatrix(res);
}