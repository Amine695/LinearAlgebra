#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include "LinearSystem.h"
#include "MatrixInv.h"
#include "Strassen.h"



int main(void)
{
    //int k = 2;
    //int n = pow(2,k);
    int n = 4;
    int p = 293;
    //clock_t start,end;
    //double time_elapsed;
    srand(time(NULL));
    printf("************** Algorithmic Project ***************\n");
    //sleep(1);
    printf("************** LU decomposition ***************\n");
    //printf("\nEnter the dimension of the matrix: ");
    //scanf("%d",&n);

     // Initialisation
    matrix A, P, L, U;
    L = init_matrix(n);
    P = init_matrix(n);
    U = init_matrix(n);
/*

    // Load A3 into A
    //A = loadMatrix(A3,n);
    A = generate(n,p);
    
    // Computation and time measurements
    start = clock();
    LU(A,L,U,P,n,p);
    end = clock();
    time_elapsed = ((double)end - start) / CLOCKS_PER_SEC;
    
    // Show results
    printf("A matrix:\n\n");
    displayMatrix(A,n);
    printf("\nL matrix:\n\n");
    displayMatrix(L,n);
    printf("\nU matrix:\n\n");
    displayMatrix(U,n);
    printf("\nP matrix:\n\n");
    displayMatrix(P,n);

     // Check if PA == LU
    correction(A,L,U,P,n,p);
    printf("Time duration : %f s\n",time_elapsed);

    sleep(1);
    printf("\n\n************** Linear system solving ***************\n");

    int B[n]; // solutions
    int Z[n]; // LZ = B
    int X[n]; // UX = Z
    generateSolution(B,n,p);
    InitVector(Z,X,n);
    

    // Computation and time measurements
    start = clock();
    LinearSystem(L,U,Z,B,X,n,p);
    end = clock();
    time_elapsed = ((double)end - start) / CLOCKS_PER_SEC;

    // Show results
    DisplaySolutions(B,X,Z,n);

    // Check if LZ = B and UX = Z
    correctionLS(L,U,Z,B,X,n,p);
    printf("Time duration : %f s\n",time_elapsed);

    sleep(1);
    printf("\n\n************** Inverse matrix **********************\n");

    // The inverse matrix 
    int s[n+1][n+1]; 

    // Computation and time measurements
    start = clock();
    InverseMatrix(A,s,n,p);
    end = clock();

    // Show results
    DisplayInv(s,n);
    time_elapsed = ((double)end - start) / CLOCKS_PER_SEC;


    // Check if A * InvA = I
    correctionInv(A,s,n,p);  

    sleep(1); */
    printf("\n\n************** Strassen Algorithm **********************\n");  
    int A5[][4] = { { 3, 4, 7, 2},
                    { 1, 9, 7, 6},
                    { 2, 4, 1, 9},
                    { 10, 4, 12, 11}};


    //displayStrassen("test",t,2,2);

    matrix As = loadMatrix(A5,n);
    displayStrassen("As",As,n);
    matrix inv = InvStrassen(As,n,p);
    //displayStrassen("A",As,n);
    displayStrassen("Inv:",inv,n); 


    freeMatrix(As);
    freeMatrix(inv);
    //freeMatrix(A);
    freeMatrix(L);
    freeMatrix(U);
    freeMatrix(P);


    return 0;
}