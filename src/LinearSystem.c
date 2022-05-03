#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include "LinearSystem.h"
#include "Tools.h"

/** 
 * @file file that solve a linear system using the LU decomposition AX = B. 
*/


/** @brief
 * Step 1: Solve LZ = B for Z 
 * Step 2: Solve UX = Z for X
 * Step 3: Check 
*/

/**
 *  @brief Compute the B matrix and Z matrix using forward and backward permutations.
*/
void LinearSystem(matrix L, matrix U, int Z[], int B[], int X[],DIM,prime)
{
    
    int sum;
    int tmp = 0;
    
    // Finding Z forward substitution
    foreach(i,0,n) 
    {
        sum = 0;
        foreach(j,0,i){
            sum += (long)(L[i][j]*Z[j]) % p;
        }
        tmp = sub(B[i],sum,p);
        Z[i]= (long) ( tmp * inv(L[i][i],p)) % p; 
        if(Z[i] < 0) Z[i] += p;
    }
    
    // Finding X  backward substitution
    forneg(i,n-1,0)
    {
        sum = 0;
        forneg(k,n-1,i){
            sum += (long)((U[i][k]*X[k]) % p);
        }
        tmp = sub(Z[i],sum,p);
        X[i]= (long) ( tmp * inv(U[i][i],p)) % p;

    }

}

/**
 * @brief Function that generates a vector of size n with random values between 0 and p-1
*/
 void generateSolution(int sol[],DIM,prime)
{
    for_(i,n)
        sol[i] = rand() % p;
    
}

/**
 * @brief Function that initialize the Z and X vectors to 0
*/
void InitVector(int Z[], int X[], DIM)
{
    for_(i,n)
    {
        Z[i] = 0;
        X[i] = 0;
    }
}

/**
 *  @brief Display the X and Z vectors, with X the solutions of the linear system
*/
void DisplaySolutions(int B[],int X[],int Z[],DIM)
{
    printf("1 :Solution vector \n");
    for_(i,n)
        printf("B%d = %d\n",i,B[i]);

    printf("\n2 : Solve LZ = B for Z\n\n");
    foreach(i,0,n)
        printf("Z%d = %d\n",i+1,Z[i]);
    

    printf("\n3 : Solve UX = Z for X\n\n");
    foreach(i,0,n)
        printf("X%d = %d\n",i+1,X[i]);
    
    
}

/**
 *  @brief Function verifying that we have the equality PA = LU
*/
void correctionLS(matrix L, matrix U, int Z[], int B[], int X[],DIM,prime)
{
    bool ok = true;
    int tmp[n];
    int tmp2[n];

    for_(i,n)
    {
        tmp[i] = 0;
        tmp2[i] = 0;
    }

    for_(i,n){
        for_(k,n){
            tmp[i]  += (long) ((L[i][k] * Z[k]) % p); // LZ = B
            tmp2[i] += (long) ((U[i][k] * X[k]) % p); // UX = Z
            if(tmp[i] > p) tmp[i] -= p;
            if(tmp2[i] > p) tmp2[i] -= p;

        }

    }
    printf("\nChecking :\n");
    printf("\nLZ\n");
    for_(i,n)
        printf("%3.1d ",tmp[i]);
   
    printf("\nUX\n");
    for_(i,n)    
        printf("%3.1d ",tmp2[i]);

    for_(i,n)
        if((tmp[i] != B[i]) && (tmp2[i] != Z[i]))
            ok = false;
        
    if(ok)
        printf("\n\nLZ == B : OK !\nUX = Z : OK!\n");
    else
        printf("\n\nLZ == B : KO !\nUX = Z : KO!\n");


}

/**
 *  @brief Function that run the linear system solving from main
*/
void RunLS(matrix L, matrix U, int B[], int Z[], int X[], DIM, prime)
{
    printf("\n\n************** Linear system solving ***************\n");
    double start,end,time_elapsed;

    // creation of a random solution vector modulo p
    generateSolution(B,n,p);
    // Initialization of Z and X with 0
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
    printf("Time duration with n = %d : %f s\n",n,time_elapsed);
}
