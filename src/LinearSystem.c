#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include "LinearSystem.h"
#include "Tools.h"

/** 
 * Solve a linear system using the LU decomposition AX = B. 
 * Step 1: retreive the L and U matrices computed before
 * Step 2: Solve LZ = B for Z 
 * Step 3: Solve UX = Z for X
*/

/**
 * Compute the B matrix and Z matrix using forward permutations.
*/
void LinearSystem(matrix L, matrix U, int Z[], int B[], int X[],DIM,prime)
{
    // Finding Z 
    double sum;
    foreach(i,0,n)
    {
        sum = 0;
        foreach(j,0,n){
            sum += (long)((L[i][j]*Z[j]) % p);
        }
        Z[i]= (long) ((B[i]-sum) * modInverse(L[i][i],p)) % p; // division modulaire
    }

    // Finding X 
    forneg(i,n-1,0)
    {
        sum = 0;
        forneg(k,n-1,i)
            sum += (long)((U[i][k]*X[k]));
        X[i]= (long) ((Z[i]-sum) * modInverse(U[i][i],p)) % p;
    }

}

/**
 * @brief Function who generates a vector of size n with random values between 0 and p-1
 * @return sol : solution vector
*/
 void generateSolution(int sol[],DIM,prime)
{
    for_(i,n)
        sol[i] = rand() % p;
    
}


void InitVector(int Z[], int X[], DIM)
{
    for_(i,n)
    {
        Z[i] = 0;
        X[i] = 0;
    }
}

/**
 * Display the X and Z vectors, with X the solutions of the linear system
*/
void DisplaySolutions(int B[],int X[],int Z[],DIM)
{

    for_(i,n)
        printf("B%d = %d\n",i,B[i]);

    printf("\n1 : Solve LZ = B for Z\n\n");
    foreach(i,0,n)
        printf("Z%d = %d\n",i+1,Z[i]);
    

    printf("\n2 : Solve UX = Z for X\n\n");
    foreach(i,0,n)
        printf("X%d = %d\n",i+1,X[i]);
    
    
}

/**
 * Function verifying that we have the equality PA = LU
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

    for_(i,n)
        for_(k,n){
            tmp[i]  += (long) ((L[i][k] * Z[k]) % p); // LZ = B
            tmp2[i] += (long) ((U[i][k] * X[k]) % p); // UX = Z
        }

    for_(i,n)
        if((tmp[i] != B[i]) && (tmp2[i] != Z[i]))
            ok = false;
        
    if(ok)
        printf("\nLZ == B : OK !\nUX = Z : OK!\n");
    else
        printf("LZ == B : KO !\nUX = Z : KO!\n");


}