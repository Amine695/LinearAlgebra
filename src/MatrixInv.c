#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "MatrixInv.h"
#include "Tools.h"

/**
 * Function similar to LU(), we use x, the pivot to make arithmetics operations on A
*/
void InvLU(matrix A,DIM,prime)
{
    int x;
    for_(k,n-1)
    {
        foreach(j,k+1,n)
        {
            x= (long)  ((A[j][k]) * modInverse(A[k][k],p)) % p;
            foreach(i,k,n)
            {  
                A[j][i]=  (long) (((A[j][i]-x) * A[k][i]) % p);
            }
            A[j][k]=x;
        }
    }
}



/**
 * Compute the inverse of the A  matrix using the InvLU function
*/
void InverseMatrix(matrix A,int s[][SIZE],DIM, prime)
{
    int x,y[n],d[n];
    
    InvLU(A,n,p);

    for_(i,n)
        for_(j,n)
            s[i][j] = 0;

    for_(m,n)
    {
        d[0] = 0;d[1] = 0;d[2] = 0;
        d[m] = 1;

        for_(i,n)
        {
            x = 0;
            foreach(j,0,i)
                x += (long) ((A[i][j] * y[j]) % p);
            y[i] = d[i]-x;
        }

        forneg(i,n-1,0)
        {
            x = 0;
            foreach(j,i+1,n)
                x += (long) ((A[i][j] * s[j][m]) % p);
            s[i][m]= (long) ((y[i]-x) * modInverse(A[i][i],p)) % p;
            
        }
    }
    
}


void DisplayInv(int s[][SIZE], DIM)
{
    printf("\nInverse matrix : \n");
    
    for_(i,n)
    {   
        for_(j,n)
        {
            printf("%7.1d   ",s[i][j]);
        }
        printf("\n");
    }
}


/**
 * Check if A.A^-1 = Identity matrix
*/
void correctionInv(matrix A, int s[][SIZE], DIM,prime)
{
    bool ok = false;
    int tmp[n+1][n+1];
    for_(i,n)
        for_(j,n)
            tmp[i][j] = 0;
    for_(i,n)
        for_(j,n)
            for_(k,n)
                tmp[i][j] += (long) ((A[i][k] * s[k][j]) % p);
    
    printf("\n");

    for_(i,n)
        for_(j,n){
            tmp[i][j] = fabs(tmp[i][j]);
            if((i == j && tmp[i][j] == 1) || (i != j && tmp[i][j] == 0))
                ok = true;
        }
        

    if(ok)
        printf("\nA.InvA = Identity : OK !\n");
    else
    {
        printf("\nA.InvA != Identity : KO !\n");
        for_(i,n){
            for_(j,n)
                printf("%d  ",tmp[i][j]);
            printf("\n");
        }
    }


}