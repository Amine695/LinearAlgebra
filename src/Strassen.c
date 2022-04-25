#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "Strassen.h"


/**
 * @file Implementation of the Strassen's inversion algorithm
*/




matrix matmul(matrix a, matrix b, DIM , prime)
{
    matrix res = init_matrix(n);
    int tmp = 0;
    for ( int i = 0; i < n; i++ )
    {
        for ( int j = 0; j < n; j++ )
        {
            
            for (int k = 0; k < n; k++) {
                tmp = ((long) a[i][k] * b[k][j]) % p;
                if(tmp < 0) tmp += p;
                res[i][j] = (res[i][j] + tmp) % p; 
                //printf("a = %d ,b = %d , res = :%d\n",a[i][k],b[k][j],(a[i][k] * b[k][j]) % p);
            }
            if(res[i][j] < 0)
                printf("matmul\n");
        }

    }
    return res;
 
}


/**
 * @brief Substract two matrices 
*/
matrix SubstractMatrix(matrix A, matrix B,DIM,prime)
{
    matrix res = init_matrix(n);
    int tmp = 0;
    for_(i,n)
        for_(j,n){
            tmp = (A[i][j] - B[i][j]) % p;
            if(tmp < 0) tmp += p;
            res[i][j] =  tmp;

        }
    
    return res;

}

/**
 * @brief Change the sign of the matrix
*/
void ChangeSign(matrix A,DIM,prime)
{
    for_(i,n)
        for_(j,n)
            if(A[i][j])
                A[i][j] = (-1 * A[i][j]) + p;
            

}




/**
 * @brief Decompose the matrix in sub-matrix
*/
matrix SubMatrix(matrix A,DIM, int d)
{
    
    if( d < 1 || d > 4)
    {
        printf("Error value !\n");
        exit(EXIT_FAILURE);
    }
    int cpt = 0;
    int k = n/2;
    matrix res = init_matrix(k);

    for_(i,n)
        for_(j,n)
        {
            if(d == 1 && i < k && j < k)
            {
                res[i][j] = A[i][j];
                cpt++;

            }                
            if(d == 2 && i < k && j >= k)
            {
                res[i][j-k] = A[i][j];
                cpt++;

            }  
        
            if(d == 3 && i >= k && j < k)
            {
                res[i-k][j] = A[i][j];
                cpt++;

            }  
        
            if(d == 4 && i >= k && j >= k) 
            {
                res[i-k][j-k] = A[i][j];
                cpt++;
            }  

            if(cpt >= 2 * k)
                return res;
        }

    return res;
    
}

/**
 * @brief Build the matrix A from its sub-matrices
*/
matrix assembly(matrix A11, matrix A12, matrix A21, matrix A22,DIM)
{
    int  k = n/2;
    matrix A = init_matrix(n);
   

    for(int i = 0; i < k; i++)
        for(int j = 0; j < k; j++)
            A[i][j] = A11[i][j];

    for(int i = 0; i < k; i++)
        for(int j = k; j < n; j++)
            A[i][j] = A12[i][j-k];

    for(int i = k; i < n; i++)
        for(int j = 0; j < k; j++)
            A[i][j] = A21[i-k][j];

    for(int i = k; i < n; i++)
        for(int j = k; j < n; j++)
            A[i][j] = A22[i-k][j-k];
        
    //displayStrassen("A12",A12,k,n);
    return A;
    

}

void displayStrassen(char title[],matrix m, DIM)
{
    
    printf("%s\n",title);
    printf("\n");
    for_(i,n)
    {   
        for_(j,n)
        {
            printf("%4.1d     ",m[i][j]);
        }
        printf("\n");
    }
    
}

matrix InvStrassen(matrix A,DIM,prime)
{
    //printf("n = %d\n",n);
    if(n == 1)
    {
        matrix a = init_matrix(n);
        int value = A[0][0];
        a[0][0] =  modInverse(value,p); 
        return a;
    }
  
    int k = n/2;

    //split into four
    matrix A11 = SubMatrix(A,n,1);  
    matrix A12 = SubMatrix(A,n,2);  
    matrix A21 = SubMatrix(A,n,3);  
    matrix A22 = SubMatrix(A,n,4);   
    
    // Strassen steps
    matrix R1  = InvStrassen( A11,k ,p);
    matrix R2  = matmul( A21, R1, k, p);
    matrix R3  = matmul( R1, A12 ,k ,p);
    matrix R4  = matmul( A21, R3 ,k ,p); 
    matrix R5  = SubstractMatrix( R4, A22,k, p);
    matrix R6  = InvStrassen( R5, k, p);
    matrix X12 = matmul( R3, R6 , k, p);
    matrix X21 = matmul( R6, R2 , k, p);
    matrix R7  = matmul( R3, X21, k, p);
    matrix X11 = SubstractMatrix( R1, R7, k, p);
    ChangeSign( R6,k ,p);
    
    //displayStrassen("A22",R4,k,n);
    matrix tmp = assembly(X11,X12,X21,R6,n);

    freeMatrix(R1);
    freeMatrix(R2);
    freeMatrix(R3);
    freeMatrix(R4);
    freeMatrix(R5);
    freeMatrix(R6);
    freeMatrix(R7);
    freeMatrix(X12);
    freeMatrix(X21);
    freeMatrix(X11);
    freeMatrix(A11);
    freeMatrix(A12);
    freeMatrix(A21);
    freeMatrix(A22);

    return tmp;  
   
}



