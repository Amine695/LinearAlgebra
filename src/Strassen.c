/**
 * @file Strassen.c
 * @author Amine Berbagui
 * @brief File that compute all Strassen's algorithms
 * @date 2022-05-05
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "Strassen.h"
#include "MatrixInv.h"



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
 * @brief Add two matrices
*/
matrix AddMatrix(matrix A, matrix B,DIM,prime)
{
    matrix res = init_matrix(n);
    int tmp = 0;
    for_(i,n)
        for_(j,n){
            tmp = (A[i][j] + B[i][j]) % p;
            if(tmp > p) tmp -= p;
            res[i][j] =  tmp;

        }
    
    return res;
}

/**
 * @brief Change the sign of the matrix
 * @return res : the matrix of opposite sign
*/
matrix ChangeSign(matrix A,DIM,prime)
{
    matrix res = init_matrix(n);
    for_(i,n)
        for_(j,n)
            if(A[i][j])
                res[i][j] = (-1 * A[i][j]) + p;
    return res;

}


/**
 * @brief Decompose the matrix in sub-matrix
 * @return res the sub-matrix
*/
matrix SubMatrix(matrix A,DIM, int d)
{
    // exclude error
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
            // Compute A11
            if(d == 1 && i < k && j < k)
            {
                res[i][j] = A[i][j];
                cpt++;

            }             
            // Compute A12
            if(d == 2 && i < k && j >= k)
            {
                res[i][j-k] = A[i][j];
                cpt++;

            }  
            // Compute A21
            if(d == 3 && i >= k && j < k)
            {
                res[i-k][j] = A[i][j];
                cpt++;

            }  
            // Compute A22
            if(d == 4 && i >= k && j >= k) 
            {
                res[i-k][j-k] = A[i][j];
                cpt++;
            }  

        }

    return res;
    
}

/**
 * @brief Build the matrix A from its sub-matrices
 * @return res the complete matrix
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
        
    return A;
    

}

/**
 * @brief print the matrix to stdout
*/
void displayStrassen(char title[],matrix m, DIM)
{
    
    printf("%s\n",title);
    printf("\n");
    for_(i,n)
    {   
        for_(j,n)
        {
            printf("%7.1d     ",m[i][j]);
        }
        printf("\n");
    }
    
}

/**
 * @brief Compute the Strassen's inversion algorithm
 * @return res the Inverse of A
*/
matrix InvStrassen(matrix A,DIM,prime)
{
    
    // if n == 1, we return the modular inverse of the only coefficient
    if(n == 1)
    {
        matrix a = init_matrix(n);
        int value = A[0][0];
        a[0][0] =  inv(value,p); 
        return a;
    }
    
    // split the size by 2
    int k = n/2;

    //split A into four sub-matrices
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
    matrix X22 = ChangeSign( R6,k ,p);
    
    // Assemble the complete matrix
    matrix res = assembly(X11,X12,X21,X22,n);

    // free all the matrices
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
    freeMatrix(X22);
    freeMatrix(A11);
    freeMatrix(A12);
    freeMatrix(A21);
    freeMatrix(A22);

    return res;  
   
}

/**
 * @brief Compute the Strassen's multiplication algorithm
 * @return res the result of the product of A and B
*/
matrix StrassenMult(matrix A, matrix B, DIM,prime)
{
    // minimum size below which the naive product is applied
    if(n <= 32)
    {
        matrix res = matmul(A,B,n,p);
        return res;
    }
    else
    {
        matrix res = init_matrix(n);        
        
        //Divide the matrix A and B into 4 sub-matrices
        matrix a11 = SubMatrix(A, n, 1);
        matrix a12 = SubMatrix(A, n, 2);
        matrix a21 = SubMatrix(A, n, 3);
        matrix a22 = SubMatrix(A, n, 4);	
        matrix b11 = SubMatrix(B, n, 1);
        matrix b12 = SubMatrix(B, n, 2);
        matrix b21 = SubMatrix(B, n, 3);
        matrix b22 = SubMatrix(B, n, 4);


        // Algorithm steps
        matrix m1 = StrassenMult(AddMatrix(a11,a22,n/2,p),AddMatrix(b11,b22,n/2,p),n/2,p);
        matrix m2 = StrassenMult(AddMatrix(a21,a22,n/2,p),b11,n/2,p);
        matrix m3 = StrassenMult(a11,SubstractMatrix(b12,b22,n/2,p) ,n/2,p);
        matrix m4 = StrassenMult(a22,SubstractMatrix(b21,b11,n/2,p) ,n/2,p);
        matrix m5 = StrassenMult(AddMatrix(a11,a12,n/2,p),b22,n/2,p);
        matrix m6 = StrassenMult(SubstractMatrix(a21,a11,n/2,p),AddMatrix(b11,b12,n/2,p),n/2,p);
        matrix m7 = StrassenMult(SubstractMatrix(a12,a22,n/2,p),AddMatrix(b21,b22,n/2,p),n/2,p);

        // build the sub-matrices 
        matrix c11 = AddMatrix(SubstractMatrix(AddMatrix(m1,m4,n/2,p),m5,n/2,p),m7,n/2,p);
        matrix c12 = AddMatrix(m3,m5,n/2,p);
        matrix c21 = AddMatrix(m2,m4,n/2,p);
        matrix c22 = AddMatrix(SubstractMatrix(AddMatrix(m1,m3,n/2,p),m2,n/2,p),m6,n/2,p);

        // Assemble the complete matrix
        res = assembly(c11,c12,c21,c22,n);

        // memory free
        freeMatrix(a11);
        freeMatrix(a12);
        freeMatrix(a21);
        freeMatrix(a22);
        freeMatrix(b11);
        freeMatrix(b12);
        freeMatrix(b21);
        freeMatrix(b22);
        freeMatrix(m1);
        freeMatrix(m2);
        freeMatrix(m3);
        freeMatrix(m4);
        freeMatrix(m5);
        freeMatrix(m6);
        freeMatrix(m7);
        freeMatrix(c11);
        freeMatrix(c12);
        freeMatrix(c21);
        freeMatrix(c22);
        
        return res;
    }

}

/**
 * @brief Compute the Strassen's inversion algorithm with the Strassen's multiplication algorithm
 * @return res the Inverse of A
*/
matrix InvStrassenMultStrassen(matrix A,matrix I,matrix L,matrix U,matrix P,matrix Z,DIM,prime)
{

    // if n == 1, we return the modular inverse of the only coefficient
    /* if(n == 1)
    {
        matrix a = init_matrix(n);
        int value = A[0][0];
        a[0][0] =  inv(value,p); 
        return a;
    } */

    // apply the naive inverse instead of the modular inverse
    if(n <= 64)
    {
        matrix res;
        res = InverseMatrix(I,L,U,P,Z,n,p);
        return res;
    }
  
    // split the size by 2
    int k = n/2;

    //split A into four sub-matrices
    matrix A11 = SubMatrix(A,n,1);  
    matrix A12 = SubMatrix(A,n,2);  
    matrix A21 = SubMatrix(A,n,3);  
    matrix A22 = SubMatrix(A,n,4);   
    
    // Strassen steps
    matrix R1  = InvStrassenMultStrassen( A11,I,L,U,P,Z,k ,p);
    matrix R2  = StrassenMult( A21, R1, k, p);   
    matrix R3  = StrassenMult( R1, A12 ,k ,p);   
    matrix R4  = StrassenMult( A21, R3 ,k ,p);    
    matrix R5  = SubstractMatrix( R4, A22,k, p);
    matrix R6  = InvStrassenMultStrassen( R5,I,L,U,P,Z, k, p);
    matrix X12 = StrassenMult( R3, R6 , k, p);  
    matrix X21 = StrassenMult( R6, R2 , k, p);   
    matrix R7  = StrassenMult( R3, X21, k, p);   
    matrix X11 = SubstractMatrix( R1, R7, k, p);
    matrix X22 = ChangeSign( R6,k ,p);
    
    // Assemble the complete matrix
    matrix res = assembly(X11,X12,X21,X22,n);

    // free all the matrices
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

    return res;  
   
}

/**
 * @brief Run the Strassen's inverse algorithm, not for benchmarks
*/
void RunStrassenInv(matrix A,DIM,prime)
{
    printf("\n\n************** Strassen's Inversion Algorithm **********************\n");  
    double start,end,time_elapsed;
    displayStrassen("A",A,n); 
    start = clock();
    matrix inv = InvStrassen(A,n,p);
    end = clock();
    time_elapsed = ((double)end - start) / CLOCKS_PER_SEC;
    displayStrassen("Inverse matrix:",inv,n); 

    printf("Strassen's inverse time duration with n = %d : %f s\n",n,time_elapsed);
    freeMatrix(inv);
  
}

/**
 * @brief Run the Strassen's multiplication algorithm, not for benchmarks
*/
void RunStrassenMult(matrix A,matrix B,DIM,prime)
{
    printf("\n\n************** Strassen's multiplication Algorithm **********************\n");  
    displayMatrix("A matrix\n",A,n);
    displayMatrix("\nB matrix\n",B,n);
    double start,end,time_elapsed;
    start = clock();
    matrix res = StrassenMult(A,B,n,p);
    end = clock();
    time_elapsed = ((double)end - start) / CLOCKS_PER_SEC;

    displayStrassen("\nMatrix product :",res,n);
    printf("\nTime duration with n = %d : %f s\n",n,time_elapsed);

    freeMatrix(res);

}


