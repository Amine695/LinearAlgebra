/**
 * @file LU.c
 * @author Amine Berbagui
 * @brief File that compute the LU decomposition
 * @date 2022-05-05
 * @copyright Copyright (c) 2022
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "LU.h"
#include "Tools.h"




/**
 *  @brief initialize the square matrix to 0
*/
void null_matrix(matrix m, DIM)
{
    for_(i,n)
        for_(j,n)
            m[i][j] = 0;
}



/**
 * @brief allocate memory for the matrix and initialize it
 * @return m : matrix n*n initialized
*/
matrix init_matrix(DIM)
{
    matrix m = malloc(sizeof(int*) * n);
    m[0] = malloc(sizeof(int) * n * n);

    for_(i,n)
        m[i] = m[0] + n * i;
    null_matrix(m,n);

    return m;
}

/**
 * @brief Copy a matrix to a new one
*/
matrix Copy(matrix m,DIM)
{
    matrix C = init_matrix(n);

    for_(i,n)
        for_(j,n)
            C[i][j] = m[i][j];
    return C;
}

/**
 *  @brief load a 2D array into a new square matrix
 * @return m : the new matrix
*/
matrix loadMatrix(void *s, DIM)
{
    matrix m = init_matrix(n);
    for_(i,n)
        for_(j,n)
            m[i][j] = ((int(*)[n]) s)[i][j]; // s is a pointer to an array of n int
    return m;
}

/**
 *  @brief Generate a n*n matrix of random integers number
 * We can adjust the random range for performances mesurements 
 * @return A : the random matrix
*/
matrix generate(DIM,prime)
{
    matrix A = init_matrix(n);
    for_(i,n) 
    	for_(j,n) 
    		A[i][j] = rand() % p;
    
    return A;
}

/**
 *  @brief free the matrix and all elements within it
*/
void freeMatrix(matrix m)
{   
    free(m[0]);
    free(m);
}


/**
 *  @brief Display the matrix m to stdout 
*/
void displayMatrix(char title[],matrix m,DIM)
{
    printf("\n%s \n",title);
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
 * @brief Function that compute the matrices product
*/
matrix matmul(matrix a, matrix b, DIM , prime)
{
    matrix res = init_matrix(n);
    int tmp;
    for ( int i = 0; i < n; i++ )
    {
        for ( int j = 0; j < n; j++ )
        {
            tmp = 0;    
            for (int k = 0; k < n; k++) {
                tmp = (long) ( a[i][k] * b[k][j]) % p;
                if(tmp < 0) tmp += p;
                res[i][j] = (res[i][j] + tmp) % p; 
                
            }
        }

    }
    return res;
 
}

/**
 * @brief Function that create an identity matrix
 * @return res
*/
matrix Identity(DIM)
{
    matrix identity = init_matrix(n);
    for_(i,n)
        for_(j,n)
        {
            if(i == j)
                identity[i][j] = 1;
            else
                identity[i][j] = 0;
        }
    return identity;
}



/**
 *  @brief Compute the permutation (swapping) of rows using the pivoting  process 
*/
void pivot(matrix a, matrix p, DIM)
{
    for_(i,n)
        for_(j,n)
            p[i][j] = (i == j); // p is an identity matrix
            
    
    for_(i,n)
    {
        int max_j = i;
        foreach(j,i,n)
            if(a[j][i] > a[max_j][i])
                max_j = j;
        
        if(max_j != i)
            for_(k,n)
            {
                int tmp;
                tmp = p[i][k];
                p[i][k] = p[max_j][k];
                p[max_j][k] = tmp;
            }
    }
}

/**
 *  @brief Compute the LU decomposition of P*A 
*/
void LU(matrix A, matrix L, matrix U,matrix P, DIM,prime)
{
    
    int tmp = 0,tmp2 = 0;
    for_(i,n)
        L[i][i] = 1; // identity matrix
    
    for_(i,n)
        for_(j,n)
        {
            int s;
            // on calcul U
            if(j <= i)
            {
                s = 0;
                foreach(k, 0, j)
                    s+= (long) ((L[j][k] * U[k][i]));
                tmp2 = s % p;
                tmp = sub(A[j][i],tmp2,p);  
                U[j][i] = sub(A[j][i],tmp2,p);

                // readjustment in case we exceed p or it goes below 0
                if(U[j][i] < 0) U[j][i] += p;
                if(U[j][i] > p) U[j][i] -= p;
            }   
                
            // On calcul L
            if(j >= i)
            {
                s = 0;
                foreach(k, 0, i)
                    s+= (long) ((L[j][k] * U[k][i]));
                tmp2 = s % p;
                tmp = sub(A[j][i],tmp2,p);
                L[j][i] = (long) (tmp * inv(U[i][i],p)) % p;// modular inverse

                // readjustment in case we exceed p or it goes below 0
                if(L[j][i] < 0) L[j][i] += p;
                if(L[j][i] > p) L[j][i] -= p;
            }
        }


}


/**
 *  @brief Check if the LU decomposition is correct by making the product of L and U
*/
void correction(matrix A, matrix L, matrix U, matrix P, DIM, prime)
{
    matrix PA = init_matrix(n);
    matrix LU = init_matrix(n);
    bool ok = false;

    // Compute P.A and L.U
    for_(i,n)
        for_(j,n)
            for_(k,n){
                PA[i][j] += (long) ((P[i][k] * A[k][j]) % p);
                LU[i][j] += (long) ((L[i][k] * U[k][j]) % p);

                if(LU[i][j] < 0) LU[i][j] += p; //readjustment in case L[i][j] < 0
                if(LU[i][j] > p) LU[i][j] -= p; //readjustment in case we exceed p
            }

    for_(i,n)
        for_(j,n)
        {
            if(PA[i][j] == LU[i][j])
                ok = true; 
            else
                ok = false;
        }

    if(ok)
        printf("\nPA == LU : OK !\n"); 
    else
        printf("\nPA != LU : KO !\n"); 

    freeMatrix(PA);
    freeMatrix(LU);

}


/**
 *  @brief Run the LU decomposition to make the main function cleaner
 *  It uses the pivot function to compute the permutation matrix
*/
void RunLU(matrix A, matrix P, matrix L, matrix U, DIM, prime)
{
    printf("************** LU decomposition ***************\n");
    // calcul the permutation matrix
    pivot(A,P,n);
    // compute P*A 
    matrix PA = matmul(P,A,n,p);
    // Computation and time measurements
    double start,end,time_elapsed;
    start = clock();
    LU(PA,L,U,P,n,p); // Here PA = P*A so it computes PA = LU
    end = clock();
    time_elapsed = ((double)end - start) / CLOCKS_PER_SEC;
    
    // Show results
    displayMatrix("A matrix",A,n);
    displayMatrix("L matrix",L,n);
    displayMatrix("U matrix",U,n);
    displayMatrix("P matrix",P,n);

     // Check if PA == LU
    correction(A,L,U,P,n,p);
    printf("Time duration with n = %d : %f s\n",n,time_elapsed);

    freeMatrix(PA);
}