#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "LU.h"
#include "Tools.h"




/**
 * initialize the square matrix
*/
void null_matrix(matrix m, DIM)
{
    for_(i,n)
        for_(j,n)
            m[i][j] = 0;
}



/**
 * allocate memory for the matrix and initialize it
 * @return m : matrix n*n initialized
*/
matrix init_matrix(DIM)
{
    matrix m = malloc(sizeof(int*) * n);
    m[0] = malloc(sizeof(int) * n * n);

    for_(i,n)
        m[i] = m[0] + n * i;
        //m[i] = malloc(sizeof(int) * n); //each rows
    null_matrix(m,n);

    return m;
}



/**
 * load the values given in parameter into a new square matrix
 * @return m : the new matrix
*/
matrix loadMatrix(void *s, DIM)
{
    matrix m = init_matrix(n);
    for_(i,n)
        for_(j,n)
            m[i][j] = ((int(*)[n]) s)[i][j]; // s is a pointer to an array of n doubles
    return m;
}

/**
 * Generate a n*n matrix of random integers number
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
 * free the matrix and all elements within it
*/
void freeMatrix(matrix m)
{   
    free(m[0]);
    free(m);
}


/**
 * Display the matrix m to stdout
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
 * Compute the matrix product between a and b
 * @return : matrix c , result of the product
*/
matrix mat_mult(matrix a, matrix b, DIM, prime)
{
    matrix c = init_matrix(n);

    for_(i,n)
        for_(j,n)
            for_(k,n)
                c[i][j] += (long) ((a[i][k] * b[k][j]) % p) ;
    return c;
}


/**
 * Compute the permutation (swapping) of matrices using the pivoting  process 
*/
void pivot(matrix a, matrix p, DIM)
{
    for_(i,n)
        for_(j,n)
            p[i][j] = (i == j); // only the diagonal coefficients
    
    for_(i,n)
    {
        int max_j = i;
        foreach(j,i,n)
            if(abs(a[j][i]) > abs(a[max_j][i]))
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
 * Compute the LU decomposition using the pivot function
*/
void LU(matrix A, matrix L, matrix U, matrix P, DIM,prime)
{
    null_matrix(L,n);
    null_matrix(U,n);
    pivot(A,P,n);


    matrix Abis = mat_mult(P,A,n,p);

    int tmp = 0;
    for_(i,n)
        L[i][i] = 1;
    
    for_(i,n)
        for_(j,n)
        {
            int s;
            if(j <= i)
            {
                s = 0;
                foreach(k, 0, j)
                    s+= (long) ((L[j][k] * U[k][i]) % p);
                tmp = sub(A[j][i],s,p);
                U[j][i] = sub(Abis[j][i],s,p);
                if(U[j][i] < 0) U[j][i] += p;
                if(U[j][i] > p) U[j][i] -= p;
            }
            if(j > i)
            {
                s = 0;
                foreach(k, 0, i)
                    s+= (long) ((L[j][k] * U[k][i]) % p);
                tmp = sub(Abis[j][i],s,p);
                L[j][i] = (long) (tmp * inv(U[i][i],p)) % p;// inversion modulaire
                if(L[j][i] < 0) L[j][i] += p;
                if(L[j][i] > p) L[j][i] -= p;
            }
        }

    freeMatrix(Abis);

}


/**
 * Check if the LU decomposition is correct
*/
void correction(matrix A, matrix L, matrix U, matrix P, DIM, prime)
{
    matrix PA = init_matrix(n);
    matrix LU = init_matrix(n);
    bool ok = false;
    for_(i,n)
        for_(j,n)
            for_(k,n){
                PA[i][j] += (long) ((P[i][k] * A[k][j]) % p);
                LU[i][j] += (long) ((L[i][k] * U[k][j]) % p);

                if(LU[i][j] < 0) LU[i][j] += p;
                if(LU[i][j] > p) LU[i][j] -= p;
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


void RunLU(matrix A, matrix P, matrix L, matrix U, DIM, prime)
{
    printf("************** LU decomposition ***************\n");

    // Computation and time measurements
    double start,end,time_elapsed;
    start = clock();
    LU(A,L,U,P,n,p);
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
}