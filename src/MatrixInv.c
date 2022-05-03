#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "MatrixInv.h"
#include "Tools.h"

/**
 * @brief Function that compute the Z matrix for each row of the matrix
 * @return result 
*/
int computeZ(matrix L, matrix I,matrix Z,int col, int row,DIM,prime) 
{
    int sum = 0;
    int tmp = 0;
    for(int i = 0; i < n; i++) 
        if(i != row) 
            tmp += L[row][i] * Z[i][col]; // product row of L by column of Z
    sum = tmp % p; // modulo p

    int result = sub(I[row][col],sum,p); // substract sum to the identity matrix
    result = (long) (result * inv(L[row][row],p)) % p; // compute the modular inverse

    return result;
}

/**
 * @brief Function that compute the inverse matrix 
 * @return result 
*/
int computeInverse(matrix inverse, matrix U,matrix Z, int col, int row, DIM,prime) 
{
    int sum = 0;
    int tmp = 0;
    for(int i = 0; i < n; i++) 
        if(i != row)
            tmp += U[row][i] * inverse[i][col]; // product row of U by column of inverse matrix
        
    sum = (tmp % p); // result mod p
    int result = sub(Z[row][col],sum,p); // substract sum to the Z vector
    result = (long) (result *  inv(U[row][row],p)) % p; // compute the modular inverse

    return result;
}

/**
 * @brief Function in which we call the two functions above to build the Z vector and the inverse matrix.
 * @return res : the inverse matrix. 
*/
matrix InverseMatrix(matrix I,matrix L, matrix U,matrix P, matrix Z, DIM,prime) 
{
    matrix inverse = init_matrix(n);
    // compute z
    for(int col = 0; col < n; col++) {
        for(int row = 0; row < n; row++) {
            Z[row][col] = computeZ(L,I,Z,col, row,n,p);
        }
    }

    // compute inverse
    for(int col = 0; col < n; col++) {
        for(int row = n - 1; row >= 0; row--) {
            inverse[row][col] = computeInverse(inverse,U,Z,col,row,n,p);
        }
    }
   
    return inverse;
}




/**
 * @brief Check if A.A^-1 = Identity matrix
*/
void correctionInv(matrix A, matrix s, DIM,prime)
{
    bool ok = false;

    // matrix to check that A * InvA = identity
    matrix V = init_matrix(n);
    
    // product of A . s with s the inverse matrix
    for_(i,n)
        for_(j,n){
            for_(k,n){
                V[i][j] += (long) ((A[i][k] * s[k][j]));
                
            }
            V[i][j] = V[i][j] % p;
        }
    printf("\n");
    displayMatrix("A x Inv A",V,n);
       
    // check that V is the identity matrix
    for_(i,n)
        for_(j,n)
        {
            if((i == j && V[i][j] == 1) || (i != j && V[i][j] == 0))
                ok = true;
        }

    
    if(ok)
    {
        printf("\nA.InvA = Identity : OK !\n");
    }
    else
        printf("\nA.InvA != Identity : KO !\n");
    freeMatrix(V);    
}


/**
 * @brief Function that launch the Naive Inverse called in the main
*/
void RunNaiveInverse(matrix A,matrix L,matrix U,matrix P,DIM,prime)
{
    printf("\n\n************** Naive Inverse matrix **********************\n");

    // Inv : the inverse matrix 
    matrix Inv ;
    matrix Z = init_matrix(n);
    matrix identity = Identity(n);

    // Computation and time measurements
    double start,end,time_elapsed;
    start = clock();
    Inv = InverseMatrix(identity,L,U,P,Z,n,p);
    end = clock();
    
    // Show results
    displayMatrix("\nInverse matrix",Inv,n);
    printf("\n");
 
    time_elapsed = ((double)end - start) / CLOCKS_PER_SEC;
    printf("Naive Inverse      time duration with n = %d : %f s\n",n,time_elapsed);
   
    // Check if A * InvA = I
    correctionInv(A,Inv,n,p);
    freeMatrix(Z);
    freeMatrix(identity);
    freeMatrix(Inv);
    
}
