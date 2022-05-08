/**
 * @file main.c
 * @author Amine Berbagui
 * @brief file that contains the main function and launch the program.
 * @date 2022-05-05
 * @copyright Copyright (c) 2022
 * 
 */


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "LinearSystem.h"
#include "MatrixInv.h"
#include "Strassen.h"
#include "Benchmarks.h"


/**
 * @brief Function that launch the project when the user input n and p.
*/
void LaunchProject(matrix A,matrix Abis,matrix L,matrix U, matrix P, matrix S1, matrix S2,
            int Z[], int X[], int B[], DIM, prime)
{
    printf("************** Algorithmic Project ***************\n\n");
    //LU decomposition
    RunLU(A,P,L,U,n,p);
    // Linear system solving
    RunLS(L,U,B,Z,X,n,p);
    // Naive Inverse 
    RunNaiveInverse(A,L,U,P,n,p);
    // Strassen Inverse 
    RunStrassenInv(Abis,n,p);
    // Strassen multiplication
    RunStrassenMult(S1,S2,n,p); 

}

/**
 * @brief Function that launch the benchmarks wanted with real time execution
 * There are 4 possible options.
*/
void LaunchBenchmarks(matrix A,matrix L,matrix U,matrix P,matrix Abis,matrix S1, matrix S2,
                        DIM, prime,FILE *fp1, FILE *fp2,FILE *fp3,FILE *fp4,FILE *fp5,int choice)
{
    switch(choice)
    {
        case 1:
            BenchmarkNaiveInverse(A,L,U,P,n,p,fp1); 
            BenchmarkStrassenInv(Abis,n,p,fp2);
            break;

        case 2:
            BenchmarkNaiveInverse(A,L,U,P,n,p,fp1); 
            BenchmarkStrassenInv(Abis,n,p,fp2);
            BenchmarkStrassenInvAndMult(Abis,n,p,fp5);
            break;

        case 3:
            BenchmarksNaiveMul(S1,S2,n,p,fp3);
            BenchmarkStrassenMult(S1,S2,n,p,fp4);
            break;
        default:
            break;
    }

}

/**
 * @brief Display the menu for the benchmark
 * @return choice : the choice chosen by the user
 */
int ChooseBenchmark()
{
    printf("*************** Welcome to the benchmark menu!****************\n");
    printf("Press 1 : Naive matrix inversion VS Strassen's inversion using naive multiplication\n");
    printf("Press 2 : Naive matrix inversion VS Strassen's inversion using naive multiplication VS Strassen's inversion using Strassen's multiplication\n");
    printf("Press 3 : Naive multiplication   VS Strassen's multiplication\n");
    int choice;
    if(scanf("%d", &choice) != 0 && (choice - 1)*(choice - 3) > 0)
    {
        printf("Error, the input does not match the possible options\n");
        exit(1);
    }
    return choice;
}


/**
 * @brief Function that take the arguments and check that all is good
*/
void GetArguments(int argc, char *argv[], int *n, int *p)
{
    if(argc != 5 && argc != 2)
    {
        printf("Error : argument missing required\n");
        exit(1);
    }
    
    if(strcmp(argv[1],"--n")==0)
    {
        *n = atoi(argv[2]);
        if(*n % 2 != 0)
        {
            printf("Error : n must be a power of 2\n");
            exit(1);
        }

    }
    else
    {
        printf("Error : size argument incorrect\n");
        exit(-1);
    }
        
    if(strcmp(argv[3],"--p")==0)
    {
        *p = atoi(argv[4]);
        if(isPrime(p) != 0)
        {
            printf("Error : p must be a prime number\n");
            exit(-1);
        }
    }
        
}


/**
 * @brief Main function that launch the program
*/
int main(int argc, char* argv[])
{
    int n,p,choice;
    srand(time(NULL));

    // if the user type benchmarks, we run LaunchBenchmarks
    if(strcmp(argv[1],"benchmarks") == 0)
    {
        // Which functions to benchmark
        choice = ChooseBenchmark();

        // maximum size of 2^k
        int max;      
        int k = 1;
        // a prime number fixed
        int p = 1151; 

        // To export data
        FILE * fp1;
        FILE * fp2;
        FILE * fp3;
        FILE * fp4;
        FILE * fp5;

        // we grab the argument from user
        max = atoi(argv[2]);
        
        // path location 
        char NaiveInvPath[] = "benchmarks/NaiveInversion.txt";
        char StrassenInvPath[] = "benchmarks/StrassenInversion.txt";
        char NaiveMultPath[] = "benchmarks/NaiveMultiplication.txt";
        char StrassenMultPath[] = "benchmarks/StrassenMultiplication.txt";
        char StrassenInvMultPath[] = "benchmarks/StrassenInversionMultiplication.txt";
        
        // open files
        fp1 = fopen(NaiveInvPath,"w+");
        fp2 = fopen(StrassenInvPath,"w+");
        fp3 = fopen(NaiveMultPath,"w+");
        fp4 = fopen(StrassenMultPath,"w+");
        fp5 = fopen(StrassenInvMultPath,"w+");

        fprintf(fp1,"Size\t\tTime(s)");
        fprintf(fp2,"Size\t\tTime(s)");
        fprintf(fp3,"Size\t\tTime(s)");
        fprintf(fp4,"Size\t\tTime(s)");
        fprintf(fp5,"Size\t\tTime(s)");
        
        // benchmarks loop : we start from 2^1 to 2^max
        while(k <= max)
        {
            // new size
            int n = pow(2,k);

            // initialize size matrices
            matrix A,L,U,P,Abis,S1,S2;
            L = init_matrix(n);
            P = init_matrix(n);
            U = init_matrix(n);

            // generate a random matrix,Abis for Strassen Inverse,S1 and S2 for the Strassen multiplication
            A = generate(n,p);
            S1 = generate(n,p);
            S2 = generate(n,p);
            Abis = Copy(A,n);
            
            // Launch the benchmarks
            LaunchBenchmarks(A,L,U,P,Abis,S1,S2,n,p,fp1,fp2,fp3,fp4,fp5,choice);
            printf("\n");
            // increment size counter
            k++;

            // Memory free
            freeMatrix(Abis);
            freeMatrix(A);
            freeMatrix(S1);
            freeMatrix(S2);
            freeMatrix(L);
            freeMatrix(U);
            freeMatrix(P);
        }
        printf("\nData successfully exported in the benchmarks directory! \n");
        printf("Generating plot...\n");

        // close files
        fclose(fp1);
        fclose(fp2);
        fclose(fp3);
        fclose(fp4);
        fclose(fp5);

        // run python script for plot 
        // !! WARNING : you might have a different python interpreter path !!
        // !! if so, please change it at the top of the python script !!
        if (system("./plots/plot.py") == -1)
            exit(1);
        printf("Plot generated, take a look in the plot directory! \n");
    }
    // otherwise, we run LaunchProject
    else
    {
        // we take arguments from user 
        GetArguments(argc,argv,&n,&p);

        // same procedure as in the if
        matrix A,Abis, P, L, U,S1,S2;
        L = init_matrix(n);
        P = init_matrix(n);
        U = init_matrix(n);
        Abis = init_matrix(n);
        A = generate(n,p);
        S1 = generate(n,p);
        S2 = generate(n,p);
        Abis = Copy(A,n);
        
        // for the linear system solving
        int B[n]; // Solutions vector
        int Z[n]; // LZ = B
        int X[n]; // UX = Z
        
        // Launch all the algorithms
        LaunchProject(A,Abis,L,U,P,S1,S2,Z,X,B,n,p);
        printf("\n");

        // Memory free
        freeMatrix(Abis);
        freeMatrix(A);
        freeMatrix(L);
        freeMatrix(U);
        freeMatrix(P);
        freeMatrix(S1);
        freeMatrix(S2);


    }

    return 0;
}