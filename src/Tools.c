#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Tools.h"


/**
 *  @file  This file contains all the needed functions for the project 
*/


/**
 * @brief Compute the modular inverse
 * @return the modular inverse of x mod p
*/
int inv(int x, int p)
{
    if(x == 1) return 1;

    if(x < 0) x = abs(x);

    if(x == 0)
        x += 1;


    int r0 = x, u0 = 1;
    int r1 = p, u1 = 0;
    int q, tmp;

    while(r1 != 0)
    {
        q = r0 / r1;
        
        tmp = u1;
        u1 = u0 - q * u1;
        u0 = tmp;

        tmp = r1;
        r1 = r0 % r1;
        r0 = tmp;
    }
    
    if(r0 != 1)
    {
        printf("Inverse of %d mod %d doesn't exist\n",x,p);
        exit(1);
        return -1;
    }

    if(u0 < 0) u0 += p;
    return u0;
}




/**
 * @brief compute modular addition
*/
int add(int x, int y, int p)
{
    int a = x + y;
    if(a >= p) return a - p;
    if(a < 0) return a + p;
    return a;
}

/**
 * @brief compute modular substraction
*/
int sub(int x, int y, int p)
{
    int s = x - y;
    if(s >= p) return s - p;
    if(s < 0) return s + p;
    return s;
}

/**
 * @brief Check if a number is a prime number
*/
int isPrime(int *n)
{
    // Corner case
    if (*n <= 1)
        return -1;
 
    // Check from 2 to square root of n
    for (int i = 2; i <= sqrt(*n); i++)
        if (*n % i == 0)
            return -1;
 
    return 0;
}

