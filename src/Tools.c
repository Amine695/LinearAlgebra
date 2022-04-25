#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "Tools.h"


/* This file contains all the needed functions for the project */


/**
 * Extended Euclid algorithm
*/
int gcdExtended(int a, int b, int *x, int *y)
{
    // Base Case
    if (a == 0)
    {
        *x = 0;
        *y = 1;
        return b;
    }
    if(a < 0)
        a = abs(a);
 
    int x1, y1; // To store results of recursive call
    int gcd = gcdExtended(b%a, a, &x1, &y1);

    // Update x and y using results of recursive
    // call
    *x = y1 - (b/a) * x1;
    *y = x1;
 
    return gcd;
}


int modInverse(int a, int m)
{
    int x, y;
    int g = gcdExtended(a, m, &x, &y);
    if (g != 1) {
        printf("\nInverse doesn't exist, gcd = %d != 1\n",g);
        return g;
    }
    
    
    // m is added to handle negative x
    int res = (x % m + m) % m;
    return res;

}

