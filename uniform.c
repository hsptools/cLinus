#include <stdio.h>
#include "linus.h"
double uniform(double start, double end, double *s)
{
    double diff, ranval;
    diff = end-start; 
    ranval = random1(s);
    return (start+ranval*diff);
}

