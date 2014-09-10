#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define M1   4294967087
#define M2   4294944443

void gen_seeds(double *s)
{
    srand (time(NULL));
    s[0] = (double)(rand()%((M1-1) +1)+1);
    s[1] = (double)(rand()%((M1-1) +1)+1);
    s[2] = (double)(rand()%((M1-1) +1)+1);
    s[3] = (double)(rand()%((M2-1) +1)+1);
    s[4] = (double)(rand()%((M2-1) +1)+1);
    s[5] = (double)(rand()%((M2-1) +1)+1);
    
}
