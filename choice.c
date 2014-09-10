#include <stdio.h>
#include <string.h>
#include "linus.h"
int choice(int size, double *s)
{
    double ranval;
    int index;
    ranval = random1(s);
    //printf("Ranval %lf\n",ranval);
    index = (int)(size * ranval);
    //strcpy(res,list[index]);
    //printf("Index %d\n",index);
    return index;
}

