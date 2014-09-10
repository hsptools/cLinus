#include <stdio.h>
#include <math.h>
#include "linus.h"

double _distance(Atom3d *a,double *b)
{
    //printf("%lf %lf %lf %lf %lf %lf\n",a->x,a->y,a->z,b[0],b[1],b[2]);
    return sqrt( pow( (a->x -b[0]) , 2) + pow( (a->y -b[1]) , 2) + pow( (a->z -b[2]) , 2) );
}

