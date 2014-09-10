#include <stdio.h>
#include <math.h>
//include "Atom3d.h"
#include "linus.h"
double distance(Atom3d *a1, Atom3d *a2)
{
    double d;
    //printf ("X1 %lf Y1 %lf Z1 %lf X2 %lf Y2 %lf Z2 %lf\n", a1->x,a1->y,a1->z,a2->x,a2->y,a2->z);
    double dx = a1->x - a2->x;
    double dy = a1->y - a2->y;
    double dz = a1->z - a2->z;
    d = sqrt(dx*dx + dy*dy + dz*dz);
    
    //printf("Distance parents %lf \n", d);  
    return d;
    //printf("Saved Distance %s %d %s %d %lf\n",a1->name,a1->resnum,a2->name,a2->resnum, *d);   
}
