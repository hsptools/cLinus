#include <stdio.h>
#include <math.h>
#include <string.h>
//#include "Atom3d.h"
#include "linus.h"


#define RADIANS_TO_DEGREES 57.29577951308232e0
#define DEGREES_TO_RADIANS 0.017453292519943295e0


double angle(Atom3d *a1, Atom3d *a2, Atom3d *a3)
{
    double ang;
    //printf("Before Angle parents %s %s %lf %lf\n", a2->name,a3->name, a2->x, a3->x);
    double x1 = a1->x - a2->x;
    double y1 = a1->y - a2->y;
    double z1 = a1->z - a2->z;

    double x2 = a3->x - a2->x;
    double y2 = a3->y - a2->y;
    double z2 = a3->z - a2->z;

    double xx = y1*z2 - y2*z1;
    double yy = z1*x2 - z2*x1;
    double zz = x1*y2 - x2*y1;

    double cross = sqrt(xx*xx + yy*yy + zz*zz);
    double dot = x1*x2 + y1*y2 + z1*z2;

    ang = fabs(atan2(cross, dot)) * RADIANS_TO_DEGREES;

    return ang;

    //printf("After Angle parents %s %s\n", a2->name,a3->name);

    //a1->sp_angle = ang;
   //printf("Angle %lf %s %d\n", a1->sp_angle, a1->name, a1->resnum);
}



