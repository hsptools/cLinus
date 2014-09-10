#include <stdio.h>
#include <math.h>
#include <string.h>
//#include "Atom3d.h"
//#include "Atom3d.h"
#include "linus.h"

#define RADIANS_TO_DEGREES 57.29577951308232e0
#define DEGREES_TO_RADIANS 0.017453292519943295e0

double torsion(Atom3d *a1, Atom3d *a2, Atom3d *a3, Atom3d *a4)
{
    double tor;

    double a = a1->x - a2->x;
    double b = a1->y - a2->y;
    double c = a1->z - a2->z;

    double d = a3->x - a2->x;
    double e = a3->y - a2->y;
    double f = a3->z - a2->z;

    double g = a4->x - a3->x;
    double h = a4->y - a3->y;
    double i = a4->z - a3->z;

    double ax = b*f - e*c;
    double ay = c*d - f*a;
    double az = a*e - b*d;

    double bx = h*f - e*i;
    double by = i*d - f*g;
    double bz = g*e - h*d;

    double xx = ay*bz - by*az;
    double yy = az*bx - bz*ax;
    double zz = ax*by - bx*ay;

    double cross = sqrt(xx*xx + yy*yy + zz*zz);
    double dot = ax*bx + ay*by + az*bz;

    tor = fabs(atan2(cross, dot)) * RADIANS_TO_DEGREES;

    if ((ax*(e*bz-f*by) + ay*(f*bx-d*bz) +  az*(d*by-e*bx)) >=0.0e0)
        tor = -tor;

    return tor;
    //printf("Main %lf %lf %lf\n",a1->x, a1->y, a1->z);
    //printf("First %lf %lf %lf\n",a2->x, a2->y, a2->z);
    //printf("second %lf %lf %lf\n",a3->x, a3->y, a3->z);
    //printf("third %lf %lf %lf\n",a4->x, a4->y, a4->z);

    //a1->tp_torsion = ang;
    //printf("Torsion %lf for atom %s num %d\n", a1->tp_torsion, a1->name, a1->resnum);
}

