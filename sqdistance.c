#include <stdio.h>
//#include "Atom3d.h"
#include "linus.h"

double sqdistance(Atom3d *a1, Atom3d *a2)
{
    double sqdist;
    double dx = a1->x - a2->x;
    double dy = a1->y - a2->y;
    double dz = a1->z - a2->z;
    //printf("%lf\n", dx*dx + dy*dy + dz*dz);
    sqdist = (dx*dx + dy*dy + dz*dz);
    return sqdist;
}
