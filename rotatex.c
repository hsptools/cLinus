#include <stdio.h>
#include <math.h>
#include <string.h>
//#include "Atom3d.h"
#include "linus.h"

#define RADIANS_TO_DEGREES 57.29577951308232e0
#define DEGREES_TO_RADIANS 0.017453292519943295e0

void rotatex(Atom3d *atom, double theta)
{
    double y, z;
    double trad = theta * DEGREES_TO_RADIANS;
    double sint = sin(trad);
    double cost = cos(trad);

    y = atom->y*cost - atom->z*sint;
    z = atom->y*sint + atom->z*cost;
    atom->y = y;
    atom->z = z;

    //printf("%lf %lf\n", atom->y,atom->z);
}
