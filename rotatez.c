#include <stdio.h>
#include <math.h>
#include <string.h>
//#include "Atom3d.h"
#include "linus.h"

#define RADIANS_TO_DEGREES 57.29577951308232e0
#define DEGREES_TO_RADIANS 0.017453292519943295e0

void rotatez(Atom3d *atom, double theta)
{
    double x, y;
    double trad = theta * DEGREES_TO_RADIANS;
    double sint = sin(trad);
    double cost = cos(trad);

    x = atom->x*cost - atom->y*sint;
    y = atom->y*cost + atom->x*sint;
    atom->y = y;
    atom->x = x;
}

