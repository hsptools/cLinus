#include <stdio.h>
#include "Atom3d.h"

void get_coords(Atom3d *self, double * coords)
{
    coords[0] = self->x;
    coords[1] = self->y;
    coords[2] = self->z;    
}

