#include <stdio.h>
#include <math.h>
//#include "Atom3d.h"
#include "linus.h"
void geocenter(Atom3d *atoms,int start, int end, double *x, double *y, double *z )
{
    int i,n;
    double xc=0.0e0, yc=0.0e0, zc=0.0e0;

    for (i=start; i<end; i++)
    {
        Atom3d *atom = &atoms[i];
        xc += atom->x;
        yc += atom->y;
        zc += atom->z;
    }

    *x = xc/end;
    *y = yc/end;
    *z = zc/end;
}

double radius_of_gyration(Atom3d *atoms, int start, int end)
{


  /*
        o *radius_of_gyration* - compute the radius of gyration of
        the molecule

        """Radius of gyration of molecule or a subset of atoms in
        the molecule.  With no arguments specified the raidus of
        gyration of the entire molecule is computed.  If the
        radius of gyration of a contiguous subset of the atoms
        is required then the indices of the first and (last + 1)
        atoms should be supplied.

        Arguments

            o *start* - integer - index of the first atom

            o *end* - integer - index of last atom

        Result

            The radius of gyration as a floating point number
    */
    double rg;

    double xc, yc, zc;
    int i;

    geocenter(atoms, start, end, &xc, &yc, &zc);
    for(i=start; i<end; i++)
    {
        Atom3d *atom = &atoms[i];
        double dx = atom->x - xc;
        double dy = atom->y - yc;
        double dz = atom->z - zc;
        rg += (dx*dx + dy*dy + dz*dz);
    }

    rg /= end;
    rg = sqrt(rg);

    return rg;
    //printf("Radius of gyration %lf",*rg);
}

