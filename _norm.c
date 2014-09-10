#include <stdio.h>
#include <math.h>
//#include "Atom3d.h"
#include "linus.h"


void _norm(double x1, double y1, double z1, double x2, double y2, double z2, double * x_norm, double *y_norm, double *z_norm)
{
    //printf ("X1 %lf Y1 %lf Z1 %lf X2 %lf Y2 %lf Z2 %lf\n", x1,y1,z1,x2,y2,z2);
    double dx = x1 - x2;
    double dy = y1 - y2;
    double dz = z1 - z2;
    double d = 1/sqrt(dx*dx + dy*dy + dz*dz);
    *x_norm = dx*d;
    *y_norm = dy*d;
    *z_norm = dz*d;

    //printf("Norm %lf %lf %lf\n",*x_norm,*y_norm,*z_norm );
}
