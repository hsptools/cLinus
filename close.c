#include <stdio.h>
#include <math.h>
//#include "Atom3d.h"
#include "linus.h"


int close(Atom3d *a1, Atom3d *a2, double tol)
{
    double dx, dy, dz;

    if (tol < 0.0e0)
    {
        tol= 0.90*(a1->radius) + 0.90*(a2->radius);
        if (tol<=1.0e-6) return 0;
    }
    dx = fabs(a1->x - a2->x);
    dy = fabs(a1->y - a2->y);
    dz = fabs(a1->z - a2->z);
  
    //printf("%s %s %lf %lf %lf %lf\n",a1->name,a2->name,dx,dy,dz,tol);

    if(dx >= tol) return 0;
    if(dy >= tol) return 0;
    if(dz >= tol) return 0;
    if((dx*dx + dy*dy + dz*dz) > tol*tol) return 0;
    return 1;
    
    
    
    
    
   /* d=d*d;
    //printf(" d %lf ",d);
    dx = (a1->x - a2->x);
    dx=dx*dx;
    //printf(" dx %lf ",dx);
    if (dx >= d) return 0;
    dy = (a1->y - a2->y);
    dy=dy*dy;
    //printf(" dy %lf ",dy);
    if (dy >= d) return 0;
    dz = (a1->z - a2->z);
    dz=dz*dz;
    //printf(" dz %lf ",dz);
    if (dz >= d) return 0;
    //printf(" dxyz %lf ",dx+dy+dz);

    if ((dx + dy + dz) < d)
        return 1;
    else
        return 0;
        
    */    
}

