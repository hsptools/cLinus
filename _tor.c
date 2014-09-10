#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "linus.h"
#define RADIANS_TO_DEGREES 57.29577951308232e0


/*
    """Calculates the torsion angle between four atoms.
       This function is based on the fact that the cross product of
       two vectors will produce a vector normal to the plane defined
       by those two vectors.  Then, a torsion angle can be defined by
       determining the angle between those two normal vectors.
    """
*/

double _tor(Atom3d *a1, Atom3d *a2, Atom3d *a3, Atom3d *a4)
{
    double v1p1[3] = {a1->x - a2->x, a1->y - a2->y, a1->z - a2->z};
    
    double v2p1[3] = {a3->x - a2->x, a3->y - a2->y, a3->z - a2->z};
    
    double v1p2[3] = {a2->x - a3->x, a2->y - a3->y, a2->z - a3->z};
    
    double v2p2[3] = {a4->x - a3->x, a4->y - a3->y, a4->z - a3->z};
    
    double np1[3] = {   v1p1[1] * v2p1[2] - v1p1[2] * v2p1[1],
                        v1p1[2] * v2p1[0] - v1p1[0] * v2p1[2],
                        v1p1[0] * v2p1[1] - v1p1[1] * v2p1[0]};
                        
    double np2[3] = {   v1p2[1] * v2p2[2] - v1p2[2] * v2p2[1],
                        v1p2[2] * v2p2[0] - v1p2[0] * v2p2[2],
                        v1p2[0] * v2p2[1] - v1p2[1] * v2p2[0]};
   
    double dp = 0.0;
    double np1mag = 0.0;
    double np2mag = 0.0;
    double stp = 0.0;
    
    int i;
    for(i=0;i<3;i++)
    {
        dp += np1[i] * np2[i];
        np1mag = np1mag + np1[i]*np1[i];
        np2mag = np2mag + np2[i]*np2[i];
        stp = stp + np2[i]*v1p1[i];
    }
    
    np1mag = sqrt(np1mag);
    np2mag = sqrt(np2mag);
    
    double sign;
    if(stp > 0.0) sign = 1.0;
    else sign = -1.0;
    
    double cost = dp/np1mag/np2mag;
    
    //printf("_tor %lf\n",sign * acos(cost) * RADIANS_TO_DEGREES);
    return sign * acos(cost) * RADIANS_TO_DEGREES;
}

