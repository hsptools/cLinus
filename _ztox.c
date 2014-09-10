#include <stdio.h>
#include <math.h>
#include <string.h>
//#include "Atom3d.h"
#include "linus.h"


#define DEGRAD 0.017453292519943295e0


void _ztox(double *coordsfp, double *coordssp, double *coordstp, double *zmvalu, double *atcoor)
{

   // Convert internal coordinates to cartesian coordinates
   //
   //    Arguments
   //
   //       o *coords* - list of type 'f' and 3 elements
   //
   //       o *zmvalu* - list of type 'f' and 3 elements, where each
   //       element is respecively the distance, angle and torsion angle of
   //       the atom to it's parent, second parent and third parent.
   //
   //       o *zmindx* - list of type 'i' and 3 elements, where each
   //       element is respectively the integer index of the atoms' first,
   //       second and third parents.
   //
   //    Returns
   //
   //       atcoor - the converted x,y,z coords

 
   double fpd, spa, tpt;
   double sina, cosa, sint, cost;
   double fpx, fpy, fpz;
   double spx, spy, spz;
   double tpx, tpy, tpz;
   double u1x, u1y, u1z;
   double u2x, u2y, u2z;
   double u3x, u3y, u3z;
   double u4x, u4y, u4z;
   double sine, cosine;

   fpd = zmvalu[0];
   spa = zmvalu[1];
   tpt = zmvalu[2];

   spa = spa * DEGRAD;
   tpt = tpt * DEGRAD;

   sina = sin(spa);
   cosa = -cos(spa);

   sint = sina * sin(tpt);
   cost = sina * cos(tpt);

   fpx = coordsfp[0];
   fpy = coordsfp[1];
   fpz = coordsfp[2];

   spx = coordssp[0];
   spy = coordssp[1];
   spz = coordssp[2];

   tpx = coordstp[0];
   tpy = coordstp[1];
   tpz = coordstp[2];

   //u2x = x_norm(fpx, fpy, fpz, spx, spy, spz);
   //u2y = y_norm(fpx, fpy, fpz, spx, spy, spz);
   //u2z = z_norm(fpx, fpy, fpz, spx, spy, spz);

   //u1x = x_norm(spx, spy, spz, tpx, tpy, tpz);
   //u1y = y_norm(spx, spy, spz, tpx, tpy, tpz);
   //u1z = z_norm(spx, spy, spz, tpx, tpy, tpz);

   _norm(fpx, fpy, fpz, spx, spy, spz,&u2x,&u2y,&u2z);

   _norm(spx, spy, spz, tpx, tpy, tpz,&u1x,&u1y,&u1z);


   cosine = u1x*u2x + u1y*u2y + u1z*u2z;
   //printf("Cosin %lf %lf %lf \n",cosine, );

   if (fabs(cosine) < 1.0)
   {
      sine = 1.0/sqrt(1.0 - (cosine*cosine));
      //printf("Sin %lf %lf \n",sine,cosine*cosine,1.0 - (cosine*cosine));
   }
   else
      sine = 1.0/sqrt((cosine*cosine) - 1.0);
   //printf("Sin %lf\n",sine);

   u3x = sine * (u1y*u2z - u1z*u2y);
   u3y = sine * (u1z*u2x - u1x*u2z);
   u3z = sine * (u1x*u2y - u1y*u2x);

   u4x = cost * (u3y*u2z - u3z*u2y);
   u4y = cost * (u3z*u2x - u3x*u2z);
   u4z = cost * (u3x*u2y - u3y*u2x);

   //printf("%lf %lf %lf %lf %lf %lf \n",u3x,u3y,u3z,u4x,u4y,u4z);

   atcoor[0] = fpx + fpd*(u2x*cosa + u3x*sint + u4x);
   atcoor[1] = fpy + fpd*(u2y*cosa + u3y*sint + u4y);
   atcoor[2] = fpz + fpd*(u2z*cosa + u3z*sint + u4z);

   //return(atcoor);
}
