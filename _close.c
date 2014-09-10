#include <stdio.h>
#include <math.h>


int _close(double *coords1, double *coords2, double hsd)
{
   double x1, y1, z1, x2, y2, z2, dist;

   x1 = coords1[0];
   y1 = coords1[1];
   z1 = coords1[2];
   x2 = coords2[0];
   y2 = coords2[1];
   z2 = coords2[2];

   dist = sqrt(pow((x1 - x2), 2) + pow((y1 - y2), 2) + pow((z1 - z2), 2));
   //printf("distance %lf\n",dist);
   
   if (dist < hsd)
      return(1);

   return(0);
}
