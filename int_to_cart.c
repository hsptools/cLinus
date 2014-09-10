#include <stdio.h>
#include <math.h>
#include <string.h>
//#include "Atom3d.h"
#include "linus.h"


#define RADIANS_TO_DEGREES 57.29577951308232e0
#define DEGREES_TO_RADIANS 0.017453292519943295e0

// takes atoms a1 and its parents a2,a3,a4
void int_to_cart(Atom3d *a1, Atom3d *a2,Atom3d *a3,Atom3d *a4)
{

    double dist;
    double rangl, rtors, sina, cosa, sint, cost;
    double u1x, u1y, u1z;
    double u2x, u2y, u2z;
    double u3x, u3y, u3z;
    double u4x, u4y, u4z;
    double sine, cosine, d;
    
    if (a1->fp_ind == -1)
    {
        printf("No first parent!\n");
        return ;
    }

    if (a1->sp_ind == -1)
    {
        printf("No second parent!\n");
        return;
    }

    if (a1->tp_ind == -1)
    {
        printf("No Third parent!\n");
        return;
    }

    rangl = DEGREES_TO_RADIANS * (a1->sp_angle);
    rtors = DEGREES_TO_RADIANS * (a1->tp_torsion);
    dist = a1->fp_distance;

    //printf("d %lf a %lf t %lf\n",a1->fp_distance,a1->sp_angle,a1->tp_torsion);
    
    sina = sin(rangl); cosa = -cos(rangl);
    sint = sina * sin(rtors); cost = sina * cos(rtors);

    u1x = a3->x - a4->x;
    u1y = a3->y - a4->y;
    u1z = a3->z - a4->z;
    d = 1.0e0 / sqrt(u1x*u1x + u1y*u1y + u1z*u1z);
    u1x *= d;
    u1y *= d;
    u1z *= d;

    //printf("u1x %lf u1y %lf u1z %lf\n",u1x,u1y,u1z);

    u2x = a2->x - a3->x;
    u2y = a2->y - a3->y;
    u2z = a2->z - a3->z;
    d = 1.0e0 / sqrt(u2x*u2x + u2y*u2y + u2z*u2z);
    u2x *= d;
    u2y *= d;
    u2z *= d;

    cosine = u1x*u2x + u1y*u2y + u1z*u2z;

    if (fabs(cosine) < 1.0)
        sine = 1/sqrt(1.0 - cosine*cosine);
    else
        sine = 1/sqrt(cosine*cosine - 1.0);

    u3x = sine * (u1y*u2z - u1z*u2y);
    u3y = sine * (u1z*u2x - u1x*u2z);
    u3z = sine * (u1x*u2y - u1y*u2x);

    u4x = cost * (u3y*u2z - u3z*u2y);
    u4y = cost * (u3z*u2x - u3x*u2z);
    u4z = cost * (u3x*u2y - u3y*u2x);

    //printf("u1x %lf u1y %lf u1z %lf u2x %lf u2y %lf u2z %lf u3x %lf u3y %lf u3z %lf u4x %lf u4y %lf u4z %lf\n",u1x,u1y,u1z,u2x,u2y,u2z,u3x,u3y,u3z,u4x,u4y,u4z);

    //printf("Old: X %lf Y %lf Z %lf a2x %lf a2y %lf a2z %lf  %lf %lf %lf \n",a1->x,a1->y,a1->z,a2->x,a2->y,a2->z,dist*(u2x*cosa + u4x + u3x*sint),dist*(u2y*cosa + u4y + u3y*sint),dist*(u2z*cosa + u4z + u3z*sint));
    a1->x = a2->x + dist*(u2x*cosa + u4x + u3x*sint);
    a1->y = a2->y + dist*(u2y*cosa + u4y + u3y*sint);
    a1->z = a2->z + dist*(u2z*cosa + u4z + u3z*sint);
    //printf("New: X %lf Y %lf Z %lf\n",a1->x,a1->y,a1->z);
}

void ztox(Atom3d *atoms, int start, int end) // start=None, end=None)
{
   // Save internal coordinates following a conformation change
   int i;
   for (i = start; i < end; i++)
   {
        //printf("Ztox %s %d %d %d\n",atoms[i].name,atoms[i].fp_ind,atoms[i].sp_ind,atoms[i].tp_ind);
        if(atoms[i].fp_ind > -1 && atoms[i].sp_ind > -1 && atoms[i].tp_ind > -1 )
        {
            //printf("Searching int to cart %s\n",atoms[i].name,&atoms[atoms[i].fp_ind].name, &atoms[atoms[i].sp_ind].name, &atoms[atoms[i].tp_ind].name);
            int_to_cart(&atoms[i],&atoms[atoms[i].fp_ind], &atoms[atoms[i].sp_ind], &atoms[atoms[i].tp_ind] );
        }    

   }
}


void ztox_waters(Atom3d *water,Atom3d *atoms, int start, int end) // start=None, end=None)
{
   // Save internal coordinates following a conformation change
   int i;
   for (i = start; i < end; i++)
   {
        //printf("Ztox %s %d %d %d\n",atoms[i].name,atoms[i].fp_ind,atoms[i].sp_ind,atoms[i].tp_ind);
        if(water[i].fp_ind > -1 && water[i].sp_ind > -1 && water[i].tp_ind > -1 )
        {
            //printf("Searching int to cart %s\n",atoms[i].name,&atoms[atoms[i].fp_ind].name, &atoms[atoms[i].sp_ind].name, &atoms[atoms[i].tp_ind].name);
            int_to_cart(&water[i],&atoms[water[i].fp_ind], &atoms[water[i].sp_ind], &atoms[water[i].tp_ind] );
        }    

   }
}

