#include <stdio.h>
#include <math.h>
#include <string.h>
//#include "Atom3d.h"
#include "linus.h"


#define RADIANS_TO_DEGREES 57.29577951308232e0
#define DEGREES_TO_RADIANS 0.017453292519943295e0

void bond_rotate(Atom3d *batom1, Atom3d *batom2, Atom3d *atoms, double theta, int n)
{
    int i;
    double trad = theta * DEGREES_TO_RADIANS;
    double sinth = sin(trad);
    double costh = cos(trad);
    double dist= distance(batom1, batom2);
    double ox = batom1->x;
    double oy = batom1->y;
    double oz = batom1->z;

    //printf("ox %lf oy %lf oz %lf dist %lf\n",ox,oy,oz,dist);

    // direction cosines of the bond axis

    double lambda = (batom2->x - ox) / dist;
    double mu = (batom2->y - oy) / dist;
    double nu = (batom2->z - oz) / dist;

    //printf("lambda %lf mu %lf nu %lf\n",lambda, mu, nu);

    // construct the rotation matrix

    double zip = 1.0e0 - costh;
    double rot00=costh+lambda*lambda*zip;
    double rot01=lambda*mu*zip-nu*sinth;
    double rot02=lambda*nu*zip+mu*sinth;
    double rot10=lambda*mu*zip+nu*sinth;
    double rot11=costh+mu*mu*zip;
    double rot12=mu*nu*zip-lambda*sinth;
    double rot20=lambda*nu*zip-mu*sinth;
    double rot21=mu*nu*zip+lambda*sinth;
    double rot22=costh+nu*nu*zip;


    //printf("rot00 %lf rot22 %lf\n", rot00, rot22);

    // apply rotation

    for (i=0; i<n; i++)
    {

        Atom3d *atom = &atoms[i];

        //printf("Old X %lf Y %lf Z %lf",atom->x,atom->y,atom->z);
        double dx = atom->x - ox;
        double dy = atom->y - oy;
        double dz = atom->z - oz;
        atom->x = ox + dx*rot00 + dy*rot01 + dz*rot02;
        atom->y = oy + dx*rot10 + dy*rot11 + dz*rot12;
        atom->z = oz + dx*rot20 + dy*rot21 + dz*rot22;
        //printf("New X %lf Y %lf Z %lf\n",atom->x,atom->y,atom->z);
    }
}

