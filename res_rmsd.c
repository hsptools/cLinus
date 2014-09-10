#include <stdio.h>
#include <math.h>
#include "linus.h"


double res_rmsd(Atom3d *a,Atom3d *b,Atom3d *c,Atom3d *d,double *e,double *f,double *g,double *h)
{
    int N=4;
    double sum2 = 0.0;
    double dist,rmsd;
    
    
    dist = _distance(a,e);
    sum2 = sum2 + (dist *dist);
    //printf("dist %lf , sum2 %lf ",dist, sum2);
    
    dist = _distance(b,f);
    sum2 = sum2 + (dist *dist);
    //printf("dist %lf , sum2 %lf ",dist, sum2);
    
    
    dist = _distance(c,g);
    sum2 = sum2 + (dist *dist);
    //printf("dist %lf , sum2 %lf ",dist, sum2);
    
    
    dist = _distance(d,h);
    sum2 = sum2 + (dist *dist);
    //printf("dist %lf , sum2 %lf\n",dist, sum2);
    
    
    
    rmsd = sqrt((int)(sum2/N));
    return rmsd;
    

}


