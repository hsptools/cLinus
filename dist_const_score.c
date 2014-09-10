#include <stdio.h>
#include <math.h>
#include "linus.h"

double dist_const_score(LinusProtein *p,DcList *dclist, int size)
{
    double edist = 0.0;
    int i;
    for(i=0;i<size;i++)
    {
        double dist = distance(&p->atoms[dclist[i].a1],&p->atoms[dclist[i].a2]);
        double dmax = dclist[i].d + dclist[i].dmax;
        double dmin = dclist[i].d - dclist[i].dmin;
        if(dist > dmax)
        {
            edist+= sqrt(dist-dmax);
        }
        else if(dist < dmax && dist > dmin)
        {
            edist += 0.0;
        }
        else if(dist < dmin)
        {
            edist += sqrt(dmin-dist);
        }
    }
    return edist;
}

