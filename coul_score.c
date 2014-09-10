#include <stdio.h>
#include <math.h>
#include "linus.h"

double coul_score(LinusProtein *p,AtomList *qlist, int size, double scale)
{
    double ecoul = 0.0;
    double D = 4.0;
    int i;
    double V;
    for(i=0;i<size;i++)
    {
        double r = distance(&p->atoms[qlist[i].a1],&p->atoms[qlist[i].a2]);
        V = 331.0 * (qlist[i].value/(D*r));
        ecoul += V;
        
    }
    ecoul *= scale;
    return ecoul;
}

