#include <stdio.h>
#include <math.h>
#include "linus.h"

/*
#   
     #
    # V = Eps,i,j[(Sigma,i,j/ri,j)**12]
    #
    # where Sigma = (Rmin/2)**(1/6)
    # and Rmin/2 is the vdw radius
    #
    # epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
    # but we will use esp = 0.1 for all atom types
    #
    # Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
    # where Rmin/2 is the atom radius
    #
    # Sigma is precalculated and is in vlist (see construct.py)
*/

double soft_debump_score(LinusProtein *p, AtomList *vlist, int size)
{
    double esdebmp = 0.0;
    double eps = 1.0;
    int i;
    double V;
    for(i=0;i<size;i++)
    {
        if(close(&p->atoms[vlist[i].a1], &p->atoms[vlist[i].a2],-1.0))
        {
            double r = distance(&p->atoms[vlist[i].a1], &p->atoms[vlist[i].a2]);
            double Sigma = vlist[i].value;
            V = eps * (pow(Sigma/r,12));
        }
        else
        {
            V=0.0;
        }
        
        esdebmp += V;
    }    
    return esdebmp;         
}

