#include <stdio.h>
#include <math.h>
//#include "Atom3d.h"
//#include "LinusProtein.h"
//#include "AtomList.h"
#include "linus.h"
//int make_vdw_list(LinusProtein * p, AtomList * vlist, int size);

/*
    #
    # V(Lennard-Jones) = 4Eps,i,j[(Sigma,i,j/ri,j)**12 - (Sigma,i,j/ri,j)**6]
    #
    # epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
    # but we will use esp = 0.15 for all atom types
    #
    # Sigma/2 = Rmin/2**(1/6)
    # where Rmin/2 is the atom radius
    #
    # Sigma is precalculated and is in vlist (see construct.py)
*/
    
double LJ_score(LinusProtein *p, AtomList *vlist, int size, double scale)
{
    //scale = 1.0
    double LJ = 0.0;
    double eps = 0.15;
    double foureps = 4.0 * eps;
    
    int i;
    double V;
    
    for(i=0;i<size;i++)
    {
        double r = distance(&p->atoms[vlist[i].a1], &p->atoms[vlist[i].a2]);
        double Sigma = vlist[i].value;
        
        V = foureps*(pow(Sigma/r,12)- pow(Sigma/r,6));
        //printf("%d %s %d %s %lf %lf %lf\n",p->atoms[vlist[i].a1].resnum,p->atoms[vlist[i].a1].name,
        //p->atoms[vlist[i].a2].resnum,p->atoms[vlist[i].a2].name,r,Sigma,V);
        
        LJ += V;
        //printf("V %lf LJ %lf\n",V,LJ);
    }
    LJ *= scale;
    return LJ;
}

