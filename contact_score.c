#include <stdio.h>
//#include "Atom3d.h"
#include "linus.h"

/*

Compute contact score for a molecule

Arguments

    o *clist* - object - the object returned by the *make_contact_list*
    function in the *construct* module

    o *cprobe* - float - the distance over which the contact energy
    is to be scaled to zero.

Result

    The contact energy as a floating point number

*/

double contact_score(LinusProtein *p,AtomList *clist, int size,double cprobe)
{

    double eborpc, dx, dy, dz, dist;
    double econ = 0.0e0;
    int i;
    eborpc = 1.0e0/cprobe;

    for(i=0;i<size;i++)
    {
        double dmin = clist[i].value;
        double dmax = dmin + cprobe;
        if(close(&p->atoms[clist[i].a1],&p->atoms[clist[i].a2],dmax))
        {
            double d12=distance(&p->atoms[clist[i].a1],&p->atoms[clist[i].a2]);
            if(d12 < dmin)
                econ += clist[i].econ;
            else
                econ += clist[i].econ*(dmax-d12)*eborpc;
        }


    }

   return -econ;
}

/*
print_contact_score

Print contact energy details to a file

Arguments

    o *mol* - instance of LinusMol - molecule whose contact score
    is to be computed and listed to file

    o *clist* - object - the object returned by the *make_contact_list*
    function in the *construct* module

    o *cprobe* - float - the distance over which the contact energy
    is to be scaled to zero.

    o *file* - string or file object - name of a file or a file
    previously opened for writing.

Result

    The contact energy of the molecule

*/
double print_contact_score(LinusProtein *p, AtomList *clist, int size,double cprobe, char *filename)
{
    FILE *file = fopen(filename,"w+");

    fprintf(file,"%s","\n------ Cotact Score Listing-----\n");

    double eborpc = 1.0/cprobe;
    double econ =0.0;
    int i;

    for(i=0;i<size;i++)
    {
        double dmin = clist[i].value;
        double dmax = dmin + cprobe;
        double etemp;
        if(close(&p->atoms[clist[i].a1],&p->atoms[clist[i].a2],dmax))
        {
            double d12=distance(&p->atoms[clist[i].a1],&p->atoms[clist[i].a2]);
            if(d12 < dmin)
                etemp = clist[i].econ;
            else
                etemp = clist[i].econ*(dmax-d12)*eborpc;

            fprintf(file,"%3d %s %s %3d %s %s Distance = %lf Score = %lf\n",
                    p->atoms[clist[i].a1].resnum,p->res_names[p->atoms[clist[i].a1].resnum], p->atoms[clist[i].a1].name,
                    p->atoms[clist[i].a2].resnum,p->res_names[p->atoms[clist[i].a2].resnum], p->atoms[clist[i].a2].name,
                    d12,etemp);
            econ += etemp;

        }

    }
    fprintf(file,"\nTotal Contact Score = %g\n",-econ);
    fclose(file);
    return -econ;
}

