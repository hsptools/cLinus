#include <stdio.h>
#include <math.h>
//#include "Atom3d.h"
#include "linus.h"


int _ishbonded_loos(LinusProtein *p, Atom3d *atom, AtomList *hbtab_loos,int size, int *numint_loos)
{
    int i;
    for (i=0;i<size;i++)
    {
        Atom3d d1 = p->atoms[hbtab_loos[i].a1];
        Atom3d acp = p->atoms[hbtab_loos[i].a2];
        double dist = hbtab_loos[i].value;
        if(atom->resnum == d1.resnum && !strcmp(p->res_names[atom->resnum],p->res_names[d1.resnum]) && !strcmp(atom->name,d1.name))
        {
            *numint_loos=*numint_loos+1;
            return 1;
        }
        else if(atom->resnum == acp.resnum && !strcmp(p->res_names[atom->resnum],p->res_names[acp.resnum]) && !strcmp(atom->name,acp.name))
        {
            *numint_loos=*numint_loos+1;
            return 1;
        }
        
    }
    return 0;

}

