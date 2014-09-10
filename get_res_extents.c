#include <stdio.h>
//#include "Atom3d.h"
#include "linus.h"
#include <string.h>

void get_res_extents(LinusProtein *p, int *minres, int *maxres)
{
    if(!strcmp(p->res_names[0],"ACE"))
    {
        *minres =1;
    }
    else
    {
        *minres = 0;
    }
    if(!strcmp(p->res_names[p->num_res-1],"NME"))
    {
        *maxres = p->num_res-1;
    }
    else
    {
        *maxres = p->num_res;
    }

}
