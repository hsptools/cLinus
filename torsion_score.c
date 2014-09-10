#include <stdio.h>
#include <string.h>
//#include "Atom3d.h"
#include <stdlib.h>
//#include "LinusProtein.h"
#include "linus.h"

double torsion_score(LinusProtein *p, double tp)
{
    //Default  tp = 1.0;

    double etot=0.0;


    int i;

    for (i=0;i<p->num_res;i++)
    {
        //printf("i %d\n",i);
        //printf("Atom %s %d in res %s\n",p->phi_atoms[i].name,p->phi_atoms[i].resnum, p->res_names[i]);
        if(p->atoms[p->phi_atoms[i]].tp_torsion > 0.0 && p->atoms[p->phi_atoms[i]].tp_torsion <150.0)
        {
            if(!strcmp(p->res_names[i],"GLY"))
            {
                etot = etot -tp;
            }
            else if(!strcmp(p->res_names[i],"ASN"))
            {
                etot = etot;
            }
            else
            {
                //printf("Before tp is %lf\n", *etot);
                 etot = etot + tp;
                //printf("After etot is %lf\n",*etot);
            }

        }
    }
    return etot;
    //printf("Final etot is %lf\n",*etot);
    //tor_score = etot;
    //printf("Tor score is %lf\n",*tor_score);
}

