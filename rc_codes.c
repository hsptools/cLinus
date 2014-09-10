#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "linusSecStr.h"
#include "linus.h"

//rc_codes returns mesostate codes for protein
void rc_codes(LinusProtein *p, char **codes)
{
    int i;
    for(i=0;i<p->num_res;i++)
    {
        strcpy(codes[i],"??");
        //printf("Codes %s\n",codes[i]);
    }
    
    for(i=2;i<p->num_res-1;i++)
    {
        double phi = p->atoms[p->phi_atoms[i]].tp_torsion;
        double psi = p->atoms[p->psi_atoms[i]].tp_torsion -180.0;
        if (psi > 180.0)
        {
            psi -= 360.0;
        }
        else if(psi < -180)
        {
            psi += 360;
        }
        double ome = p->atoms[p->omega_atoms[i]].tp_torsion;
        res_rc(phi,psi,ome,codes[i]);
    }

}
