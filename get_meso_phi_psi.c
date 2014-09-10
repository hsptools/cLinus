#include <stdio.h>
#include <string.h>
#include "linus.h"
#include "meso.h"

void get_meso_phi_psi(char **mesolist, double *philist, double *psilist,double *omelist,double *s)
{
    int i,j;
    double phi1,phi2,psi1,psi2,phi,psi,ome;
    
    //printf("Searching meso\n");
    for(i=0;i<4;i++)
    {
        for(j=0;j<145;j++)
        {
            if(!strcmp(mesolist[i],meso_name[j]))
            {
                phi1 = meso_phi1[j];
                phi2 = meso_phi2[j];
                psi1 = meso_psi1[j];
                psi2 = meso_psi2[j];
                break;
            }
        }
        
        phi = uniform(phi1,phi2,s);
        //printf("Phi %lf Phi 1 %lf Phi2 %lf\n",phi,phi1,phi2);
        if(phi > 180.0)
            phi = phi - 360.0;
        else if (phi < - 180.0)
            phi = phi + 360.0;

        psi = uniform(psi1,psi2,s);
        //printf("psi %lf psi 1 %lf psi2 %lf\n",psi,psi1,psi2);
        psi = psi+180.0;
        if(psi > 180.0)
            psi = psi - 360.0;
        else if (psi < - 180.0)
            psi = psi + 360.0;
            
        ome = uniform(175.0,185.0,s);
        //printf("%s phi %lf (%lf,%lf) psi %lf (%lf,%lf) ome %lf\n",mesolist[i],phi,phi1,phi2,psi,psi1,psi2,ome);
        philist[i+1] = phi;
        psilist[i+1] = psi;
        omelist[i+1] = ome;
    }
}

