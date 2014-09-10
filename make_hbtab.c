#include <stdio.h>
#include <math.h>
#include "linus.h"
#define DEGREES_TO_RADIANS 0.017453292519943295e0

// for create_poly_ala_tet_meso1


int make_hbtab(LinusProtein *p, HbList *hblist, int size, int loos_size,AtomList *hbtab_loos)
{
    int nacp,i,j;
    double dist;
    int num_loos=0;
    
    for(i=0; i<size; i++)
    {
        Atom3d d1 = p->atoms[hblist[i].donor];
        Atom3d da1 = p->atoms[hblist[i].da1];
        Atom3d da2 = p->atoms[hblist[i].da2];
        nacp = hblist[i].acps_size;

        for (j=0; j<nacp; j++)
        {
            Atom3d acp = p->atoms[hblist[i].acps[j].acp];
            Atom3d acp1 = p->atoms[hblist[i].acps[j].acp1];
            double hbd = hblist[i].acps[j].hbdist;
            double hbm = hblist[i].acps[j].hbdmax;
            double hbs = hblist[i].acps[j].hbene;
            //printf("%s %s %s %s %s\n",d1.name,da1.name,da2.name,acp.name,acp1.name);
            // first three keep oxygen in +/- 10 deg cone of NH
            // Equivalent to angle >= 170 in Taylor and Kennard
            // Last angle keeps N in +/- 70 deg cone of CO
            // Equivalent to phi >= 20 in Taylor and Kennard
           
            if (close(&d1, &acp, hbm))
            {
                //printf("Close\n");
                if( fabs(torsion(&acp,&d1,&da1,&da2)) > 170.0 &&
                    angle(&da1,&d1,&acp) > 109.0 &&            
                    angle(&da2,&d1,&acp) > 109.0 &&
                    angle(&d1,&acp,&acp1) > 110.0 )
               
                {
                    //printf("Hbonded\n");
                    dist = distance(&d1, &acp);
                    if(loos_size)
                    {
                        hbtab_loos[num_loos].a1 = hblist[i].donor; 
                        hbtab_loos[num_loos].a2 = hblist[i].acps[j].acp;
                        hbtab_loos[num_loos].value = dist;
                    }
                    num_loos++;  
                    
                }
                
            }
        }
    }

    return num_loos; 
}    

