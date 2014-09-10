#include <stdio.h>
#include <math.h>
#include <string.h>
//#include "Atom3d.h"
#include "linus.h"
#define DEGREES_TO_RADIANS 0.017453292519943295e0


int make_hbtab_loos(LinusProtein *p, HbList *hblist, int size, int loos_size, AtomList *hbtab_loos)
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
            
            if (!(strcmp(d1.name, " OG ") && strcmp(d1.name, " OG1") &&
            strcmp(d1.name, " OH ") && strcmp(d1.name, " NH1") &&
            strcmp(d1.name, " NZ ") && strcmp(d1.name, " NH2") &&
            strcmp(d1.name, " ND1") && strcmp(d1.name, " ND2") &&
            strcmp(d1.name, " OD2") && strcmp(d1.name, " SG ") &&
            strcmp(d1.name, " NE2") && strcmp(d1.name, " OE2") &&
            strcmp(d1.name, " NE1") && strcmp(d1.name, " OE1") &&
            strcmp(d1.name, " OD1")))
            {
                if (close(&d1, &acp, hbm))
                {
                    if (!(angle(&d1, &acp, &acp1) > 90.0e0))
                    {
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

            else if(!(strcmp(acp.name," OG ") && strcmp(acp.name," OG1")) && (strcmp(acp.name," OH ") && 
                strcmp(acp.name," ND1") && strcmp(acp.name," ND2") && strcmp(acp.name," OD2") && 
                strcmp(acp.name," SG ") && strcmp(acp.name," NE2") && strcmp(acp.name," OE2") && 
                strcmp(acp.name," NE1") && strcmp(acp.name," OE1") && strcmp(acp.name," OD1"))) 
            {
                if (close(&d1, &acp, hbm))
                {
                    if (!(angle(&d1, &acp, &acp1) > 90.0e0))
                    {
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
            else            
            {
                if (close(&d1, &acp, hbm))
                {
                    //Orig 
                    if(fabs(torsion(&acp,&d1,&da1,&da2)) > 130.0 && angle(&da1,&d1,&acp) >69.0 &&
                    angle(&da2,&d1,&acp) >69.0 && angle(&d1,&acp,&acp1) >90.0)
                    
                    // from dr. panasik script
                    //if( fabs(torsion(&acp,&d1,&da1,&da2)) > 170.0 &&
                    //angle(&da1,&d1,&acp) > 109.0 &&            
                    //angle(&da2,&d1,&acp) > 109.0 &&
                    //angle(&d1,&acp,&acp1) > 110.0 )
                    
                    {
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
    }    

    return num_loos;
}

