#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "linus.h"

void set_weights_all(Linus_sim *sim,char *move_type,double weight)
{
    int i;
    double *w;
    //printf("Type is %s %lf\n",move_type,weight);
    
    if(!strcmp(move_type,"helix"))
        w = sim->protein.res_helix_wt;
    else if(!strcmp(move_type,"strand"))
        w = sim->protein.res_strand_wt;
    else if(!strcmp(move_type,"coil"))
        w = sim->protein.res_coil_wt;
    else if(!strcmp(move_type,"PII"))
        w = sim->protein.res_PII_wt;
    else  if(!strcmp(move_type,"turn1"))
        w = sim->protein.res_turn1_wt;
    else if(!strcmp(move_type,"turn2"))
        w = sim->protein.res_turn2_wt;
        

    
    for(i=0;i<sim->protein.num_res;i++)
    {
        //printf("%d\n",i); 
        w[i] = weight;
    }
   
  
}


void set_weights_int(Linus_sim *sim,char *move_type,double weight, int *res,int size)
{
    double *w;
    
    if(!strcmp(move_type,"helix"))
        w = sim->protein.res_helix_wt;
    else if(!strcmp(move_type,"strand"))
        w = sim->protein.res_strand_wt;
    else if(!strcmp(move_type,"coil"))
        w = sim->protein.res_coil_wt;
    else if(!strcmp(move_type,"PII"))
        w = sim->protein.res_PII_wt;
    else  if(!strcmp(move_type,"turn1"))
        w = sim->protein.res_turn1_wt;
    else if(!strcmp(move_type,"turn2"))
        w = sim->protein.res_turn2_wt;
        
    int i;
    
    for(i=0;i<size;i++)
    {
            w[res[i]] = weight;
    }
  
}

void set_weights_char(Linus_sim *sim,char *move_type,double weight, char *res,int size)
{
    double *w;
    
    if(!strcmp(move_type,"helix"))
        w = sim->protein.res_helix_wt;
    else if(!strcmp(move_type,"strand"))
        w = sim->protein.res_strand_wt;
    else if(!strcmp(move_type,"coil"))
        w = sim->protein.res_coil_wt;
    else if(!strcmp(move_type,"PII"))
        w = sim->protein.res_PII_wt;
    else  if(!strcmp(move_type,"turn1"))
        w = sim->protein.res_turn1_wt;
    else if(!strcmp(move_type,"turn2"))
        w = sim->protein.res_turn2_wt;
        
    int i,j;
    
    for(i=0;i<sim->protein.num_res;i++)
    {
        for(j=0;j<size;j++)
        {
            if(!strcmp(sim->protein.res_names[i],"res[j]"))
                w[i] = weight;
        }    
    }
  
}






