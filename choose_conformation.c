#include <stdio.h>
#include <string.h>
#include "linusSetConf.h"
#include "linus.h"
#include "Chidict.h"
#include "meso.h"

/*
  """Choose a conformation given the index of a residue

    Arguments

        o *p* - instance of linusMol - the protein whose conformation is
        to be changed

        o *res_id* - integer - index of residue

    Description

        If mesostates have been specified for the residue then one of
        the allowed mesostates is chosen.  Otherwise one of helix,
        strand, turn1, turn2 and coil is chosen in accordance with the
        sampling weights of the residue.

    Result

        return a tuple of 2 elements, the first being a single character
        representation ('C' for mesostate choice, 'H' for helix, 'E' for
        strand, '1' for turn1, '2' for turn2 and 'C' for coil) of the
        move chosen and the second item being the function corresponding
        to the move chosen
    """    

*/


void choose_conformation(LinusProtein *p, int res_id, char **meso_key, int *meso_val,int meso_size,double *s, char *state, char *func,
    int trial, int use_pms_sampling, int use_msd_focus, int use_local_mov, int use_ms_sampling, int use_msd_sampling)
{
    if(trial < 0) trial =0;
    if(use_pms_sampling < 0) use_pms_sampling =0;
    if(use_msd_focus < 0) use_msd_focus =0;
    if(use_local_mov < 0) use_local_mov =0;
    if(use_ms_sampling < 0) use_ms_sampling =0;
    if(use_msd_sampling < 0) use_msd_sampling =0;
    
    if(use_local_mov == 1 && res_id <= (p->num_res-6) && use_ms_sampling == 0)
    {
        strcpy(state,"C");
        strcpy(func,"local_move");
    }
     else if(use_local_mov == 1 && res_id <= (p->num_res-6) && use_ms_sampling == 1)
    {
        strcpy(state,"C");
        strcpy(func,"local_ms_move");
    }    
    else if(use_local_mov == 1 && res_id > (p->num_res-6))
    {
        if(use_ms_sampling ==1 && strcmp(p->mesostates[res_id],"--"))
        {
            strcpy(state,"C");
            strcpy(func,"set_mesostate_conformation");
        }
        else
        {
            strcpy(state,"C");
            strcpy(func,"set_coil");
        }
    }    
    
    if(trial > 1 && use_msd_focus == 1)
    {
        if(p->res_f_sd[res_id] < 5.0 && p->res_f_sd[res_id] > 0.0001 && p->res_y_sd[res_id] < 5.0 && p->res_y_sd[res_id] > 0.0001 )
        {
            strcpy(state,"C");
            strcpy(func,"set_msd_conformation");
        }
    
    }
    
    if(use_msd_sampling ==1 && use_pms_sampling ==0)
    {
         if(p->res_f_mean[res_id-1] < 999.0 && p->res_f_mean[res_id] < 999.0 && p->res_f_mean[res_id+1] < 999.0)
         {
            strcpy(state,"C");
            strcpy(func,"set_msd_conformation");    
         }
         else
         {
            strcpy(state,"C");
            strcpy(func,"set_coil");
         }
    }
    
    int i,j;
    double sum=0.0;
    
    if(use_msd_sampling ==1 && use_pms_sampling ==1)
    {
        if(p->res_f_mean[res_id-1] < 999.0 && p->res_f_mean[res_id] < 999.0 && p->res_f_mean[res_id+1] < 999.0)
        {
            strcpy(state,"C");
            strcpy(func,"set_msd_conformation");    
        }
        else
        {
             double r, rsum;
             int pos,found;
             
            // get pmeso_states
            // residue 1 has res_id 1 but
            //  pmeso 1 has index 0
            
            // get mesostate for res_id - 1
         
            for(i=0;i<144;i++)
            {
                sum += p->pmeso_wt[res_id-2][i];
            }
           
            r = uniform(1.0e-9, sum,s);
            pos =0;
            found = 0;
            rsum = 0.0;
            for(i=0;i<144;i++)
            {
                if(!found)
                {
                    rsum = rsum + p->pmeso_wt[res_id-2][i];
                    if(r <= rsum)
                    {
                        for(j=0;j<meso_size;j++)
                        {
                            if(meso_val[j] == pos)
                            {
                                strcpy(p->mesostates[res_id-1],meso_key[j]);
                                found ==1;
                                break;
                            }
                        }
                    }
                    pos +=1;
               }
               
           }
           
            // get mesostate for res_id + 1
         
            for(i=0;i<144;i++)
            {
                sum += p->pmeso_wt[res_id-2][i];
            }
           
            r = uniform(0.0001, sum,s);
            pos =0;
            found = 0;
            rsum = 0.0;
            for(i=0;i<144;i++)
            {
                if(!found)
                {
                    rsum = rsum + p->pmeso_wt[res_id+0][i];
                    if(r <= rsum)
                    {
                        for(j=0;j<meso_size;j++)
                        {
                            if(meso_val[j] == pos)
                            {
                                strcpy(p->mesostates[res_id+1],meso_key[j]);
                                found ==1;
                                break;
                            }
                        }
                    }
                    pos +=1;
               }
               
           }
           
            // get mesostate for res_id
         
            for(i=0;i<144;i++)
            {
                sum += p->pmeso_wt[res_id-1][i];
            }
           
            r = uniform(0.0001, sum,s);
            pos =0;
            found = 0;
            rsum = 0.0;
            for(i=0;i<144;i++)
            {
                if(!found)
                {
                    rsum = rsum + p->pmeso_wt[res_id-2][i];
                    if(r <= rsum)
                    {
                        for(j=0;j<meso_size;j++)
                        {
                            if(meso_val[j] == pos)
                            {
                                strcpy(p->mesostates[res_id],meso_key[j]);
                                strcpy(state,"C");
                                strcpy(func,"set_mesostate_conformation");
                            }
                        }
                    }
                    pos +=1;
               }
               
           }
          
        }
    }
       
    if(use_pms_sampling ==1 && use_msd_sampling ==0)
    {
         
            double r, rsum;
             int pos,found;
             
            // get pmeso_states
            // residue 1 has res_id 1 but
            //  pmeso 1 has index 0
            
            // get mesostate for res_id - 1
         
            for(i=0;i<144;i++)
            {
                sum += p->pmeso_wt[res_id-2][i];
            }
           
            r = uniform(1.0e-9, sum,s);
            pos =0;
            found = 0;
            rsum = 0.0;
            for(i=0;i<144;i++)
            {
                if(!found)
                {
                    rsum = rsum + p->pmeso_wt[res_id-2][i];
                    if(r <= rsum)
                    {
                        for(j=0;j<meso_size;j++)
                        {
                            if(meso_val[j] == pos)
                            {
                                strcpy(p->mesostates[res_id-1],meso_key[j]);
                                strcpy(p->mesostates[res_id+1],meso_key[j]);
                                found ==1;
                                break;
                            }
                        }
                    }
                    pos +=1;
               }
               
           }
           
            // get mesostate for res_id + 1
         
            for(i=0;i<144;i++)
            {
                sum += p->pmeso_wt[res_id-2][i];
            }
           
            r = uniform(0.0001, sum,s);
            pos =0;
            found = 0;
            rsum = 0.0;
            for(i=0;i<144;i++)
            {
                if(!found)
                {
                    rsum = rsum + p->pmeso_wt[res_id+0][i];
                    if(r <= rsum)
                    {
                        for(j=0;j<meso_size;j++)
                        {
                            if(meso_val[j] == pos)
                            {
                                strcpy(p->mesostates[res_id+1],meso_key[j]);
                                found ==1;
                                break;
                            }
                        }
                    }
                    pos +=1;
               }
               
           }
           
            // get mesostate for res_id
         
            for(i=0;i<144;i++)
            {
                sum += p->pmeso_wt[res_id-1][i];
            }
           
            r = uniform(0.0001, sum,s);
            pos =0;
            found = 0;
            rsum = 0.0;
            for(i=0;i<144;i++)
            {
                if(!found)
                {
                    rsum = rsum + p->pmeso_wt[res_id-2][i];
                    if(r <= rsum)
                    {
                        for(j=0;j<meso_size;j++)
                        {
                            if(meso_val[j] == pos)
                            {
                                strcpy(p->mesostates[res_id],meso_key[j]);
                                strcpy(state,"C");
                                strcpy(func,"set_mesostate_conformation");
                            }
                        }
                    }
                    pos +=1;
               }
               
           }
         
         
    }
    
    if(use_ms_sampling ==1 && strcmp(p->mesostates[res_id],"--"))
    {
        strcpy(state,"C");
        strcpy(func,"set_mesostate_conformation");
    }
    
    else
    {
        double hw = p->res_helix_wt[res_id];
        double sw = p->res_strand_wt[res_id];
        double t1w = p->res_turn1_wt[res_id];
        double t2w = p->res_turn2_wt[res_id];
        double cw = p->res_coil_wt[res_id];
        double pw = p->res_PII_wt[res_id];
        double r = uniform(0.0,pw,s);
        if(r<=hw)
        {
            strcpy(state,"H");
            strcpy(func,"set_helix");
        }
        else if(r<=sw)
        {
            strcpy(state,"E");
            strcpy(func,"set_strand");
        }
        else if(r<=t1w)
        {
            strcpy(state,"1");
            strcpy(func,"set_turn");
        }
        else if(r<=t2w)
        {
            strcpy(state,"2");
            strcpy(func,"set_turn");
        }
        else if(r<=cw)
        {
            strcpy(state,"C");
            strcpy(func,"set_coil");
        }
        else if(r<=pw)
        {
            strcpy(state,"P");
            strcpy(func,"set_PII");
        }
        
    }
          
    
}

















































