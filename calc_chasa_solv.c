#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linus.h"
//#include "CntmultList.h"


/*
"""Return the CHASA solvation energy for use in simulations

    Arguments

      o *mol* - instance of linusMol - molecule whose accessible surface
        area is to be calculated
      o *hblist* - list of potential donors and acceptors
        In order for this list to be correct for the chasa_solv the 
        following hbond parameters must be set in the command file:

#   sim.set_chasa_parameters(use_chasa=1)
#   # if chasa is used must use following
#   sim.set_hbond_parameters(use_hbond=1,
#                            hbond_score_short=?.?,
#                            hbond_score_long=?.?,
#                            hbond_winmax=NUMRES,
#                            use_sidechain_hbond=1)

    """

*/

double calc_chasa_solv(LinusProtein *p, HbList *hblist, int hblist_size)
{
      
    int minres,maxres;
    get_res_extents(p,&minres,&maxres);
    
    int use_data  =1;
    double data[p->num_atoms];
    
    
    FILE *file;
    
    int numint_loos, numvirt_loos, bbtot,numbb;
    CntmultList solv_list[p->num_res];
    double *wat_list;
    
    int wlist = mk_virt_bb_cntmult_loosHbd(p,hblist,hblist_size,0,0,file,0,1.25,&numint_loos,&numvirt_loos,&bbtot,&numbb,wat_list,0,solv_list);
    wat_list = (double *) malloc(wlist*sizeof(double));
    numint_loos= numvirt_loos=bbtot=numbb=0;
    mk_virt_bb_cntmult_loosHbd(p,hblist,hblist_size,0,0,file,0,1.25,&numint_loos,&numvirt_loos,&bbtot,&numbb,wat_list,wlist,solv_list);
    //printf("numint_loos %d numvirt_loos %d bbtot %d numbb %d numext %d \n", numint_loos, numvirt_loos,bbtot, numbb,wlist);
    int i;
    
    int flags[p->num_atoms];
    make_asa_list(p,-1,-1,-1,flags);
    double probe = 1.4;
    int ndiv = 3;
    double ext_radius = 1.4;
    double tot_asa = asa_eval(p->atoms, p->num_atoms, data, wat_list, flags, probe,
                        ext_radius, use_data, wlist, ndiv);
    //printf("Asa %lf\n",tot_asa);


    double p_solv_nrg = 0.0;
    double ap_solv_nrg = 0.0;
    double Gamma_p = 3.0/5.0;
    double Gamma_hb_oxy = 2.0;
    double Gamma_ap = 0.03;
    double CHASA = 0.0;
    double non_hbd_score = 5.0;
  
    for(i=minres;i<maxres;i++)
    {
        double occ = 0.0;
        int start = p->res_fai[i];
        int end = p->res_fai[i+1];
        int j;
        for (j=start;j<end;j++)
        {
            
            //printf("Atom %s %d\n",p->atoms[j].name, p->atoms[j].resnum);
            if(!strcmp(p->atoms[j].name," N  "))
            {
                if (solv_list[i].num_N >0)
                {
                    p_solv_nrg = p_solv_nrg - (Gamma_p *(double)(solv_list[i].num_N));
                    //printf("N %d %lf\n",solv_list[i].num_N,p_solv_nrg);
                } 
                //if not solvated add penalty                
                
            }
            
            else if (!strcmp(p->atoms[j].name," O  "))
            {
                if (solv_list[i].num_O > 0)
                {
                    if(solv_list[i].flag_O ==0)
                    {
                        p_solv_nrg = p_solv_nrg - (Gamma_p *(double)(solv_list[i].num_O));
                        //printf("O1 %d %lf\n",solv_list[i].num_O,p_solv_nrg);
                    }
                    else if(solv_list[i].flag_O > 0)
                    {
                        p_solv_nrg = p_solv_nrg - (Gamma_hb_oxy);
                        //printf("O2 - %lf",p_solv_nrg);
                    }
                } 
            }
            
            else if(strchr(p->atoms[j].name,'C') != NULL)
            {
                //printf("C atom %s %lf\n",p->atoms[j].name,data[j]);
                ap_solv_nrg = ap_solv_nrg + (Gamma_ap * data[j]);
                
            }
           
        }
    
       //printf("ap %lf p %lf\n",ap_solv_nrg, p_solv_nrg);
    }
    
      double  tot_solv_nrg = ap_solv_nrg + p_solv_nrg;
      
    //printf("tot_solv_nrg %lf\n", tot_solv_nrg);  
    return tot_solv_nrg;
    
 
}


