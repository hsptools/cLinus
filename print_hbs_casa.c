#include <stdio.h>
#include <string.h>
#include "linus.h"
//#include "CntmultList.h"


/*
"""Prints PDB file containing:
           o number of solvating waters on backbone groups - occupancy column
           o CASA - Bfactor column
           o number of unsatisfied backbone Hbond D and A - in TER record

    Arguments

      o *filename* - filename or open file - location to which ASA
        information will be written

      o *mol* - instance of linusMol - molecule whose accessible surface
        area is to be calculated

      o *flags* - object - the object returned by the *make_asa_list*
        function in the **construct** module.

      o *probe* - float - the water probe size to use for accessible surface

      o *ndiv* - integer - a number from 1 to 5 representing the number of
        sampling points used in the area calculation.  From least (1) to
        most (5) accurate, this is either 60, 240, 960, 3840, or 15,350
        points.

      o *ext_atoms* - list of coordinates - a list of tuples, each
        tuple containing (x, y, z) of an additional atom to include
        in the ASA calculation

      o *ext_radius* - float - the radius of the atoms in ext_atoms.

      o *title* - string - An optional title for the PDB output
    """
*/
void print_hbs_casa(LinusProtein *p, int * flags, double probe, int ndiv,char *filename)
{
    if(probe < 0) probe = 1.4;
    if(ndiv < 0) ndiv = 3;
    double ext_radius = 1.4;
    
    int minres,maxres;
    get_res_extents(p,&minres,&maxres);
    
    CntmultList solv_list[p->num_res];
    int i;
    for(i=minres;i<maxres+1;i++)
    {
        solv_list[i].res_id_N = i;
        solv_list[i].is_N = 1;
        solv_list[i].num_N = 0;
        solv_list[i].flag_N = 0;
        solv_list[i].res_id_O = i;
        solv_list[i].is_O = 2;
        solv_list[i].num_O = 0;
        solv_list[i].flag_O = 0;
        //printf("%d\n",solv_list[i].num_N);
    }
    
    
    
    int use_ext = 0;
    double *ext_coords;
    
    int use_data  =1;
    double data[p->num_atoms];
    //printf("Calling asa_eval %d\n",p->num_atoms);
  
    
    
    double tot_asa = asa_eval(p->atoms, p->num_atoms, data, ext_coords, flags, probe,
                        ext_radius, use_data, use_ext, ndiv);
                       // printf("Asa %lf\n",tot_asa);
                        
    FILE *file = fopen (filename,"w+");                    
    
    double p_solv_nrg = 0.0;
    double ap_solv_nrg = 0.0;
    double Gamma_p = 3.0/5.0;
    double Gamma_hb_oxy = 2.0;
    double Gamma_ap = 0.03;
    double CHASA = 0.0;
    
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
                    p_solv_nrg = p_solv_nrg - (Gamma_p *(solv_list[i].num_N));
                } 
                //if not solvated add penalty                
                else if(solv_list[i].num_N < 0)
                {
                     p_solv_nrg = p_solv_nrg + 1.0;
                }
                fprintf(file,"ATOM  %5d %4s %3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",j+1,p->atoms[j].name,p->res_names[i],p->res_num[i],p->atoms[j].x,p->atoms[j].y,p->atoms[j].z,solv_list[i].num_N,data[j]);
            }
            
            else if (!strcmp(p->atoms[j].name," O  "))
            {
                if (solv_list[i].num_O > 0)
                {
                    if(solv_list[i].flag_O ==0)
                    {
                        p_solv_nrg = p_solv_nrg - (Gamma_p *(solv_list[i].num_O));
                    }
                    else if(solv_list[i].flag_O > 0)
                    {
                        p_solv_nrg = p_solv_nrg - (Gamma_hb_oxy);
                    }
                } 
                              
                //if not solvated add penalty                
                else if(solv_list[i].num_O < 0)
                {
                     p_solv_nrg = p_solv_nrg + 1.0;
                }
                fprintf(file,"ATOM  %5d %-4s %-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",j+1,p->atoms[j].name,p->res_names[i],p->res_num[i],p->atoms[j].x,p->atoms[j].y,p->atoms[j].z,solv_list[i].num_O,data[j]);
            }
            else
            {
                 fprintf(file,"ATOM  %5d %-4s %-3s  %4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",j+1,p->atoms[j].name,p->res_names[i],p->res_num[i],p->atoms[j].x,p->atoms[j].y,p->atoms[j].z,occ,data[j]);
            }
            if(strchr(p->atoms[j].name,'C') != NULL)
            {
                ap_solv_nrg = ap_solv_nrg + (Gamma_ap * data[j]);
                CHASA = CHASA + data[j];
            }
           
        }
    
    
    
    }
    double  tot_solv_nrg = ap_solv_nrg + p_solv_nrg;
    fprintf(file,"ER %8.3f %8.3f\n",CHASA,tot_solv_nrg);    
    fclose(file);
    
    

}


