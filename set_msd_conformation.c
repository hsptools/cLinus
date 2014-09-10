#include <stdio.h>
#include <string.h>
#include "linusSetConf.h"
#include "linus.h"
#include "Chidict.h"
#include "meso.h"

//set the conformation of a specified segement to helical conformation

/*
 """Set the conformation of a segment in accordance with the
    accumulated phi/psi mean and SD for each residue in the segment.

    Arguments

        o *p* - instance of linusMol - the protein whose conformation is
        to be changed

        o *start* - integer - index of the first residue in the segment

        o *end* - integer - index of the last residue in the segment -
        *note* - only the conformation of residues start to end -1 are
        changed.

        o *generator* - object - a random number generator

    Description

        Chi-values are chosen at random.  
    """    


*/


void set_msd_conformation(LinusProtein *p, int start, int end, double *s)
{
    int *phi_atoms = p->phi_atoms;
    int *psi_atoms = p->psi_atoms;
    int *ome_atoms = p->omega_atoms;
    int *chi_atoms = p->chi_atoms;
    double *fm = p->res_f_mean;
    double *fsd = p->res_f_sd;
    double *ym = p->res_y_mean;
    double *ysd = p->res_y_sd;
    
    int index=-1;
    char ms[3]="";
    double phi1,phi2,psi1,psi2,phi,psi,ome=0.0;
    int res_ind,ran_ind,pos,atom_num;
    double chi;
   
    
    int i,j,k;
    for(i=start;i<end;i++)
    {
        phi = uniform(fm[i]-fsd[i],fm[i]+fsd[i],s);
        //printf("Phi %lf Phi 1 %lf Phi2 %lf\n",phi,phi1,phi2);
        if(phi > 180.0)
            phi = phi - 360.0;
        else if (phi < - 180.0)
            phi = phi + 360.0;

        psi = uniform(ym[i]-ysd[i],ym[i]+ysd[i],s);
        //printf("psi %lf psi 1 %lf psi2 %lf\n",psi,psi1,psi2);
        psi = psi+180.0;
        if(psi > 180.0)
            psi = psi - 360.0;
        else if (psi < - 180.0)
            psi = psi + 360.0;
            
        ome = uniform(175.0,185.0,s);
        //printf("Ome %lf Ome 1 %lf Ome2 %lf\n",ome,175.0,185.0);
        if(ome > 180.0)
            ome = ome - 360.0;
        
                
        p->atoms[phi_atoms[i]].tp_torsion = phi;
        p->atoms[psi_atoms[i]].tp_torsion = psi;
        p->atoms[ome_atoms[i]].tp_torsion = ome;
        p->atoms[phi_atoms[i]].sp_angle = uniform(111.0, 113.0,s);    
        printf("%s %s %d Phi %lf Psi 1 %lf Ome %lf\n",p->atoms[i].name,p->res_names[i],p->atoms[i].resnum,phi,psi,ome);
    }  // for range i 
    
    
    int curr = -1;
    
    
    for(k=0;k<p->num_chi_atoms;k++)
    {
        if(p->atoms[chi_atoms[k]].resnum != curr)
        {
            curr = p->atoms[chi_atoms[k]].resnum;
            res_ind = get_chi_res_index(p->res_names[curr]);
            atom_num =0;
            if(res_ind > -1)
                ran_ind = choice(chi_size[res_ind],s);
        }
        //printf("Res id %d ran id %d\n",res_ind,ran_ind);
        if(res_ind > -1)
        {
            //ran_ind = choice(0,chi_start[res_ind]);
            pos = (ran_ind * chi_size[res_ind]) + chi_start[res_ind]+(2*atom_num);
            //printf("%s %s %d %d %lf %lf\n",p->res_names[curr],p->chi_atoms[k].name,pos,pos+1,chi_value[pos],chi_value[pos+1]);
            chi = uniform(chi_value[pos],chi_value[pos+1],s);
            //printf("Before Chi %lf ",chi);  
            atom_num+=1;
            if(chi > 180.0)
                chi = chi -360.0;
            else if (chi < -180.0)
                chi = chi + 360.0;
            p->atoms[chi_atoms[k]].tp_torsion = chi; 
            //printf("Chi %lf\n",chi);    
        }
    }
    
    
   
}


