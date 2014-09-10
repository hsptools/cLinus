#include <stdio.h>
#include <string.h>
#include "linusSetConf.h"
#include "linus.h"
#include "Chidict.h"
#include "meso.h"

//set the conformation of a specified segement to helical conformation

/*
 """Set a segment to  COIL-PII-COIL conformation

    Arguments

        o *p* - instance of linusMol - the protein whose conformation is
          to be changed
          
        o *start* - integer - index of the first residue of the segment
          whose conformation is to be changed

        o *end* - integer - index of the last residue of the segment whose
          conformation is to be changed. **Note** only the conformation
          of residues *start* to *end-1* will be changed

        o *generator* - object - a random number generator
        

    Description

        The central residue is set to PII, the 1st and 3rd to COIL.
        Chi-values for the residues
        are chosen at random.  This function modifies the Z-matrix of
        the molecule.  Cartesian coordinates are *not* calculated.
    """    
*/


void set_PII(LinusProtein *p, int start, int end, double *s)
{
    int *phi_atoms = p->phi_atoms;
    int *psi_atoms = p->psi_atoms;
    int *ome_atoms = p->omega_atoms;
    int *chi_atoms = p->chi_atoms;
    double *pw = p->res_PII_wt;
    
    //printf("%s %s %s %s %d",phi_atoms[0].name,psi_atoms[0].name,omega_atoms[0].name,chi_atoms[0].name,hw[0]);
    
    int i,j,k;
    double phi1,phi2,psi1,psi2,phi,psi,ome;
    char ms[3];
    int index;
    int res_ind,ran_ind,pos,atom_num;
    double chi;
    
    // First residue in moveset
    //printf("\n\nFirst residue\n");
    i= start;
    
    if(!strcmp(p->res_names[i],"PRO"))
    {
        index = choice(18,s);
        strcpy(ms,PRO_COIL_MS[index]);
        //printf("Meso %s\n",ms);
           
        for(j=0;j<145;j++)
        {
            if(!strcmp(ms,meso_name[j]))
            {
                phi1 = meso_phi1[j];
                phi2 = meso_phi2[j];
                psi1 = meso_psi1[j];
                psi2 = meso_psi2[j];
                //printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                break;
            }
        }
    }    
    else if(!strcmp(p->res_names[i],"GLY"))
    {
        index = choice(93,s);
        strcpy(ms,GLY_COIL_MS[index]);
        //printf("Meso %s\n",ms);
           
        for(j=0;j<145;j++)
        {
            if(!strcmp(ms,meso_name[j]))
            {
                phi1 = meso_phi1[j];
                phi2 = meso_phi2[j];
                psi1 = meso_psi1[j];
                psi2 = meso_psi2[j];
                //printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                break;
            }
        }
    }
    else
    {
        index = choice(94,s);
        strcpy(ms,COIL_MS[index]);
        //printf("Meso %s\n",ms);
           
        for(j=0;j<145;j++)
        {
            if(!strcmp(ms,meso_name[j]))
            {
                phi1 = meso_phi1[j];
                phi2 = meso_phi2[j];
                psi1 = meso_psi1[j];
                psi2 = meso_psi2[j];
                //printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                break;
            }
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
    //printf("Ome %lf Ome 1 %lf Ome2 %lf\n",ome,175.0,185.0);
    if(ome > 180.0)
        ome = ome - 360.0;
    
        
    /*
    try:
        phi_atoms[i].tp_torsion = phi
    except TypeError:
        pass
    else:
        phi_atoms[i].sp_angle = uniform(108.0, 111.0)
    
    */
    
     
    p->atoms[phi_atoms[i]].tp_torsion = phi;
    p->atoms[psi_atoms[i]].tp_torsion = psi;
    p->atoms[ome_atoms[i]].tp_torsion = ome;
    p->atoms[phi_atoms[i]].sp_angle = uniform(108.0, 111.0,s);
    //printf("%s %s %d Phi %lf Psi 1 %lf Ome %lf\n",p->atoms[i].name,p->res_names[i],p->atoms[i].resnum,phi,psi,ome);
       
    //int curr = -1;
    
    
    for(k=0;k<p->num_chi_atoms;k++)
    {
        if(p->atoms[chi_atoms[k]].resnum == i)
        {
            res_ind = get_chi_res_index(p->res_names[i]);
            atom_num =0;
            if(res_ind > -1)
            {
                ran_ind = choice(chi_size[res_ind],s);
                pos = (ran_ind * chi_size[res_ind]) + chi_start[res_ind]+(2*atom_num);
                chi = uniform(chi_value[pos],chi_value[pos+1],s);
                atom_num+=1;
                if(chi > 180.0)
                    chi = chi -360.0;
                else if (chi < -180.0)
                    chi = chi + 360.0;
                p->atoms[chi_atoms[k]].tp_torsion = chi; 
                printf("Chi %lf\n",chi);    
                
            }    
        }
    }    
            
  
    // Second residue in moveset
    //printf("\n\nSecond residue\n");
    i= start+1;
    index = choice(2,s);
    strcpy(ms,PII_MS[index]);
    //printf("Meso %s\n",ms);
           
    for(j=0;j<145;j++)
    {
        if(!strcmp(ms,meso_name[j]))
        {
            phi1 = meso_phi1[j];
            phi2 = meso_phi2[j];
            psi1 = meso_psi1[j];
            psi2 = meso_psi2[j];
            //printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
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
    //printf("Ome %lf Ome 1 %lf Ome2 %lf\n",ome,175.0,185.0);
    if(ome > 180.0)
        ome = ome - 360.0;
     
    p->atoms[phi_atoms[i]].tp_torsion = phi;
    p->atoms[phi_atoms[i]].sp_angle = uniform(108.0, 111.0,s);
    p->atoms[psi_atoms[i]].tp_torsion = psi;
    p->atoms[ome_atoms[i]].tp_torsion = ome;
    //printf("%s %s %d Phi %lf Psi 1 %lf Ome %lf\n",p->atoms[i].name,p->res_names[i],p->atoms[i].resnum,phi,psi,ome);
       
   
    for(k=0;k<p->num_chi_atoms;k++)
    {
        if(p->atoms[chi_atoms[k]].resnum == i)
        {
            res_ind = get_chi_res_index(p->res_names[i]);
            atom_num =0;
            if(res_ind > -1)
            {
                ran_ind = choice(chi_size[res_ind],s);
                pos = (ran_ind * chi_size[res_ind]) + chi_start[res_ind]+(2*atom_num);
                chi = uniform(chi_value[pos],chi_value[pos+1],s);
                atom_num+=1;
                if(chi > 180.0)
                    chi = chi -360.0;
                else if (chi < -180.0)
                    chi = chi + 360.0;
                p->atoms[chi_atoms[k]].tp_torsion = chi; 
                printf("Chi %lf\n",chi);    
                
            }    
        }
    }
  
  
  // Third residue in moveset
    
    //printf("\n\nThird residue\n");    
    
    i= start +2;
        
    if(!strcmp(p->res_names[i],"PRO"))
    {
        index = choice(18,s);
        strcpy(ms,PRO_COIL_MS[index]);
        //printf("Meso %s\n",ms);
           
        for(j=0;j<145;j++)
        {
            if(!strcmp(ms,meso_name[j]))
            {
                phi1 = meso_phi1[j];
                phi2 = meso_phi2[j];
                psi1 = meso_psi1[j];
                psi2 = meso_psi2[j];
                //printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                break;
            }
        }
    }    
    else if(!strcmp(p->res_names[i],"GLY"))
    {
        index = choice(93,s);
        strcpy(ms,GLY_COIL_MS[index]);
        //printf("Meso %s\n",ms);
           
        for(j=0;j<145;j++)
        {
            if(!strcmp(ms,meso_name[j]))
            {
                phi1 = meso_phi1[j];
                phi2 = meso_phi2[j];
                psi1 = meso_psi1[j];
                psi2 = meso_psi2[j];
                //printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                break;
            }
        }
    }
    else
    {
        index = choice(94,s);
        strcpy(ms,COIL_MS[index]);
        //printf("Meso %s\n",ms);
           
        for(j=0;j<145;j++)
        {
            if(!strcmp(ms,meso_name[j]))
            {
                phi1 = meso_phi1[j];
                phi2 = meso_phi2[j];
                psi1 = meso_psi1[j];
                psi2 = meso_psi2[j];
                //printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                break;
            }
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
    //printf("Ome %lf Ome 1 %lf Ome2 %lf\n",ome,175.0,185.0);
    if(ome > 180.0)
        ome = ome - 360.0;
    
      
    p->atoms[phi_atoms[i]].tp_torsion = phi;
    p->atoms[phi_atoms[i]].sp_angle = uniform(108.0, 111.0,s);
    p->atoms[psi_atoms[i]].tp_torsion = psi;
    p->atoms[ome_atoms[i]].tp_torsion = ome;
    //printf("%s %s %d Phi %lf Psi 1 %lf Ome %lf\n",p->atoms[i].name,p->res_names[i],p->atoms[i].resnum,phi,psi,ome);
       
    
    for(k=0;k<p->num_chi_atoms;k++)
    {
        if(p->atoms[chi_atoms[k]].resnum == i)
        {
            res_ind = get_chi_res_index(p->res_names[i]);
            atom_num =0;
            if(res_ind > -1)
            {
                ran_ind = choice(chi_size[res_ind],s);
                pos = (ran_ind * chi_size[res_ind]) + chi_start[res_ind]+(2*atom_num);
                chi = uniform(chi_value[pos],chi_value[pos+1],s);
                atom_num+=1;
                if(chi > 180.0)
                    chi = chi -360.0;
                else if (chi < -180.0)
                    chi = chi + 360.0;
                p->atoms[chi_atoms[k]].tp_torsion = chi; 
                printf("Chi %lf\n",chi);    
                
            }    
        }
    }
  
   
}
