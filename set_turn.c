#include <stdio.h>
#include <string.h>
#include "linusSetConf.h"
#include "linus.h"
#include "Chidict.h"
#include "meso.h"

//set the conformation of a specified segement to helical conformation

/*
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

        For each residue in the specified segment the conformation is
        set to a strand conformation, i.e. phi in the range -150 to -120,
        psi in the range 130 to 160 and omega in the range 180 +/- 5.
        The conformation of a residue is changed if and only if it
        has a non-zero strand weight.  Chi-values for the residues
        are chosen at random.  This function modifies the Z-matrix of
        the molecule.  Cartesian coordinates are *not* calculated.


*/


void set_turn(LinusProtein *p, int res1, int res2, double *s)
{

    double phi1_1,phi2_1,psi1_1,psi2_1,phi1,psi1,ome1;
    double phi1_2,phi2_2,psi1_2,psi2_2,phi2,psi2,ome2;
    double w1,w2;
    
    char ms1[3],ms2[3];
    int index;
    char *r1name = p->res_names[res1];
    char *r2name = p->res_names[res2];
    
    int res_ind,ran_ind,pos,atom_num;
    double chi;
    int j,k;

    if(strcmp("PRO",r1name) && strcmp("PRO",r2name))
    {
         index = choice(54,s);
         strcpy(ms1,TURN_MS[index][0]); 
         strcpy(ms2,TURN_MS[index][1]);
       
        for(j=0;j<145;j++)
        {
            if(!strcmp(ms1,meso_name[j]))
            {
                phi1_1 = meso_phi1[j];
                phi2_1 = meso_phi2[j];
                psi1_1 = meso_psi1[j];
                psi2_1 = meso_psi2[j];
                ////printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                break;
            }
        }

        for(j=0;j<145;j++)
        {
            if(!strcmp(ms2,meso_name[j]))
            {
                phi1_2 = meso_phi1[j];
                phi2_2 = meso_phi2[j];
                psi1_2 = meso_psi1[j];
                psi2_2 = meso_psi2[j];
                ////printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                break;
            }
        }
        
        w1=w2= uniform(175.0,185.0,s);
        if(w1>180.0) 
            w1 = w1-360.0;
        if(w2>180.0) 
            w2 = w2-360.0;
            
    }
    
    
    else if(!strcmp("PRO",r1name))
    {
        if(!strcmp("PRO",r2name))
        {
             index = choice(4,s);
             strcpy(ms1,PRO_PRO_TURN[index][0]); 
             strcpy(ms2,PRO_PRO_TURN[index][1]);
           
            for(j=0;j<145;j++)
            {
                if(!strcmp(ms1,meso_name[j]))
                {
                    phi1_1 = meso_phi1[j];
                    phi2_1 = meso_phi2[j];
                    psi1_1 = meso_psi1[j];
                    psi2_1 = meso_psi2[j];
                    ////printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                    break;
                }
            }

            for(j=0;j<145;j++)
            {
                if(!strcmp(ms2,meso_name[j]))
                {
                    phi1_2 = meso_phi1[j];
                    phi2_2 = meso_phi2[j];
                    psi1_2 = meso_psi1[j];
                    psi2_2 = meso_psi2[j];
                    ////printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                    break;
                }
            }
            
            w1=w2= uniform(175.0,185.0,s);
            if(w1>180.0) 
                w1 = w1-360.0;
            if(w2>180.0) 
                w2 = w2-360.0;
        }    // r2 PRO
        
        
        else if(!strcmp("GLY",r2name))
        {
             index = choice(3,s);
             strcpy(ms1,PRO_GLY_TURN[index][0]); 
             strcpy(ms2,PRO_GLY_TURN[index][1]);
           
            for(j=0;j<145;j++)
            {
                if(!strcmp(ms1,meso_name[j]))
                {
                    phi1_1 = meso_phi1[j];
                    phi2_1 = meso_phi2[j];
                    psi1_1 = meso_psi1[j];
                    psi2_1 = meso_psi2[j];
                    ////printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                    break;
                }
            }

            for(j=0;j<145;j++)
            {
                if(!strcmp(ms2,meso_name[j]))
                {
                    phi1_2 = meso_phi1[j];
                    phi2_2 = meso_phi2[j];
                    psi1_2 = meso_psi1[j];
                    psi2_2 = meso_psi2[j];
                    ////printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                    break;
                }
            }
            
            w1=w2= uniform(175.0,185.0,s);
            if(w1>180.0) 
                w1 = w1-360.0;
            if(w2>180.0) 
                w2 = w2-360.0;
        }    // r2 GLY
        
        else
        {
             index = choice(15,s);
             strcpy(ms1,PRO_X_TURN[index][0]); 
             strcpy(ms2,PRO_X_TURN[index][1]);
           
            for(j=0;j<145;j++)
            {
                if(!strcmp(ms1,meso_name[j]))
                {
                    phi1_1 = meso_phi1[j];
                    phi2_1 = meso_phi2[j];
                    psi1_1 = meso_psi1[j];
                    psi2_1 = meso_psi2[j];
                    ////printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                    break;
                }
            }

            for(j=0;j<145;j++)
            {
                if(!strcmp(ms2,meso_name[j]))
                {
                    phi1_2 = meso_phi1[j];
                    phi2_2 = meso_phi2[j];
                    psi1_2 = meso_psi1[j];
                    psi2_2 = meso_psi2[j];
                    ////printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                    break;
                }
            }
            
            w1=w2= uniform(175.0,185.0,s);
            if(w1>180.0) 
                w1 = w1-360.0;
            if(w2>180.0) 
                w2 = w2-360.0;
        }    // r2 X
        
        
        
        
            
    }
    
    else if(!strcmp("PRO",r2name))
    {
        if(!strcmp("GLY",r1name))
        {
             index = choice(2,s);
             strcpy(ms1,GLY_PRO_TURN[index][0]); 
             strcpy(ms2,GLY_PRO_TURN[index][1]);
           
            for(j=0;j<145;j++)
            {
                if(!strcmp(ms1,meso_name[j]))
                {
                    phi1_1 = meso_phi1[j];
                    phi2_1 = meso_phi2[j];
                    psi1_1 = meso_psi1[j];
                    psi2_1 = meso_psi2[j];
                    ////printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                    break;
                }
            }

            for(j=0;j<145;j++)
            {
                if(!strcmp(ms2,meso_name[j]))
                {
                    phi1_2 = meso_phi1[j];
                    phi2_2 = meso_phi2[j];
                    psi1_2 = meso_psi1[j];
                    psi2_2 = meso_psi2[j];
                    ////printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                    break;
                }
            }
            
            w1=w2= uniform(175.0,185.0,s);
            if(w1>180.0) 
                w1 = w1-360.0;
            if(w2>180.0) 
                w2 = w2-360.0;
        }    // r1 GLY
        
        else
        {
             index = choice(54,s);
             strcpy(ms1,TURN_MS[index][0]); // maybe X_PRO_TURN Check
             strcpy(ms2,TURN_MS[index][1]);
           
            for(j=0;j<145;j++)
            {
                if(!strcmp(ms1,meso_name[j]))
                {
                    phi1_1 = meso_phi1[j];
                    phi2_1 = meso_phi2[j];
                    psi1_1 = meso_psi1[j];
                    psi2_1 = meso_psi2[j];
                    ////printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                    break;
                }
            }

            for(j=0;j<145;j++)
            {
                if(!strcmp(ms2,meso_name[j]))
                {
                    phi1_2 = meso_phi1[j];
                    phi2_2 = meso_phi2[j];
                    psi1_2 = meso_psi1[j];
                    psi2_2 = meso_psi2[j];
                    ////printf("%lf %lf %lf %lf\n",phi1,phi2,psi1,psi2);
                    break;
                }
            }
            
            w1=w2= uniform(175.0,185.0,s);
            if(w1>180.0) 
                w1 = w1-360.0;
            if(w2>180.0) 
                w2 = w2-360.0;
        }    // r2 PRO
    
    }
    
    phi1 = uniform(phi1_1,phi2_1,s);
    ////printf("Phi %lf Phi 1 %lf Phi2 %lf\n",phi,phi1,phi2);
    if(phi1 > 180.0)
        phi1 = phi1 - 360.0;
    else if (phi1 < - 180.0)
        phi1 = phi1 + 360.0;

    psi1 = uniform(psi1_1,psi2_1,s);
    ////printf("psi %lf psi 1 %lf psi2 %lf\n",psi,psi1,psi2);
    psi1 = psi1+180.0;
    if(psi1 > 180.0)
        psi1 = psi1 - 360.0;
    else if (psi1 < - 180.0)
        psi1 = psi1 + 360.0;

    
    phi2 = uniform(phi1_2,phi2_2,s);
    ////printf("Phi %lf Phi 1 %lf Phi2 %lf\n",phi,phi1,phi2);
    if(phi2 > 180.0)
        phi2 = phi2 - 360.0;
    else if (phi2 < - 180.0)
        phi2 = phi2 + 360.0;

    psi2 = uniform(psi1_2,psi2_2,s);
    ////printf("psi %lf psi 1 %lf psi2 %lf\n",psi,psi1,psi2);
    psi2 = psi2+180.0;
    if(psi2 > 180.0)
        psi2 = psi2 - 360.0;
    else if (psi2 < - 180.0)
        psi2 = psi2 + 360.0;


    if(p->res_turn1_wt[res1] >= 1.0e-9)
    {
        p->atoms[p->phi_atoms[res1]].tp_torsion = phi1;
        p->atoms[p->psi_atoms[res1]].tp_torsion = psi1;
        p->atoms[p->omega_atoms[res1]].tp_torsion = w1;
        p->atoms[p->phi_atoms[res1]].sp_angle = uniform(111.0, 114.0,s);
        //printf("%s %s %d Phi %lf Psi 1 %lf Ome %lf\n",p->atoms[res1].name,p->res_names[res1],p->atoms[res1].resnum,phi1,psi1,ome1);
        
        for(k=0;k<p->num_chi_atoms;k++)
        {
            if(p->atoms[p->chi_atoms[k]].resnum == res1)
            {
                res_ind = get_chi_res_index(p->res_names[res1]);
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
                    p->atoms[p->chi_atoms[k]].tp_torsion = chi; 
                    //printf("Chi %lf\n",chi);    
                    
                }    
            }
        }    
    }
    
    if(p->res_turn2_wt[res2] >= 1.0e-9)
    {
        p->atoms[p->phi_atoms[res1]].tp_torsion = phi1;
        p->atoms[p->psi_atoms[res1]].tp_torsion = psi1;
        p->atoms[p->omega_atoms[res1]].tp_torsion = w1;
        p->atoms[p->phi_atoms[res1]].sp_angle = uniform(111.0, 114.0,s);
        //printf("%s %s %d Phi %lf Psi 1 %lf Ome %lf\n",p->atoms[res2].name,p->res_names[res2],p->atoms[res2].resnum,phi2,psi2,ome2);
        
        for(k=0;k<p->num_chi_atoms;k++)
        {
            if(p->atoms[p->chi_atoms[k]].resnum== res2)
            {
                res_ind = get_chi_res_index(p->res_names[res2]);
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
                    p->atoms[p->chi_atoms[k]].tp_torsion = chi; 
                    //printf("Chi %lf\n",chi);    
                    
                }    
            }
        }    
    }

}

