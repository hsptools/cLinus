#include <stdio.h>
#include <string.h>
#include "linusSetConf.h"
#include "linus.h"
#include "Chidict.h"
#include "meso.h"

int local_move(LinusProtein *p, int start, int end, double *s,double *rmsd)
{

     // res is the middle of three residues to be moved in typical Linus speak
    // the psi and phi of res-1 (start) is randomly moved by 'vary' amount below
    // then psi/phi/ome of start,start+1,start+2,start+3 are minimized to make the CA
    // of start+4,start+5,start+6,start+7 have an rmsd of less than ftol compared
    // to what they were before the psi/phi move of res-1

    // resnums passed are:
    // start = res - 1
    // end = res + 2
    // and we want res-1,res,res+1,res+2
    
    
    int res = start+1;
    Atom3d *atoms = p->atoms;
    
//   print protein.res_pdb_number[start], protein.res_pdb_number[end]

    // get the first atom indices of within moveset
    // for restoring below
    
    int *fai = p->res_fai;
    int fai_start,fai_end;
    
    if(fai[start] > 0) 
        fai_start = fai[start];
    else
        fai_start = 0;
    
    // here want fai of residue after moveset as downstream bracket
    
    if(fai[start+4] < p->num_atoms) 
        fai_end = fai[start+4];
    else
        fai_start = p->num_atoms;
    
    
    // get CA coords of downstream residues before change
    // these are the target coords for rmsd calculation
    Atom3d t1 = atoms[fai[start+4]];
    double e[3] = {t1.x,t1.y,t1.z};
    Atom3d t2 = atoms[fai[start+5]]; 
    double f[3] = {t2.x,t2.y,t2.z};
    Atom3d t3 = atoms[fai[start+6]];
    double g[3] = {t3.x,t3.y,t3.z};
    Atom3d t4 = atoms[fai[start+7]];
    double h[3] = {t4.x,t4.y,t4.z};
    //printf("%s %s %s %s\n",t1.name,t2.name,t3.name,t4.name);
                                                                             
    // try to find valid angle set "finish" number of attempts
    int count = 0;                                                                
    int finish = 10000;                                                             
    int success = 0;

    
    // vary psi/phi angle a random amount +/- vary                        
    double vary = 30.0;                                                           
    // vary omega angle a random amount +/- varyo                        
    double varyo = 5.0;                                                          
              
    double posneg;
              
    while (count < finish)
    {    
    
    //===========================================================================================================================
        //get psi, phi, omega values for residue 0 (start) of moveset
        
        
        // randomly pick plus or minus for varying the psi angles from current
        if(uniform(0.0,1.0,s) > 0.5)                                                   
            posneg = 1.0;                                                     
        else                                                               
            posneg = -1.0;                                                    
                                                                             
        // first vary psi                                                     
        double psi0 = p->atoms[p->psi_atoms[start]].tp_torsion + (posneg * uniform(0.0,vary,s));                                  

        // make correction in case it went out of bounds
        if (psi0 > 180.0)                                                        
            psi0 = psi0 - 360.0;                                                  
        else if (psi0 < -180.0)                                                     
            psi0 = psi0 + 360.0;                                                  

        // change internal coords
        p->atoms[p->psi_atoms[start]].tp_torsion = psi0;
                                                                             
        
        
        
        // randomly pick plus or minus for varying the phi angles from current
        if (uniform(0.0,1.0,s) > 0.5)                                                   
            posneg = 1.0;                                                     
        else                                                               
            posneg = -1.0;                                                    
                                                                             
        double phi0 = p->atoms[p->phi_atoms[start]].tp_torsion + (posneg * uniform(0.0,vary,s));                                  

        // make correction in case it went out of bounds
        if (phi0 > 180.0)
            phi0 = phi0 - 360.0;
        else if (phi0 < -180.0)
            phi0 = phi0 + 360.0;

        // change internal coords
        p->atoms[p->phi_atoms[start]].tp_torsion = phi0; 
        
        
        
        // randomly pick plus or minus for varying the omega angles from current
        if (uniform(0.0,1.0,s) > 0.5)
            posneg = 1.0;
        else
            posneg = -1.0;

        // omega_atoms refer to z-matrix parents so add 1
        double ome0 = p->atoms[p->omega_atoms[start+1]].tp_torsion + (posneg * uniform(0.0,varyo,s));
                                                                             
        // change internal coords
        p->atoms[p->omega_atoms[start+1]].tp_torsion = ome0;
    
    
    //===========================================================================================================================
        //get psi, phi, omega values for moveset index1 (start+1) of moveset
        
        // randomly pick plus or minus for varying the psi angles from current
        if(uniform(0.0,1.0,s) > 0.5)                                                   
            posneg = 1.0;                                                     
        else                                                               
            posneg = -1.0;                                                    
                                                                             
        // first vary psi                                                     
        double psi1 = p->atoms[p->psi_atoms[start+1]].tp_torsion + (posneg * uniform(0.0,vary,s));                                  

        // make correction in case it went out of bounds
        if (psi1 > 180.0)                                                        
            psi1 = psi1 - 360.0;                                                  
        else if (psi1 < -180.0)                                                     
            psi1 = psi1 + 360.0;                                                  

        // change internal coords
        p->atoms[p->psi_atoms[start+1]].tp_torsion = psi1;
                                                                             
        
        
        
        // randomly pick plus or minus for varying the phi angles from current
        if (uniform(0.0,1.0,s) > 0.5)                                                   
            posneg = 1.0;                                                     
        else                                                               
            posneg = -1.0;                                                    
                                                                             
        double phi1 = p->atoms[p->phi_atoms[start+1]].tp_torsion + (posneg * uniform(0.0,vary,s));                                  

        // make correction in case it went out of bounds
        if (phi1 > 180.0)
            phi1 = phi1 - 360.0;
        else if (phi1 < -180.0)
            phi1 = phi1 + 360.0;

        // change internal coords
        p->atoms[p->phi_atoms[start+1]].tp_torsion = phi1; 
        
        
        
        // randomly pick plus or minus for varying the omega angles from current
        if (uniform(0.0,1.0,s) > 0.5)
            posneg = 1.0;
        else
            posneg = -1.0;

        // omega_atoms refer to z-matrix parents so add 1
        double ome1 = p->atoms[p->omega_atoms[start+2]].tp_torsion + (posneg * uniform(0.0,varyo,s));
                                                                             
        // change internal coords
        p->atoms[p->omega_atoms[start+2]].tp_torsion = ome1;
    
    //===========================================================================================================================
        //get psi, phi, omega values for moveset index1 (start+2) of moveset
        
               // randomly pick plus or minus for varying the psi angles from current
        if(uniform(0.0,1.0,s) > 0.5)                                                   
            posneg = 1.0;                                                     
        else                                                               
            posneg = -1.0;                                                    
                                                                             
        // first vary psi                                                     
        double psi2 = p->atoms[p->psi_atoms[start+2]].tp_torsion + (posneg * uniform(0.0,vary,s));                                  

        // make correction in case it went out of bounds
        if (psi2 > 180.0)                                                        
            psi2 = psi2 - 360.0;                                                  
        else if (psi2 < -180.0)                                                     
            psi2 = psi2 + 360.0;                                                  

        // change internal coords
        p->atoms[p->psi_atoms[start+2]].tp_torsion = psi2;
                                                                             
        
        
        
        // randomly pick plus or minus for varying the phi angles from current
        if (uniform(0.0,1.0,s) > 0.5)                                                   
            posneg = 1.0;                                                     
        else                                                               
            posneg = -1.0;                                                    
                                                                             
        double phi2 = p->atoms[p->phi_atoms[start+2]].tp_torsion + (posneg * uniform(0.0,vary,s));                                  

        // make correction in case it went out of bounds
        if (phi2 > 180.0)
            phi2 = phi2 - 360.0;
        else if (phi2 < -180.0)
            phi2 = phi2 + 360.0;

        // change internal coords
        p->atoms[p->phi_atoms[start+2]].tp_torsion = phi2; 
        
        
        
        // randomly pick plus or minus for varying the omega angles from current
        if (uniform(0.0,1.0,s) > 0.5)
            posneg = 1.0;
        else
            posneg = -1.0;

        // omega_atoms refer to z-matrix parents so add 1
        double ome2 = p->atoms[p->omega_atoms[start+3]].tp_torsion + (posneg * uniform(0.0,varyo,s));
                                                                             
        // change internal coords
        p->atoms[p->omega_atoms[start+3]].tp_torsion = ome2;
    
        
        //===========================================================================================================================
        //get psi, phi, omega values for moveset index1 (start+3) of moveset
        
        // randomly pick plus or minus for varying the psi angles from current
        if(uniform(0.0,1.0,s) > 0.5)                                                   
            posneg = 1.0;                                                     
        else                                                               
            posneg = -1.0;                                                    
                                                                             
        // first vary psi                                                     
        double psi3 = p->atoms[p->psi_atoms[start+3]].tp_torsion + (posneg * uniform(0.0,vary,s));                                  

        // make correction in case it went out of bounds
        if (psi3 > 180.0)                                                        
            psi3 = psi3 - 360.0;                                                  
        else if (psi3 < -180.0)                                                     
            psi3 = psi3 + 360.0;                                                  

        // change internal coords
        p->atoms[p->psi_atoms[start+3]].tp_torsion = psi3;
                                                                             
        
        
        
        // randomly pick plus or minus for varying the phi angles from current
        if (uniform(0.0,1.0,s) > 0.5)                                                   
            posneg = 1.0;                                                     
        else                                                               
            posneg = -1.0;                                                    
                                                                             
        double phi3 = p->atoms[p->phi_atoms[start+3]].tp_torsion + (posneg * uniform(0.0,vary,s));                                  

        // make correction in case it went out of bounds
        if (phi3 > 180.0)
            phi3 = phi3 - 360.0;
        else if (phi3 < -180.0)
            phi3 = phi3 + 360.0;

        // change internal coords
        p->atoms[p->phi_atoms[start+3]].tp_torsion = phi3; 
        
        
        
        // randomly pick plus or minus for varying the omega angles from current
        if (uniform(0.0,1.0,s) > 0.5)
            posneg = 1.0;
        else
            posneg = -1.0;

        // omega_atoms refer to z-matrix parents so add 1
        double ome3 = p->atoms[p->omega_atoms[start+4]].tp_torsion + (posneg * uniform(0.0,varyo,s));
                                                                             
        // change internal coords
        p->atoms[p->omega_atoms[start+4]].tp_torsion = ome3; 
        
        
        
        
        // update the coords after above changes
        
        ztox(atoms,0,p->num_atoms);

        // pass the 'CA' atoms of start+4,start+5,start+6,start+7       
        // for comparison with original postions
                                 
        Atom3d a = atoms[fai[start+4]];                                      
        Atom3d b = atoms[fai[start+5]];                                     
        Atom3d c = atoms[fai[start+6]];                                      
        Atom3d d = atoms[fai[start+7]];                                      
        
        //printf("%s %s %s %s\n",a.name,b.name,c.name,d.name);
        
        // calculate rmsd
        
        *rmsd = res_rmsd(&a,&b,&c,&d,e,f,g,h);  
        //printf("Local rmsd %lf\n",*rmsd);        
                                                                             
        double ftol = 0.5;

        // if rmsd is below tolerance exit
        if (*rmsd <= ftol)                                                   
        {    
            count = count + 1;
            success = 1;
            break;
        }        
                                                                             
        // otherwise increment counter and try again
        count = count + 1;                                                    
        restore_coords_range(atoms,0,p->num_atoms);
        
    }
    if(success)
    {
        return 1;
    }
    else
    {
        restore_coords_range(atoms,0,p->num_atoms);
        ztox(atoms,0,p->num_atoms);
        return 0;    
    }
        
}

