#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "Atom3d.h"
#include "linus.h"


int mk_amid_solv(LinusProtein *p, int start, int end, Atom3d * donor, Atom3d *sp, Atom3d *tp, 
                int *numvirt, int tryall, FILE *watpdb, int watpdb_write,double ext_rad, double *wat_coord)
{
   int ctr =0; // keep count of water coords
    Atom3d **nzmindx = malloc(3*sizeof(Atom3d *));
    int i;
    for(i=0;i<3;i++)
    {
        nzmindx[i] = malloc(sizeof(Atom3d *));
    }
    nzmindx[0] = donor; // O atom
    nzmindx[1] = sp;  // C atom
    nzmindx[2] = tp;  // CA atom
   
    
    //coords of the first parent, second and third   
    double coordsfp[3];
    get_coords(donor,coordsfp);
   
    
    double coordssp[3];
    get_coords(sp,coordssp);
   
    
    double coordstp[3];
    get_coords(tp,coordstp);
   

    //define the distance, angle, and torsion for the 180 deg "on top"
    double nzmvalu[3]={2.95,180.0,180.0};
    
    double nwat2_crds[3];
    _ztox(coordsfp,coordssp,coordstp,nzmvalu,nwat2_crds);
   
    
    // heck for bumps of N hbonded virtual O with all other atoms
    
    int isbump = 0;
    int gotone = 0;
    double scale = 0.9;
    double hsd;
    for(i=0;i<p->num_atoms;i++)
    {
        if( strcmp(p->res_names[p->atoms[i].resnum],"NME") && strcmp(p->res_names[p->atoms[i].resnum],"ACE"))
        {
            double atm_crds[3];
            get_coords(&p->atoms[i],atm_crds);
            double atm_rad = scale * p->atoms[i].radius;
            if(!strcmp(p->atoms[i].name, " H  "))
            {
                hsd = 1.94;
            }
            else
            {
                hsd = (ext_rad + atm_rad);
            }
            if(_close(nwat2_crds,atm_crds,hsd))
            {
                isbump = 1;
                break;
            }
        
        }
    
    }
    
    // if no bump write out virtual water
    if(!isbump)
    {
        *numvirt += 1;
        
        if(watpdb_write && (tryall || !gotone))
        {
            fprintf(watpdb,"ATOM    900  O   HOH   911    %8.3f%8.3f%8.3f%6.2f%6.2f\n",nwat2_crds[0],nwat2_crds[1],nwat2_crds[2],0.0,50.0);
        }
        if(gotone == 0)
        {
   
            wat_coord[ctr++]=nwat2_crds[0];
            wat_coord[ctr++]=nwat2_crds[1];
            wat_coord[ctr++]=nwat2_crds[2];
            gotone =1;
        }
    }
    else if(watpdb_write && tryall)
    {
                fprintf(watpdb,"ATOM    900  O   HOH   911    %8.3f%8.3f%8.3f%6.2f%6.2f\n",nwat2_crds[0],nwat2_crds[1],nwat2_crds[2],0.0,0.0);
    }
    
    
    
   
 
    return ctr;
}

