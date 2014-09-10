#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "Atom3d.h"
#include "linus.h"


int mk_nh_solv(LinusProtein *p, int start, int end, Atom3d * donor, Atom3d *sp, Atom3d *tp, 
                int *numvirt, int tryall, FILE *watpdb, int watpdb_write,double ext_rad, double *wat_coord)
{
    //printf("Entered NH\n");
  
    //coords of the first parent, second and third   
    double coordsfp[3];
    get_coords(donor,coordsfp);
        
    double coordssp[3];
    get_coords(sp,coordssp);
       
    double coordstp[3];
    get_coords(tp,coordstp);
    
    //define the distance, angle, and torsion for the 180 deg "on top"
    double nzmvalu[3]={2.95,119.0,180.0};
    
    double nwat2_crds[3];
    _ztox(coordsfp,coordssp,coordstp,nzmvalu,nwat2_crds);
        
    // check for bumps of N hbonded virtual O with all other atoms
   
     int i;
    int ctr = 0;   
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
            fprintf(watpdb,"ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n",nwat2_crds[0],nwat2_crds[1],nwat2_crds[2],0.0,50.0);
        }
        if(gotone == 0)
        {
           
            wat_coord[ctr]=nwat2_crds[0];
            ctr++;
            wat_coord[ctr]=nwat2_crds[1];
            ctr++;
            wat_coord[ctr]=nwat2_crds[2];
            ctr++;
            gotone =1;
        }
    }
    else if(watpdb_write && tryall)
    {
                fprintf(watpdb,"ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n",nwat2_crds[0],nwat2_crds[1],nwat2_crds[2],0.0,0.0);
    }
    
    
    
    // define the distance, angle, and torsion for C side
    double nzmvalu1[3]={2.95,109.0,180.0};
    
    double nwat1_crds[3];
    _ztox(coordsfp,coordssp,coordstp,nzmvalu1,nwat1_crds);
        
    //    #check for bumps of N hbonded virtual O with all other atoms
    isbump =0;
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
            if(_close(nwat1_crds,atm_crds,hsd))
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
            fprintf(watpdb,"ATOM    900  O   HOH   910   %8.3f%8.3f%8.3f%6.2f%6.2f\n",nwat1_crds[0],nwat1_crds[1],nwat1_crds[2],0.0,50.0);
        }
        if(gotone == 0)
        {
           wat_coord[ctr]=nwat1_crds[0];
            ctr++;
            wat_coord[ctr]=nwat1_crds[1];
            ctr++;
            wat_coord[ctr]=nwat1_crds[2];
            ctr++;
            gotone =1;
        }
    }
    else if(watpdb_write && tryall)
    {
                fprintf(watpdb,"ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n",nwat1_crds[0],nwat1_crds[1],nwat1_crds[2],0.0,0.0);
    }
    
      // define the distance, angle, and torsion for the "right side"
    double nzmvalu2[3]={2.95,129.0,180.0};
    
    double nwat3_crds[3];
    _ztox(coordsfp,coordssp,coordstp,nzmvalu2,nwat3_crds);
    //printf("nwat3 %lf %lf %lf\n",nwat3_crds[0],nwat3_crds[1],nwat3_crds[2]);
    
    isbump =0;
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
            if(_close(nwat3_crds,atm_crds,hsd))
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
            fprintf(watpdb,"ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n",nwat3_crds[0],nwat3_crds[1],nwat3_crds[2],0.0,50.0);
        }
        if(gotone == 0)
        {
            wat_coord[ctr]=nwat3_crds[0];
            ctr++;
            wat_coord[ctr]=nwat3_crds[1];
            ctr++;
            wat_coord[ctr]=nwat3_crds[2];
            ctr++;
            gotone =1;
        }
    }
    else if(watpdb_write && tryall)
    {
                fprintf(watpdb,"ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n",nwat3_crds[0],nwat3_crds[1],nwat3_crds[2],0.0,0.0);
    }
    
      // define the distance, angle, and torsion for the left side
    double nzmvalu3[3]={2.95,118.5,170.0};
    
    double nwat4_crds[3];
    _ztox(coordsfp,coordssp,coordstp,nzmvalu3,nwat4_crds);
    //printf("nwat4 %lf %lf %lf\n",nwat4_crds[0],nwat4_crds[1],nwat4_crds[2]);
    
    isbump =0;
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
            if(_close(nwat4_crds,atm_crds,hsd))
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
            fprintf(watpdb,"ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n",nwat4_crds[0],nwat4_crds[1],nwat4_crds[2],0.0,50.0);
        }
        if(gotone == 0)
        {
            wat_coord[ctr]=nwat4_crds[0];
            ctr++;
            wat_coord[ctr]=nwat4_crds[1];
            ctr++;
            wat_coord[ctr]=nwat4_crds[2];
            ctr++;
            gotone =1;
        }
    }
    else if(watpdb_write && tryall)
    {
                fprintf(watpdb,"ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n",nwat4_crds[0],nwat4_crds[1],nwat4_crds[2],0.0,0.0);
    }
    
         // define the distance, angle, and torsion for the "far lone pair" sol5
    double nzmvalu4[3]={2.95,118.5,-170.0};
    
    double nwat5_crds[3];
    _ztox(coordsfp,coordssp,coordstp,nzmvalu4,nwat5_crds);
    //printf("nwat5 %lf %lf %lf\n",nwat5_crds[0],nwat5_crds[1],nwat5_crds[2]);
    
    isbump =0;
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
            if(_close(nwat5_crds,atm_crds,hsd))
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
            fprintf(watpdb,"ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n",nwat5_crds[0],nwat5_crds[1],nwat5_crds[2],0.0,50.0);
        }
        if(gotone == 0)
        {
            wat_coord[ctr]=nwat5_crds[0];
            ctr++;
            wat_coord[ctr]=nwat5_crds[1];
            ctr++;
            wat_coord[ctr]=nwat5_crds[2];
            ctr++;
            gotone =1;
        }
    }
    else if(watpdb_write && tryall)
    {
                fprintf(watpdb,"ATOM    900  O   HOH   910    %8.3f%8.3f%8.3f%6.2f%6.2f\n",nwat5_crds[0],nwat5_crds[1],nwat5_crds[2],0.0,0.0);
    }
    return ctr;
}

