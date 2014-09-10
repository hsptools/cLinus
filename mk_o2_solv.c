#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "Atom3d.h"
#include "linus.h"


int mk_o2_solv(LinusProtein *p, int start, int end, Atom3d * acp, Atom3d *sp, Atom3d *tp, 
                int *numvirt, int tryall, FILE *watpdb, int watpdb_write,double ext_rad, double *wat_coord)
{
    //printf("entered o2\n");
    int i;
    int ctr =0; // keep count of water coords
        
    //coords of the first parent, second and third   
    double coordsfp[3];
    get_coords(acp,coordsfp);
        
    double coordssp[3];
    get_coords(sp,coordssp);
        
    double coordstp[3];
    get_coords(tp,coordstp);
    
    //define the distance, angle, and torsion for the 180 deg "on top"
    double ozmvalu[3]={2.95,180.0,0.0};
    
    double owat2_crds[3];
    _ztox(coordsfp,coordssp,coordstp,ozmvalu,owat2_crds);
    //printf("Owat2 %lf %lf %lf\n",owat2_crds[0],owat2_crds[1],owat2_crds[2]);
    
    // Check for bumps of virt water with other atoms
    
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
            if(_close(owat2_crds,atm_crds,hsd))
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
            fprintf(watpdb,"ATOM    900  O   HOH   921    %8.3f%8.3f%8.3f%6.2f%6.2f\n",owat2_crds[0],owat2_crds[1],owat2_crds[2],0.0,50.0);
        }
        if(gotone == 0)
        {
            wat_coord[ctr]=owat2_crds[0];
            ctr++;
            wat_coord[ctr]=owat2_crds[1];
            ctr++;
            wat_coord[ctr]=owat2_crds[2];
            ctr++;
            gotone =1;
           
           gotone =1;
        }
    }
    else if(watpdb_write && tryall)
    {
                fprintf(watpdb,"ATOM    900  O   HOH   921    %8.3f%8.3f%8.3f%6.2f%6.2f\n",owat2_crds[0],owat2_crds[1],owat2_crds[2],0.0,0.0);
    }
    
    
    
    // define the distance, angle, and torsion for the "lower side" sol2
    double ozmvalu1[3]={2.95,130.0,90.0};
    
    double owat1_crds[3];
    _ztox(coordsfp,coordssp,coordstp,ozmvalu1,owat1_crds);
    //printf("Owat1 %lf %lf %lf\n",owat1_crds[0],owat1_crds[1],owat1_crds[2]);
    
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
            if(_close(owat1_crds,atm_crds,hsd))
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
            fprintf(watpdb,"ATOM    900  O   HOH   922    %8.3f%8.3f%8.3f%6.2f%6.2f\n",owat1_crds[0],owat1_crds[1],owat1_crds[2],0.0,50.0);
        }
        if(gotone == 0)
        {
           wat_coord[ctr]=owat1_crds[0];
            ctr++;
            wat_coord[ctr]=owat1_crds[1];
            ctr++;
            wat_coord[ctr]=owat1_crds[2];
            ctr++;
            gotone =1;
        }
    }
    else if(watpdb_write && tryall)
    {
                fprintf(watpdb,"ATOM    900  O   HOH   922    %8.3f%8.3f%8.3f%6.2f%6.2f\n",owat1_crds[0],owat1_crds[1],owat1_crds[2],0.0,0.0);
    }
    
      // define the distance, angle, and torsion for the "upper side" sol3
    double ozmvalu2[3]={2.95,130.0,-90.0};
    
    double owat3_crds[3];
    _ztox(coordsfp,coordssp,coordstp,ozmvalu2,owat3_crds);
    //printf("Owat3 %lf %lf %lf\n",owat3_crds[0],owat3_crds[1],owat3_crds[2]);
    
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
            if(_close(owat3_crds,atm_crds,hsd))
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
            fprintf(watpdb,"ATOM    900  O   HOH   923    %8.3f%8.3f%8.3f%6.2f%6.2f\n",owat3_crds[0],owat3_crds[1],owat3_crds[2],0.0,50.0);
        }
        if(gotone == 0)
        {
            wat_coord[ctr]=owat3_crds[0];
            ctr++;
            wat_coord[ctr]=owat3_crds[1];
            ctr++;
            wat_coord[ctr]=owat3_crds[2];
            ctr++;
            gotone =1;
        }
    }
    else if(watpdb_write && tryall)
    {
                fprintf(watpdb,"ATOM    900  O   HOH   923    %8.3f%8.3f%8.3f%6.2f%6.2f\n",owat3_crds[0],owat3_crds[1],owat3_crds[2],0.0,0.0);
    }
    
      // define the distance, angle, and torsion for the "near lone pair" sol4
    double ozmvalu3[3]={2.95,130.0,0.0};
    
    double owat4_crds[3];
    _ztox(coordsfp,coordssp,coordstp,ozmvalu3,owat4_crds);
    //printf("Owat4 %lf %lf %lf\n",owat4_crds[0],owat4_crds[1],owat4_crds[2]);
    
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
            if(_close(owat4_crds,atm_crds,hsd))
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
            fprintf(watpdb,"ATOM    900  O   HOH   924    %8.3f%8.3f%8.3f%6.2f%6.2f\n",owat4_crds[0],owat4_crds[1],owat4_crds[2],0.0,50.0);
        }
        if(gotone == 0)
        {
            wat_coord[ctr]=owat4_crds[0];
            ctr++;
            wat_coord[ctr]=owat4_crds[1];
            ctr++;
            wat_coord[ctr]=owat4_crds[2];
            ctr++;
            gotone =1;
        }
    }
    else if(watpdb_write && tryall)
    {
                fprintf(watpdb,"ATOM    900  O   HOH   924    %8.3f%8.3f%8.3f%6.2f%6.2f\n",owat4_crds[0],owat4_crds[1],owat4_crds[2],0.0,0.0);
    }
    
         // define the distance, angle, and torsion for the "far lone pair" sol5
    double ozmvalu4[3]={2.95,-130.0,0.0};
    
    double owat5_crds[3];
    _ztox(coordsfp,coordssp,coordstp,ozmvalu4,owat5_crds);
    //printf("owat5 %lf %lf %lf\n",owat5_crds[0],owat5_crds[1],owat5_crds[2]);
    
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
            if(_close(owat5_crds,atm_crds,hsd))
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
            fprintf(watpdb,"ATOM    900  O   HOH   925    %8.3f%8.3f%8.3f%6.2f%6.2f\n",owat5_crds[0],owat5_crds[1],owat5_crds[2],0.0,50.0);
        }
        if(gotone == 0)
        {
            wat_coord[ctr]=owat5_crds[0];
            ctr++;
            wat_coord[ctr]=owat5_crds[1];
            ctr++;
            wat_coord[ctr]=owat5_crds[2];
            ctr++;
            gotone =1;
        }
    }
    else if(watpdb_write && tryall)
    {
                fprintf(watpdb,"ATOM    900  O   HOH   925    %8.3f%8.3f%8.3f%6.2f%6.2f\n",owat5_crds[0],owat5_crds[1],owat5_crds[2],0.0,0.0);
    }
    
    return ctr;
}

