#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linus.h"


char *sc_res[34] ={ "ARG","ARG","ARG","ARG","ASN","ASN","ASN","ASN","ASP","ASP","ASP","ASP","GLN","GLN","GLN","GLN","GLU","GLU","GLU","GLU","HIS","HIS","LYS","LYS","LYS","SER","SER","SER","THR","THR","THR","TYR","TYR","TRP"};
char *sc_fp[34] ={" NH1"," NH1"," NH2"," NH2"," OD1"," OD1"," ND2"," ND2"," OD1"," OD1"," OD2"," OD2"," OE1"," OE1"," NE2"," NE2"," OE1"," OE1"," OE2"," OE2"," ND1"," NE2"," NZ "," NZ "," NZ "," OG "," OG "," OG "," OG1"," OG1"," OG1"," OH "," OH "," NE1"};
char *sc_sp[34]={" CZ ", " CZ ", " CZ ", " CZ ", " CG ", " CG ", " CG ", " CG ", " CG ", " CG ", " CG ", " CG ", " CD ", " CD ", " CD ", " CD ", " CD ", " CD ", " CD ", " CD ", " CG ", " CD2", " CE ", " CE ", " CE ", " CB ", " CB ", " CB ", " CB ", " CB ", " CB ", " CZ ", " CZ ", " CD1"};
char *sc_tp[34] = {" NE ", " NE ", " NE ", " NE ", " CB ", " CB ", " CB ", " CB ", " CB ", " CB ", " CB ", " CB ", " CG ", " CG ", " CG ", " CG ", " CG ", " CG ", " CG ", " CG ", " CE1", " CE1", " CD ", " CD ", " CD ", " CA ", " CA ", " CA ", " CA ", " CA ", " CA ", " CE1", " CE1", " CG "};
double fpd[34] = {2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8,2.8};
double spa[34] = {120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,120,111.58};
double tpt[34] ={180,0,180,0,180,0,180,0,180,0,180,0,180,0,180,0,180,0,180,0,180,180,180,60,-60,180,60,-60,180,60,-60,180,0,180};


int add_waters(LinusProtein *p)
{
    int i,j,num_wat=0;
    for(i=0;i<p->num_res;i++)
    {
        for(j=0;j<34;j++)
        {
            if(!strcmp(p->res_names[i],sc_res[j]))
                num_wat++;
        }
    }
    
    p->waters = malloc(num_wat * sizeof(Atom3d));
    
    /*
    for(i=0;i<num_wat;i++)
    {
        waters[i].first_parent = malloc(sizeof(Atom3d));
        waters[i].second_parent = malloc(sizeof(Atom3d));
        waters[i].third_parent = malloc(sizeof(Atom3d));
    }
    
    */
    //printf("Done allocating\n");       
    int ctr=0;
    for(i=0;i<p->num_res;i++)
    {
        for(j=0;j<34;j++)
        {
            if(!strcmp(p->res_names[i],sc_res[j]))
            {
                //printf("%s %s %s %s %s\n",p->res_names[i],sc_res[j],sc_fp[j],sc_sp[j],sc_tp[j]);
                //get_atom_with_name(waters[ctr].first_parent,p->atoms,p->res_fai[i],p->res_fai[i+1],sc_fp[j]); 
                //get_atom_with_name(waters[ctr].second_parent,p->atoms,p->res_fai[i],p->res_fai[i+1],sc_sp[j]); 
                //get_atom_with_name(waters[ctr].third_parent,p->atoms,p->res_fai[i],p->res_fai[i+1],sc_tp[j]); 
                p->waters[ctr].fp_ind = get_atom_ind_with_name(p->atoms,p->res_fai[i],p->res_fai[i+1],sc_fp[j]);
                //printf("Fp ind %d %d %d\n", p->waters[ctr].fp_ind,p->res_fai[i],p->res_fai[i+1]);
                p->waters[ctr].sp_ind = get_atom_ind_with_name(p->atoms,p->res_fai[i],p->res_fai[i+1],sc_sp[j]);
                //printf("Sp ind %d %d %d \n", p->waters[ctr].sp_ind,p->res_fai[i],p->res_fai[i+1]);
                p->waters[ctr].tp_ind = get_atom_ind_with_name(p->atoms,p->res_fai[i],p->res_fai[i+1],sc_tp[j]);
                //printf("Tp ind %d %d %d\n", p->waters[ctr].tp_ind,p->res_fai[i],p->res_fai[i+1]);
                p->waters[ctr].fp_distance = fpd[j];
                p->waters[ctr].sp_angle = spa[j];
                p->waters[ctr].tp_torsion = tpt[j];
                p->waters[ctr].resnum = i;
                ctr++;
            }    
        }
    }
    //printf("Done %d\n",num_wat);
    ztox_waters(p->waters,p->atoms,0,num_wat);
    
    
return num_wat;
}

