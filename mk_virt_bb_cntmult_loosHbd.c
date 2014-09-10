#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "Atom3d.h"
#include "linus.h"
//#include "CntmultList.h"
/*
""" Calculates hydrogen bond satisfaction.

Internal hydrogen bonded backbone N and O are identified using criteria
extracted from the PDB.  (See Kortemme et al., JMB, 2003) Sidechains are
included as possible donors and acceptors.

Each backbone N not hydrogen bonded will be solvated if sterically possible.
This is accomplished by attempting to place an oxygen (hard sphere radius,
ext_rad Ang.) proximate to the nitrogen at hydrogen bonding distance and
orientation. Five attempts are made at different postions within the cone
of approach. The number of successful placements is tracked.

Each backbone O whether or not it is internally hydrogen bonded is then
solvated, if possible, as described above for the N.

The module returns:

numint, numvirt, bbtot, numbb,wat_list, cntmult_list 

where,

numint         = number of internal hydrogen bonds
numvirt        = number of possible solvating waters (each backbone N
                 and O could have up to 5)
bbtot          = number of total backbone N and O (minus 1 because first
                 N is not looked at)
numbb          = number of backbone N and O that are either internally
                 hydrogen bonded or can be solvated
wat_list       = x,y,z coordinate list of solvating oxygens suitable for
                 passing to CHASA calculation
cntmult_list   = list of backbone N and O with two variables: number of
                 solvating waters (-1 if not solvated) and if an O whether
                 or not it is both hydrogen bonded and solvated. The
                 format is:

                 ([i,1,0,0],[i,2,0,0])

                 where first index is residue ID, second is N (1) or O (2), 
                 third is for number of solvating waters (or -1) and fourth
                 is flag for solvated O which is both internally hydrogen
                 bonded and solvated.
"""

  """
    Places oxygen in hydrogen bonding orientation to all donors and acceptors.
    Goes thru polypeptide and calls,
		mk_o2_solv, mk_nh_solv, mk_oh_solv, mk_amid_solv
    as appropriate.
    """
*/

int mk_virt_bb_cntmult_loosHbd(LinusProtein *p,HbList *hblist,int hblist_size, int res, int tryall,FILE *watpdb,int watpdb_write,double ext_rad,
int *numint_loos, int *numvirt_loos, int *bbtot, int *numbb, double *wat_list, int wat_list_size, CntmultList *cntmult_list)


{
    int ctr_wat_list=0;

    // res=0,tryall=0,watpdb=None,ext_rad=1.25
    if(res < 0) res = 0;
    if(tryall < 0) tryall =0;
    if(ext_rad <0.0) ext_rad = 1.25;
    if(watpdb_write < 0) watpdb_write= 0;

    int minres,maxres;
    get_res_extents(p,&minres,&maxres);
    //printf("Min %d Max %d\n",minres,maxres);

    //printf("hbtab started\n");
    AtomList *hbtab_loos;
    int hbtab_size = make_hbtab_loos(p,hblist,hblist_size,0,hbtab_loos);
    hbtab_loos=malloc(hbtab_size*sizeof(AtomList));
    make_hbtab_loos(p,hblist,hblist_size,hbtab_size,hbtab_loos);
    
    //printf("hbtab finished\n");
    int i;
    for(i=minres;i<maxres+1;i++)
    {
        cntmult_list[i].res_id_N = i;
        cntmult_list[i].is_N = 1;
        cntmult_list[i].num_N = 0;
        cntmult_list[i].flag_N = 0;
        cntmult_list[i].res_id_O = i;
        cntmult_list[i].is_O = 2;
        cntmult_list[i].num_O = 0;
        cntmult_list[i].flag_O = 0;
        //printf("i %d fai %d\n",p->res_fai[i]);
    }
    //printf("cntmultlist finished\n");
    *bbtot = 0;
    *numbb = 0;
    int numposs=0;
    *numint_loos =0;
    int numint_sc_loos =0;
    *numvirt_loos =0;
    int numvirt_sc_loos =0;

    // now see if only local score is requested

    if(res>0)
    {
        minres = res - 1;
        maxres = res + 2;
    }

    //printf("Min %d Max %d\n",minres,maxres);
    
    // backbone first

       
    for(i=minres;i<maxres;i++)
    {
        //printf("Getting start end\n");
        int start = p->res_fai[i];
        int end = p->res_fai[i+1];
        
        //printf("Getting atoms %d %d\n",start,end);
        Atom3d acp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," O  ")];
        Atom3d sp=   p->atoms[get_atom_ind_with_name(p->atoms,start,end," C  ")];
        Atom3d tp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CA ")];
        //get_atom_with_name(&acp,p->atoms,start,end," O  ");
        //get_atom_with_name(&sp,p->atoms,start,end," C  ");
        //get_atom_with_name(&tp,p->atoms,start,end," CA ");
        
        //printf("%s %s %s\n",acp.name,sp.name,tp.name);
        
        *bbtot+=1;
        numposs+=1;
        
        int hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,numint_loos);
        //printf("Is Hbond %d\n",hbond);
        
        int old_numvirt_loos = *numvirt_loos;
        
        double wcoord[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
   
             
        mk_o2_solv(p,start,end,&acp,&sp,&tp,numvirt_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
        
        //printf("\n mk_o2_solv: Numvirt is %d wcoord %lf %lf %lf\n",*numvirt_loos, wcoord[0],wcoord[1],wcoord[2]);

        
         
        if(wcoord[0] != 0)
        { 
            if(wat_list_size >0)    
            {
                wat_list[ctr_wat_list] = wcoord[0];
                ctr_wat_list++;
                wat_list[ctr_wat_list] = wcoord[1];
                ctr_wat_list++;
                wat_list[ctr_wat_list] = wcoord[2];
                ctr_wat_list++;
            }
            else
                ctr_wat_list+=3;
        }   

        if(hbond)
        {
            *numbb +=1;
            if(*numvirt_loos > old_numvirt_loos)
            {
                // if hbonded and solvent accessible just add one water
                cntmult_list[i].num_O = 1;
                //# but flag this water as special for linusScore
                cntmult_list[i].flag_O = 1;
            }
        }
        
        if(!hbond)
        {
            if(*numvirt_loos > old_numvirt_loos)
            {
                *numbb +=1;
                int num_solv = *numvirt_loos - old_numvirt_loos;
                cntmult_list[i].num_O = num_solv;
            }
            else
            {
                cntmult_list[i].num_O = -1;
            }
         }
          
        
   
   }// backbone

  
    
    
    //construct a virtual water along the NH axis at 2.95 from the N
    for(i=minres;i<maxres;i++)
    {
        if(strcmp(p->res_names[i],"PRO"))
        {
            if(!(i==0 && minres ==0))
            {
                int start = p->res_fai[i];
                int end = p->res_fai[i+1];
                int prev = p->res_fai[i-1];
                
                Atom3d donor = p->atoms[get_atom_ind_with_name(p->atoms,start,end," N  ")];
                Atom3d sp=   p->atoms[get_atom_ind_with_name(p->atoms,start,end," CA ")];
                Atom3d tp = p->atoms[get_atom_ind_with_name(p->atoms,prev,start," C  ")];
                //get_atom_with_name(&donor,p->atoms,start,end," N  ");
                //get_atom_with_name(&sp,p->atoms,start,end," CA ");
                //get_atom_with_name(&tp,p->atoms,prev,start," C  ");
                
                *bbtot+=1;
                numposs+=1;
                
              
                
                int hbond = _ishbonded_loos(p,&donor,hbtab_loos,hbtab_size,numint_loos);
                if(hbond)
                {
                    *numbb +=1;
                }
                if(!hbond)
                {
                    int old_numvirt_loos = *numvirt_loos;
                    double wcoord[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
                    
                    
                    mk_nh_solv(p,start,end,&donor,&sp,&tp,numvirt_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                    //printf("\nSize of NH is %d & Numvirt is %d\n",wcoord_size,*numvirt_loos);
                   // printf("\n mk_nh_solv: Numvirt is %d wcoord %lf %lf %lf\n",*numvirt_loos, wcoord[0],wcoord[1],wcoord[2]);
                    if(wcoord[0] != 0)
                    { 
                        if(wat_list_size >0)    
                        {
                            wat_list[ctr_wat_list] = wcoord[0];
                            ctr_wat_list++;
                            wat_list[ctr_wat_list] = wcoord[1];
                            ctr_wat_list++;
                            wat_list[ctr_wat_list] = wcoord[2];
                            ctr_wat_list++;
                        }
                        else
                            ctr_wat_list+=3;
                    }
                    if(*numvirt_loos > old_numvirt_loos)
                    {
                        *numbb +=1;
                        int num_solv = *numvirt_loos - old_numvirt_loos;
                        cntmult_list[i].num_N = num_solv;
                    }
                    else
                    {
                        cntmult_list[i].num_N = -1;
                    }                    
        
                }
                
            }
        }
    }
     

    
    
    // now for sidechains
    
    for(i=minres;i<maxres;i++)
    {
        if(!strcmp(p->res_names[i],"ARG"))
        {
            int start = p->res_fai[i];
            int end = p->res_fai[i+1];
            double wcoord[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            //construct virtual waters at 2.95 from the NH1
            
           Atom3d acp = p->atoms[get_atom_ind_with_name(&p->atoms,start,end," NH1")];
           Atom3d sp = p->atoms[get_atom_ind_with_name(&p->atoms,start,end," CZ ")];
           Atom3d tp = p->atoms[get_atom_ind_with_name(&p->atoms,start,end," NE ")];
            
            numposs+=1;
            int hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,&numint_sc_loos);
            
            if(!hbond)
            {
                  
                mk_oh_solv(p,start,end,&acp,&sp,&tp,&numvirt_sc_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                if(wcoord[0] != 0)
                { 
                    if(wat_list_size >0)    
                    {
                        wat_list[ctr_wat_list++] = wcoord[0];
                        wat_list[ctr_wat_list++] = wcoord[1];
                        wat_list[ctr_wat_list++] = wcoord[2];
                    }
                    else
                        ctr_wat_list+=3;
                }    
            }
            //construct virtual waters at 2.95 from the NH2
           acp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," NH2")];
           sp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CZ ")];
           tp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," NE ")];
            
            numposs+=1;
            hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,&numint_sc_loos);
            
            if(!hbond)
            {
                  
                mk_oh_solv(p,start,end,&acp,&sp,&tp,&numvirt_sc_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                if(wcoord[0] != 0)
                { 
                    if(wat_list_size >0)    
                    {
                        wat_list[ctr_wat_list++] = wcoord[0];
                        wat_list[ctr_wat_list++] = wcoord[1];
                        wat_list[ctr_wat_list++] = wcoord[2];
                    }
                    else
                        ctr_wat_list+=3;
                }    
            }
        } //ARG
        
        
        if(!strcmp(p->res_names[i],"GLU"))
        {
            int start = p->res_fai[i];
            int end = p->res_fai[i+1];
            double wcoord[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            //construct virtual waters at 2.95 from the O
            
           Atom3d acp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," OE1")];
           Atom3d sp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CD ")];
           Atom3d tp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CG ")];
            
            numposs+=1;
            int hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,&numint_sc_loos);
            
            if(!hbond)
            {
                  
                mk_o2_solv(p,start,end,&acp,&sp,&tp,&numvirt_sc_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                if(wcoord[0] != 0)
                { 
                    if(wat_list_size >0)    
                    {
                        wat_list[ctr_wat_list++] = wcoord[0];
                        wat_list[ctr_wat_list++] = wcoord[1];
                        wat_list[ctr_wat_list++] = wcoord[2];
                    }
                    else
                        ctr_wat_list+=3;
                }    
            }
            //construct virtual waters at 2.95 from OE2
           acp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," OE2")];
           sp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CD ")];
           tp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CG ")];
            
            numposs+=1;
            hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,&numint_sc_loos);
            
            if(!hbond)
            {
                  
                mk_o2_solv(p,start,end,&acp,&sp,&tp,&numvirt_sc_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                if(wcoord[0] != 0)
                { 
                    if(wat_list_size >0)    
                    {
                        wat_list[ctr_wat_list++] = wcoord[0];
                        wat_list[ctr_wat_list++] = wcoord[1];
                        wat_list[ctr_wat_list++] = wcoord[2];
                    }
                    else
                        ctr_wat_list+=3;
                }    
            }
        } //GLU
        
        
        if(!strcmp(p->res_names[i],"ASN"))
        {
            int start = p->res_fai[i];
            int end = p->res_fai[i+1];
            double wcoord[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            //construct virtual waters at 2.95 from the OD1
            
           Atom3d acp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," OD1")];
           Atom3d sp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CG ")];
           Atom3d tp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CB ")];
            
            numposs+=1;
            int hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,&numint_sc_loos);
            
            if(!hbond)
            {
                  
                mk_o2_solv(p,start,end,&acp,&sp,&tp,&numvirt_sc_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                if(wcoord[0] != 0)
                { 
                    if(wat_list_size >0)    
                    {
                        wat_list[ctr_wat_list++] = wcoord[0];
                        wat_list[ctr_wat_list++] = wcoord[1];
                        wat_list[ctr_wat_list++] = wcoord[2];
                    }
                    else
                        ctr_wat_list+=3;
                }    
            }
            //construct virtual waters at 2.95 from ND2
           acp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," ND2")];
           sp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CG ")];
           tp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CB ")];
            
            numposs+=1;
            hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,&numint_sc_loos);
            
            if(!hbond)
            {
                  
                mk_amid_solv(p,start,end,&acp,&sp,&tp,&numvirt_sc_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                if(wcoord[0] != 0)
                { 
                    if(wat_list_size >0)    
                    {
                        wat_list[ctr_wat_list++] = wcoord[0];
                        wat_list[ctr_wat_list++] = wcoord[1];
                        wat_list[ctr_wat_list++] = wcoord[2];
                    }
                    else
                        ctr_wat_list+=3;
                }    
            }
        } //ASN
        
        
        
        if(!strcmp(p->res_names[i],"ASP"))
        {
            int start = p->res_fai[i];
            int end = p->res_fai[i+1];
            double wcoord[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            //construct virtual waters at 2.95 from the OD1
            
           Atom3d acp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," OD1")];
           Atom3d sp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CG ")];
           Atom3d tp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CB ")];
            
            numposs+=1;
            int hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,&numint_sc_loos);
            
            if(!hbond)
            {
                  
                mk_o2_solv(p,start,end,&acp,&sp,&tp,&numvirt_sc_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                if(wcoord[0] != 0)
                { 
                    if(wat_list_size >0)    
                    {
                        wat_list[ctr_wat_list++] = wcoord[0];
                        wat_list[ctr_wat_list++] = wcoord[1];
                        wat_list[ctr_wat_list++] = wcoord[2];
                    }
                    else
                        ctr_wat_list+=3;
                }    
            }
            //construct virtual waters at 2.95 from OD2
           acp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," OD2")];
           sp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CG ")];
           tp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CB ")];
            
            numposs+=1;
            hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,&numint_sc_loos);
            
            if(!hbond)
            {
                  
                mk_o2_solv(p,start,end,&acp,&sp,&tp,&numvirt_sc_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                if(wcoord[0] != 0)
                { 
                    if(wat_list_size >0)    
                    {
                        wat_list[ctr_wat_list++] = wcoord[0];
                        wat_list[ctr_wat_list++] = wcoord[1];
                        wat_list[ctr_wat_list++] = wcoord[2];
                    }
                    else
                        ctr_wat_list+=3;
                }    
            }
        } //ASP
        
        
        if(!strcmp(p->res_names[i],"GLN"))
        {
            int start = p->res_fai[i];
            int end = p->res_fai[i+1];
            double wcoord[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            //construct virtual waters at 2.95 from the OE1
            
           Atom3d acp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," OE1")];
           Atom3d sp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CD ")];
           Atom3d tp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CG ")];
            
            numposs+=1;
            int hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,&numint_sc_loos);
            
            if(!hbond)
            {
                  
                mk_o2_solv(p,start,end,&acp,&sp,&tp,&numvirt_sc_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                if(wcoord[0] != 0)
                { 
                    if(wat_list_size >0)    
                    {
                        wat_list[ctr_wat_list++] = wcoord[0];
                        wat_list[ctr_wat_list++] = wcoord[1];
                        wat_list[ctr_wat_list++] = wcoord[2];
                    }
                    else
                        ctr_wat_list+=3;
                }    
            }
            //construct virtual waters at 2.95 from NE2
           acp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," NE2")];
           sp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CD ")];
           tp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CG ")];
            
            numposs+=1;
            hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,&numint_sc_loos);
            
            if(!hbond)
            {
                  
                mk_amid_solv(p,start,end,&acp,&sp,&tp,&numvirt_sc_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                if(wcoord[0] != 0)
                { 
                    if(wat_list_size >0)    
                    {
                        wat_list[ctr_wat_list++] = wcoord[0];
                        wat_list[ctr_wat_list++] = wcoord[1];
                        wat_list[ctr_wat_list++] = wcoord[2];
                    }
                    else
                        ctr_wat_list+=3;
                }    
            }
        } //GLN
        
        if(!strcmp(p->res_names[i],"HIS"))
        {
            int start = p->res_fai[i];
            int end = p->res_fai[i+1];
            double wcoord[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            //construct virtual waters at 2.95 from the NE2
            
           Atom3d acp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," NE2")];
           Atom3d sp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CD2")];
           Atom3d tp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CG ")];
            
            numposs+=1;
            int hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,&numint_sc_loos);
            
            if(!hbond)
            {
                  
                mk_nh_solv(p,start,end,&acp,&sp,&tp,&numvirt_sc_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                if(wcoord[0] != 0)
                { 
                    if(wat_list_size >0)    
                    {
                        wat_list[ctr_wat_list++] = wcoord[0];
                        wat_list[ctr_wat_list++] = wcoord[1];
                        wat_list[ctr_wat_list++] = wcoord[2];
                    }
                    else
                        ctr_wat_list+=3;
                }    
            }
            //construct virtual waters at 2.95 from ND1
            acp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," ND1")];
            sp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CG ")];
            tp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CB ")];
            
            numposs+=1;
            hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,&numint_sc_loos);
            
            if(!hbond)
            {
                  
                mk_nh_solv(p,start,end,&acp,&sp,&tp,&numvirt_sc_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                if(wcoord[0] != 0)
                { 
                    if(wat_list_size >0)    
                    {
                        wat_list[ctr_wat_list++] = wcoord[0];
                        wat_list[ctr_wat_list++] = wcoord[1];
                        wat_list[ctr_wat_list++] = wcoord[2];
                    }
                    else
                        ctr_wat_list+=3;
                }    
            }
        } //HIS
        
        if(!strcmp(p->res_names[i],"LYS"))
        {
            int start = p->res_fai[i];
            int end = p->res_fai[i+1];
            double wcoord[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            //construct virtual waters at 2.95 from the NZ
            
           Atom3d acp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," NZ ")];
           Atom3d sp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CE ")];
           Atom3d tp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CD ")];
            
            numposs+=1;
            int hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,&numint_sc_loos);
            
            if(!hbond)
            {
                  
                mk_oh_solv(p,start,end,&acp,&sp,&tp,&numvirt_sc_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                if(wcoord[0] != 0)
                { 
                    if(wat_list_size >0)    
                    {
                        wat_list[ctr_wat_list++] = wcoord[0];
                        wat_list[ctr_wat_list++] = wcoord[1];
                        wat_list[ctr_wat_list++] = wcoord[2];
                    }
                    else
                        ctr_wat_list+=3;
                }    
            }
            
        } //LYS
        
        if(!strcmp(p->res_names[i],"SER"))
        {
            int start = p->res_fai[i];
            int end = p->res_fai[i+1];
            double wcoord[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            //construct virtual waters at 2.95 from the OG
            
           Atom3d acp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," OG ")];
           Atom3d sp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CB ")];
           Atom3d tp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CA ")];
            
            numposs+=1;
            int hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,&numint_sc_loos);
            
            if(!hbond)
            {
                  
                mk_oh_solv(p,start,end,&acp,&sp,&tp,&numvirt_sc_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                if(wcoord[0] != 0)
                { 
                    if(wat_list_size >0)    
                    {
                        wat_list[ctr_wat_list++] = wcoord[0];
                        wat_list[ctr_wat_list++] = wcoord[1];
                        wat_list[ctr_wat_list++] = wcoord[2];
                    }
                    else
                        ctr_wat_list+=3;
                }    
            }
            
        } //SER
        
         if(!strcmp(p->res_names[i],"THR"))
        {
            int start = p->res_fai[i];
            int end = p->res_fai[i+1];
            double wcoord[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            //construct virtual waters at 2.95 from the OG1
            
           Atom3d acp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," OG1")];
           Atom3d sp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CB ")];
           Atom3d tp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CA ")];
            
            numposs+=1;
            int hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,&numint_sc_loos);
            
            if(!hbond)
            {
                  
                mk_oh_solv(p,start,end,&acp,&sp,&tp,&numvirt_sc_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                if(wcoord[0] != 0)
                { 
                    if(wat_list_size >0)    
                    {
                        wat_list[ctr_wat_list++] = wcoord[0];
                        wat_list[ctr_wat_list++] = wcoord[1];
                        wat_list[ctr_wat_list++] = wcoord[2];
                    }
                    else
                        ctr_wat_list+=3;
                }    
            }
            
        } //THR
        
         if(!strcmp(p->res_names[i],"TYR"))
        {
            int start = p->res_fai[i];
            int end = p->res_fai[i+1];
            double wcoord[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
            //construct virtual waters at 2.95 from the OH
            
           Atom3d acp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," OH ")];
           Atom3d sp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CZ ")];
           Atom3d tp = p->atoms[get_atom_ind_with_name(p->atoms,start,end," CE1")];
            
            numposs+=1;
            int hbond = _ishbonded_loos(p,&acp,hbtab_loos,hbtab_size,&numint_sc_loos);
            
            if(!hbond)
            {
                  
                mk_oh_solv(p,start,end,&acp,&sp,&tp,&numvirt_sc_loos,tryall, watpdb,watpdb_write,ext_rad,wcoord);
                if(wcoord[0] != 0)
                { 
                    if(wat_list_size >0)    
                    {
                        wat_list[ctr_wat_list++] = wcoord[0];
                        wat_list[ctr_wat_list++] = wcoord[1];
                        wat_list[ctr_wat_list++] = wcoord[2];
                    }
                    else
                        ctr_wat_list+=3;
                }    
            }
            
        } //TYR
        
    }    

    int total_loos = *numint_loos + *numvirt_loos;
    int bb_nonbd_loos = *bbtot - *numbb;
    
    //printf("numint_loos %d numvirt_loos %d bbtot %d numbb %d\n", *numint_loos, *numvirt_loos,*bbtot, *numbb);
    /*
    printf("Size of wlist %d\n",ctr_wat_list);
    if(wat_list_size >0)
    {
        int x;
        for(x=0;x<ctr_wat_list-2;x+=3)
        {
            printf(" %lf %lf %lf \n",wat_list[x],wat_list[x+1],wat_list[x+2]);
        } 
    
    }
    printf("\nCntList\n");
    
    for(i=minres;i<maxres+1;i++)
    {
        printf("%d %d %d %d %d %d %d %d\n",cntmult_list[i].res_id_N,cntmult_list[i].is_N,
        cntmult_list[i].num_N,cntmult_list[i].flag_N,cntmult_list[i].res_id_O,cntmult_list[i].is_O,
        cntmult_list[i].num_O,cntmult_list[i].flag_O);
        
    }
    */
    return ctr_wat_list;
}//mk_virt_bb_cntmult_loosHbd    

