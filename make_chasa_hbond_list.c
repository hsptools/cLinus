#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
//#include "Atom3d.h"
#include "linus.h"


/*
    'use_hbond': 0,
    'use_sidechain_hbond': 0,
    'hbond_distance': 3.0,
    'hbond_probe': 1.5,
    'hbond_score_short': 0.5,
    'hbond_score_long': 1.5,
    'hbond_torsion': 130.0,
    'sidechain_hbond_distance': 3.0,
    'sidechain_hbond_score': 0.5,
    'sidechain_hbond_torsion': 130.0,
    'hbond_winmin': 2,
    'hbond_winmax': 6,
    'use_hbs_local' : 0

*/

int make_chasa_hbond_list(LinusProtein * p, HbList * hblist,
                     int wmin, int wmax,
                     double hbdist, double hbprobe,
                     double hbtor,  double shbdist,
                     double shbtor, int *size, int *max_acp_size)
{

    double hbdmax = hbdist + hbprobe;
    double shbdmax = shbdist + hbprobe;

    wmin =2;
    wmax = p->num_res -2;
    
    int j;

    double ene = 0.0;
    
    int ctr_hb =0;
    int ctr_acp=0;

    int minres,maxres;
    get_res_extents(p,&minres,&maxres);

    int i;
    for(i=minres; i < maxres;i++)
    {
        //printf("Res %s \n",p->res_names[i]);
        int start = p->res_fai[i];
        int end  = p->res_fai[i+1];
        if(strcmp(p->res_names[i],"PRO"))
        {

        
            // NH as donor

            if(!(i==0 && minres ==0))
            {
                ctr_acp=0;
                for(j=wmin;j<wmax;j++)
                {
                    int pos;

                    pos=i-j;
                    if(pos >= minres)
                    {
                        
                            if(*size>0)
                            {
                                hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = hbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = hbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = hbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                                
                            }
                            ctr_acp++;
                        
                    }//pos >=minres

                    pos = i+j;
                    if(pos<maxres)
                    {
                      
                            if(*size>0)
                            {
                                hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = hbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = hbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = hbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;

                            }
                            ctr_acp++;
                        
                    }//pos < maxres
                }//for wmin < j <wmax   #1



                // look for sidechain acceptors with wmin = 0, so GLU can Hbond itself,etc
                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {
                                                   int scp[p->num_atoms];
                            int ctr_sc = get_scacceptor(p->atoms,p->res_fai[pos],p->res_fai[pos+1],p->res_names[pos],scp);
                            int m;
                            for(m=0;m<ctr_sc;m++)
                            {
                                //printf("SCatoms %s\n",scp[m].name);
                                if(*size>0)
                                {
                                    hblist[ctr_hb].acps[ctr_acp].acp = scp[m];
                                    hblist[ctr_hb].acps[ctr_acp].acp1 = p->atoms[scp[m]].fp_ind;
                                    hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                    hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                    hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                    hblist[ctr_hb].acps[ctr_acp].hbene = ene;

                                }
                                ctr_acp++;
                            }
                        
                    }

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                            int scp[p->num_atoms];
                            int ctr_sc = get_scacceptor(p->atoms,p->res_fai[pos],p->res_fai[pos+1],p->res_names[pos],scp);
                            int m;
                            for(m=0;m<ctr_sc;m++)
                            {
                                //printf("SCatoms %s\n",scp[m].name);
                                if(*size>0)
                                {
                                    hblist[ctr_hb].acps[ctr_acp].acp = scp[m];
                                    hblist[ctr_hb].acps[ctr_acp].acp1 = p->atoms[scp[m]].fp_ind;
                                    hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                    hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                    hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                    hblist[ctr_hb].acps[ctr_acp].hbene = ene;

                                }
                                ctr_acp++;

                            }
                        

                    }

                }//for wmin < j <wmax   #2



                if(*size>0)
                {
                
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," N  ");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CA ");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,p->res_fai[i-1],start," C  ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;

            }//i!=0, minres!=0

            //for GLU OE1 as the donor

            if(!strcmp(p->res_names[i],"GLU"))
            {
                //printf("Reached GLU OE1\n");
                ctr_acp=0;


                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {
                       
                        
                            //printf("OE1 i-j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                                hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    } //if pos i-j

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                        
                            //printf("OE1 i+j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                                hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    }//if pos i+j

                }//for

                if(*size>0)
                {
                    
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," OE1 ");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CD ");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,start,end," CG ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                    
                    //printf("Size of hblist %d is %d\n",ctr_hb,ctr_acp);
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;


            } // GLU OE1



            //for GLU OE2 as the donor

            if(!strcmp(p->res_names[i],"GLU"))
            {
                //printf("Reached GLU OE2\n");
                ctr_acp=0;


                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {
                       
                        
                            //printf("OE2 i-j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                                hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    } //if pos i-j

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                        
                            //printf("OE2 i+j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                                hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    }//if pos i+j

                }//for

                if(*size>0)
                {
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," OE2");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CD ");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,start,end," CG ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                    //printf("Size of hblist %d is %d\n",ctr_hb,ctr_acp);
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;


            } // GLU OE2



              //for ARG NH1 as the donor

            if(!strcmp(p->res_names[i],"ARG"))
            {
                //printf("Reached ARG NH1\n");
                ctr_acp=0;


                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {
                       
                        
                            //printf("OE2 i-j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                                hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    } //if pos i-j

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                        
                            //printf("NH1 i+j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                                hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    }//if pos i+j

                }//for

                if(*size>0)
                {
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," NH1");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CZ ");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,start,end," NE ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                    //printf("Size of hblist %d is %d\n",ctr_hb,ctr_acp);
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;


            } // ARG NH1


           //for ARG NH2 as the donor

            if(!strcmp(p->res_names[i],"ARG"))
            {
                //printf("Reached ARG NH2\n");
                ctr_acp=0;


                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {
                       
                        
                            //printf("OE2 i-j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                                hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    } //if pos i-j

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                        
                            //printf("NH1 i+j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    }//if pos i+j

                }//for

                if(*size>0)
                {
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," NH2");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CZ ");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,start,end," NE ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                    //printf("Size of hblist %d is %d\n",ctr_hb,ctr_acp);
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;


            } // ARG NH2

           //for ASN OD1 as the donor

            if(!strcmp(p->res_names[i],"ASN"))
            {
                //printf("Reached ARG OD1\n");
                ctr_acp=0;


                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {

                            //printf("OE2 i-j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    } //if pos i-j

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                        
                            //printf("NH1 i+j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                              hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    }//if pos i+j

                }//for

                if(*size>0)
                {
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," OD1");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CG ");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,start,end," CB ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                    ////printf("Size of hblist %d is %d\n",ctr_hb,ctr_acp);
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;


            } // ASN OD1

             //for ASN ND2/ND1 as the donor

            if(!strcmp(p->res_names[i],"ASN"))
            {
                //printf("Reached ASN ND2/1\n");
                ctr_acp=0;


                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {
                       
                        
                            //printf("OE2 i-j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    } //if pos i-j

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                        
                            //printf("NH1 i+j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                                hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    }//if pos i+j

                }//for

                if(*size>0)
                {
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," ND2");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CG ");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,start,end," CB ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                    ////printf("Size of hblist %d is %d\n",ctr_hb,ctr_acp);
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;


            } // ASN ND2/ND1


          //for ASP OD1 as the donor

            if(!strcmp(p->res_names[i],"ASP"))
            {
                //printf("Reached ARG NH1\n");
                ctr_acp=0;


                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {
                       
                        
                            //printf("OE2 i-j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    } //if pos i-j

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                        
                            //printf("NH1 i+j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                             hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    }//if pos i+j

                }//for

                if(*size>0)
                {
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," OD1");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CG ");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,start,end," CB ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                    ////printf("Size of hblist %d is %d\n",ctr_hb,ctr_acp);
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;


            } // ASP OD1

             //for ASP OD2 as the donor

            if(!strcmp(p->res_names[i],"ASP"))
            {
                //printf("Reached ASP OD2\n");
                ctr_acp=0;


                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {
                       
                        
                            //printf("OE2 i-j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    } //ifpos i-j

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                        
                            //printf("NH1 i+j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    }//if pos i+j

                }//for

                if(*size>0)
                {
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," OD2");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CG ");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,start,end," CB ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                    ////printf("Size of hblist %d is %d\n",ctr_hb,ctr_acp);
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;


            } // ASP OD2


             //for GLN OE1 as the donor

            if(!strcmp(p->res_names[i],"GLN"))
            {
                //printf("Reached GLN OE1\n");
                ctr_acp=0;


                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {
                       
                        
                            //printf("OE2 i-j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                                hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    } //if pos i-j

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                        
                            //printf("NH1 i+j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    }//if pos i+j

                }//for

                if(*size>0)
                {
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," OE1");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CD ");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,start,end," CG ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                    ////printf("Size of hblist %d is %d\n",ctr_hb,ctr_acp);
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;


            } // GLN OE1




             //for GLN NE2 as the donor

            if(!strcmp(p->res_names[i],"GLN"))
            {
                //printf("Reached GLN NE2\n");
                ctr_acp=0;


                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {
                       
                        
                            //printf("OE2 i-j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                                hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    } //if pos i-j

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                        
                           //printf("NH1 i+j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    }//if pos i+j

                }//for

                if(*size>0)
                {
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," NE2");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CD ");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,start,end," CG ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                    ////printf("Size of hblist %d is %d\n",ctr_hb,ctr_acp);
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;


            } // GLN NE2



             //for HIS NE2 as the donor

            if(!strcmp(p->res_names[i],"HIS"))
            {
                //printf("Reached HIS NE2\n");
                ctr_acp=0;


                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {
                       
                        
                            //printf("OE2 i-j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    } //if pos i-j

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                        
                            //printf("NH1 i+j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    }//if pos i+j

                }//for

                if(*size>0)
                {
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," NE2");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CD2");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,start,end," CG ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                    ////printf("Size of hblist %d is %d\n",ctr_hb,ctr_acp);
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;


            } // HIS NE2


            //for HIS ND1 as the donor

            if(!strcmp(p->res_names[i],"HIS"))
            {
                //printf("Reached HIS ND1\n");
                ctr_acp=0;


                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {
                       
                        
                            //printf("OE2 i-j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    } //if pos i-j

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                        
                            //printf("NH1 i+j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    }//if pos i+j

                }//for

                if(*size>0)
                {
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," ND1");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CG ");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,start,end," CB ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                    ////printf("Size of hblist %d is %d\n",ctr_hb,ctr_acp);
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;


            } // HIS ND1


            //for LYS NZ as the donor

            if(!strcmp(p->res_names[i],"LYS"))
            {
                //printf("Reached LYS NZ\n");
                ctr_acp=0;


                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {
                       
                        
                            //printf("OE2 i-j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    } //if pos i-j

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                        
                            //printf("NH1 i+j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    }//if pos i+j

                }//for

                if(*size>0)
                {
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," NZ ");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CE ");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,start,end," CD ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                    ////printf("Size of hblist %d is %d\n",ctr_hb,ctr_acp);
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;


            } // LYS NZ


            //for SER OG as the donor

            if(!strcmp(p->res_names[i],"SER"))
            {
                //printf("SER OG\n");
                ctr_acp=0;


                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {
                       
                        
                            //printf("OE2 i-j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                                hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    } //if pos i-j

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                        
                            //printf("NH1 i+j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    }//if pos i+j

                }//for

                if(*size>0)
                {
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," OG ");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CB ");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,start,end," CA ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                    ////printf("Size of hblist %d is %d\n",ctr_hb,ctr_acp);
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;


            } // SER OG

            //for THR OG1 as the donor

            if(!strcmp(p->res_names[i],"THR"))
            {
                //printf("THR OG1\n");
                ctr_acp=0;


                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {
                       
                        
                            //printf("OE2 i-j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                                hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    } //if pos i-j

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                        //printf("NH1 i+j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                                hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    }//if pos i+j

                }//for

                if(*size>0)
                {
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," OG1");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CB ");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,start,end," CA ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                    ////printf("Size of hblist %d is %d\n",ctr_hb,ctr_acp);
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;


            } // THR OG1

            //for TYR OH as the donor

            if(!strcmp(p->res_names[i],"TYR"))
            {
                //printf("TYR OH\n");
                ctr_acp=0;


                for(j=wmin;j<wmax;j++)
                {
                    int pos;
                    pos = i-j;
                    if(pos>=minres)
                    {
                       
                        
                            //printf("OE2 i-j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    } //if pos i-j

                    pos = i+j;
                    if(pos < maxres)
                    {
                       
                        
                            //printf("NH1 i+j ctr_hb %d ctr_acp %d\n",ctr_hb,ctr_acp);
                            if(*size>0)
                            {
                               hblist[ctr_hb].acps[ctr_acp].acp = get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," O  ");
                                hblist[ctr_hb].acps[ctr_acp].acp1 =  get_atom_ind_with_name(p->atoms,p->res_fai[pos],p->res_fai[pos+1]," C  ");
                                hblist[ctr_hb].acps[ctr_acp].hbdist = shbdist;
                                hblist[ctr_hb].acps[ctr_acp].hbdmax = shbdmax;
                                hblist[ctr_hb].acps[ctr_acp].hbtor = shbtor;
                                hblist[ctr_hb].acps[ctr_acp].hbene = ene;
                            }
                            ctr_acp++;
                        
                    }//if pos i+j

                }//for

                if(*size>0)
                {
                    hblist[ctr_hb].donor = get_atom_ind_with_name(p->atoms,start,end," OH ");
                    hblist[ctr_hb].da1 = get_atom_ind_with_name(p->atoms,start,end," CZ ");
                    hblist[ctr_hb].da2 = get_atom_ind_with_name(p->atoms,start,end," CG ");
                    hblist[ctr_hb].acps_size  = ctr_acp;
                    ////printf("Size of hblist %d is %d\n",ctr_hb,ctr_acp);
                }
                if(ctr_acp > *max_acp_size) *max_acp_size = ctr_acp;
                ctr_hb++;


            } // TYR OH




        }//res != "PRO"

    }//for minres < i < maxres

    *size = ctr_hb;

    //printf ("Final size of hblist is %d\n",*size);
    return ctr_hb;
}//make_hbond_list

