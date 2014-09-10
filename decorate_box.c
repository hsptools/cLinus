#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "linus.h"
#include "Chidict.h"

char buff[BUFSIZ];

void decorate_box(char *pdbfile, char *allowed, char *out)
{
    
    struct timeval tstart; 
    // Get start time and print
    int rc;
    rc = gettimeofday(&tstart, NULL);
    if (rc == 0)
    {
        //printf("gettimeofday() successful.\n");
        //printf("Start %u.%06u\n", tstart.tv_sec, tstart.tv_usec);
    }
    
    int i;
    
    LinusProtein prot;
    strcpy(prot.filename, pdbfile);
    protein_from_pdb(pdbfile,&prot);
    //printf("Parse protein \n");  
    int num_res = prot.num_res;
    Hbdpar hbparms;
    default_hbdpar(&hbparms);
    hbparms.use_hbond = 1;
    hbparms.use_sidechain_hbond = 1;
    
    int wmin = 1; //2
    int wmax = num_res - 1; // or 1
   
    double BETA = 2.0;

    int jx;
    /*
     for(jx=0;jx<prot.num_atoms;jx++)
     {
       printf("%s res %d  tp %lf x %lf y %lf z %lf\n",prot.atoms[jx].name,prot.atoms[jx].resnum,prot.atoms[jx].tp_torsion,prot.atoms[jx].x,prot.atoms[jx].y,prot.atoms[jx].z);
                            
                         
      }
    */

    HbList *hblist;
    int max_acp_size=0,size=0;
    int num_hb = make_hbond_list(&prot,hblist,wmin,wmax,hbparms.hbond_distance,hbparms.hbond_probe,hbparms.hbond_torsion,
               hbparms.hbond_score_short,hbparms.hbond_score_long,hbparms.sidechain_hbond_distance,hbparms.sidechain_hbond_torsion,
               hbparms.sidechain_hbond_score,hbparms.use_hbond,hbparms.use_sidechain_hbond, &size,&max_acp_size);

    hblist = malloc(num_hb*sizeof(HbList));
       
    for(i=0;i<size;i++)
    {
        hblist[i].acps = malloc(max_acp_size*sizeof(Acps));
    }   
    make_hbond_list(&prot,hblist,wmin,wmax,hbparms.hbond_distance,hbparms.hbond_probe,hbparms.hbond_torsion,
    hbparms.hbond_score_short,hbparms.hbond_score_long,hbparms.sidechain_hbond_distance,hbparms.sidechain_hbond_torsion,
    hbparms.sidechain_hbond_score,hbparms.use_hbond,hbparms.use_sidechain_hbond,&num_hb,&max_acp_size);
    
    //printf("Hbond list \n");  
    
    AtomList *vlist;
    int num_v = make_vdw_list(&prot,vlist,0);
    vlist=malloc(num_v*sizeof(AtomList));
    make_vdw_list(&prot,vlist,num_v);
    /*
    int x;
    for(x=0;x<num_v;x++)
    {
        printf("%s %d %s %d %lf %d\n",prot.atoms[vlist[x].a1].name, prot.atoms[vlist[x].a1].resnum, prot.atoms[vlist[x].a2].name, prot.atoms[vlist[x].a2].resnum,vlist[x].value, x);    
    }
    */
    
    //printf("VDW list \n");  
    char out1[300];
    strcpy(out1,out);
    //strcat(out1,allowed);
    strcat(out1,"_successful_decorations");
    
    //char out2[300];
    //strcpy(out2,out);
    //strcat(out2,allowed);
    //strcat(out2,"_successful_decorations_summary");
    
    FILE *infile = fopen(allowed, "r");
    FILE *ofile = fopen(out1, "w+");
    //FILE *ofile2 = fopen(out2, "w+");
    
    /*
    printf("Side chain rotamer set size %s %d %s %d %s %d %s %d. Total possible combinations: %d\n",
    prot.res_names[2],get_chi_res_combinations(prot.res_names[2]),
    prot.res_names[3],get_chi_res_combinations(prot.res_names[3]),
    prot.res_names[4],get_chi_res_combinations(prot.res_names[4]),
    prot.res_names[5],get_chi_res_combinations(prot.res_names[5]),
    get_chi_res_combinations(prot.res_names[2]) * get_chi_res_combinations(prot.res_names[3]) *
    get_chi_res_combinations(prot.res_names[4]) * get_chi_res_combinations(prot.res_names[5])
    );
    */
   
    int accpt = 0;
    int num =0;
    double all_energy_for_set = 0.0;
    double all_Henergy_for_set = 0.0;
    int total_attmp = 0;
    
    // for writing outfile lines in batches rather than 1 at a time (for speed)
    
    int outnum =0;
    //outlist;
    
    
    // for debugging
    int failedintra = 0;
    int bump = 0;
    int set = 1;
    int startnum = 0;
    char *line2 = "default";
    
    int ljnum = 0;
    int no_bump =0;
    
    double s[6];
    gen_seeds(s);
    
  
    int i1,i2,i3,i4;
    int j1,j2,j3,j4;
    
    int s1 = get_chi_res_start(prot.res_names[1]);
    int s2 = get_chi_res_start(prot.res_names[2]);
    int s3 = get_chi_res_start(prot.res_names[3]);
    int s4 = get_chi_res_start(prot.res_names[4]);
    
    int e1 = get_chi_res_end(prot.res_names[1]);
    int e2 = get_chi_res_end(prot.res_names[2]);
    int e3 = get_chi_res_end(prot.res_names[3]);
    int e4 = get_chi_res_end(prot.res_names[4]);
    
    int res1_size = 2 * get_chi_res_size_2(prot.res_names[1]);
    int res2_size = 2 * get_chi_res_size_2(prot.res_names[2]);
    int res3_size = 2 * get_chi_res_size_2(prot.res_names[3]);
    int res4_size = 2 * get_chi_res_size_2(prot.res_names[4]);
    
    //printf("1 %d %d %d \n",s1,e1,res1_size);
    //printf("2 %d %d %d \n",s2,e2,res2_size);
    //printf("3 %d %d %d\n",s3,e3,res3_size);
    //printf("4 %d %d %d\n",s4,e4,res4_size);
    
    double phi1,phi2,phi3,phi4,psi1,psi2,psi3,psi4,ome1,ome2,ome3,ome4,energy,DH1,DH2,DH3;
    int nhb,identnum,hexname;
    
    while(fgets(buff, sizeof buff, infile) != NULL)
    {
    
        if(sscanf(buff,"%lf %lf %lf  %lf %lf %lf  %lf %lf %lf  %lf %lf %lf  %lf  %d  %lf %lf %lf %d %d \n",
           &phi1, &psi1, &ome1, &phi2, &psi2, &ome2, &phi3, &psi3, &ome3, &phi4, &psi4, &ome4, &energy, &hexname, &DH1, &DH2, &DH3, &nhb, &identnum) == 19)
        {
            
            //printf("DH1 %lf DH2 %lf DH3 %lf nhb %d identnum %d\n",DH1,DH2,DH3,nhb,identnum);
            prot.atoms[prot.phi_atoms[1]].tp_torsion =  phi1;
            prot.atoms[prot.psi_atoms[1]].tp_torsion = set_sy(psi1);
            prot.atoms[prot.omega_atoms[1]].tp_torsion = ome1;
            
            prot.atoms[prot.phi_atoms[2]].tp_torsion =  phi2;
            prot.atoms[prot.psi_atoms[2]].tp_torsion = set_sy(psi2);
            prot.atoms[prot.omega_atoms[2]].tp_torsion = ome2;
            
            prot.atoms[prot.phi_atoms[3]].tp_torsion =  phi3;
            prot.atoms[prot.psi_atoms[3]].tp_torsion = set_sy(psi3);
            prot.atoms[prot.omega_atoms[3]].tp_torsion = ome3;
            
            prot.atoms[prot.phi_atoms[4]].tp_torsion =  phi4;
            prot.atoms[prot.psi_atoms[4]].tp_torsion = set_sy(psi4);
            prot.atoms[prot.omega_atoms[4]].tp_torsion = ome4;
            
                      
            int chi_acpt = 0;
            double total_e = 0.0;
            double total_hbonde = 0.0;
            double lowest_e = 100000.0;
          
            num  +=1;
            int curr_res =-1; 
            int curr_atm = -1;
            int hbs;
            
            for(i1= s1;i1<=e1;i1+=res1_size)
            {     
                
                for(i2= s2;i2<= e2;i2+=res2_size)
                {
                    
                    for(i3= s3;i3<= e3;i3+=res3_size)
                    {
                        
                        for(i4= s4;i4<= e4;i4+=res4_size)
                        {      
                            
                            for(j1=0;j1<prot.num_chi_atoms;j1++)
                            {
                                //printf("%s %s %d\n",prot.res_names[prot.atoms[prot.chi_atoms[j1]].resnum],prot.atoms[prot.chi_atoms[j1]].name,prot.chi_atoms[j1]);
                                if(curr_res != prot.atoms[prot.chi_atoms[j1]].resnum)
                                {
                                    curr_res = prot.atoms[prot.chi_atoms[j1]].resnum;
                                    curr_atm = 0;
                                }
                             
                                if(curr_res == 1)
                                {
                                    //printf("1 Res %s Atom %s %lf %lf\n", prot.res_names[prot.atoms[prot.chi_atoms[j1]].resnum],prot.atoms[prot.chi_atoms[j1]].name,chi_value[(curr_atm*2)+i1],chi_value[(curr_atm*2)+i1+1]);
                                    double chi = uniform(chi_value[(curr_atm*2)+i1],chi_value[(curr_atm*2)+i1+1],s);
                                    if(chi > 180.0) chi = chi - 360.0;
                                    else if(chi < - 180.0) chi = chi + 360.0;
                                    //prot.atoms[prot.chi_atoms[j1]].tp_torsion = chi;
                                    //printf ("Chi 1 %lf\n",chi);
                                    prot.atoms[prot.chi_atoms[j1]].tp_torsion = chi;
                                    curr_atm++;
                                }
                                if(curr_res == 2)
                                {
                                    //printf("2 Res %s Atom %s %lf %lf\n", prot.res_names[prot.atoms[prot.chi_atoms[j1]].resnum],prot.atoms[prot.chi_atoms[j1]].name,chi_value[(curr_atm*2)+i2],chi_value[(curr_atm*2)+i2+1]);
                                    double chi = uniform(chi_value[(curr_atm*2)+i2],chi_value[(curr_atm*2)+i2+1],s);
                                    if(chi > 180.0) chi = chi - 360.0;
                                    else if(chi < - 180.0) chi = chi + 360.0;
                                    //prot.atoms[prot.chi_atoms[j1]].tp_torsion = chi;
                                    //printf (" Chi 2 %lf\n",chi);
                                    prot.atoms[prot.chi_atoms[j1]].tp_torsion = chi;
                                    curr_atm++;
                                }
                                if(curr_res == 3)
                                {
                                    //printf("3 Res %s Atom %s %lf %lf\n", prot.res_names[prot.atoms[prot.chi_atoms[j1]].resnum],prot.atoms[prot.chi_atoms[j1]].name,chi_value[(curr_atm*2)+i3],chi_value[(curr_atm*2)+i3+1]);
                                    double chi = uniform(chi_value[(curr_atm*2)+i3],chi_value[(curr_atm*2)+i3+1],s);
                                    if(chi > 180.0) chi = chi - 360.0;
                                    else if(chi < - 180.0) chi = chi + 360.0;
                                    //prot.atoms[prot.chi_atoms[j1]].tp_torsion = chi;
                                    prot.atoms[prot.chi_atoms[j1]].tp_torsion = chi;
                                    //printf (" Chi 3 %lf\n",chi);
                                    curr_atm++;
                                }
                                if(curr_res == 4)
                                {
                                    //printf("4 Res %s Atom %s %lf %lf\n", prot.res_names[prot.atoms[prot.chi_atoms[j1]].resnum],prot.atoms[prot.chi_atoms[j1]].name,chi_value[(curr_atm*2)+i4],chi_value[(curr_atm*2)+i4+1]);
                                    double chi = uniform(chi_value[(curr_atm*2)+i4],chi_value[(curr_atm*2)+i4+1],s);
                                    if(chi > 180.0) chi = chi - 360.0;
                                    else if(chi < - 180.0) chi = chi + 360.0;
                                    //prot.atoms[prot.chi_atoms[j1]].tp_torsion = chi;
                                    prot.atoms[prot.chi_atoms[j1]].tp_torsion = chi;
                                    //printf (" Chi 4 %lf\n",chi);
                                    curr_atm++;
                                }
                               
                            }//chiatoms
                          
                            /*
                            prot.atoms[prot.chi_atoms[0]].tp_torsion =   -51.852886483;
                            prot.atoms[prot.chi_atoms[1]].tp_torsion = 72.8314042178;
                            prot.atoms[prot.chi_atoms[2]].tp_torsion =    -173.789881915;
                            prot.atoms[prot.chi_atoms[3]].tp_torsion =   -173.938292315 ;
                            prot.atoms[prot.chi_atoms[4]].tp_torsion =     66.764297861 ;
                            prot.atoms[prot.chi_atoms[5]].tp_torsion =     -70.0729974944;
                            prot.atoms[prot.chi_atoms[6]].tp_torsion = -62.3021384154 ;
                            prot.atoms[prot.chi_atoms[7]].tp_torsion = -37.6958239053;
                            */
                            total_attmp+=1;    
                            ztox(prot.atoms, 0,prot.num_atoms);
                           /*
                           for(j2=0;j2<prot.num_atoms;j2++)
                            {
                                printf("%s res %d  tp %lf x %lf y %lf z %lf\n",prot.atoms[j2].name,prot.atoms[j2].resnum,prot.atoms[j2].tp_torsion,prot.atoms[j2].x,prot.atoms[j2].y,prot.atoms[j2].z);
                            
                         
                            }*/
                         //                           
                            
                            
                            int int_ok = 1;
                            
                            for (j2=0;j2<prot.num_res;j2++)
                            {
                                if(bump_check_intra(prot.atoms,prot.res_fai,j2))
                                {
                                    int_ok = 0;
                                    break;
                                }
                            }
                                  
                            //int h = bumpcheck(prot.atoms, prot.res_fai, 1, num_res-2, 2, num_res-1);
                            
                            
                            //printf("Bump Check %d\n", h);
                            if(int_ok && !bumpcheck(prot.atoms, prot.res_fai, 1, num_res-2, 2, num_res-1))
                            {
                                //printf("Entered here\n");
                                FILE *file;
                                int numint,numvirt,bbtot,numbb;
                                double *wat_list;
                                int ctr_wat_list =0;
                                CntmultList cntmult_list[num_res];
                                int wlist = mk_virt_bb_cntmult_loosHbd(&prot,hblist,num_hb,0,0,file,0,1.25,&numint,&numvirt,&bbtot,&numbb,wat_list,0,cntmult_list);
                                wat_list = (double *) malloc(wlist*sizeof(double));
                                mk_virt_bb_cntmult_loosHbd(&prot,hblist,num_hb,0,0,file,0,1.25,&numint,&numvirt,&bbtot,&numbb,wat_list,wlist,cntmult_list);
                                //printf("Finished mk virt\n");
                      
                                hbs = numbb;
                                //printf("HBS %d\n",hbs);
                                
                                if(hbs > 7)
                                {
                                    double hbondE = hbond_score(&prot,hblist, num_hb);
                                    double solvE = calc_chasa_solv2(&prot,numint,numvirt,bbtot,numbb,wat_list,wlist,cntmult_list);
                                    //double calc_chasa_solv2(LinusProtein *p, int numint_loos, int numvirt_loos, int btot, int numbb, double *wat_list, int wlist,CntmultList *solv_list);
                                    //double calc_chasa_solv2(LinusProtein *p, int numint_loos, int numvirt_loos, int btot, int numbb, double *wat_list, int wlist,CntmultList *solv_list)
                                    
                                    double LJE = LJ_score(&prot,vlist,num_v,1.0);
                                    double energyV = hbondE + solvE + LJE;
                                    double energyH = exp(-1*hbondE*BETA); //exp(-1*hbondE*BETA);
                                    double energy = exp(-1*energyV*BETA); //exp(-1*energyV*BETA);
                                    if(lowest_e < energy) lowest_e = energy;
                                    total_hbonde += energyH;
                                    total_e += energy;
                                    accpt += 1;
                                    chi_acpt += 1;
                                    //printf("%lf %lf %lf %lf %lf %lf %d %d\n",hbondE, solvE, LJE, energyV, energyH,energy, accpt, chi_acpt);
                                }
                                
                            }  
                        
                        }//i4
                    
                    }//i3   
                
                }//i2      
            
            }//i1
            
            //all_energy_for_set += total_e;
            //all_Henergy_for_set += total_hbonde;
             
            if(chi_acpt > 0)
            {
                fprintf(ofile,"%8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %4.3f  %3d  %6.1f %6.1f %6.1f %4d %4d %13d %5d %26.3f %8.3f %26.3f\n",
                        phi1, psi1, ome1, phi2, psi2, ome2, phi3, psi3, ome3, phi4, psi4, ome4, energy, hexname, DH1, DH2, DH3,nhb,hbs,identnum, chi_acpt, total_e, total_hbonde,lowest_e);
            } 
            
        }//if
    
    } // while   
    
    
    //printf("Number of backbone struct = %d\n", num);
    //printf("Total attempts to decorate = %d\n",total_attmp);
   // printf("Total successes of decoration = %d\n",accpt);
    //printf("Sum of all energy of complete set = %lf\n", all_energy_for_set);
   // fprintf(ofile2,"num backbones: %d  \ntotal attempts to decorate: %d  \nsuccess decorations %d \nTotal Energy: %lf \nTotal Hbond Energy: %lf\n",
                   // num,total_attmp,accpt, all_energy_for_set, all_Henergy_for_set);
                    
    
    //fclose(ofile2);    
    
    struct timeval tend; 
    // Get start time and print
    rc = gettimeofday(&tend, NULL);
    if (rc == 0)
    {
        //printf("gettimeofday() successful.\n");
        //printf("End %u.%06u\n", tend.tv_sec, tend.tv_usec);
    } 
    long elapsed = (tend.tv_sec-tstart.tv_sec)*1000000LL + tend.tv_usec-tstart.tv_usec;
    double seconds = (double) elapsed/1000000.0;
    fprintf(ofile,"Elapsed %d Seconds %lf\n",elapsed,seconds);
    
    fclose(infile);
    fclose(ofile);
    
}






