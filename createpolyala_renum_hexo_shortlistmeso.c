#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "linus.h"
#include "linusSecStr.h"

#define PHI_MS_SIZE 9
#define PSI_MS_SIZE 12
//#define COIL_MS_SIZE 2
#define COIL_MS_SIZE 40
#define GLY_COIL_MS_SIZE 93
#define MESO_LIST_SIZE 4
#define MESO_SIZE 145


// Command line: ./<binary> <filename> <meso1>
void create_hexo_mesoshort(char *pdbfile,char *meso1)
{
     struct timeval tstart; 
    // Get start time and print
    int rc;
    rc = gettimeofday(&tstart, NULL);
    if (rc == 0)
    {
        //printf("gettimeofday() successful.\n");
        printf("Start %u.%06u\n", tstart.tv_sec, tstart.tv_usec);
    }
    
    

   //creates poly ala box by iterating through mesostates for a hexamer: goal 10 billion structures
   char *phi_MS[PHI_MS_SIZE] = {"Aa", "Ba", "Ca", "Da", "Ea", 
                                "Fd", "Hb", "Ia", "Jd"};

   char *psi_MS[PSI_MS_SIZE] = {"Aa", "Bb", "Bc", "Bd", "Be", 
                                "Bf", "Bg", "Bh", "Bi", "Bj", 
                                "Bk", "Bl"};
    
   //char *COIL_MS[2] = {"Aa", "Aj"};
                                
   char *COIL_MS[COIL_MS_SIZE] = {"Aj", "Ak", "Al", 
                                  "Bj", "Bk", "Bl",
                                  "Ca", "Cd", "Ce", "Cf", "Cg","Cj", "Ck", "Cl", 
                                  "Da", "Dd", "De", "Df", "Dg","Di", "Dj", "Dk", "Dl", 
                                  "Ed", "Ee", "Ef", "Eg", "Ej", "Ek", "El", 
                                  "Fd", "Fe", "Ff", "Fj", "Fk", "Fl", 
                                  "Ig", "Ih", "Jf", "Jg"};  // 40 total*/

   char *GLY_COIL_MS[GLY_COIL_MS_SIZE] = {"Aa", "Ab", "Ak", "Al", "Ba", 
                                          "Bb", "Bf", "Bg", "Bh", "Bj",
                                          "Bk", "Bl", "Ca", "Cb", "Cc", 
                                          "Cd", "Ce", "Cf", "Cg", "Ch",
                                          "Cj", "Ck", "Cl", "Da", "Db", 
                                          "Dc", "Dd", "De", "Df", "Dg",
                                          "Dh", "Di", "Dj", "Dk", "Dl", 
                                          "Ea", "Eb", "Ed", "Ee", "Ef",
                                          "Eg", "Ej", "Ek", "El", "Fd", 
                                          "Fe", "Ff", "Fk", "Hb", "Hc",
                                          "Hg", "Hh", "Hi", "Ia", "Ib", 
                                          "Ic", "Ie", "If", "Ig", "Ih",
                                          "Ii", "Ik", "Il", "Ja", "Jb", 
                                          "Jc", "Jd", "Je", "Jf", "Jg",
                                          "Jh", "Ji", "Jj", "Jk", "Jl", 
                                          "Ka", "Kb", "Kc", "Ke", "Kf",
                                          "Kg", "Kh", "Ki", "Kj", "Kk", 
                                          "Kl", "La", "Lb", "Le", "Lf",
                                          "Lg", "Lk", "Ll"};

    
    
   int i, j, k, ms2, ms3, ms4;
    
   char *meso2, *meso3, *meso4;
   
   LinusProtein prot;
   strcpy(prot.filename, pdbfile);

   protein_from_pdb(pdbfile,&prot);
   
   int num_res = prot.num_res;
   int num_atom = prot.num_atoms;
   Atom3d *atoms = prot.atoms;
   
   int wmin = 1; //2
   int wmax = num_res - 1;
   Hbdpar hbparms;
   default_hbdpar(&hbparms);
   hbparms.use_hbond = 1;
   hbparms.use_sidechain_hbond = 0;
   
   
   // Make Hbond list
   
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
    
    
    
   /* 
    for(i=0;i<num_hb;i++)
    {
        for(j=0;j<hblist[i].acps_size;j++)
        {
            printf("Hblist %s %s %s %s %s\n",prot.atoms[hblist[i].donor].name,prot.atoms[hblist[i].da1].name,prot.atoms[hblist[i].da2].name,
            prot.atoms[hblist[i].acps[j].acp].name, prot.atoms[hblist[i].acps[j].acp1].name);
        }
    } 
    */
       
    int num_sterically_allowed = 0;
    int num_steric_hbs = 0;
  
    char out1[40] = "all_meso/hbs_";
    strcat(out1,meso1);
    strcat(out1,".pdb");
    //printf("Out %s\n",out);
    FILE *file1 = fopen(out1,"w+");
   
    char out2[40] = "all_meso/hexo_";
    strcat(out2,meso1);
    strcat(out2,".txt");
    //printf("Out %s\n",out);
    FILE *file2 = fopen(out2,"w+");
 
    int attmp = 0; // monte carlo attempt number
   
    double philist[MESO_LIST_SIZE+1];
    philist[0]=-140.0;
    double psilist[MESO_LIST_SIZE+1]; 
    psilist[0]=140.0;
    double omelist[MESO_LIST_SIZE+1];
    omelist[0]=180.0;
     //printf("Meso 1 %s argv %s",meso1,argv[2]);
    
    double s[6];
    gen_seeds(s);
    
    for (ms2 = 0; ms2 < COIL_MS_SIZE; ms2++) 
    {
        // Print  Monte Carlo information
        //printf("%s %d %d %d\n", COIL_MS[ms2],attmp, num_sterically_allowed, num_steric_hbs);
        
        for (ms3 = 0; ms3 < COIL_MS_SIZE; ms3++)
        {
            for (ms4 = 0; ms4 < COIL_MS_SIZE; ms4++) 
            {
                //printf("%d %d %d\n",ms2,ms3,ms4); 
                // do ten in every set of mesostates
                for (j = 0; j < 10; j++)  //change j<10
                {
                     //printf("Allocating meso\n");
                    char **meso_list = (char **)malloc(4 * sizeof(char *));
                    for (i = 0; i < 4; i++)
                    {
                        meso_list[i] = (char *)malloc(sizeof(char)*3);
                        
                    } 
                   
                    strcpy(meso_list[0],meso1); 
                    strcpy(meso_list[1],COIL_MS[ms2]);
                    strcpy(meso_list[2],COIL_MS[ms3]);
                    strcpy(meso_list[3],COIL_MS[ms4]);
                    get_meso_phi_psi(meso_list,philist,psilist,omelist,s);
                    //printf("Finished get meso\n");
                                
                    for (k = 1; k < 5; k++) 
                    {
                        atoms[prot.phi_atoms[k]].tp_torsion = philist[k];
                        atoms[prot.phi_atoms[k]].sp_angle = uniform(111.0,113.0,s);
                        atoms[prot.psi_atoms[k]].tp_torsion = psilist[k];
                        atoms[prot.omega_atoms[k]].tp_torsion = omelist[k];
                    } // end k loop
                    
                    
                    
                    ztox(atoms, 0,num_atom);
                    /*
                    printf("Coordinates\n");
                    for(i=0;i<num_atom;i++)
                    {
                        printf("%lf %lf %lf\n",atoms[i].x,atoms[i].y,atoms[i].z);
                    }*/
                    
                    attmp += 1;
                                       
                    

                    // since this is poly ala chain, no need to do bumpcheck intra
                     //bump_check(atoms, fai, 1, numres-2, 2, numres-1):
                    if (!bumpcheck(atoms, prot.res_fai, 1, num_res-2, 2, num_res-1)) 
                    {
                        //printf("No Bump %d\n",attmp);
                        //printf("attmp %d ",attmp);
                        //printf("phitp1 %lf phisp1 %lf psi1 %lf ome1 %lf ",philist[1],atoms[prot.phi_atoms[1]].sp_angle,psilist[1],omelist[1]);
                        //printf("phitp2 %lf phisp2 %lf psi2 %lf ome2 %lf ",philist[2],atoms[prot.phi_atoms[2]].sp_angle,psilist[2],omelist[2]);
                        //printf("phitp3 %lf phisp3 %lf psi3 %lf ome3 %lf ",philist[3],atoms[prot.phi_atoms[3]].sp_angle,psilist[3],omelist[3]);
                        //printf("phitp4 %lf phisp4 %lf psi4 %lf ome4 %lf ",philist[4],atoms[prot.phi_atoms[4]].sp_angle,psilist[4],omelist[4]);
                        num_sterically_allowed = num_sterically_allowed + 1;
                        //print_conf(format, &prot,file2, num_sterically_allowed);
                        //fprintf(file2,"Testing");
                        //printf("steric %d ",num_sterically_allowed);
                       
                        FILE *file;
                        int numint,numvirt,bbtot,numbb;
                        double *wat_list;
                        int ctr_wat_list =0;
                        CntmultList cntmult_list[num_res];
                        int wlist = mk_virt_bb_cntmult_loosHbd(&prot,hblist,num_hb,0,0,file,0,1.25,&numint,&numvirt,&bbtot,&numbb,wat_list,0,cntmult_list);
                        wat_list = (double *) malloc(wlist*sizeof(double));
                        mk_virt_bb_cntmult_loosHbd(&prot,hblist,num_hb,0,0,file,0,1.25,&numint,&numvirt,&bbtot,&numbb,wat_list,wlist,cntmult_list);
                        //printf("Finished mk virt\n");
                      
                        int hbs = numbb;
                        //printf("hbs %d ",hbs);
    






                       
                            num_steric_hbs+=1;
                            
                            double DH1,DH2,DH3;
                            int a[6];
                            int m;
                            for(m=0;m < num_res;m++)
                            {
                                a[m] = get_atom_ind_with_name(atoms,prot.res_fai[m],prot.res_fai[m+1]," CA ");
                            }
                            DH1 = _tor(&atoms[a[0]],&atoms[a[1]],&atoms[a[2]],&atoms[a[3]]);
                            DH2 = _tor(&atoms[a[1]],&atoms[a[2]],&atoms[a[3]],&atoms[a[4]]);
                            DH3 = _tor(&atoms[a[2]],&atoms[a[3]],&atoms[a[4]],&atoms[a[5]]);
                            
                            double phi1=prot.atoms[prot.phi_atoms[1]].tp_torsion;
                            double phi2=prot.atoms[prot.phi_atoms[2]].tp_torsion;
                            double phi3=prot.atoms[prot.phi_atoms[3]].tp_torsion;
                            double phi4=prot.atoms[prot.phi_atoms[4]].tp_torsion;
                               
                            double psi1=out_sy(prot.atoms[prot.psi_atoms[1]].tp_torsion);
                            double psi2=out_sy(prot.atoms[prot.psi_atoms[2]].tp_torsion);
                            double psi3=out_sy(prot.atoms[prot.psi_atoms[3]].tp_torsion);
                            double psi4=out_sy(prot.atoms[prot.psi_atoms[4]].tp_torsion);
                            
                            double ome1=prot.atoms[prot.omega_atoms[1]].tp_torsion;
                            double ome2=prot.atoms[prot.omega_atoms[2]].tp_torsion;
                            double ome3=prot.atoms[prot.omega_atoms[3]].tp_torsion;
                            double ome4=prot.atoms[prot.omega_atoms[4]].tp_torsion;
                            
                            int hexname = get_hexname(DH1,DH2,DH3);
   
                             char *fmt1 = "%8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %3d\n";
                            //printf(fmt,phi1, psi1, ome1, phi2, psi2, ome2, phi3, psi3, ome3,  phi4, psi4, ome4, num_steric_hbs);
                            fprintf(file2,fmt1,phi1, psi1, ome1, phi2, psi2, ome2, phi3, psi3, ome3, phi4, psi4, ome4, hexname);




                         if(hbs > 7)
                        {
                            num_steric_hbs+=1;



                         
                            double energy = 0.0;
                            
                            char *fmt = "%8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f  %3d  %6.1f %6.1f %6.1f %d %13d \n";   
                            //printf(fmt,phi1, psi1, ome1, phi2, psi2, ome2, phi3, psi3, ome3,  phi4, psi4, ome4, num_steric_hbs);    
                            fprintf(file1,fmt,phi1, psi1, ome1, phi2, psi2, ome2, phi3, psi3, ome3, phi4, psi4, ome4, energy, hexname, DH1, DH2, DH3, hbs, num_steric_hbs);
                        }
                        //printf("num_hbs_steric %d\n",num_steric_hbs);
                        
                    }
                    //else
                    //{
                        //printf("Bump %d\n",attmp);
                    //}
                   
                         
                } // end j loop
            } // end ms4 loop
        } // end ms3 loop
        
    } // end ms2 loop
    

   fclose(file1);
    //fclose(file2);
   printf("Finished attmp %d num_sterically_allowed %d num_steric_hbs %d\n",attmp,num_sterically_allowed,num_steric_hbs);
   // Get end time             
   struct timeval tend; 
    // Get start time and print
    rc = gettimeofday(&tend, NULL);
    if (rc == 0)
    {
        //printf("gettimeofday() successful.\n");
        printf("End %u.%06u\n", tend.tv_sec, tend.tv_usec);
    } 
    long elapsed = (tend.tv_sec-tstart.tv_sec)*1000000LL + tend.tv_usec-tstart.tv_usec;
    double seconds = (double) elapsed/1000000.0;
   printf("Elapsed %d Seconds %lf\n",elapsed,seconds);

}
