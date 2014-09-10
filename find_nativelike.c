#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "linus.h"

// char *protlist[54] = {
//"VVME","VMEL","MELE","ELED","LEDC","EDCA","DVAF","VAFK","DADY","ADYA","DYAL","YALL","VAKK","AKKD","KKDV","KDVK","DVKV","VKVL","DRIR","RIRR","IRRM","RRMT","DWVS","WVSM","LEIN","EINE","INEF", 
//"VMME","MMEI","MEID","EIDD","IDDC","DDCA","MTAF","TAFK","DADV","ADVA","DVAL","VALL","VASR","ASRN","SRNI","RNIK","NIKV","IKVL","SSIE","SIEK","IEKL","EKLF","KWTT","WTTM","LSID","SIDA","IDAF"};
               

void run_nativelike( char *name,double nat_phi1, double nat_psi1, double nat_ome1,
                                double nat_phi2, double nat_psi2, double nat_ome2,
                                double nat_phi3, double nat_psi3, double nat_ome3,
                                double nat_phi4, double nat_psi4, double nat_ome4, 
                                char * allowed, int rangep)
{
    
    char *fmt12 = "%11d %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f  %3d  %6.1f %6.1f %6.1f \n";
   
    char out1[200];
    strcpy(out1,name); 
    strcat(out1,"_nativelike_structures_by_phi_psi_box");
    FILE *ofile = fopen(out1,"w+");
    
    char out2[200];
    strcpy(out2,name); 
    strcat(out2,"identmatches");
    FILE *ofile2 = fopen(out2,"w+");
    
    char out3[200];
    strcpy(out3,name); 
    strcat(out3,"_native_box");
    FILE *ofile3 = fopen(out3,"w+");
    
    
    fprintf(ofile3,fmt12,nat_phi1, nat_psi1, nat_ome1, nat_phi2, nat_psi2, nat_ome2, nat_phi3, nat_psi3, nat_ome3, nat_phi4, nat_psi4, nat_ome4, 0, 0, 0, 0, 0);
    
    int num_found = 0;
    int total_struct = 0;
    printf("Matching: %s %lf %lf %lf %lf %lf %lf %lf %lf\n",name,nat_phi1,nat_psi1,nat_phi2,nat_psi2,nat_phi3,nat_psi3,nat_phi4,nat_psi4);
    
    
    double phi1,phi2,phi3,phi4,psi1,psi2,psi3,psi4,ome1,ome2,ome3,ome4,total_e, total_hbonde,lowest_e;
    int num_steric_hbs,chi_acpt;
    
    FILE *infile = fopen(allowed,"r");
    char buff[BUFSIZ];
    while(fgets(buff, sizeof buff, infile) != NULL)
    {
        if(sscanf(buff,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf\n",
                          &phi1, &psi1, &ome1, &phi2, &psi2, &ome2, &phi3, &psi3, &ome3, &phi4, &psi4, &ome4, &num_steric_hbs, 
                         &chi_acpt, &total_e, &total_hbonde,&lowest_e))
        {
            total_struct +=1;
            if(fabs(180.0-fabs(180.0-fabs(phi1-nat_phi1))) <= rangep)
            {
                if(fabs(180.0-fabs(180.0-fabs(psi1-nat_psi1))) <= rangep)
                {
                    if(fabs(180.0-fabs(180.0-fabs(phi2-nat_phi2))) <= rangep)
                    {
                        if(fabs(180.0-fabs(180.0-fabs(psi2-nat_psi2))) <= rangep)
                        {
                            if(fabs(180.0-fabs(180.0-fabs(phi3-nat_phi3))) <= rangep)
                            {
                                if(fabs(180.0-fabs(180.0-fabs(psi3-nat_psi3))) <= rangep)
                                {
                                    if(fabs(180.0-fabs(180.0-fabs(phi4-nat_phi4))) <= rangep)
                                    {
                                        if(fabs(180.0-fabs(180.0-fabs(psi4-nat_psi4))) <= rangep)
                                        {
                                            num_found+=1;
                                            fprintf(ofile,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d %lf %lf %lf\n",
                                            phi1, psi1, ome1, phi2, psi2, ome2, phi3, psi3, ome3, phi4, psi4, ome4, num_steric_hbs, 
                                            chi_acpt, total_e, total_hbonde,lowest_e);        
                                        }
                                    }
                                }
                            }  
                        }
                    }
                }
            }
        
        
        }                 
    }
    printf("For %s Matches found = %d out of %d\n",name,num_found,total_struct);
    printf("Finished with %s\n", name);
    
}               
               

void find_nativelike(char *tets_name, char *allowed)
{

    int native = 0;
    int range = 0;
    int r1 = 0;
    int r2 = 0;
    char *afilelist[10]; //allowed file list
    
   
    char *name;
    double nat_phi1, nat_psi1, nat_ome1;
    double nat_phi2, nat_psi2, nat_ome2;
    double nat_phi3, nat_psi3, nat_ome3;
    double nat_phi4, nat_psi4, nat_ome4;
    
    //takes files generated from processphipsilisttotetfile (4 letter aa code for tet, followed by phi, psis
    
    FILE *tets_file = fopen(tets_name,"r");
    
    char buff[BUFSIZ];
    
    while(fgets(buff, sizeof buff, tets_file) != NULL)
    {
    
        if(sscanf(buff,"%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d\n",
           name, &nat_phi1, &nat_psi1, &nat_ome1, &nat_phi2, &nat_psi2, &nat_ome2,&nat_phi3, &nat_psi3, &nat_ome3,&nat_phi4, &nat_psi4, &nat_ome4) == 13)
        {
            int rangep = 30;
            int i;
            for(i=0;i<54;i++)
            {
                run_nativelike(name, nat_phi1, nat_psi1, nat_ome1, 
                                     nat_phi2, nat_psi2, nat_ome2, 
                                     nat_phi3, nat_psi3, nat_ome3, 
                                     nat_phi4, nat_psi4, nat_ome4, 
                                     allowed,rangep);
            
            }
        }
    }        
    
}



