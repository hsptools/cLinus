#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "linus.h"

// char *protlist[54] = {
//"VVME","VMEL","MELE","ELED","LEDC","EDCA","DVAF","VAFK","DADY","ADYA","DYAL","YALL","VAKK","AKKD","KKDV","KDVK","DVKV","VKVL","DRIR","RIRR","IRRM","RRMT","DWVS","WVSM","LEIN","EINE","INEF", 
//"VMME","MMEI","MEID","EIDD","IDDC","DDCA","MTAF","TAFK","DADV","ADVA","DVAL","VALL","VASR","ASRN","SRNI","RNIK","NIKV","IKVL","SSIE","SIEK","IEKL","EKLF","KWTT","WTTM","LSID","SIDA","IDAF"};
               

void find_native_like( double nat_phi1, double nat_psi1, double nat_ome1,
                       double nat_phi2, double nat_psi2, double nat_ome2,
                       double nat_phi3, double nat_psi3, double nat_ome3,
                       double nat_phi4, double nat_psi4, double nat_ome4, 
                       double rangep,char *allowed, char *prefix)
{
    if(rangep<0) rangep = 30.0;
    
    char *fmt = "%8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f  %3d  %6.1f %6.1f %6.1f %4d %4d %13d %5d %8.3f %8.3f %8.3f \n";
   
    char out1[200];
    strcpy(out1,prefix); 
    strcat(out1,"_nativelike");
    FILE *ofile = fopen(out1,"w+");
    
    /*
    char out2[200];
    strcpy(out2,prefix); 
    strcat(out2,"identmatches");
    FILE *ofile2 = fopen(out2,"w+");
    
    
    char out3[200];
    strcpy(out3,name); 
    strcat(out3,"_native_box");
    FILE *ofile3 = fopen(out3,"w+");
    */
    
    //fprintf(ofile3,fmt12,nat_phi1, nat_psi1, nat_ome1, nat_phi2, nat_psi2, nat_ome2, nat_phi3, nat_psi3, nat_ome3, nat_phi4, nat_psi4, nat_ome4, 0, 0, 0, 0, 0);
    
    int num_found = 0;
    int total_struct = 0;
    //printf("Matching: %s %lf %lf %lf %lf %lf %lf %lf %lf\n",name,nat_phi1,nat_psi1,nat_phi2,nat_psi2,nat_phi3,nat_psi3,nat_phi4,nat_psi4);
    
    
    double phi1,phi2,phi3,phi4,psi1,psi2,psi3,psi4,ome1,ome2,ome3,ome4,total_e, total_hbonde,lowest_e,energy,DH1, DH2, DH3;
    int identnum,chi_acpt,hexname,nhb1,nhb2;
    
    FILE *infile = fopen(allowed,"r");
    char buff[BUFSIZ];
    while(fgets(buff, sizeof buff, infile) != NULL)
    {
        if(sscanf(buff,"%lf %lf %lf  %lf %lf %lf  %lf %lf %lf  %lf %lf %lf  %lf  %d  %lf %lf %lf %d %d %d %d %lf %lf %lf\n",
                          &phi1, &psi1, &ome1, &phi2, &psi2, &ome2, &phi3, &psi3, &ome3, &phi4, &psi4, &ome4, &energy, &hexname, &DH1, &DH2, &DH3, &nhb1,&nhb2,&identnum, 
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
                                            fprintf(ofile,"%8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f %8.3f %8.3f  %8.3f  %3d  %6.1f %6.1f %6.1f %4d %4d %13d %5d %8.3f %8.3f %8.3f \n",
                                            phi1, psi1, ome1, phi2, psi2, ome2, phi3, psi3, ome3, phi4, psi4, ome4, energy, hexname, DH1, DH2, DH3,nhb1,nhb2,identnum, chi_acpt, total_e, total_hbonde,lowest_e);      
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
    printf("Matches found = %d out of %d\n",num_found,total_struct);
        
}               
               
