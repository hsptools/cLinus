#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <regex.h>
#include "linusSecStr.h"
#include "linus.h"


void rc_ss(LinusProtein *p, char **sst)
{
    //int nres = p->num_res;
    
    int nres = 19;
    int i,j,k;
    
    
    char **codes;
    char code_string[3*nres];
    
    
    codes = (char **)malloc(nres * sizeof(char *));
    for (i = 0; i < nres; i++)
    {
        codes[i] = (char *)malloc(sizeof(char)*3);
    }
    
    rc_codes(p,codes);
    
    //                0    1    2    3     4   5    6    7    8    9    10    11  12   13   14   15   16   17   18
    //char *codes[19]={"De","De","De","De","De","De","Ib","Dg","Bj","Bk","Bl","Ek","Bj","Bj","Bj","De","De","De","De"};
    
    //printf("%s\n",codes[18]);
    // char sst[prot->num_res]; in calling function
    
    for(i=0;i<nres;i++)
    {  
        strcpy(sst[i],"C");
    }
    
    for(i=0;i<nres;i++)
    {
        if(is_PII(codes[i]))
            strcpy(sst[i],"P");        
    }
    
    for(i=0;i<nres-1;i++)
    {
        //Turn dictionary comprises of two mesostates
        char code[5];
        strcpy(code,codes[i]);
        strcat(code,codes[i+1]);
        //printf("Code %s\n",code);
        if(is_turn(code))
        {
            strcpy(sst[i],"1");
            strcpy(sst[i+1],"2");
        }    
        
    }
    
    //find helices
    int start,end,tmp_s,tmp_e,len;
    len =0;
    
    for(i=0;i<nres;i++)
    {
        //printf("%s i %d \n",codes[i],i);
        tmp_e =0;
        if(is_helix(codes[i]))
        {
            tmp_s = i;
            len+=1;
            //printf("  temp start %d len %d \n",tmp_s,len); 
            for(j=i+1;j<nres;j++)
            {
                //printf("  j %d new code %s ",j,codes[j]);                 
                if(is_helix(codes[j]))
                {
                    len++;
                    //printf("  len %d ",len);
                    i=j;
                    //printf("(i = %d)\n",i);
                    if(j==nres -1 && len >=5)
                    {
                         tmp_e=end = j;
                         start = tmp_s;
                         //printf("\nHelix - start %d end %d\n",start,end);
                    }
                }
                else
                {
                    
                     
                     if(len >=5)
                     {
                         tmp_e = end = j-1;
                         start = tmp_s;
                         //printf("\nHelix - start %d end %d\n",start,end);
                     }
                     i=j;
                     len =0;
                     //printf("Next iteration for i %d \n\n",i);
                     break;
                }
            }
        
        }
        if(tmp_e > 0)
        {
            //printf("\nHelix - start %d end %d\n",start,end);
            for(k=start;k<=end;k++)
            {
               
                strcpy(sst[k],"H");
            }
        }  
    
    }
    
    
     //find strands
    start=end=tmp_s=tmp_e=len=0;
        
    for(i=0;i<nres;i++)
    {
    
        tmp_e =0;
        if(is_strand(codes[i]))
        {
            tmp_s = i;
            len+=1;
            //printf("  temp start %d len %d \n",tmp_s,len); 
            for(j=i+1;j<nres;j++)
            {
                //printf("  j %d new code %s ",j,codes[j]);                 
                if(is_strand(codes[j]))
                {
                    len++;
                    //printf("  len %d ",len);
                    i=j;
                    //printf("(i = %d)\n",i);
                    if(j==nres -1 && len >=3)
                    {
                         tmp_e=end = j;
                         start = tmp_s;
                         //printf("\nHelix - start %d end %d\n",start,end);
                    }
                }
                else
                {
                    
                     
                     if(len >=3)
                     {
                         tmp_e = end = j-1;
                         start = tmp_s;
                         //printf("\nHelix - start %d end %d\n",start,end);
                     }
                     i=j;
                     len =0;
                     //printf("Next iteration for i %d \n\n",i);
                     break;
                }
            }
        
        }
        if(tmp_e > 0)
        {
            //printf("\nStrand - start %d end %d\n",start,end);
            for(k=start;k<=end;k++)
            {
               
                strcpy(sst[k],"E");
            }
        }  
    
    }
    /*
    for(i=0;i<nres;i++)
    {
        printf("%s",sst[i]);
    }
    printf("\n");
    */

}

