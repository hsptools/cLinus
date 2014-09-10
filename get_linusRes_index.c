#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "LinusRes.h"

int get_linusRes_index(char *rname, char *aname)
{
    int start, end; // get residues range from linusRes
    
    get_res_range(rname,&start,&end); // get indexes for linusRes for a given residue
    
    //printf("Start %d End %d for Res %s Atom %s\n",start,end,rname,aname);
         
    char res_atm_name1[8]; // concatenate residue and atom name to search in linusRes
    strcpy(res_atm_name1,rname);
    strcat(res_atm_name1,aname);
     
    //printf("-%s-%s-%s-\n",rname,aname,res_atm_name1);    
    int i;
    int idx =-1;    
    for(i=start; i< end;i++) // Search elements for a particular residue
    {
        if(!strcmp(res_atm_name1,res_atm_name_arr[i]))
        {
            idx = i;
        }
    }
    //printf("Index: %d\n",idx);
    return idx;
         

}



