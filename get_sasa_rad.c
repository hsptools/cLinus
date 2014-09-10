#include <stdio.h>
#include <string.h>
#include "LinusRes.h"
//#include "Atom3d.h"
#include "linus.h"


double get_sasa_rad(char *rname, char *aname)
{
    int i,idx;
    
    char res_atm_name1[8]; // concatenate residue and atom name to search in linusRes
    strcpy(res_atm_name1,rname);
    strcat(res_atm_name1,aname);
    
    for (i =0; i < 167; i++)
    {
        if( !strcmp(sasa_names_arr[i],res_atm_name1) )
        {
            idx =i;
        }
    }
    //printf("Index %d Value %lf\n",idx, sasa_values_arr[idx]);
    return sasa_values_arr[idx];
}


