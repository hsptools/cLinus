#include <stdio.h>
#include <string.h>
#include "Chidict.h"

int get_chi_res_index(char *res)
{
    int i;
    int idx = -1;    
    for(i=0; i< 17;i++) // 17 is the number of residues that have Chi values
    {
        if(!strcmp(res,chi_names[i]))
        {
            idx = i;
        }
    }
    //printf("Index: %d\n",idx);
    return idx;

}

