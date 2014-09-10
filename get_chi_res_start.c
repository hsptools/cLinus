#include <stdio.h>
#include <string.h>
#include "Chidict.h"

int get_chi_res_start(char *res)
{
    int i;
    int idx = -1;    
    for(i=0; i< 20;i++) // 17 is the number of residues that have Chi values
    {
        if(!strcmp(res,chi_names[i]))
        {
            idx = i;
        }
    }
    //printf("Orig Chi size: %d\n",chi_size[idx]/2);
    return chi_start[idx];

}
