#include <stdio.h>
#include "LinusRes.h"

void get_atom_chrg(char *res_name, char * atm_name, double *chrg)
{
    int idx;
    idx = get_linusRes_index(res_name,atm_name);
    
    //printf("Index %d Value %lf\n",idx, charge_arr[idx]);
    *chrg = charge_arr[idx];
}
