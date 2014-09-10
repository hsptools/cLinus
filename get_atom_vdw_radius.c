#include <stdio.h>
#include "LinusRes.h"

void get_atom_vdw_radius(char *res_name, char * atm_name, double *vdw_rad)
{
    //printf("%s\n",name);
    int idx;
    idx = get_linusRes_index(res_name,atm_name);
    
    //printf("Index %d Value %lf\n",idx, sasa_values[idx]);
    *vdw_rad = vdw_radius_arr[idx];
}

