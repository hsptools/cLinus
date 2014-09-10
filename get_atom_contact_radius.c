#include <stdio.h>
#include "LinusRes.h"

double get_atom_contact_radius(char *res_name, char * atm_name, double *contact_rad)
{
    //printf("%s\n",name);
    int idx;
    idx = get_linusRes_index(res_name,atm_name);
    
    //printf("Index %d Value %lf\n",idx, sasa_values[idx]);
    *contact_rad = contact_radius_arr[idx];
}
