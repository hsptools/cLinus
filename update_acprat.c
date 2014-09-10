#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "linus.h"

/*
"""Update phi/psi means and S.D.

        Arguments

           o *fmsd* - phi mean and stddev
           o *ymsd* - psi mean and stddev

        """
*/

void update_acprat(Linus_sim *sim,double *acprat)
{
    int i=0;
    
    for(i=0;i<sim->protein.num_res;i++)
    {
        sim->protein.res_acprat[i] =acprat[i];
        
    }
}




