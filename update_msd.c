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

void update_msd(Linus_sim *sim,double *fm, double *fsd,double *ym, double *ysd)
{
    int i=0;
    
    for(i=0;i<sim->protein.num_res;i++)
    {
        sim->protein.res_f_mean[i] =fm[i];
        sim->protein.res_f_sd[i] =fsd[i];
        sim->protein.res_y_mean[i] =ym[i];
        sim->protein.res_y_sd[i] =ysd[i];
    }
}

