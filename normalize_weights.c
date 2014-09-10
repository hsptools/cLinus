#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "linus.h"

//"Normalize each residue's sampling weights so that it sums to 1.0

void normalize_weights(Linus_sim *sim)
{
    double *hw = sim->protein.res_helix_wt;
    double *sw = sim->protein.res_strand_wt;
    double *t1w = sim->protein.res_turn1_wt;
    double *t2w = sim->protein.res_turn2_wt;
    double *cw = sim->protein.res_coil_wt;
    double *pw = sim->protein.res_PII_wt;

    int i;
    double norm;
    for(i=0;i< sim->protein.num_res;i++)
    {
        if((hw[i] + sw[i] + t1w[1] + t2w[i] + cw[i]) > 0)
        {
            norm =  1.0 / (hw[i] + sw[i] + t1w[i] + t2w[i] + cw[i]);
            hw[i] = hw[i] * norm;
            sw[i] = sw[i] * norm;
            cw[i] = cw[i] * norm;
            pw[i] = pw[i] * norm;
            t1w[i] = t1w[i] * norm;
            t2w[i] = t2w[i] * norm;
            printf("%lf %lf %lf %lf %lf %lf\n", hw[i],sw[i], cw[i],pw[i], t1w[i],t2w[i]);
        }    
    
    }
    
}

