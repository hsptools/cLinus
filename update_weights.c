#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "linus.h"

/*
"""Update sampling weights from a list of weights

        Arguments

           o *W* - list of lists - each list item is a list of length
           6 containing, respectively, the helix, strand, turn1, turn2,
           coil and PII weights.

        """   

*/


void update_weights(Linus_sim *sim, double *w)
{
    double *hw = sim->protein.res_helix_wt;
    double *sw = sim->protein.res_strand_wt;
    double *t1w = sim->protein.res_turn1_wt;
    double *t2w = sim->protein.res_turn2_wt;
    double *cw = sim->protein.res_coil_wt;
    double *pw = sim->protein.res_PII_wt;
    
    int i;
    for(i=0;i<sim->protein.num_res;i++)
    {
        hw[i] = w[(i*6)+0];
        sw[i] = w[(i*6)+1];
        t1w[i] = w[(i*6)+2];
        t2w[i] = w[(i*6)+3];
        cw[i] = w[(i*6)+4];
        pw[i] = w[(i*6)+5];
    }    
    
}




