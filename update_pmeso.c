#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "linus.h"

//subtract 2 because no pmeso_weights for ACE and NME

void update_pmeso(Linus_sim *sim, double **w)
{
    int i;
    for(i=0;i<sim->protein.num_res -2;i++)
    {
        sim->protein.pmeso_wt[i]= w[i];
    }    
    
}




