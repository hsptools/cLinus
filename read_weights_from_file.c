#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "linus.h"

/*
 """Read sampling weights from a file

        Arguments

            o *filename* - string - name of file containing the weights

            o *skip* - integer - number of leading lines to skip in file

        Description

            Weights for every residue should be specified.  Look at
            a SWT file output from the simulation for an example of
            a valid input file

        """    
*/

void read_weights_from_file(Linus_sim *sim,char *filename, int skip)
{
    //printf("%s\n",filename);
    int num = sim->protein.num_res;
    double hw[num],sw[num],t1w[num],t2w[num],cw[num],pw[num];
    char buff[BUFSIZ];
    int ctr =0,size =0;
    FILE *file = fopen(filename,"r");
    while (fgets(buff, sizeof buff, file) != NULL)
    {
        
        if(ctr > skip)
        {
            if (sscanf(buff, "%lf %lf %lf %lf %lf %lf \n",&hw[size],&sw[size],&t1w[size],&t2w[size],&cw[size],&pw[size]) == 6)
            {
                sim->protein.res_helix_wt[size] = hw[size];
                sim->protein.res_strand_wt[size] = sw[size];
                sim->protein.res_turn1_wt[size] = t1w[size];
                sim->protein.res_turn2_wt[size] = t2w[size];
                sim->protein.res_coil_wt[size] = cw[size];
                sim->protein.res_PII_wt[size] = pw[size];
                size++;
            }
        }    
           
        ctr+=1;
        
    }
    
    if(sim->protein.num_res > size+1 )
    {
        printf("File does not contain weight for all residues!");
        exit(1);
    }
    
      
}

