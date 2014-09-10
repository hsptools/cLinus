#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include "Atom3d.h"
#include "linus.h"

double compaction_score(LinusProtein *p, double ce)
{
    double comp_score;
    // ce=1.0
    double rg=radius_of_gyration(p->atoms,0,p->num_atoms);

    if(rg > ce)
    {
        comp_score = pow(rg-ce,2.0);
    }
    else
    {
        comp_score = 0.0;
    }

    return comp_score;
}

void print_compaction_score(LinusProtein *p, double cp, char *filename, double *score)
{
    FILE *file = fopen(filename,"w+");
    //fscanf(file, "%d %s %s %d %s %s Minimum = %lf Actual = %lf\n",i,rn[i],atom->name,j,rn[j],jatm->name, 0.90*(atom->radius+jatm->radius),dist);
    fprintf(file,"%s","\n------ Compaction Score -----\n");
    double rg = radius_of_gyration(p->atoms,0,p->num_atoms);
    fprintf(file,"Radius of Gyration = %lf\n",rg);
    fprintf(file,"\nCompaction Score = %lf\n",exp(-rg*rg));
    *score = rg*cp;
    fclose(file);
}

