#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include "Atom3d.h"
#include "linus.h"

void print_all_bumps(LinusProtein *prot, char *filename)
{
    int nr = prot->num_res;
    Atom3d *atoms = prot->atoms;
    int *fai = prot->res_fai;
    char **rn = prot->res_names;
    int nb =0;
    FILE *file = fopen(filename,"w+");

    int i,j,ia,ja;

    for (i=0;i<nr-1;i++)
    {
        int i1 = fai[i];
        int i2 = fai[i+1];

        for (j=i+1;j<nr;j++)
        {
            int j1 = fai[j];
            int j2 = fai[j+1];

            if(!close(&atoms[i1+1],&atoms[j1+1],15.0))
                continue;
            for (ia=i1; ia<i2;ia++)
            {
                Atom3d *atom = &atoms[ia];

                for (ja = j1; ja< j2; ja++)
                {
                    Atom3d * jatm = &atoms[ja];

                    if(close(atom,jatm,-1.0))
                    {
                        if ((j-i) > 1)
                        {
                            if (!(strcmp(atom->name, " H  ") || strcmp(jatm->name, " O  ")) ||
                              !(strcmp(atom->name, " O  ") ||strcmp(jatm->name, " H  ")))
                            {
                                double dist=distance(atom,jatm);
                                double ovrlap = (atom->radius + jatm->radius - dist);
                                if (ovrlap < 0.55)
                                  continue;
                            }
                            else
                            {
                                nb++;
                                double dist= distance(atom,jatm);
                                fprintf(file, "%d %s %s %d %s %s Minimum = %lf Actual = %lf\n",i,rn[i],atom->name,j,rn[j],jatm->name, 0.90*(atom->radius+jatm->radius),dist);
                            }
                        } //if ((j-i) > 1)

                        else if ((atom->bsep1 + jatm->bsep0) > 3)
                        {
                            if (!(strcmp(atom->name, " H  ") || strcmp(jatm->name, " O  ")) ||
                                !(strcmp(atom->name, " H  ") || strcmp(jatm->name, " N  ")) ||
                                !(strcmp(atom->name, " N  ") || strcmp(jatm->name, " H  ")) ||
                                !(strcmp(atom->name, " O  ") || strcmp(jatm->name, " H  ")))
                            {
                                double dist=distance(atom,jatm);
                                double ovrlap = (atom->radius + jatm->radius - dist);
                                if (ovrlap < 0.55)
                                    continue;
                            }
                            else
                            {
                                double dist= distance(atom,jatm);
                                fprintf(file, "%d %s %s %d %s %s Minimum = %lf Actual = %lf\n",i,rn[i],atom->name,j,rn[j],jatm->name, 0.90*(atom->radius+jatm->radius),dist);
                            }

                        }
                    }

                }

            }
        }

    }


    fclose(file);
}

